#!/bin/bash
###############################################################################
# convert_dosage_vcf_to_pgen.sh
#
# Backend: Convert imputed VCF to PGEN, then filter to keep-list subjects.
#
# Strategy (avoids VCF sample ID format issues):
#   1. Convert VCF to PGEN (all samples, no --keep)
#   2. Check the resulting .psam to see how PLINK2 named the samples
#   3. Build a keeplist that matches the .psam format
#   4. Filter the PGEN to keep only those subjects
#
# Arguments:
#   $1 = Directory containing imputed VCF files
#   $2 = Output directory for PGEN files
#   $3 = Keep list file (FID IID from .fam, space-separated)
#   $4 = (optional) .fam file for sex update
#   $5 = (optional) genome build: "hg19" (default) or "hg38" — used for chrX --split-par
###############################################################################

set -uo pipefail

vcf_dir="$1"
pgen_dir="$2"
keeplist="$3"
famfile="${4:-}"        # Optional: .fam file for sex update
genome_build="${5:-hg19}"  # Optional: hg19 (default) or hg38 for --split-par

echo "  Converting VCF to PGEN format..."
echo "  Input:  ${vcf_dir}"
echo "  Output: ${pgen_dir}"
echo "  Keep:   ${keeplist}"
[[ -n "${famfile}" ]] && echo "  .fam:   ${famfile}"
echo "  Genome build: ${genome_build} (used for chrX --split-par)"
echo ""

#=============================================================================
# Check PLINK2 version — must be recent (>=v2.00a3 or stable v2.0)
# Old alpha builds (v2.00a2LM, Jan 2019) have a "Malformed .pgen" bug
# when converting large VCFs.
#=============================================================================
PLINK2_CMD="plink2"
plink2_version=$("${PLINK2_CMD}" --version 2>/dev/null | head -1 || echo "unknown")
echo "  PLINK2 version: ${plink2_version}"

# Check if version is too old (2019 or earlier alpha)
needs_update=0
if echo "${plink2_version}" | grep -qiE "2\.00a[12]|201[0-9]|unknown"; then
    needs_update=1
fi

if [[ ${needs_update} -eq 1 ]]; then
    echo ""
    echo "  ╔═══════════════════════════════════════════════════════════════╗"
    echo "  ║  WARNING: Your PLINK2 version is too old.                    ║"
    echo "  ║  '${plink2_version}'"
    echo "  ║  has known bugs with VCF-to-PGEN conversion.                 ║"
    echo "  ║  Downloading a current PLINK2 binary...                      ║"
    echo "  ╚═══════════════════════════════════════════════════════════════╝"
    echo ""

    local_plink2="${pgen_dir}/_plink2_bin"
    mkdir -p "${local_plink2}"

    if [[ ! -f "${local_plink2}/plink2" ]]; then
        echo "  Downloading PLINK2 from cog-genomics.org..."
        wget -q "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20261007.zip" \
            -O "${local_plink2}/plink2.zip" 2>&1 || \
        wget -q "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip" \
            -O "${local_plink2}/plink2.zip" 2>&1 || \
        curl -sL "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip" \
            -o "${local_plink2}/plink2.zip" 2>&1

        if [[ -f "${local_plink2}/plink2.zip" ]]; then
            unzip -o -q "${local_plink2}/plink2.zip" -d "${local_plink2}/" 2>&1
            chmod +x "${local_plink2}/plink2"
            rm -f "${local_plink2}/plink2.zip"
        fi
    fi

    if [[ -f "${local_plink2}/plink2" ]]; then
        PLINK2_CMD="${local_plink2}/plink2"
        plink2_new_ver=$("${PLINK2_CMD}" --version 2>/dev/null | head -1 || echo "unknown")
        echo "  [OK] Using downloaded PLINK2: ${plink2_new_ver}"
    else
        echo "  [ERROR] Could not download PLINK2."
        echo "          Please update PLINK2 manually: https://www.cog-genomics.org/plink/2.0/"
        echo "          Your current version has a known 'Malformed .pgen' bug."
        exit 1
    fi
    echo ""
fi

n_success=0
n_fail=0
n_skip=0
adapted_keeplist=""  # Will be built after first chromosome

#=============================================================================
# Build adapted keeplist after first successful VCF import
#=============================================================================
build_adapted_keeplist() {
    local psam_file="$1"
    local out_keep="${pgen_dir}/_adapted_keeplist.txt"

    echo "  --- Building adapted keeplist from .psam ---"

    # With --double-id, .psam has FID=IID=VCF_sample_name
    # Read the actual sample IDs from the .psam
    local all_psam_ids
    all_psam_ids=$(awk 'NR>1 && !/^#/ {print $1}' "${psam_file}")
    echo "  First .psam IDs: $(echo "${all_psam_ids}" | head -3 | tr '\n' ' ')"

    # Read keeplist
    local keep_fid keep_iid
    keep_fid=$(head -1 "${keeplist}" | awk '{print $1}')
    keep_iid=$(head -1 "${keeplist}" | awk '{print $2}')
    echo "  Keeplist first: FID='${keep_fid}' IID='${keep_iid}'"

    # Test which keeplist format matches .psam sample IDs
    local n_match_iid n_match_fid_iid n_match_fid

    n_match_iid=$(awk '{print $2}' "${keeplist}" | sort -u | \
                  comm -12 - <(echo "${all_psam_ids}" | sort -u) | wc -l)

    n_match_fid_iid=$(awk '{print $1"_"$2}' "${keeplist}" | sort -u | \
                      comm -12 - <(echo "${all_psam_ids}" | sort -u) | wc -l)

    n_match_fid=$(awk '{print $1}' "${keeplist}" | sort -u | \
                  comm -12 - <(echo "${all_psam_ids}" | sort -u) | wc -l)

    echo "  Match tests:"
    echo "    Keeplist IID:     ${n_match_iid}"
    echo "    Keeplist FID_IID: ${n_match_fid_iid}"
    echo "    Keeplist FID:     ${n_match_fid}"

    # Pick best match
    local best_n=0 best_fmt=""
    if [[ ${n_match_fid_iid} -gt ${best_n} ]]; then best_n=${n_match_fid_iid}; best_fmt="fid_iid"; fi
    if [[ ${n_match_iid} -gt ${best_n} ]]; then best_n=${n_match_iid}; best_fmt="iid"; fi
    if [[ ${n_match_fid} -gt ${best_n} ]]; then best_n=${n_match_fid}; best_fmt="fid"; fi

    if [[ ${best_n} -eq 0 ]]; then
        echo "  [ERROR] No keeplist format matched .psam IDs."
        echo "          Please reconcile .fam and VCF sample IDs."
        cp "${keeplist}" "${out_keep}"
        adapted_keeplist="${out_keep}"
        return
    fi

    echo "  [OK] Best: ${best_fmt} (${best_n} matches)"

    # Build keeplist with correct FID/IID for PLINK2 --keep
    case "${best_fmt}" in
        fid_iid) awk '{id=$1"_"$2; print id, id}' "${keeplist}" > "${out_keep}" ;;
        iid)     awk '{print $2, $2}' "${keeplist}" > "${out_keep}" ;;
        fid)     awk '{print $1, $1}' "${keeplist}" > "${out_keep}" ;;
    esac

    # Also build an ID update file to restore original .fam FID/IID
    # Format: old_FID old_IID new_FID new_IID
    local update_ids_file="${pgen_dir}/_update_ids.txt"
    case "${best_fmt}" in
        fid_iid) awk '{id=$1"_"$2; print id, id, $1, $2}' "${keeplist}" > "${update_ids_file}" ;;
        iid)     awk '{print $2, $2, $1, $2}' "${keeplist}" > "${update_ids_file}" ;;
        fid)     awk '{print $1, $1, $1, $2}' "${keeplist}" > "${update_ids_file}" ;;
    esac
    echo "  ID update file: ${update_ids_file} (maps VCF IDs -> .fam FID/IID)"
    echo "  First mappings:"
    head -3 "${update_ids_file}" | sed 's/^/    /'

    local n_keep
    n_keep=$(wc -l < "${out_keep}")
    echo "  Adapted keeplist: ${n_keep} subjects"
    echo "  First entries:"
    head -3 "${out_keep}" | sed 's/^/    /'
    echo ""

    adapted_keeplist="${out_keep}"
}

#=============================================================================
# Convert and filter one chromosome
#=============================================================================
process_chr() {
    local chr="$1"
    local vcf_file="$2"
    local sex_file="${3:-}"   # Optional: sex update file for chrX --split-par
    local final_out="${pgen_dir}/chr${chr}"

    # Skip if final output exists
    if [[ -f "${final_out}.pgen" ]]; then
        echo "  [DONE] chr${chr}: PGEN already exists, skipping"
        n_success=$((n_success + 1))
        return
    fi

    echo "  [RUN]  chr${chr}: $(basename ${vcf_file})"

    local tmp_out="${pgen_dir}/_tmp_chr${chr}"

    # Step 1: Convert VCF to PGEN (all samples)
    # --double-id sets FID=IID=VCF_sample_name, avoiding FID/IID ambiguity
    # For chrX: add --split-par to split pseudoautosomal regions (PAR1/PAR2)
    local extra_flags=()
    if [[ "${chr}" == "X" ]] || [[ "${chr}" == "23" ]]; then
        extra_flags+=(--split-par "${genome_build}")
        echo "  [INFO] chrX: adding --split-par ${genome_build}"
        # chrX with --split-par REQUIRES sex information at import time
        if [[ -n "${sex_file}" ]] && [[ -f "${sex_file}" ]]; then
            extra_flags+=(--update-sex "${sex_file}")
            echo "  [INFO] chrX: adding --update-sex ${sex_file}"
        else
            echo "  [WARN] chrX: no sex file provided — --split-par may fail"
        fi
    fi

    if ! ${PLINK2_CMD} --vcf "${vcf_file}" dosage=DS \
                --double-id \
                "${extra_flags[@]+"${extra_flags[@]}"}" \
                --make-pgen \
                --out "${tmp_out}" \
                --memory 4000 \
                --threads 2 \
                2>&1; then
        echo "  [FAIL] chr${chr}: VCF conversion failed"
        rm -f "${tmp_out}".{pgen,pvar,psam,log}
        n_fail=$((n_fail + 1))
        return
    fi

    # Step 2: Build adapted keeplist from .psam (only once)
    if [[ -z "${adapted_keeplist}" ]]; then
        build_adapted_keeplist "${tmp_out}.psam"
    fi

    # Step 3: Filter to keep-list subjects AND update IDs to match .fam
    local filtered_out="${pgen_dir}/_filtered_chr${chr}"

    if ${PLINK2_CMD} --pfile "${tmp_out}" \
              --keep "${adapted_keeplist}" \
              --make-pgen \
              --out "${filtered_out}" \
              --memory 4000 \
              --threads 2 \
              2>&1; then

        # Step 4: Update sample IDs from VCF format to .fam FID/IID
        local update_ids_file="${pgen_dir}/_update_ids.txt"
        if [[ -f "${update_ids_file}" ]]; then
            if ${PLINK2_CMD} --pfile "${filtered_out}" \
                      --update-ids "${update_ids_file}" \
                      --make-pgen \
                      --out "${final_out}" \
                      --memory 4000 \
                      --threads 2 \
                      2>&1; then
                local n_kept
                n_kept=$(awk 'NR>1 && !/^#/' "${final_out}.psam" | wc -l)
                echo "  [OK]   chr${chr}: ${n_kept} samples, IDs updated to .fam format"
                n_success=$((n_success + 1))
            else
                echo "  [WARN] chr${chr}: ID update failed, keeping VCF IDs"
                mv "${filtered_out}.pgen" "${final_out}.pgen"
                mv "${filtered_out}.pvar" "${final_out}.pvar"
                mv "${filtered_out}.psam" "${final_out}.psam"
                n_success=$((n_success + 1))
            fi
            rm -f "${filtered_out}".{pgen,pvar,psam,log}
        else
            # No update file — just rename filtered to final
            mv "${filtered_out}.pgen" "${final_out}.pgen"
            mv "${filtered_out}.pvar" "${final_out}.pvar"
            mv "${filtered_out}.psam" "${final_out}.psam"
            local n_kept
            n_kept=$(awk 'NR>1 && !/^#/' "${final_out}.psam" | wc -l)
            echo "  [OK]   chr${chr}: ${n_kept} samples kept"
            n_success=$((n_success + 1))
        fi
    else
        echo "  [FAIL] chr${chr}: sample filtering failed"
        n_fail=$((n_fail + 1))
    fi

    # Clean up temp files
    rm -f "${tmp_out}".{pgen,pvar,psam,log}
}

#=============================================================================
# Process autosomes (chr1-22)
#=============================================================================
for chr in $(seq 1 22); do
    vcf_file=""

    for pattern in \
        "${vcf_dir}/chr${chr}.dose.vcf.gz" \
        "${vcf_dir}/${chr}.dose.vcf.gz" \
        "${vcf_dir}/chr${chr}.vcf.gz"; do
        if [[ -f "${pattern}" ]]; then
            vcf_file="${pattern}"
            break
        fi
    done

    if [[ -z "${vcf_file}" ]]; then
        for candidate in "${vcf_dir}"/*.vcf.gz; do
            [[ ! -f "${candidate}" ]] && continue
            bname=$(basename "${candidate}" .vcf.gz)
            if [[ "${bname}" =~ chr${chr}$ ]] || [[ "${bname}" =~ chr0*${chr}$ ]]; then
                vcf_file="${candidate}"
                break
            fi
        done
    fi

    if [[ -z "${vcf_file}" ]]; then
        echo "  [SKIP] chr${chr}: no VCF file found"
        n_skip=$((n_skip + 1))
        continue
    fi

    process_chr "${chr}" "${vcf_file}"
done

#=============================================================================
# Process chromosome X
# chrX with --split-par requires sex information at VCF import time.
# We build an adapted sex file BEFORE calling process_chr for chrX.
#=============================================================================
vcf_file_x=""
for pattern in \
    "${vcf_dir}/chrX.dose.vcf.gz" \
    "${vcf_dir}/X.dose.vcf.gz" \
    "${vcf_dir}/chr23.dose.vcf.gz" \
    "${vcf_dir}/chrX.vcf.gz"; do
    if [[ -f "${pattern}" ]]; then
        vcf_file_x="${pattern}"
        break
    fi
done
if [[ -z "${vcf_file_x}" ]]; then
    for candidate in "${vcf_dir}"/*.vcf.gz; do
        [[ ! -f "${candidate}" ]] && continue
        bname=$(basename "${candidate}" .vcf.gz)
        if [[ "${bname}" =~ chrX$ ]] || [[ "${bname}" =~ chr23$ ]]; then
            vcf_file_x="${candidate}"
            break
        fi
    done
fi
if [[ -n "${vcf_file_x}" ]]; then
    chrx_sex_file=""

    if [[ -n "${famfile}" ]] && [[ -f "${famfile}" ]]; then
        echo ""
        echo "  --- Building sex file for chrX import ---"
        chrx_sex_file="${pgen_dir}/_chrx_update_sex.txt"
        update_ids_file="${pgen_dir}/_update_ids.txt"

        if [[ -f "${update_ids_file}" ]]; then
            # Use _update_ids.txt (built during autosome processing) to map
            # VCF sample IDs -> .fam FID/IID, then look up sex from .fam.
            # _update_ids.txt format: vcf_fid vcf_iid fam_fid fam_iid
            # .fam format:           fam_fid fam_iid father mother sex pheno
            # Output:                vcf_fid vcf_iid sex
            awk '
            NR == FNR {
                sex[$1 SUBSEP $2] = $5
                next
            }
            {
                key = $3 SUBSEP $4
                if (key in sex)
                    print $1, $2, sex[key]
                else
                    print $1, $2, "NA"
            }
            ' "${famfile}" "${update_ids_file}" > "${chrx_sex_file}"

            n_sex=$(awk '$3==1 || $3==2' "${chrx_sex_file}" | wc -l)
            n_total=$(wc -l < "${chrx_sex_file}")
            echo "  [OK] Built from _update_ids.txt: ${n_sex}/${n_total} with valid sex"
        else
            # Fallback: extract VCF sample names directly with bcftools
            # and match against .fam using the same strategies as the keeplist
            echo "  [INFO] No _update_ids.txt yet — reading VCF sample names with bcftools"
            vcf_samples_file="${pgen_dir}/_chrx_vcf_samples.txt"
            bcftools query -l "${vcf_file_x}" > "${vcf_samples_file}" 2>/dev/null

            if [[ -s "${vcf_samples_file}" ]]; then
                # Try matching VCF IDs to .fam IID, FID_IID, or FID
                python3 -I -c "
import sys

fam_by_iid = {}
fam_by_fid = {}
fam_by_fid_iid = {}

with open('${famfile}') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 5:
            fid, iid, sex = parts[0], parts[1], parts[4]
            fam_by_iid[iid] = sex
            fam_by_fid[fid] = sex
            fam_by_fid_iid[fid + '_' + iid] = sex

with open('${vcf_samples_file}') as f:
    vcf_ids = [line.strip() for line in f if line.strip()]

# Count matches for each strategy
matches = {
    'iid':     sum(1 for v in vcf_ids if v in fam_by_iid),
    'fid_iid': sum(1 for v in vcf_ids if v in fam_by_fid_iid),
    'fid':     sum(1 for v in vcf_ids if v in fam_by_fid),
}

best = max(matches, key=matches.get)
if matches[best] == 0:
    print('  [WARN] No VCF IDs matched .fam — cannot assign sex', file=sys.stderr)
    sys.exit(1)

lookup = {'iid': fam_by_iid, 'fid_iid': fam_by_fid_iid, 'fid': fam_by_fid}[best]
matched = 0

with open('${chrx_sex_file}', 'w') as out:
    for vid in vcf_ids:
        sex = lookup.get(vid, 'NA')
        if sex in ('1', '2'):
            out.write(f'{vid}\t{vid}\t{sex}\n')
            matched += 1
        else:
            out.write(f'{vid}\t{vid}\tNA\n')

print(f'  Match strategy: {best} ({matches[best]} matches)')
print(f'  Sex assigned: {matched}/{len(vcf_ids)} samples')
"
                if [[ $? -ne 0 ]]; then
                    echo "  [WARN] Could not match VCF IDs to .fam sex values"
                    chrx_sex_file=""
                else
                    echo "  [OK] Built chrX sex file from VCF sample names"
                fi
            else
                echo "  [WARN] bcftools query -l failed — cannot build sex file"
                chrx_sex_file=""
            fi
            rm -f "${vcf_samples_file}"
        fi

        # Verify the sex file has valid entries
        if [[ -n "${chrx_sex_file}" ]] && [[ -f "${chrx_sex_file}" ]]; then
            n_valid=$(awk '$3==1 || $3==2' "${chrx_sex_file}" | wc -l)
            if [[ ${n_valid} -eq 0 ]]; then
                echo "  [WARN] Sex file has 0 valid entries — chrX --split-par will likely fail"
                chrx_sex_file=""
            else
                echo "  First entries:"
                head -3 "${chrx_sex_file}" | sed 's/^/    /'
            fi
        fi
    else
        echo "  [WARN] No .fam file — chrX --split-par requires sex info and will likely fail"
    fi

    process_chr "X" "${vcf_file_x}" "${chrx_sex_file}"
fi

#=============================================================================
# Update SEX in PGEN files from the .fam
# This ALWAYS runs, even if PGENs were already converted (in case sex was missing)
#=============================================================================
echo ""
echo "  --- Updating SEX in PGEN files ---"

if [[ -n "${famfile}" ]] && [[ -f "${famfile}" ]]; then
    # Create sex update file from the .fam (FID IID SEX)
    # .fam columns: FID IID PID MID SEX PHENO
    sex_update_file="${pgen_dir}/_update_sex.txt"
    awk '{print $1, $2, $5}' "${famfile}" > "${sex_update_file}"
    n_entries=$(wc -l < "${sex_update_file}")
    echo "  Sex update file: ${sex_update_file} (${n_entries} entries)"
    echo "  First entries:"
    head -3 "${sex_update_file}" | sed 's/^/    /'

    # Check if sex values look valid (should be 1 or 2)
    n_male=$(awk '$3==1' "${sex_update_file}" | wc -l)
    n_female=$(awk '$3==2' "${sex_update_file}" | wc -l)
    n_other=$(awk '$3!=1 && $3!=2' "${sex_update_file}" | wc -l)
    echo "  .fam sex: ${n_male} male (1), ${n_female} female (2), ${n_other} other/missing"

    if [[ ${n_male} -eq 0 ]] && [[ ${n_female} -eq 0 ]]; then
        echo "  [WARN] No valid sex values in .fam (column 5). Skipping sex update."
    else
        # Directly edit .psam files (they're just tab-separated text)
        # No plink2 needed — faster and avoids --update-sex format issues
        n_sex_updated=0
        for psam_file in "${pgen_dir}"/chr*.psam; do
            [[ ! -f "${psam_file}" ]] && continue
            chr_name=$(basename "${psam_file}" .psam)

            # Read the .psam, match FID+IID to .fam sex values, write updated file
            awk -F'\t' -v OFS='\t' '
            BEGIN {
                # Read .fam file to get sex by FID+IID
                while ((getline line < "'"${famfile}"'") > 0) {
                    n = split(line, a)
                    if (n >= 5) {
                        key = a[1] "\t" a[2]
                        sex[key] = a[5]
                    }
                }
            }
            /^#/ { print; next }     # Print header as-is
            NR == 1 { print; next }  # Safety: print first line
            {
                key = $1 "\t" $2
                if (key in sex) {
                    # Print FID, IID, sex (and any remaining columns)
                    printf "%s\t%s\t%s", $1, $2, sex[key]
                    for (i = 4; i <= NF; i++) printf "\t%s", $i
                    printf "\n"
                } else {
                    # No match — keep original (print all fields)
                    printf "%s\t%s\tNA", $1, $2
                    for (i = 4; i <= NF; i++) printf "\t%s", $i
                    printf "\n"
                }
            }
            ' "${psam_file}" > "${psam_file}.tmp"

            # Verify the temp file has content before replacing
            if [[ -s "${psam_file}.tmp" ]]; then
                mv "${psam_file}.tmp" "${psam_file}"
                n_sex_updated=$((n_sex_updated + 1))
            else
                echo "  [WARN] Sex update produced empty file for ${chr_name}, keeping original"
                rm -f "${psam_file}.tmp"
            fi
        done
        echo "  Sex updated in ${n_sex_updated} .psam files"
    fi

    # Verify
    first_psam=$(ls "${pgen_dir}"/chr*.psam 2>/dev/null | head -1)
    if [[ -n "${first_psam}" ]]; then
        echo "  Verification ($(basename ${first_psam})):"
        head -5 "${first_psam}" | sed 's/^/    /'
    fi
else
    echo "  [SKIP] No .fam file provided. Sex will not be updated in .psam files."
    echo "         PLINK2 will report '0 males, 0 females' and sex-stratified GWAS may fail."
fi

#=============================================================================
# Summary
#=============================================================================
echo ""
echo "  VCF to PGEN conversion summary:"
echo "    Successful: ${n_success}"
echo "    Failed:     ${n_fail}"
echo "    Skipped:    ${n_skip}"

n_pgen=$(ls "${pgen_dir}"/*.pgen 2>/dev/null | wc -l)
echo "    PGEN files: ${n_pgen}"

if [[ ${n_fail} -gt 0 ]]; then
    echo ""
    echo "  [WARNING] ${n_fail} chromosome(s) failed."
    echo "  You can re-run — it skips already-converted chromosomes."
fi
