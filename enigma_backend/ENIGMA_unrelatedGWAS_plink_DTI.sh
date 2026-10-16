#!/bin/bash
###############################################################################
# ENIGMA_unrelatedGWAS_plink_DTI.sh
#
# Backend: Generate and optionally run PLINK2 GWAS commands.
# Called by run2 — do not run directly.
#
# Key design decisions:
#   - MAC filter = 10 (not 20)
#   - Sex-stratified models INCLUDE AffectionStatus
#   - Metrics detected dynamically from phenotype file header
###############################################################################

set -euo pipefail

echo "=== ENIGMA-DTI GWAS Backend ==="

#=============================================================================
# Detect phenotypes from header
#=============================================================================
header=$(head -1 "${phenoFILE}")

# Detect which metrics are present
declare -a METRICS=()
for m in FA MD AD RD; do
    if echo "${header}" | tr '\t' '\n' | grep -q "^${m}_"; then
        METRICS+=("${m}")
    fi
done
echo "  Detected metrics: ${METRICS[*]}"

# ROI names (bilateral averages)
ROIS="GlobalAverage ACR ALIC BCC CC CGC CGH CR CST EC FX FXST GCC IC PCR PLIC PTR RLIC SCC SCR SFO SLF SS TAP UNC"

#=============================================================================
# Determine subgroups (read directly from covar file, not .ids files)
#=============================================================================
echo "--- Determining GWAS models ---"

n_total=$(( $(wc -l < "${covarFILE}") - 1 ))

# Find column positions in the covar file
covar_header_line=$(head -1 "${covarFILE}")
sex_col=$(echo "${covar_header_line}" | tr '\t' '\n' | { grep -nx "^Sex$" || true; } | cut -d: -f1)
aff_col=$(echo "${covar_header_line}" | tr '\t' '\n' | { grep -nx "^AffectionStatus$" || true; } | cut -d: -f1)

# Count from the covar file directly
if [[ -n "${sex_col}" ]]; then
    n_males=$(awk -F'\t' -v c="${sex_col}" 'NR>1 && $c==1' "${covarFILE}" | wc -l)
    n_females=$(awk -F'\t' -v c="${sex_col}" 'NR>1 && $c==2' "${covarFILE}" | wc -l)
else
    n_males=0; n_females=0
    echo "  [WARN] No 'Sex' column found in covar file"
fi

if [[ -n "${aff_col}" ]]; then
    n_controls=$(awk -F'\t' -v c="${aff_col}" 'NR>1 && $c==1' "${covarFILE}" | wc -l)
    n_cases=$(awk -F'\t' -v c="${aff_col}" 'NR>1 && $c==2' "${covarFILE}" | wc -l)
else
    n_controls=0; n_cases=0
    echo "  [WARN] No 'AffectionStatus' column found in covar file"
fi

echo "  Total: ${n_total} | Controls: ${n_controls} | Cases: ${n_cases} | Males: ${n_males} | Females: ${n_females}"

# Generate .ids files from covar (in case they don't exist or have wrong prefix)
ids_dir="${gwas_outdir}"
if [[ -n "${sex_col}" ]]; then
    awk -F'\t' -v c="${sex_col}" 'NR>1 && $c==1 {print $1, $2}' "${covarFILE}" > "${ids_dir}/${COHORT}_males.ids"
    awk -F'\t' -v c="${sex_col}" 'NR>1 && $c==2 {print $1, $2}' "${covarFILE}" > "${ids_dir}/${COHORT}_females.ids"
fi
if [[ -n "${aff_col}" ]]; then
    awk -F'\t' -v c="${aff_col}" 'NR>1 && $c==1 {print $1, $2}' "${covarFILE}" > "${ids_dir}/${COHORT}_Conly.ids"
    awk -F'\t' -v c="${aff_col}" 'NR>1 && $c==2 {print $1, $2}' "${covarFILE}" > "${ids_dir}/${COHORT}_Aonly.ids"
fi

do_Conly=0; do_Aonly=0; do_males=0; do_females=0

[[ ${n_controls} -ge 100 ]] && [[ ${n_cases} -ge 1 ]]  && { do_Conly=1; echo "  [YES] Controls-only"; }
[[ ${n_cases} -ge 100 ]]                                 && { do_Aonly=1; echo "  [YES] Cases-only"; }
[[ ${n_males} -ge 100 ]]                                  && { do_males=1; echo "  [YES] Males-only"; }
[[ ${n_females} -ge 100 ]]                                && { do_females=1; echo "  [YES] Females-only"; }

#=============================================================================
# Parse covariate columns (from COVAR file header, not pheno)
#=============================================================================
covar_header=$(head -1 "${covarFILE}")

# Determine if AffectionStatus should be included
# (only if there are BOTH cases and controls — otherwise it's constant and causes VIF errors)
include_affection=""
if [[ ${n_cases} -gt 0 ]] && [[ ${n_controls} -gt 0 ]]; then
    include_affection="AffectionStatus"
    echo "  AffectionStatus: included (${n_controls} controls, ${n_cases} cases)"
else
    echo "  AffectionStatus: EXCLUDED (constant — all subjects same status)"
fi

covar_base="Sex Age AgeCsq AgexSex AgeCsqxSex ${include_affection}"
# Sex-stratified: NO sex or sex-interactions, but YES AffectionStatus (if applicable)
covar_nosex="Age AgeCsq ${include_affection}"
covar_noaff="Sex Age AgeCsq AgexSex AgeCsqxSex"

# Clean up whitespace from empty include_affection
covar_base=$(echo ${covar_base} | xargs)
covar_nosex=$(echo ${covar_nosex} | xargs)

pc_covars=$(echo "${covar_header}" | tr '\t' '\n' | { grep "^PC[0-9]" || true; } | tr '\n' ' ')
echo "  PCs: ${pc_covars:-<none>}"

# Detect scanner PCA covariates (scanPC1, scanPC2, ...)
scanpc_covars=$(echo "${covar_header}" | tr '\t' '\n' | { grep "^scanPC[0-9]" || true; } | tr '\n' ' ')
if [[ -n "${scanpc_covars}" ]]; then
    echo "  Scanner PCs: ${scanpc_covars}"
fi

# Identify dummy covariates (anything in covar file not standard)
standard="FID IID Sex Age AgeCsq AgexSex AgeCsqxSex AffectionStatus ${pc_covars} ${scanpc_covars}"
for m in ${METRICS[*]}; do standard="${standard} ${m}_GlobalAverage"; done

dummy_covars=""
for col in $(echo "${covar_header}" | tr '\t' '\n'); do
    is_std=0
    for s in ${standard}; do [[ "${col}" == "${s}" ]] && { is_std=1; break; }; done
    [[ ${is_std} -eq 0 ]] && dummy_covars="${dummy_covars} ${col}"
done
dummy_covars=$(echo "${dummy_covars}" | xargs)
[[ -n "${dummy_covars}" ]] && echo "  Dummies: ${dummy_covars}" || echo "  Dummies: <none>"

#=============================================================================
# Build chromosome list
#=============================================================================
chr_list=""
for c in $(seq ${chrStart} ${chrEnd}); do chr_list="${chr_list} ${c}"; done
[[ ${doChrX} -eq 1 ]] && chr_list="${chr_list} X"

#=============================================================================
# Command generator
#=============================================================================
gen_cmd() {
    local pheno=$1 covars=$2 chr=$3 outname=$4 keep=$5
    local cmd="plink2 --pfile ${pgen_dir}/chr${chr}"
    cmd="${cmd} --pheno ${phenoFILE} --pheno-name ${pheno}"
    cmd="${cmd} --covar ${covarFILE} --covar-name ${covars}"
    cmd="${cmd} --glm omit-ref hide-covar --covar-variance-standardize"
    cmd="${cmd} --mac 10 --maf 0.01"
    [[ -n "${keep}" ]] && cmd="${cmd} --keep ${keep}"
    cmd="${cmd} --out ${gwas_outdir}/sumstats/${outname}"
    echo "${cmd}"
}

#=============================================================================
# Generate commands
#=============================================================================
echo ""
echo "--- Generating GWAS commands ---"
n_gwas=0

for mf in fullcohort_wGlobal fullcohort_woGlobal Conly Aonly males females; do
    > "${gwas_outdir}/GWAS_commands_${mf}.txt"
done

for metric in ${METRICS[*]}; do
    global="${metric}_GlobalAverage"

    for roi in ${ROIS}; do
        pheno="${metric}_${roi}"
        # Verify phenotype exists in header
        if ! echo "${header}" | tr '\t' '\n' | grep -qx "${pheno}"; then
            continue
        fi

        echo "  Generating commands for: ${pheno}"

        for chr in ${chr_list}; do
            [[ ! -f "${pgen_dir}/chr${chr}.pgen" ]] && continue

            # Model 1: Full cohort WITH global
            covars="${global} ${covar_base} ${pc_covars} ${scanpc_covars} ${dummy_covars}"
            out="${COHORT}_${ancestry}_ENIGMA-DTI_fullcohort_w${metric}_${pheno}_chr${chr}_${date}_${analyst}"
            echo "$(gen_cmd "${pheno}" "$(echo ${covars}|xargs)" "${chr}" "${out}" "")" >> "${gwas_outdir}/GWAS_commands_fullcohort_wGlobal.txt"
            echo "gzip -f ${gwas_outdir}/sumstats/${out}.${pheno}.glm.linear 2>/dev/null || true" >> "${gwas_outdir}/GWAS_commands_fullcohort_wGlobal.txt"
            n_gwas=$((n_gwas+1))

            # Model 6: Full cohort WITHOUT global
            covars_wo="${covar_base} ${pc_covars} ${scanpc_covars} ${dummy_covars}"
            out_wo="${COHORT}_${ancestry}_ENIGMA-DTI_fullcohort_wo_${pheno}_chr${chr}_${date}_${analyst}"
            echo "$(gen_cmd "${pheno}" "$(echo ${covars_wo}|xargs)" "${chr}" "${out_wo}" "")" >> "${gwas_outdir}/GWAS_commands_fullcohort_woGlobal.txt"
            echo "gzip -f ${gwas_outdir}/sumstats/${out_wo}.${pheno}.glm.linear 2>/dev/null || true" >> "${gwas_outdir}/GWAS_commands_fullcohort_woGlobal.txt"
            n_gwas=$((n_gwas+1))

            # Model 2: Controls-only
            if [[ ${do_Conly} -eq 1 ]]; then
                covars_c="${global} ${covar_noaff} ${pc_covars} ${scanpc_covars} ${dummy_covars}"
                out_c="${COHORT}_${ancestry}_ENIGMA-DTI_Conly_w${metric}_${pheno}_chr${chr}_${date}_${analyst}"
                echo "$(gen_cmd "${pheno}" "$(echo ${covars_c}|xargs)" "${chr}" "${out_c}" "${ids_dir}/${COHORT}_Conly.ids")" >> "${gwas_outdir}/GWAS_commands_Conly.txt"
                echo "gzip -f ${gwas_outdir}/sumstats/${out_c}.${pheno}.glm.linear 2>/dev/null || true" >> "${gwas_outdir}/GWAS_commands_Conly.txt"
                n_gwas=$((n_gwas+1))
            fi

            # Model 3: Cases-only
            if [[ ${do_Aonly} -eq 1 ]]; then
                covars_a="${global} ${covar_noaff} ${pc_covars} ${scanpc_covars} ${dummy_covars}"
                out_a="${COHORT}_${ancestry}_ENIGMA-DTI_Aonly_w${metric}_${pheno}_chr${chr}_${date}_${analyst}"
                echo "$(gen_cmd "${pheno}" "$(echo ${covars_a}|xargs)" "${chr}" "${out_a}" "${ids_dir}/${COHORT}_Aonly.ids")" >> "${gwas_outdir}/GWAS_commands_Aonly.txt"
                echo "gzip -f ${gwas_outdir}/sumstats/${out_a}.${pheno}.glm.linear 2>/dev/null || true" >> "${gwas_outdir}/GWAS_commands_Aonly.txt"
                n_gwas=$((n_gwas+1))
            fi

            # Model 4: Males-only (includes AffectionStatus, no sex terms)
            if [[ ${do_males} -eq 1 ]]; then
                covars_m="${global} ${covar_nosex} ${pc_covars} ${scanpc_covars} ${dummy_covars}"
                out_m="${COHORT}_${ancestry}_ENIGMA-DTI_males_w${metric}_${pheno}_chr${chr}_${date}_${analyst}"
                echo "$(gen_cmd "${pheno}" "$(echo ${covars_m}|xargs)" "${chr}" "${out_m}" "${ids_dir}/${COHORT}_males.ids")" >> "${gwas_outdir}/GWAS_commands_males.txt"
                echo "gzip -f ${gwas_outdir}/sumstats/${out_m}.${pheno}.glm.linear 2>/dev/null || true" >> "${gwas_outdir}/GWAS_commands_males.txt"
                n_gwas=$((n_gwas+1))
            fi

            # Model 5: Females-only (includes AffectionStatus, no sex terms)
            if [[ ${do_females} -eq 1 ]]; then
                covars_f="${global} ${covar_nosex} ${pc_covars} ${scanpc_covars} ${dummy_covars}"
                out_f="${COHORT}_${ancestry}_ENIGMA-DTI_females_w${metric}_${pheno}_chr${chr}_${date}_${analyst}"
                echo "$(gen_cmd "${pheno}" "$(echo ${covars_f}|xargs)" "${chr}" "${out_f}" "${ids_dir}/${COHORT}_females.ids")" >> "${gwas_outdir}/GWAS_commands_females.txt"
                echo "gzip -f ${gwas_outdir}/sumstats/${out_f}.${pheno}.glm.linear 2>/dev/null || true" >> "${gwas_outdir}/GWAS_commands_females.txt"
                n_gwas=$((n_gwas+1))
            fi
        done
    done
done

echo "  Total GWAS commands: ${n_gwas}"

#=============================================================================
# Execute or report
#=============================================================================
if [[ "${mode}" == "run" ]]; then
    echo ""
    echo "--- Executing GWAS ---"
    for cmd_file in "${gwas_outdir}"/GWAS_commands_*.txt; do
        [[ ! -s "${cmd_file}" ]] && continue
        n_cmds=$(grep -c "^plink2" "${cmd_file}" 2>/dev/null || echo 0)
        echo "  Running $(basename ${cmd_file}): ${n_cmds} GWAS"
        if [[ ${Nnodes} -le 1 ]]; then
            while IFS= read -r line; do
                eval "${line}" 2>&1 || true
            done < "${cmd_file}"
        else
            if command -v parallel &>/dev/null; then
                cat "${cmd_file}" | parallel -j ${Nnodes}
            else
                cat "${cmd_file}" | xargs -P ${Nnodes} -I {} bash -c '{}'
            fi
        fi
    done
    echo "  All GWAS complete."
else
    echo ""
    echo "--- Manual mode ---"
    for cmd_file in "${gwas_outdir}"/GWAS_commands_*.txt; do
        [[ -s "${cmd_file}" ]] && echo "  $(basename ${cmd_file}): $(grep -c '^plink2' "${cmd_file}" 2>/dev/null || echo 0) commands"
    done
fi

echo ""
echo "=== GWAS Backend Complete ==="
