#!/bin/bash
###############################################################################
# run2_ENIGMA-DTI_GWAS_step2.sh
#
# ENIGMA-DTI GWAS Protocol — Step 2: Run GWAS
#
# ┌──────────────────────────────────────────────────────────────────────────┐
# │  After ALL GWAS complete, package and upload results.                   │
# │  Email enigma.DTIgenetics@gmail.com for upload credentials.             │
# └──────────────────────────────────────────────────────────────────────────┘
###############################################################################

set -uo pipefail

###############################################################################
# USER SETTINGS — Edit these to match your system and data
#
#   *** EDIT EVERYTHING IN THIS SECTION ***
#   *** DO NOT EDIT ANYTHING BELOW "END OF USER SETTINGS" ***
#
# This is the THIRD and final script (Step 2).
# Only run this AFTER the ENIGMA team has approved your run0+run1 QC results.
#
# It requires:
#   - formatted/${COHORT}_enigma_dti_gwas.pheno  (CREATED BY run1)
#   - formatted/${COHORT}_enigma_dti_gwas.covar  (CREATED BY run1)
#   - IMPUTED/pgen/chr*.{pgen,pvar,psam}          (CREATED BY run0)
#
# It produces:
#   - GWAS_output/sumstats/*.glm.linear.gz  → per-chr GWAS summary statistics
#   - GWAS_output/combined/*_allchr.txt.gz  → combined across chromosomes
#   - GWAS_output/plots/*_manhattan.png     → Manhattan/QQ plots
#   - GWAS_output/meta/*_formatted.txt.gz   → meta-analysis-ready files
###############################################################################

#==========================================================================
# SCRIPT & BASE DIRECTORIES
#==========================================================================
#-- Path to the enigma_backend/ folder (same as in run0 and run1)
run_directory=$HOME/enigma/DTI_GWAS/SCRIPTS/enigma_backend

#-- Base directory for the ENIGMA-DTI GWAS project
#   This is the parent directory containing formatted/, IMPUTED/, SCRIPTS/, etc.
#   All output goes under ${ENIGMA_DTI_GWAS_dir}/GWAS_output/
ENIGMA_DTI_GWAS_dir=$HOME/enigma/DTI_GWAS

#==========================================================================
# COHORT IDENTIFIERS
#   Must match EXACTLY what you used in run0 and run1.
#==========================================================================
COHORT="MYCOHORT"
ancestry="EUR"
analyst="AB"
date=$(date +%Y%m%d)

#==========================================================================
# INPUT FILES — All created by previous steps
#==========================================================================
#-- PGEN directory (converted genotype files)
#   >>> CREATED BY run0 (Step 0). Contains chr1.pgen, chr1.pvar, chr1.psam, etc.
#   If you set pgen_outdir="" in run0, files are in ${imputed_dir}/pgen
pgen_dir=$HOME/enigma/DTI_GWAS/IMPUTED/pgen

#-- Phenotype file for PLINK2
#   >>> CREATED BY run1 (Step 1). Tab-separated, FID IID <metric>_<ROI> columns.
#   Do NOT create or modify this file yourself.
phenoFILE=$HOME/enigma/DTI_GWAS/formatted/${COHORT}_enigma_dti_gwas.pheno

#-- Covariate file for PLINK2
#   >>> CREATED BY run1 (Step 1). Tab-separated, FID IID Sex Age AgeCsq ...
#   Contains recoded Sex (1/2), derived age terms (AgeCsq, AgexSex, AgeCsqxSex),
#   ancestry PCs, and any site dummy or scanner PCA covariates.
#   Do NOT create or modify this file yourself.
covarFILE=$HOME/enigma/DTI_GWAS/formatted/${COHORT}_enigma_dti_gwas.covar

#-- Directory containing the .pheno and .covar files above
#   >>> Same as outDir in run1 / outFolder in run0.
pheno_covar_dir=$HOME/enigma/DTI_GWAS/formatted

#==========================================================================
# COHORT TYPE
#   0 = unrelated (most common), 1 = related/family-based
#   Must match cohort_unrelated from run0 (but inverted: unrelated=1 in run0
#   corresponds to related=0 here).
#==========================================================================
related=0

#==========================================================================
# EXECUTION MODE
#   "qsub"   — submit GWAS jobs to SGE/PBS cluster (recommended for HPC)
#   "run"    — run all GWAS locally in this shell session (for small cohorts)
#   "manual" — only generate command files, do not run anything
#              (commands saved to GWAS_output/GWAS_commands_*.txt)
#==========================================================================
mode="qsub"

#==========================================================================
# QSUB / CLUSTER SETTINGS (only used when mode="qsub")
#   Adjust these for your HPC cluster.
#==========================================================================
qsub_mem="8G"              # Memory per job (e.g., "8G", "16G")
#qsub_walltime="02:00:00"  # Uncomment and set if your cluster requires walltime
qsub_queue=""               # Queue/partition name (e.g., "normal", "short")
qsub_extra=""               # Any extra qsub flags (e.g., "-l h_rt=02:00:00")
qsub_jobname="eDTI"        # Job name prefix (shows up in qstat)

#==========================================================================
# CHROMOSOMES
#   Which chromosomes to include in the GWAS.
#   chrStart/chrEnd: range of autosomes (normally 1–22)
#   doChrX: 1 = include chrX, 0 = skip chrX
#==========================================================================
chrStart=1; chrEnd=22; doChrX=1

###############################################################################
# END OF USER SETTINGS
###############################################################################

user_mode="${mode}"

LOGDIR="${ENIGMA_DTI_GWAS_dir}/formatted/logs"
mkdir -p "${LOGDIR}"
LOGFILE="${LOGDIR}/${COHORT}_run2.log"
> "${LOGFILE}"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "=================================================================="
echo " ENIGMA-DTI GWAS — Step 2"
echo " Cohort: ${COHORT} | Date: ${date} | Mode: ${user_mode}"
echo " Started: $(date)"
echo "=================================================================="
echo ""

for f in "${phenoFILE}" "${covarFILE}"; do
    [[ -f "${f}" ]] && echo "  [OK] $(basename ${f})" || { echo "  [ERROR] Not found: ${f}"; exit 1; }
done
[[ -d "${pgen_dir}" ]] && echo "  [OK] PGEN dir" || { echo "  [ERROR] ${pgen_dir}"; exit 1; }
echo ""

gwas_outdir="${ENIGMA_DTI_GWAS_dir}/GWAS_output"
mkdir -p "${gwas_outdir}"/{logs,sumstats,meta,plots,job_scripts}

export COHORT ancestry analyst date
export pgen_dir phenoFILE covarFILE pheno_covar_dir
export gwas_outdir run_directory related
export chrStart chrEnd doChrX ENIGMA_DTI_GWAS_dir
export mode="manual"
export Nnodes=1

echo "--- Generating GWAS commands ---"
chmod +x "${run_directory}/ENIGMA_unrelatedGWAS_plink_DTI.sh"
"${run_directory}/ENIGMA_unrelatedGWAS_plink_DTI.sh" 2>&1 | tee "${LOGDIR}/${COHORT}_gwas_backend.log"
echo ""

# Count commands
total_cmds=0
for cmd_file in "${gwas_outdir}"/GWAS_commands_*.txt; do
    [[ ! -s "${cmd_file}" ]] && continue
    n=$(grep -c "^plink2" "${cmd_file}" 2>/dev/null || echo 0)
    total_cmds=$((total_cmds + n))
    echo "  $(basename ${cmd_file}): ${n} GWAS"
done
echo "  Total: ${total_cmds}"
echo ""

if [[ ${total_cmds} -eq 0 ]]; then
    echo "  [ERROR] No GWAS commands generated."
    exit 1
fi

# Build chromosome list
echo "--- Building chromosome list ---"
echo "  pgen_dir: ${pgen_dir}"
echo "  chrStart: ${chrStart}, chrEnd: ${chrEnd}, doChrX: ${doChrX}"

chr_list=""
n_chrs=0
for chr in $(seq ${chrStart} ${chrEnd}); do
    if [[ -f "${pgen_dir}/chr${chr}.pgen" ]]; then
        chr_list="${chr_list} ${chr}"
        n_chrs=$((n_chrs + 1))
    fi
done
if [[ ${doChrX} -eq 1 ]] && [[ -f "${pgen_dir}/chrX.pgen" ]]; then
    chr_list="${chr_list} X"
    n_chrs=$((n_chrs + 1))
fi

echo "  Found ${n_chrs} chromosomes with PGEN files"

if [[ ${n_chrs} -eq 0 ]]; then
    echo "  [ERROR] No chr*.pgen files found in ${pgen_dir}"
    echo "  Files in pgen_dir:"
    ls "${pgen_dir}"/*.pgen 2>/dev/null | head -5 || echo "    (none)"
    ls "${pgen_dir}"/ 2>/dev/null | head -10
    exit 1
fi

chr_array_file="${gwas_outdir}/job_scripts/_chr_list.txt"
mkdir -p "$(dirname ${chr_array_file})"
echo ${chr_list} | tr ' ' '\n' > "${chr_array_file}"
echo "  Chromosome list file: ${chr_array_file}"
echo "  Chromosomes: $(cat ${chr_array_file} | tr '\n' ',' | sed 's/,$//')"
echo ""

echo "--- Execution mode: ${user_mode} ---"

#=============================================================================
# QSUB MODE
#=============================================================================
if [[ "${user_mode}" == "qsub" ]]; then

    echo "=================================================================="
    echo " Submitting array jobs via qsub"
    echo " ${n_chrs} tasks per array (one per chromosome)"
    echo "=================================================================="
    echo ""

    job_dir="${gwas_outdir}/job_scripts"
    n_arrays=0
    n_errors=0
    all_job_ids=""

    for cmd_file in "${gwas_outdir}"/GWAS_commands_*.txt; do
        [[ ! -s "${cmd_file}" ]] && continue
        model_name=$(basename "${cmd_file}" .txt | sed 's/GWAS_commands_//')
        echo "  Model: ${model_name}"

        # Extract unique phenotype names (portable — no grep -P)
        pheno_list=$(grep "pheno-name" "${cmd_file}" | \
                     sed 's/.*--pheno-name //; s/ .*//' | sort -u)

        for pheno_name in ${pheno_list}; do

            # Get one example plink2 command for this phenotype
            plink_line=$(grep "^plink2" "${cmd_file}" | grep "pheno-name ${pheno_name} " | head -1)
            [[ -z "${plink_line}" ]] && plink_line=$(grep "^plink2" "${cmd_file}" | grep "pheno-name ${pheno_name}" | head -1)
            [[ -z "${plink_line}" ]] && continue

            # Replace chromosome-specific parts with placeholder
            generic_cmd=$(echo "${plink_line}" | \
                sed "s|${pgen_dir}/chr[0-9X]*|${pgen_dir}/chrCHRNUM|g" | \
                sed 's|_chr[0-9X]*_|_chrCHRNUM_|g')

            # Get output prefix
            out_prefix=$(echo "${generic_cmd}" | sed 's/.*--out //; s/ .*//')

            # Write job script
            script="${job_dir}/${qsub_jobname}_${model_name}_${pheno_name}.sh"

            {
                echo '#!/bin/bash'
                echo "#\$ -cwd"
                echo "#\$ -N ${qsub_jobname}_${model_name}_${pheno_name}"
                echo "#\$ -l h_vmem=${qsub_mem}"
                # No time limit
                echo '#$ -o '"${gwas_outdir}/logs/"'$JOB_NAME.$TASK_ID.stdout'
                echo '#$ -e '"${gwas_outdir}/logs/"'$JOB_NAME.$TASK_ID.stderr'
                echo "#\$ -t 1-${n_chrs}"
                [[ -n "${qsub_queue}" ]] && echo "#\$ -q ${qsub_queue}"
                [[ -n "${qsub_extra}" ]] && echo "#\$ ${qsub_extra}"
                echo ""
                echo "# Clear module conflicts and activate conda environment"
                echo 'module purge 2>/dev/null || true'
                echo 'eval "$(conda shell.bash hook)" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh"'
                echo 'conda activate enigma_dti_gwas 2>/dev/null || true'
                echo ""
                echo '# Get chromosome for this array task'
                echo 'CHR=$(sed -n "${SGE_TASK_ID}p" '"${chr_array_file}"')'
                echo 'if [[ -z "${CHR}" ]]; then echo "No chr for task ${SGE_TASK_ID}"; exit 1; fi'
                echo 'echo "'"${pheno_name}"' chr${CHR} task ${SGE_TASK_ID}"'
                echo ""
                # Build output filename to check if it already exists
                out_with_chr=$(echo "${out_prefix}" | sed 's/CHRNUM/${CHR}/g')
                echo "# Skip if output already exists"
                echo 'OUTFILE='"$(echo "${out_prefix}" | sed 's/CHRNUM/${CHR}/g')"'.'"${pheno_name}"'.glm.linear.gz'
                echo 'if [[ -f "${OUTFILE}" ]]; then echo "Already done: ${OUTFILE}"; exit 0; fi'
                echo ""
                # Write plink2 command with CHRNUM replaced by ${CHR}
                echo "${generic_cmd}" | sed 's/CHRNUM/${CHR}/g'
                echo ""
                echo "# Compress output"
                echo "${out_prefix}" | sed 's/CHRNUM/${CHR}/g' | \
                    xargs -I{} echo 'gzip -f {}.'"${pheno_name}"'.glm.linear 2>/dev/null || true'
            } > "${script}"

            chmod +x "${script}"

            # Submit
            if qsub_out=$(qsub "${script}" 2>&1); then
                job_id=$(echo "${qsub_out}" | grep -o '[0-9]*' | head -1)
                echo "    ${pheno_name}: job ${job_id:-?}"
                n_arrays=$((n_arrays + 1))
                [[ -n "${job_id}" ]] && all_job_ids="${all_job_ids:-}${all_job_ids:+,}${job_id}"
            else
                echo "    [ERROR] ${pheno_name}: ${qsub_out}"
                n_errors=$((n_errors + 1))
            fi

        done
    done

    echo ""
    echo "  ┌─────────────────────────────────────────────────────────────┐"
    echo "  │ SUBMITTED: ${n_arrays} array jobs × ${n_chrs} chr = $((n_arrays * n_chrs)) tasks"
    [[ ${n_errors} -gt 0 ]] && \
    echo "  │ ERRORS: ${n_errors}"
    echo "  │ Monitor: qstat -u \$(whoami)"
    echo "  │ Scripts: ${job_dir}/"
    echo "  └─────────────────────────────────────────────────────────────┘"

    #=========================================================================
    # Postprocessing: (A) array job to combine chr files in parallel
    #                 (B) dependent plotting job
    #=========================================================================
    echo ""
    echo "--- Submitting postprocessing ---"

    group_list_file="${gwas_outdir}/job_scripts/_combine_groups.txt"
    # Build group list from GWAS command files (not sumstats — those don't exist yet)
    for cmd_file in "${gwas_outdir}"/GWAS_commands_*.txt; do
        [[ ! -s "${cmd_file}" ]] && continue
        grep "^plink2" "${cmd_file}" | sed 's/.*--out [^ ]*sumstats\///; s/ .*//' | \
            sed 's/_chr[0-9X]*_.*//g'
    done | sort -u > "${group_list_file}"
    n_groups=$(wc -l < "${group_list_file}")
    echo "  Groups to combine: ${n_groups}"

    if [[ ${n_groups} -gt 0 ]]; then
        # (A) Array job: combine chromosomes (one task per group, all parallel)
        combine_script="${job_dir}/${qsub_jobname}_combine.sh"
        {
            echo '#!/bin/bash'
            echo "#\$ -cwd"
            echo "#\$ -N ${qsub_jobname}_combine"
            echo "#\$ -l h_vmem=4G"
            echo '#$ -o '"${gwas_outdir}/logs/"'$JOB_NAME.$TASK_ID.stdout'
            echo '#$ -e '"${gwas_outdir}/logs/"'$JOB_NAME.$TASK_ID.stderr'
            echo "#\$ -t 1-${n_groups}"
            [[ -n "${qsub_queue}" ]] && echo "#\$ -q ${qsub_queue}"
            echo ""
            echo 'module purge 2>/dev/null || true'
            echo 'eval "$(conda shell.bash hook)" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh"'
            echo 'conda activate enigma_dti_gwas 2>/dev/null || true'
            echo ""
            echo "# Rebuild group list at runtime from actual output files"
            echo 'ls '"${gwas_outdir}/sumstats/"'*.glm.linear.gz 2>/dev/null | \'
            echo '    xargs -I{} basename {} | \'
            echo "    sed 's/_chr[0-9X]*_.*//g' | \\"
            echo '    sort -u > '"${group_list_file}"
            echo ""
            echo "bash ${run_directory}/ENIGMA_DTI_combine_one.sh \\"
            echo "    ${group_list_file} \\"
            echo "    ${gwas_outdir}/sumstats \\"
            echo "    ${gwas_outdir}/combined"
        } > "${combine_script}"
        chmod +x "${combine_script}"

        combine_jid=""
        if comb_out=$(qsub "${combine_script}" 2>&1); then
            combine_jid=$(echo "${comb_out}" | grep -o '[0-9]*' | head -1)
            echo "  Combine array: job ${combine_jid:-?} (${n_groups} tasks in parallel)"
        else
            echo "  [ERROR] Combine submission failed: ${comb_out}"
        fi

        # (B) Array job: generate plots (one task per combined file, all parallel)
        plots_list_file="${gwas_outdir}/job_scripts/_plot_files.txt"
        # This file gets populated AFTER combining finishes, so the plots job
        # needs to build it at runtime
        plots_script="${job_dir}/${qsub_jobname}_plots.sh"
        {
            echo '#!/bin/bash'
            echo "#\$ -cwd"
            echo "#\$ -N ${qsub_jobname}_plots"
            echo "#\$ -l h_vmem=8G"
            echo '#$ -o '"${gwas_outdir}/logs/"'$JOB_NAME.$TASK_ID.stdout'
            echo '#$ -e '"${gwas_outdir}/logs/"'$JOB_NAME.$TASK_ID.stderr'
            echo "#\$ -t 1-${n_groups}"
            [[ -n "${qsub_queue}" ]] && echo "#\$ -q ${qsub_queue}"
            echo ""
            echo 'module purge 2>/dev/null || true'
            echo 'eval "$(conda shell.bash hook)" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh"'
            echo 'conda activate enigma_dti_gwas 2>/dev/null || true'
            echo ""
            echo "# Build file list from combined directory"
            echo "PLOT_LIST=${plots_list_file}"
            echo 'ls '"${gwas_outdir}/combined/"'*_allchr.txt.gz 2>/dev/null | sort > "${PLOT_LIST}"'
            echo ""
            echo "bash ${run_directory}/ENIGMA_DTI_plot_one.sh \\"
            echo "    \${PLOT_LIST} \\"
            echo "    ${gwas_outdir}/plots \\"
            echo "    ${COHORT}"
        } > "${plots_script}"
        chmod +x "${plots_script}"

        if [[ -n "${combine_jid}" ]]; then
            if plt_out=$(qsub -hold_jid "${combine_jid}" "${plots_script}" 2>&1); then
                plt_id=$(echo "${plt_out}" | grep -o '[0-9]*' | head -1)
                echo "  Plots array: job ${plt_id:-?} (${n_groups} tasks, runs after combining)"
            else
                echo "  [WARN] Plots failed. Run manually: bash ${plots_script}"
            fi
        else
            if plt_out=$(qsub "${plots_script}" 2>&1); then
                plt_id=$(echo "${plt_out}" | grep -o '[0-9]*' | head -1)
                echo "  Plots array: job ${plt_id:-?} (ensure combining is done first)"
            else
                echo "  Run manually: bash ${plots_script}"
            fi
        fi

        # (C) Summary job: collects all per-file summaries into one CSV
        # (runs after plots finish — or user can run manually)
        echo ""
        echo "  After all jobs complete, generate summary with:"
        echo "    (echo -e "Group\tPhenotype\tN_variants\tLambda\tMin_P"; cat ${gwas_outdir}/plots/*_summary.txt | sort) > ${gwas_outdir}/${COHORT}_GWAS_summary.tsv"

        echo ""
        echo "  After completion, send to ENIGMA-DTI team:"
        echo "    ${gwas_outdir}/plots/       (Manhattan + QQ)"
        echo "    ${gwas_outdir}/combined/    (combined sumstats)"
        echo "    ${gwas_outdir}/${COHORT}_GWAS_summary.csv"
    fi


#=============================================================================
# RUN MODE
#=============================================================================
elif [[ "${user_mode}" == "run" ]]; then
    echo "--- Running locally ---"
    for cmd_file in "${gwas_outdir}"/GWAS_commands_*.txt; do
        [[ ! -s "${cmd_file}" ]] && continue
        echo "  $(basename ${cmd_file})"
        while IFS= read -r line; do
            [[ -z "${line}" ]] && continue
            eval "${line}" 2>&1 || true
        done < "${cmd_file}"
    done
    echo "  All GWAS runs complete."
    echo ""
    echo "--- Running postprocessing (combine + plots) ---"

    export COMBINED_DIR="${gwas_outdir}/combined"
    export PLOTS_DIR="${gwas_outdir}/plots"
    export COHORT_NAME="${COHORT}"

    # Combine all groups
    group_list_file="${gwas_outdir}/job_scripts/_combine_groups.txt"
    ls "${gwas_outdir}/sumstats/"*.glm.linear.gz 2>/dev/null | \
        xargs -I{} basename {} | sed 's/_chr[0-9X]*_.*//g' | sort -u > "${group_list_file}"
    n_groups=$(wc -l < "${group_list_file}")
    for i in $(seq 1 ${n_groups}); do
        SGE_TASK_ID=${i} bash "${run_directory}/ENIGMA_DTI_combine_one.sh" \
            "${group_list_file}" "${gwas_outdir}/sumstats" "${gwas_outdir}/combined"
    done

    # Generate plots
    bash "${run_directory}/ENIGMA_DTI_plots.sh" \
        "${gwas_outdir}/combined" "${gwas_outdir}/plots" "${COHORT}" \
        2>&1 | tee "${LOGDIR}/${COHORT}_plots.log"

#=============================================================================
# MANUAL MODE
#=============================================================================
else
    echo "--- Manual mode ---"
    echo "  Command files in: ${gwas_outdir}/"
    for cmd_file in "${gwas_outdir}"/GWAS_commands_*.txt; do
        [[ -s "${cmd_file}" ]] && echo "    $(basename ${cmd_file})"
    done
    echo "  Run: cat ${gwas_outdir}/GWAS_commands_*.txt | bash"
fi

echo ""
echo "=================================================================="
echo " Step 2 complete. $(date)"
echo "=================================================================="
