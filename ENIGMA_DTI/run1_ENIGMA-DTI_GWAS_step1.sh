#!/bin/bash
###############################################################################
# run1_ENIGMA-DTI_GWAS_step1.sh
#
# ENIGMA-DTI GWAS Protocol — Step 1: Covariate Formatting & QC
#
# ┌──────────────────────────────────────────────────────────────────────────┐
# │  IMPORTANT — After this step, send ALL files from:                      │
# │    formatted/logs/                                                       │
# │    formatted/QC_statistical/                                             │
# │    formatted/QC_ENIGMA/                                                  │
# │  to enigma.DTIgenetics@gmail.com and WAIT for OK before run2.           │
# └──────────────────────────────────────────────────────────────────────────┘
#
# This script:
#   1. Reformats covariates for PLINK2 (recode Sex, AffectionStatus)
#   2. Generates derived age covariates
#   3. Merges ancestry PCs and global DTI metrics
#   4. Validates that site dummy variables are included
#   5. Produces QC histograms and summary statistics
#   6. Logs min/max age and all subject counts
#
# ENIGMA-DTI Genetics Working Group
# Date: February 2026
###############################################################################

set -euo pipefail

###############################################################################
# USER SETTINGS — Edit these to match your system and data
#
#   *** EDIT EVERYTHING IN THIS SECTION ***
#   *** DO NOT EDIT ANYTHING BELOW "END OF USER SETTINGS" ***
#
# This is the SECOND script you run (Step 1).
# It requires:
#   - Your Covariates.csv              (YOU provide — same file used in run0)
#   - formatted/Phenotypes.csv         (CREATED BY run0)
#   - PCA eigenvec file                (from ENIGMA pre-imputation QC)
#   - QC4 .fam file                    (from ENIGMA pre-imputation QC)
#
# It produces (used by run2):
#   - formatted/${COHORT}_enigma_dti_gwas.pheno  → phenotype file for PLINK2
#   - formatted/${COHORT}_enigma_dti_gwas.covar  → covariate file for PLINK2
#   - formatted/QC_ENIGMA/                       → QC plots (send to ENIGMA team)
#   - formatted/logs/${COHORT}_run1.log          → log (send to ENIGMA team)
#
# After this step, send ALL files from formatted/logs/, QC_statistical/,
# and QC_ENIGMA/ to enigma.DTIgenetics@gmail.com. WAIT for OK before run2.
###############################################################################

#==========================================================================
# SCRIPT DIRECTORY & R BINARY
#   Same as in run0 — path to the enigma_backend/ folder.
#==========================================================================
run_directory=$HOME/enigma/DTI_GWAS/SCRIPTS/enigma_backend
Rbin=$(which R 2>/dev/null || echo "/usr/local/bin/R")

#==========================================================================
# INPUT FILES
#==========================================================================
#-- Covariates CSV file
#   >>> YOU provide this: the SAME Covariates.csv you used in run0.
#   Must contain at least: SubjID, Age, Sex
#   Optionally: AffectionStatus, site dummy columns
csvFILE=$HOME/enigma/DTI_GWAS/Covariates.csv

#-- Combined DTI phenotype file
#   >>> CREATED BY run0 (Step 0). Do not create this yourself.
#   This is the merged phenotype table produced by run0.
#   Location: ${outDir}/Phenotypes.csv  (where outDir = outFolder from run0)
combinedROItableFILE=$HOME/enigma/DTI_GWAS/formatted/Phenotypes.csv

#-- PCA eigenvec file (ancestry principal components)
#   >>> YOU provide this: from ENIGMA pre-imputation genotype QC.
#   This is the .eigenvec file produced by PLINK2 --pca during QC Step 4.
#   Contains FID, IID, PC1, PC2, ..., PC20 (or however many PCs you computed).
#   Naming convention: ${COHORT}_${ancestry}_${analyst}_${date}_PCACovariates.eigenvec
pcaFILE=$HOME/enigma/DTI_GWAS/${COHORT}_${ancestry}_${analyst}_${date}_PCACovariates.eigenvec

#-- QC4 .fam file (final subject list after genotype QC)
#   >>> YOU provide this: from ENIGMA pre-imputation genotype QC.
#   This is the .fam file from the QC4 step (after relatedness filtering).
#   Used to identify which genotyped subjects should be included.
#   Naming convention: ${COHORT}_${ancestry}_${analyst}_${date}_QC4.fam
localfamFILE=$HOME/enigma/DTI_GWAS/${COHORT}_${ancestry}_${analyst}_${date}_QC4.fam

#==========================================================================
# COLUMN NAME MAPPINGS
#   Tell the script how your Covariates.csv encodes age, sex, and case/control.
#   These must match the EXACT column names and values in your CSV.
#==========================================================================
#-- Age column name in Covariates.csv (numeric, in years)
ageColumnHeader="Age"

#-- Sex column name and male indicator in Covariates.csv
#   sexColumnHeader: the column name (e.g., "Sex", "sex", "Gender")
#   maleIndicator:   the value representing males (e.g., "M", "Male", "1")
#   Females are assumed to be any value that is NOT the maleIndicator.
sexColumnHeader="Sex"
maleIndicator="M"

#-- Case/control cohort settings
#   CaseControlCohort: 1 = your cohort has cases and controls, 0 = healthy only
#   If 1, you must set the column name and affected indicator below.
CaseControlCohort=0

#   affectionStatusColumnHeader: column name in Covariates.csv (e.g., "AffectionStatus")
#   affectedIndicator: value representing AFFECTED/CASE subjects (e.g., "A", "Affected", "2")
#   Controls are assumed to be any value that is NOT the affectedIndicator.
affectionStatusColumnHeader="AffectionStatus"
affectedIndicator="A"

#-- If CaseControlCohort=1, describe what the affected status represents:
#   e.g., "Major Depressive Disorder", "Schizophrenia", "Chronic Pain"
#   This will be included in QC reports and logs sent to the ENIGMA team.
affectionDescription=""

#==========================================================================
# COHORT IDENTIFIERS
#   Must match EXACTLY what you used in run0.
#==========================================================================
COHORT="MYCOHORT"
date=$(date +%Y%m%d)

#==========================================================================
# OUTPUT DIRECTORY
#   >>> SAME as outFolder in run0. This is where run0 put Phenotypes.csv
#   and where this script will write the formatted .pheno and .covar files.
#==========================================================================
outDir=$HOME/enigma/DTI_GWAS/formatted

#==========================================================================
# SITE/SCANNER CORRECTION (choose ONE method or neither)
#==========================================================================
#-- Option A: Site dummy variables in Covariates.csv (space-separated column names)
#   If your cohort has multiple scanners/sites and you added dummy-coded columns
#   to Covariates.csv, list their column names here.
#   Example: site_dummies="Dummy1_MR Dummy2_MR Dummy3_MR"
#   These columns must exist in your Covariates.csv.
site_dummies=""

#-- Option B: Scanner PCA file (alternative to site dummies)
#   If your cohort's scanner/site variability is captured by PCA, provide
#   a file with PCA scores. The first scannerPCA_n components will be added
#   as covariates (scanPC1, scanPC2, ...).
#   >>> YOU provide this: a tab-separated file with subject IDs and PCA scores.
#   Set scannerPCA_file="" if not using scanner PCA.
scannerPCA_file=""
scannerPCA_n=6           # Number of scanner PCA components to include (default: 6)
scannerPCA_idcol="eid"   # Subject ID column name in the scanner PCA file

#==========================================================================
# OPTIONAL: INTRACRANIAL VOLUME (ICV)
#   If you have ICV data (e.g., from FreeSurfer), provide the file path.
#   This is used in QC plots only — it is NOT a GWAS covariate.
#   Set icvFILE="" if not available.
#==========================================================================
icvFILE=""
icvColumnHeader="segVOL_EstimatedTotalIntraCranialVol"

###############################################################################
# END OF USER SETTINGS
###############################################################################

#=============================================================================
# Activate conda environment (ensures we use conda R/Python, not system ones)
# Without this, the system R may hit GLIBC version mismatches on HPC nodes.
#=============================================================================
if command -v conda &>/dev/null; then
    # Temporarily disable set -e AND set -u: conda's shell hook references
    # variables (CONDA_SHLVL, _CE_CONDA, etc.) that may not exist yet, and
    # set -u kills the script on any unset variable even with set +e.
    set +eu
    eval "$(conda shell.bash hook)" 2>/dev/null || \
        source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate enigma_dti_gwas 2>/dev/null || true
    set -eu
    # Update Rbin to use conda R if available
    conda_R=$(which R 2>/dev/null || true)
    if [[ -n "${conda_R}" ]] && [[ "${conda_R}" == *"enigma_dti_gwas"* ]]; then
        Rbin="${conda_R}"
    fi
fi

LOGDIR="${outDir}/logs"
mkdir -p "${LOGDIR}"
LOGFILE="${LOGDIR}/${COHORT}_run1.log"
> "${LOGFILE}"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "=================================================================="
echo " ENIGMA-DTI GWAS — Step 1: Covariate Formatting & QC"
echo " Cohort: ${COHORT} | Date: ${date}"
echo " Started: $(date)"
echo " Log: ${LOGFILE}"
echo "=================================================================="
echo ""

eName="${COHORT}"

# Validate inputs
for f in "${csvFILE}" "${combinedROItableFILE}" "${pcaFILE}" "${localfamFILE}"; do
    [[ -f "${f}" ]] && echo "  [OK] $(basename ${f})" || { echo "  [ERROR] Not found: ${f}"; exit 1; }
done

# Check AffectionStatus description
if [[ ${CaseControlCohort} -eq 1 ]]; then
    if [[ -z "${affectionDescription}" ]]; then
        echo ""
        echo "  ╔═══════════════════════════════════════════════════════════════╗"
        echo "  ║  WARNING: CaseControlCohort=1 but affectionDescription is    ║"
        echo "  ║  empty. Please set affectionDescription to describe what     ║"
        echo "  ║  the affected/case status represents.                        ║"
        echo "  ║  e.g., affectionDescription=\"Major Depressive Disorder\"      ║"
        echo "  ╚═══════════════════════════════════════════════════════════════╝"
        echo ""
    else
        echo "  Affection status: ${affectionDescription}"
        echo "    Column: ${affectionStatusColumnHeader}, Affected indicator: ${affectedIndicator}"
    fi
fi
echo ""

mkdir -p "${outDir}"

#=============================================================================
# Check site dummy variables
#=============================================================================
echo "--- Checking site dummy variables ---"
if [[ -n "${site_dummies}" ]]; then
    cov_header=$(head -1 "${csvFILE}")
    for dummy in ${site_dummies}; do
        if echo "${cov_header}" | tr ',' '\n' | grep -qx "${dummy}"; then
            echo "  [OK] Site dummy found: ${dummy}"
        else
            echo "  [WARNING] Site dummy NOT found in Covariates.csv: ${dummy}"
            echo "            Please add this column or correct the name."
        fi
    done
else
    echo "  [NOTE] No site dummy variables specified."
    echo "         If your cohort has multiple sites/scanners, you should add"
    echo "         dummy variables to Covariates.csv and list them in site_dummies."
fi
echo ""

#=============================================================================
# Reformat covariates
#=============================================================================
echo "--- Reformatting covariates ---"

(  "${Rbin}" --no-save --slave --args \
    "${csvFILE}" "${combinedROItableFILE}" "${pcaFILE}" \
    "${ageColumnHeader}" "${sexColumnHeader}" "${maleIndicator}" \
    "${CaseControlCohort}" "${affectionStatusColumnHeader}" \
    "${affectedIndicator}" "${outDir}" "${eName}" \
    "${scannerPCA_file:-}" "${scannerPCA_n:-0}" "${scannerPCA_idcol:-eid}" \
    < "${run_directory}/ENIGMA_DTI_create_formatted_files.R"  ) 2>&1 | tee "${LOGDIR}/${COHORT}_covar_format.log"

echo ""

#=============================================================================
# Log age range from formatted covariate file
#=============================================================================
covar_formatted="${outDir}/${COHORT}_enigma_dti_gwas.covar"
if [[ -f "${covar_formatted}" ]]; then
    echo "--- Age summary ---"
    # Extract Age column (4th column in standard format)
    age_col=$(head -1 "${covar_formatted}" | tr '\t' '\n' | grep -n "^Age$" | cut -d: -f1)
    if [[ -n "${age_col}" ]]; then
        ages=$(tail -n +2 "${covar_formatted}" | cut -f"${age_col}")
        min_age=$(echo "${ages}" | sort -n | head -1)
        max_age=$(echo "${ages}" | sort -n | tail -1)
        n_subj=$(echo "${ages}" | wc -l)
        echo "  N subjects: ${n_subj}"
        echo "  Min age: ${min_age}"
        echo "  Max age: ${max_age}"
    fi
    echo ""
fi

#=============================================================================
# QC plots
#=============================================================================
echo "--- Running QC plots ---"
qcOutDir="${outDir}/QC_ENIGMA"
mkdir -p "${qcOutDir}"

qc_args="--cohort ${COHORT}"
qc_args="${qc_args} --covarFILE ${covar_formatted}"
qc_args="${qc_args} --phenoFILE ${outDir}/${COHORT}_enigma_dti_gwas.pheno"
qc_args="${qc_args} --outDir ${qcOutDir}"

if [[ -n "${icvFILE}" ]] && [[ -f "${icvFILE}" ]]; then
    qc_args="${qc_args} --icvFILE ${icvFILE} --icvColumnHeader ${icvColumnHeader}"
fi

RSCRIPT_CMD="$(dirname "${Rbin}")/Rscript"
if [[ ! -f "${RSCRIPT_CMD}" ]]; then
    RSCRIPT_CMD="Rscript"  # fallback to PATH
fi

${RSCRIPT_CMD} "${run_directory}/ENIGMA_DTI_GWAS_QC_plots_stats.R" ${qc_args} \
    2>&1 | tee "${LOGDIR}/${COHORT}_qc_plots.log"

echo ""

#=============================================================================
# N < 75 warning
#=============================================================================
if [[ -f "${covar_formatted}" ]]; then
    n_final=$(( $(wc -l < "${covar_formatted}") - 1 ))
    if [[ ${n_final} -lt 75 ]]; then
        echo " ╔═══════════════════════════════════════════════════════════════╗"
        echo " ║  WARNING: Only ${n_final} subjects. Below recommended minimum 75. ║"
        echo " ║  Contact enigma.DTIgenetics@gmail.com before proceeding.     ║"
        echo " ╚═══════════════════════════════════════════════════════════════╝"
    fi
fi

echo "=================================================================="
echo " Step 1 complete. $(date)"
echo ""
echo " Output: ${outDir}/${COHORT}_enigma_dti_gwas.{pheno,covar}"
echo " QC:     ${qcOutDir}/"
echo ""
echo " ┌─────────────────────────────────────────────────────────────┐"
echo " │ Send ALL files from formatted/logs/, QC_statistical/,      │"
echo " │ and QC_ENIGMA/ to enigma.DTIgenetics@gmail.com             │"
echo " │ WAIT for OK before running run2.                           │"
echo " └─────────────────────────────────────────────────────────────┘"
echo ""
echo " Next: After receiving OK, edit and run run2_ENIGMA-DTI_GWAS_step2.sh"
echo "=================================================================="
