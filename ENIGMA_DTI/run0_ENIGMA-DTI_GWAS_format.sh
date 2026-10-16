#!/bin/bash
###############################################################################
# run0_ENIGMA-DTI_GWAS_format.sh
#
# ENIGMA-DTI GWAS Protocol — Step 0: Phenotype QC, Formatting & Genotype Prep
#
# ┌──────────────────────────────────────────────────────────────────────────┐
# │  IMPORTANT — FILES TO SEND BACK TO THE ENIGMA-DTI GENETICS TEAM        │
# │                                                                         │
# │  After running run0 and run1, please send the following to              │
# │  enigma.DTIgenetics@gmail.com for review BEFORE running the GWAS:      │
# │                                                                         │
# │    • All log files from formatted/logs/                                 │
# │    • All QC files from formatted/QC_statistical/                        │
# │    • All PDF/TXT/CSV from formatted/QC_ENIGMA/                         │
# │    • The phenotype_merge_QC_log.txt                                    │
# │    • The site_summary_statistics.csv (if multi-site)                    │
# │    • The ${COHORT}_subject_matching_report.txt                         │
# │                                                                         │
# │  Wait for confirmation before proceeding to run2.                       │
# └──────────────────────────────────────────────────────────────────────────┘
#
# This script:
#   1. Checks required software
#   2. Validates input files, detects available DTI metrics, checks column names
#   3. Harmonises subjects across all phenotype, covariate, and genotype files
#   4. Runs eDTI_outliers statistical QC and removes flagged subjects
#   5. Filters related individuals (if cohort is unrelated)
#   6. Converts imputed VCF to PGEN format
#   7. Merges DTI phenotype files into Phenotypes.csv
#   8. Produces subject matching report and site-level statistics
#
# All output is logged to ${COHORT}_run0.log
#
# ENIGMA-DTI Genetics Working Group
# Date: February 2026
###############################################################################

set -euo pipefail

###############################################################################
# USER SETTINGS — Edit these paths to match your system and data
#
#   *** EDIT EVERYTHING IN THIS SECTION ***
#   *** DO NOT EDIT ANYTHING BELOW "END OF USER SETTINGS" ***
#
# This is the FIRST script you run (Step 0).
# It requires:
#   - Your DTI phenotype CSV files (from the ENIGMA-DTI TBSS pipeline)
#   - Your covariates CSV file (you create this)
#   - Your QCed genotype files (from ENIGMA pre-imputation QC)
#   - Your imputed VCF files (from Michigan/TOPMed imputation server)
#
# It produces (used by run1 and run2):
#   - formatted/Phenotypes.csv           → used by run1
#   - IMPUTED/pgen/chr*.{pgen,pvar,psam} → used by run2
#   - formatted/QC_statistical/          → send to ENIGMA team
#   - formatted/${COHORT}_run0.log       → send to ENIGMA team
###############################################################################

#==========================================================================
# COHORT IDENTIFIERS
#   These four values are used to name output files throughout the pipeline.
#   Use the SAME values in run0, run1, and run2.
#==========================================================================
#-- Short name for your cohort (no spaces or special characters)
#   Example: "HABS_HD", "TOP", "QTIM", "HCP"
COHORT="MYCOHORT"

#-- Ancestry group (must match your genotype QC)
#   Common values: "EUR", "AFR", "AMR", "EAS", "SAS"
ancestry="EUR"

#-- Your initials (analyst who ran the pipeline)
analyst="AB"

#-- Date stamp — auto-generated, do not change
date=$(date +%Y%m%d)

#==========================================================================
# SCRIPT DIRECTORY
#   Path to the enigma_backend/ folder that came with this pipeline.
#   This folder contains the R scripts and backend shell scripts.
#   >>> YOU provide this: it's wherever you extracted the pipeline zip.
#==========================================================================
run_directory=$HOME/enigma/DTI_GWAS/SCRIPTS/enigma_backend

#==========================================================================
# R BINARY
#   Path to R (version ≥ 4.0 required).
#   If you set up the enigma_dti_gwas conda environment (recommended),
#   the default line below will find it automatically.
#   Otherwise, specify the full path. Find it with:  which R
#==========================================================================
Rbin=$(which R 2>/dev/null || echo "/usr/local/bin/R")

#==========================================================================
# DTI PHENOTYPE FILES (from the ENIGMA-DTI TBSS pipeline)
#   >>> YOU provide these: output of the ENIGMA-DTI TBSS ROI extraction.
#   Each file is a CSV with one row per subject and columns:
#     SubjID, <metric>_<ROI1>, <metric>_<ROI2>, ...
#   Example columns: SubjID, FA_ACR, FA_ALIC, FA_BCC, ...
#
#   You MUST have at least the FA file. Comment out any metric you don't have
#   by putting a # at the start of the line.
#
#   These are the combinedROItable files produced by the ENIGMA-DTI
#   protocol's ROI extraction step (TBSS → extract ROIs → average L/R).
#==========================================================================
csvFILE_FA=$HOME/enigma/DTI_GWAS/combinedROItable_FA.csv
csvFILE_MD=$HOME/enigma/DTI_GWAS/combinedROItable_MD.csv
csvFILE_AD=$HOME/enigma/DTI_GWAS/combinedROItable_AD.csv
csvFILE_RD=$HOME/enigma/DTI_GWAS/combinedROItable_RD.csv

#==========================================================================
# GENOTYPE FILE (QCed, pre-imputation)
#   >>> YOU provide this: output of ENIGMA pre-imputation genotype QC.
#   This is the PLINK file PREFIX (no .bed/.bim/.fam extension).
#   The script expects files: <prefix>.bed, <prefix>.bim, <prefix>.fam
#
#   Naming convention: ${COHORT}_${ancestry}_${analyst}_${date}_QC3
#   The "_QC3" suffix indicates it passed steps QC1–QC3 of the ENIGMA
#   genotype QC protocol.
#
#   The .fam file is used to match genotyped subjects to phenotype data.
#==========================================================================
datafile=$HOME/enigma/DTI_GWAS/${COHORT}_${ancestry}_${analyst}_${date}_QC3

#==========================================================================
# IMPUTED VCF DIRECTORY
#   >>> YOU provide this: the directory containing your imputed VCF files.
#   These come from the imputation server (Michigan, TOPMed, Sanger).
#   Expected files: chr1.dose.vcf.gz, chr2.dose.vcf.gz, ..., chrX.dose.vcf.gz
#   (Also accepts: chr1.vcf.gz naming)
#
#   This script will convert these VCFs to PLINK2 PGEN format.
#==========================================================================
imputed_dir=$HOME/enigma/DTI_GWAS/IMPUTED

#==========================================================================
# GENOME BUILD
#   The coordinate system of your imputed VCF files.
#   Used for chrX --split-par (separating PAR1/PAR2 from non-PAR X).
#   Almost all imputation servers output hg19/GRCh37 coordinates.
#   Set to "hg38" ONLY if you know your VCFs use GRCh38 coordinates.
#   If unsure, leave as "hg19".
#==========================================================================
genome_build="hg19"

#==========================================================================
# PGEN OUTPUT DIRECTORY
#   Where to write the converted PGEN files (chr1.pgen, chr1.pvar, chr1.psam, ...).
#   >>> CREATED BY THIS SCRIPT (run0). Used by run2 for GWAS.
#   Can be several GB per chromosome — ensure enough disk space.
#   If left empty (""), defaults to ${imputed_dir}/pgen
#==========================================================================
pgen_outdir=""

#==========================================================================
# OUTPUT DIRECTORY (formatted files)
#   >>> CREATED BY THIS SCRIPT (run0).
#   All formatted phenotype/covariate files and QC results go here.
#   This directory is used by run1 (Step 1) as input.
#   Key outputs created here:
#     Phenotypes.csv                    → merged DTI phenotype file (used by run1)
#     QC_statistical/                   → outlier QC results (send to ENIGMA team)
#     logs/${COHORT}_run0.log          → full log (send to ENIGMA team)
#     ${COHORT}_subject_matching_report.txt → subject overlap report
#==========================================================================
outFolder=$HOME/enigma/DTI_GWAS/formatted

#==========================================================================
# SUBJECT ID COLUMN
#   The name of the subject ID column in your DTI CSV files AND Covariates.csv.
#   Must be IDENTICAL across ALL input files.
#   Common values: "SubjID", "IID", "subjectkey"
#==========================================================================
subjID_column="SubjID"

#==========================================================================
# COVARIATES CSV FILE
#   >>> YOU provide this: a CSV file you create with one row per subject.
#   Required columns:
#     SubjID   — must match the IDs in your DTI CSV files and genotype files
#     Age      — age in years (numeric)
#     Sex      — biological sex (e.g., "M"/"F" or "Male"/"Female")
#   Optional columns:
#     AffectionStatus — case/control status (e.g., "A"/"C" or "Affected"/"Control")
#     Site/scanner dummy variables — if your cohort has multiple scanners
#   The Sex/Age coding is configured in run1 (maleIndicator, ageColumnHeader, etc.)
#==========================================================================
covariatesCSV=$HOME/enigma/DTI_GWAS/Covariates.csv

#==========================================================================
# COHORT RELATEDNESS
#   1 = cohort is UNRELATED (most common — population-based samples)
#   0 = cohort has RELATED individuals (family-based samples)
#   If 1, the script will use PLINK2 --king-cutoff to remove one of each
#   related pair before GWAS. If 0, all subjects are kept.
#==========================================================================
cohort_unrelated=1

#===========================================================================
# STATISTICAL QC SETTINGS (eDTI_outliers)
#   These settings control the automated outlier detection for DTI phenotypes.
#   The defaults below follow the ENIGMA-DTI protocol. Change only if
#   instructed by the ENIGMA-DTI genetics team.
#===========================================================================
#-- Run statistical QC? (1 = yes [recommended], 0 = skip)
run_statistical_qc=1

#-- Demographics CSV for QC (optional, for site-stratified QC)
#   If your cohort has multiple sites, provide a CSV with SubjID, site, and
#   diagnosis columns. Set to "" if not available or single-site.
demogCSV=""
dxColumn=""                # Column name for diagnosis in demogCSV (e.g., "Dx")
siteColumn=""              # Column name for site/scanner in demogCSV (e.g., "Site")

#-- Site dummy variable column names in Covariates.csv (space-separated)
#   If your cohort has multiple scanners/sites, you MUST include dummy variables
#   in Covariates.csv and list their column names here.
#   Example: site_dummies="Dummy1_MR Dummy2_MR Dummy3_MR"
#   Leave empty if single-site or using scanner PCA instead (see run1).
site_dummies=""

#-- Outlier thresholds (ENIGMA-DTI protocol defaults — do not change
#   unless instructed by the ENIGMA team)
quant1_lower=0.001; quant1_upper=0.999; perc1=0.05
quant2_lower=0.005; quant2_upper=0.995; perc2=0.10

#===========================================================================
# VISUAL QC SETTINGS (eDTI_visualQC)
#===========================================================================
#-- Run visual QC on subjects flagged by statistical QC? (1=yes, 0=skip)
#   Generates PNG images of FA registrations, skeleton overlays, and
#   warped diffusion maps so you can visually inspect outliers.
#   Requires: FSL installed and on PATH, and your TBSS directory.
run_visual_qc=1

#-- Path to your TBSS processing directory
#   >>> YOU provide this: the directory from the ENIGMA-DTI TBSS pipeline
#   that contains FA_individ/, stats/, etc.
#   Set to "" if TBSS files are not available (visual QC will be skipped).
#   Example: TBSS_dir=$HOME/enigmaDTI/TBSS/run_tbss
TBSS_dir=""

###############################################################################
# END OF USER SETTINGS — Do not edit below this line
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

#=============================================================================
# 0. Set up logging
#=============================================================================
LOGDIR="${outFolder}/logs"
mkdir -p "${LOGDIR}"
LOGFILE="${LOGDIR}/${COHORT}_run0.log"
> "${LOGFILE}"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "=================================================================="
echo " ENIGMA-DTI GWAS — Step 0"
echo " Cohort: ${COHORT} | Ancestry: ${ancestry} | Analyst: ${analyst}"
echo " Date: ${date}"
echo " Started: $(date)"
echo " Log: ${LOGFILE}"
echo "=================================================================="
echo ""
echo "┌──────────────────────────────────────────────────────────────────┐"
echo "│ REMINDER: After run0 and run1, send all QC files to             │"
echo "│ enigma.DTIgenetics@gmail.com and wait for OK before run2.       │"
echo "└──────────────────────────────────────────────────────────────────┘"
echo ""

export ENIGMA_DATE="${date}"
datafile_prefix=$(basename "${datafile}")
mkdir -p "${outFolder}" "${outFolder}/QC_statistical"

#=============================================================================
# 1. Check required software
#=============================================================================
echo "=================================================================="
echo " [1/8] Checking required software"
echo "=================================================================="

check_tool() {
    if command -v "$1" &>/dev/null; then
        echo "  [OK]      $1 — $(command -v $1)"
        return 0
    else
        echo "  [MISSING] $1"
        return 1
    fi
}

FAIL=0
check_tool plink2    || FAIL=1
check_tool "${Rbin}" || FAIL=1
check_tool vcftools  || FAIL=1
check_tool bcftools  || FAIL=1
check_tool tabix     || FAIL=1
check_tool git       || FAIL=1
check_tool pip       || check_tool pip3 || echo "  [WARN] pip not found — needed if eDTI_outliers is not installed"

if [[ ${FAIL} -eq 1 ]]; then
    echo "  [ERROR] Missing required tools. See SOP for installation links."
    exit 1
fi
echo ""

#=============================================================================
# 2. Validate input files, detect metrics, check column names
#=============================================================================
echo "=================================================================="
echo " [2/8] Validating input files & checking column names"
echo "=================================================================="

declare -a DTI_FILES=()
declare -a DTI_METRICS=()

for metric in FA MD AD RD; do
    varname="csvFILE_${metric}"
    filepath="${!varname:-}"
    if [[ -n "${filepath}" ]] && [[ -f "${filepath}" ]]; then
        n_lines=$(($(wc -l < "${filepath}") - 1))
        echo "  [OK]    ${metric}: $(basename ${filepath}) — ${n_lines} subjects"
        DTI_FILES+=("${filepath}")
        DTI_METRICS+=("${metric}")
    elif [[ -n "${filepath}" ]]; then
        echo "  [SKIP]  ${metric}: file not found (${filepath})"
    else
        echo "  [SKIP]  ${metric}: variable empty or commented out"
    fi
done

if [[ ${#DTI_FILES[@]} -eq 0 ]]; then
    echo "  [ERROR] No DTI phenotype files found. At least FA is required."
    exit 1
fi
echo "  Available metrics: ${DTI_METRICS[*]}"
echo ""

# --- Check subject ID column name across ALL files ---
echo "  Checking subject ID column '${subjID_column}' across files..."
ID_MISMATCH=0

for f in "${DTI_FILES[@]}"; do
    file_header=$(head -1 "${f}")
    if echo "${file_header}" | tr ',' '\n' | tr -d '"' | tr -d "'" | grep -qx "${subjID_column}"; then
        echo "    [OK]  $(basename ${f})"
    else
        first_col=$(echo "${file_header}" | cut -d',' -f1 | tr -d '"' | tr -d "'")
        echo "    [WARNING] $(basename ${f}): column '${subjID_column}' NOT FOUND"
        echo "              First column is '${first_col}' — is this your subject ID?"
        echo "              If so, change subjID_column=\"${first_col}\" in user settings."
        ID_MISMATCH=1
    fi
done

# Check covariates file too
if [[ -f "${covariatesCSV}" ]]; then
    cov_header=$(head -1 "${covariatesCSV}")
    if echo "${cov_header}" | tr ',' '\n' | tr -d '"' | tr -d "'" | grep -qxE "(${subjID_column}|FID|IID)"; then
        echo "    [OK]  Covariates.csv (has ${subjID_column} or FID/IID)"
    else
        first_col=$(echo "${cov_header}" | cut -d',' -f1 | tr -d '"' | tr -d "'")
        echo "    [WARNING] Covariates.csv: no '${subjID_column}', 'FID', or 'IID' column found"
        echo "              First column is '${first_col}'"
        ID_MISMATCH=1
    fi
fi

if [[ ${ID_MISMATCH} -eq 1 ]]; then
    echo ""
    echo "  ╔═══════════════════════════════════════════════════════════════╗"
    echo "  ║  WARNING: Subject ID column names do not match across files. ║"
    echo "  ║  The script will try to use the first column as fallback,    ║"
    echo "  ║  but this may cause mismatches. Please verify your files     ║"
    echo "  ║  have a consistent subject ID column name.                   ║"
    echo "  ╚═══════════════════════════════════════════════════════════════╝"
fi
echo ""

# Check genotype files
FAIL=0
for ext in bed bim fam; do
    if [[ -f "${datafile}.${ext}" ]]; then
        echo "  [OK]    $(basename ${datafile}).${ext}"
    else
        echo "  [ERROR] Not found: ${datafile}.${ext}"
        FAIL=1
    fi
done

if [[ -d "${imputed_dir}" ]]; then
    n_vcf=$(ls "${imputed_dir}"/*.vcf.gz 2>/dev/null | wc -l || echo 0)
    echo "  [OK]    Imputed dir: ${imputed_dir} (${n_vcf} VCF files)"
else
    echo "  [ERROR] Imputed dir not found: ${imputed_dir}"
    FAIL=1
fi

if [[ -f "${covariatesCSV}" ]]; then
    n_cov=$(($(wc -l < "${covariatesCSV}") - 1))
    echo "  [OK]    Covariates.csv — ${n_cov} subjects"
else
    echo "  [ERROR] Covariates.csv not found: ${covariatesCSV}"
    FAIL=1
fi

if [[ -n "${demogCSV}" ]] && [[ ! -f "${demogCSV}" ]]; then
    echo "  [WARN]  Demographics CSV not found: ${demogCSV} — QC will run without dx/site"
    demogCSV=""
fi

if [[ ${FAIL} -eq 1 ]]; then
    echo "  [ERROR] Missing required files."
    exit 1
fi

n_geno=$(wc -l < "${datafile}.fam")
echo ""
echo "  Genotype .fam subjects: ${n_geno}"
echo ""

#=============================================================================
# 3. Harmonise subjects across phenotype, covariate, and genotype files
#=============================================================================
echo "=================================================================="
echo " [3/8] Harmonising subjects across all files"
echo "=================================================================="

MATCH_REPORT="${outFolder}/QC_statistical/${COHORT}_subject_matching_report.txt"

export DTI_FILES_STR=$(IFS=';'; echo "${DTI_FILES[*]}")
export COVAR_CSV="${covariatesCSV}"
export FAM_FILE="${datafile}.fam"
export SUBJID_COL="${subjID_column}"
export MATCH_REPORT_FILE="${MATCH_REPORT}"
export COHORT

"${Rbin}" --no-save --slave 2>&1 <<'RSCRIPT' | tee "${LOGDIR}/${COHORT}_subject_harmonise.log"
dti_files_str <- Sys.getenv("DTI_FILES_STR")
covar_csv     <- Sys.getenv("COVAR_CSV")
fam_file      <- Sys.getenv("FAM_FILE")
id_col        <- Sys.getenv("SUBJID_COL")
report_file   <- Sys.getenv("MATCH_REPORT_FILE")
cohort_name   <- Sys.getenv("COHORT")

dti_files <- unlist(strsplit(dti_files_str, ";"))

cat("=== Subject Harmonisation ===\n\n")

#---------------------------------------------------------------------
# ID normalization: strip common prefixes so "sub-100206" matches "100206"
# Handles: sub-, SUB-, Sub-, sub_, SUB_, Sub_
#---------------------------------------------------------------------
normalize_id <- function(ids) {
  gsub("^[Ss][Uu][Bb][-_]", "", trimws(as.character(ids)))
}

# Track whether normalization changed any IDs (for logging)
id_was_normalized <- FALSE

# Read DTI files and find IDs, with fallback to first column
dti_ids_list     <- list()  # original IDs
dti_ids_norm     <- list()  # normalized IDs
dti_norm_map     <- list()  # normalized -> original mapping per file
id_col_used      <- list()
for (f in dti_files) {
  df <- read.csv(f, header=TRUE, stringsAsFactors=FALSE,
                  na.strings=c("NA","na","","N/A","-999"))
  col <- id_col
  if (!col %in% names(df)) {
    col <- names(df)[1]
    cat(sprintf("  [NOTE] '%s' not in %s, using '%s'\n", id_col, basename(f), col))
  }
  ids <- as.character(df[[col]])
  ids <- ids[!is.na(ids) & ids != ""]
  norm_ids <- normalize_id(ids)
  if (any(norm_ids != ids)) {
    cat(sprintf("  [NOTE] Stripped ID prefixes in %s (e.g., '%s' -> '%s')\n",
                basename(f), ids[norm_ids != ids][1], norm_ids[norm_ids != ids][1]))
    id_was_normalized <- TRUE
  }
  metric <- gsub(".*combinedROItable_|\\.csv","",basename(f))
  dti_ids_list[[metric]] <- ids
  dti_ids_norm[[metric]] <- norm_ids
  # Build normalized -> original map (first occurrence wins)
  nmap <- setNames(ids, norm_ids)
  dti_norm_map[[metric]] <- nmap
  id_col_used[[metric]] <- col
  cat(sprintf("  %s: %d subjects (ID col: %s)\n", metric, length(ids), col))
}

common_dti_norm <- Reduce(intersect, dti_ids_norm)
cat(sprintf("  Common across all DTI files: %d\n", length(common_dti_norm)))

# Covariates
cov <- read.csv(covar_csv, header=TRUE, stringsAsFactors=FALSE,
                 na.strings=c("NA","na","","N/A","-999"))
if ("IID" %in% names(cov)) {
  cov_ids_raw <- as.character(cov[["IID"]])
} else if (id_col %in% names(cov)) {
  cov_ids_raw <- as.character(cov[[id_col]])
} else {
  cov_ids_raw <- as.character(cov[[1]])
  cat(sprintf("  [NOTE] Using first column '%s' for covariate IDs\n", names(cov)[1]))
}
cov_ids_raw <- cov_ids_raw[!is.na(cov_ids_raw) & cov_ids_raw != ""]
cov_ids_norm <- normalize_id(cov_ids_raw)
if (any(cov_ids_norm != cov_ids_raw)) {
  cat(sprintf("  [NOTE] Stripped ID prefixes in Covariates.csv (e.g., '%s' -> '%s')\n",
              cov_ids_raw[cov_ids_norm != cov_ids_raw][1],
              cov_ids_norm[cov_ids_norm != cov_ids_raw][1]))
  id_was_normalized <- TRUE
}
cov_norm_map <- setNames(cov_ids_raw, cov_ids_norm)
cat(sprintf("  Covariate subjects: %d\n", length(cov_ids_raw)))

# Genotype
fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)
fam_iids <- as.character(fam[,2])
fam_fids <- as.character(fam[,1])
fam_fid_iid <- paste(fam_iids, sep="_")
fam_iids_norm <- normalize_id(fam_iids)
fam_fids_norm <- normalize_id(fam_fids)
if (any(fam_iids_norm != fam_iids)) {
  cat(sprintf("  [NOTE] Stripped ID prefixes in .fam IIDs (e.g., '%s' -> '%s')\n",
              fam_iids[fam_iids_norm != fam_iids][1],
              fam_iids_norm[fam_iids_norm != fam_iids][1]))
  id_was_normalized <- TRUE
}
cat(sprintf("  Genotype subjects: %d\n", length(fam_iids)))

if (id_was_normalized) {
  cat("\n  ┌─────────────────────────────────────────────────────────────┐\n")
  cat("  │ ID NORMALIZATION ACTIVE: prefixes like 'sub-' were stripped │\n")
  cat("  │ for matching. Output IDs will use the .fam file's format.   │\n")
  cat("  └─────────────────────────────────────────────────────────────┘\n")
}

# Overlaps (all on normalized IDs)
dti_and_cov_norm <- intersect(common_dti_norm, cov_ids_norm)
dti_only_norm    <- setdiff(common_dti_norm, cov_ids_norm)
cov_only_norm    <- setdiff(cov_ids_norm, common_dti_norm)
cat(sprintf("\n  DTI ∩ Covariates: %d\n", length(dti_and_cov_norm)))
if (length(dti_only_norm)>0) cat(sprintf("  DTI only (no covariate): %d\n", length(dti_only_norm)))
if (length(cov_only_norm)>0) cat(sprintf("  Covariate only (no DTI): %d\n", length(cov_only_norm)))

# Match to genotype using normalized IDs
in_geno_norm <- dti_and_cov_norm[dti_and_cov_norm %in% fam_iids_norm |
                                  dti_and_cov_norm %in% fam_fids_norm]
not_in_geno_norm <- setdiff(dti_and_cov_norm, in_geno_norm)
cat(sprintf("  DTI ∩ Covariates ∩ Genotype: %d\n", length(in_geno_norm)))
if (length(not_in_geno_norm)>0) cat(sprintf("  DTI+Cov but no genotype: %d\n", length(not_in_geno_norm)))

# Map matched normalized IDs back to .fam's original IIDs (for PLINK compatibility)
fam_norm_to_iid <- setNames(fam_iids, fam_iids_norm)
fam_norm_to_fid <- setNames(fam_fids, fam_fids_norm)
in_geno_fam_iid <- character(length(in_geno_norm))
for (i in seq_along(in_geno_norm)) {
  nid <- in_geno_norm[i]
  if (nid %in% fam_iids_norm) {
    in_geno_fam_iid[i] <- fam_norm_to_iid[nid]
  } else {
    in_geno_fam_iid[i] <- fam_norm_to_fid[nid]
  }
}

int_file <- sub("_report.txt","_intersection_IDs.txt", report_file)
writeLines(in_geno_fam_iid, int_file)

# Age
age_col <- NULL
for (ac in c("Age","age","AGE")) {
  if (ac %in% names(cov)) { age_col <- ac; break }
}
if (!is.null(age_col)) {
  ages <- as.numeric(cov[[age_col]])
  ages <- ages[!is.na(ages)]
  cat(sprintf("\n  Age: min=%.1f, max=%.1f, mean=%.1f\n", min(ages), max(ages), mean(ages)))
}

# Write report
sink(report_file)
cat(sprintf("ENIGMA-DTI GWAS Subject Matching Report\nDate: %s\nCohort: %s\n\n", Sys.time(), cohort_name))
if (id_was_normalized) cat("NOTE: ID prefixes (e.g., 'sub-') were stripped for matching.\n\n")
cat("=== SUBJECT COUNTS ===\n")
for (m in names(dti_ids_list)) cat(sprintf("  %s phenotype: %d (ID col: %s)\n", m, length(dti_ids_list[[m]]), id_col_used[[m]]))
cat(sprintf("  Common across DTI: %d\n  Covariates: %d\n  Genotype .fam: %d\n\n",
            length(common_dti_norm), length(cov_ids_raw), length(fam_iids)))
cat("=== OVERLAP ===\n")
cat(sprintf("  Phenotype w/ matching covariate: %d / %d\n", length(dti_and_cov_norm), length(common_dti_norm)))
cat(sprintf("  Phenotype w/ matching genotype:  %d / %d\n", length(in_geno_norm), length(common_dti_norm)))
cat(sprintf("  Covariate w/ matching phenotype: %d / %d\n", length(dti_and_cov_norm), length(cov_ids_raw)))
cat(sprintf("  Covariate w/ matching genotype:  %d / %d\n",
            length(intersect(cov_ids_norm, c(fam_iids_norm, fam_fids_norm))), length(cov_ids_raw)))
cat(sprintf("\n  ANALYSIS-READY (all 3 match): %d\n", length(in_geno_norm)))
if (!is.null(age_col)) cat(sprintf("\n=== AGE ===\n  Min: %.1f\n  Max: %.1f\n  Mean: %.1f\n", min(ages), max(ages), mean(ages)))
if (length(dti_only_norm)>0) { cat(sprintf("\n=== DTI w/o covariates (%d) ===\n",length(dti_only_norm))); for(id in dti_only_norm) cat(sprintf("  %s\n",id)) }
if (length(not_in_geno_norm)>0) { cat(sprintf("\n=== DTI+Cov w/o genotype (%d) ===\n",length(not_in_geno_norm))); for(id in not_in_geno_norm) cat(sprintf("  %s\n",id)) }
sink()
cat(sprintf("\n  Report: %s\n", report_file))
RSCRIPT

echo ""

#=============================================================================
# 4. ENIGMA-DTI Statistical QC (eDTI_outliers) — BEFORE relatedness filter
#=============================================================================
echo "=================================================================="
echo " [4/8] ENIGMA-DTI Statistical QC (eDTI_outliers)"
echo "=================================================================="

qc_stat_dir="${outFolder}/QC_statistical"
qc_prefix="${qc_stat_dir}/${COHORT}_eDTI_statQC"
OUTLIER_IDS_FILE="${qc_stat_dir}/${COHORT}_statistical_outlier_IDs.txt"

if [[ ${run_statistical_qc} -eq 1 ]]; then

    # --- Install eDTI_outliers if needed ---
    # Strategy: try to use it directly first; if not found, pip install
    # without creating a full conda environment (avoids long env creation).

    if command -v eDTI_outliers &>/dev/null; then
        echo "  [OK] eDTI_outliers already available on PATH."
    else
        echo "  eDTI_outliers not found on PATH. Attempting installation..."
        echo ""

        # Clone the repo if not already present
        edti_qc_dir="${outFolder}/_eDTI_additional_QC_install"
        if [[ ! -d "${edti_qc_dir}/eDTI_additional_QC" ]]; then
            echo "  Cloning eDTI_additional_QC repository..."
            mkdir -p "${edti_qc_dir}"
            git clone https://github.com/lizhaddad/eDTI_additional_QC.git \
                "${edti_qc_dir}/eDTI_additional_QC" 2>&1 | tail -3
        else
            echo "  Repository already cloned."
        fi

        # Try pip install directly (no conda env needed — just Python + deps)
        echo "  Installing via pip (user-level)..."
        PIP_CMD=""
        if command -v pip3 &>/dev/null; then
            PIP_CMD="pip3"
        elif command -v pip &>/dev/null; then
            PIP_CMD="pip"
        fi

        if [[ -n "${PIP_CMD}" ]]; then
            # Install required dependencies first (the package's setup.py
            # may not pull these automatically with --user)
            echo "  Installing dependencies: pandas, numpy, openpyxl..."
            ${PIP_CMD} install --user pandas numpy openpyxl xlsxwriter \
                2>&1 | tee "${LOGDIR}/${COHORT}_pip_deps.log" | tail -3

            echo "  Installing eDTI_outliers package..."
            ${PIP_CMD} install --user "${edti_qc_dir}/eDTI_additional_QC" \
                2>&1 | tee "${LOGDIR}/${COHORT}_pip_install.log" | tail -5

            # Add user site-packages bin to PATH
            USER_BASE=$(python3 -m site --user-base 2>/dev/null || python -m site --user-base 2>/dev/null || echo "")
            if [[ -n "${USER_BASE}" ]] && [[ -d "${USER_BASE}/bin" ]]; then
                export PATH="${USER_BASE}/bin:${PATH}"
            fi
        fi

        # If pip install didn't work, try conda as fallback
        if ! command -v eDTI_outliers &>/dev/null; then
            echo "  pip install did not put eDTI_outliers on PATH."
            echo "  Trying conda environment as fallback..."

            if command -v conda &>/dev/null; then
                # Source conda
                CONDA_EXE=$(which conda)
                CONDA_BASE=$(dirname $(dirname "${CONDA_EXE}"))
                if [[ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
                    source "${CONDA_BASE}/etc/profile.d/conda.sh"
                fi

                if ! conda env list 2>/dev/null | grep -q "ENIGMA_DTI_env"; then
                    echo "  Creating ENIGMA_DTI_env (this may take a few minutes)..."
                    conda env create -f "${edti_qc_dir}/eDTI_additional_QC/environment.yaml" \
                        2>&1 | tee "${LOGDIR}/${COHORT}_conda_create.log" | tail -5
                fi
                conda activate ENIGMA_DTI_env 2>/dev/null || true
                pip install "${edti_qc_dir}/eDTI_additional_QC" 2>&1 | tail -3
            fi
        fi

        # Final check
        if command -v eDTI_outliers &>/dev/null; then
            echo "  [OK] eDTI_outliers installed successfully."
        else
            echo "  [ERROR] Could not install eDTI_outliers."
            echo "          You can install it manually:"
            echo "            git clone https://github.com/lizhaddad/eDTI_additional_QC.git"
            echo "            pip install ./eDTI_additional_QC"
            echo "          Then re-run this script."
            echo ""
            echo "  Continuing without statistical QC..."
            run_statistical_qc=0
        fi
    fi
fi

if [[ ${run_statistical_qc} -eq 1 ]]; then

    # Create DTI inputs file
    DTI_INPUTS_FILE="${qc_stat_dir}/${COHORT}_DTI_inputs.txt"
    > "${DTI_INPUTS_FILE}"
    for f in "${DTI_FILES[@]}"; do
        echo "${f}" >> "${DTI_INPUTS_FILE}"
    done
    echo "  DTI inputs: ${DTI_INPUTS_FILE}"
    cat "${DTI_INPUTS_FILE}" | sed 's/^/    /'

    # Build command as array (avoids string/eval quoting issues)
    # NOTE: eDTI_outliers defaults are already what we want:
    #   --quant1 0.001 0.999  --perc1 0.05
    #   --quant2 0.005 0.995  --perc2 0.10
    # Only pass these flags if the user changed them from defaults.

    declare -a QC_ARGS=()
    QC_ARGS+=(--subjID "${subjID_column}")
    QC_ARGS+=(--DTIinputs "${DTI_INPUTS_FILE}")
    QC_ARGS+=(--output "${qc_prefix}")

    # Only pass threshold args if non-default (avoids argparse type bugs)
    if [[ "${quant1_lower}" != "0.001" ]] || [[ "${quant1_upper}" != "0.999" ]]; then
        QC_ARGS+=(--quant1 "${quant1_lower}" "${quant1_upper}")
    fi
    if [[ "${perc1}" != "0.05" ]]; then
        QC_ARGS+=(--perc1 "${perc1}")
    fi
    if [[ "${quant2_lower}" != "0.005" ]] || [[ "${quant2_upper}" != "0.995" ]]; then
        QC_ARGS+=(--quant2 "${quant2_lower}" "${quant2_upper}")
    fi
    if [[ "${perc2}" != "0.10" ]]; then
        QC_ARGS+=(--perc2 "${perc2}")
    fi

    if [[ -n "${demogCSV}" ]] && [[ -f "${demogCSV}" ]]; then
        QC_ARGS+=(--demogCSV "${demogCSV}")
        [[ -n "${dxColumn}" ]]   && QC_ARGS+=(--dx "${dxColumn}")
        [[ -n "${siteColumn}" ]] && QC_ARGS+=(--site "${siteColumn}")
    fi

    echo ""
    echo "  Running: eDTI_outliers ${QC_ARGS[*]}"
    echo ""

    # Fix known bug: numeric subject IDs cause "expected str, got numpy.int64"
    # in eDTI_outliers .join() calls. We patch the installed source file with
    # Python (handles \r\n natively) and clear __pycache__ so the fix takes effect.
    python3 -c "
import os, sys, shutil
try:
    from eDTI_outliers.scripts import cli
    f = cli.__file__
    with open(f) as fh:
        src = fh.read()
    old = '.join(outlier_subjects_array)'
    new = '.join([str(x) for x in outlier_subjects_array])'
    if old in src:
        with open(f, 'w') as fh:
            fh.write(src.replace(old, new))
        pycache = os.path.join(os.path.dirname(f), '__pycache__')
        if os.path.isdir(pycache):
            shutil.rmtree(pycache)
        print(f'  [FIX] Patched numpy.int64 bug in {f}')
    else:
        print(f'  [OK] {f} already patched')
except Exception as e:
    print(f'  [WARN] Could not patch eDTI_outliers: {e}', file=sys.stderr)
"

    eDTI_outliers "${QC_ARGS[@]}" 2>&1 | tee "${LOGDIR}/${COHORT}_eDTI_outliers.log"
    echo ""

    # ---------------------------------------------------------------
    # TWO-PASS OUTLIER REMOVAL
    #
    # Pass 1: Use the SUMMARY file (default eDTI_outliers flagging).
    #         If <5% of total sample is flagged, remove all of them.
    #
    # Pass 2: If >=5% is flagged, the default criteria are too
    #         aggressive for this dataset. Tighten to only remove
    #         subjects meeting EITHER of:
    #           A) "Summary (Across Column Total)" = 2
    #           B) Any value in "5.0% of measures >0.001 or <0.999"
    #         (union of those two columns from the byCriteria CSV)
    # ---------------------------------------------------------------

    FLAGGED_ALL_FILE="${qc_stat_dir}/${COHORT}_all_flagged_IDs.txt"
    FLAGGED_REVIEW_FILE="${qc_stat_dir}/${COHORT}_flagged_review_only_IDs.txt"

    # --- Pass 1: Get ALL flagged from SUMMARY ---
    SUMMARY_FILE=$(find "${qc_stat_dir}" -name "*SUMMARY*subjects_flagged*" -type f 2>/dev/null | head -1)
    if [[ -n "${SUMMARY_FILE}" ]] && [[ -f "${SUMMARY_FILE}" ]]; then
        echo "  Pass 1: Parsing default flagged subjects from SUMMARY..."
        grep -v "^$\|^The following\|^We recommend\|flagged\|outlier\|criteria\|inspect\|subject\|categoriz\|default\|specified" \
            "${SUMMARY_FILE}" 2>/dev/null | \
            awk 'NF==1 {print $1}' | \
            sed 's/^"//;s/"$//' | \
            sort -u > "${FLAGGED_ALL_FILE}" || true
    else
        echo "  [INFO] No SUMMARY file found."
        touch "${FLAGGED_ALL_FILE}"
    fi
    sed -i '/^$/d' "${FLAGGED_ALL_FILE}" 2>/dev/null || true
    n_flagged_all=$(wc -l < "${FLAGGED_ALL_FILE}" 2>/dev/null || echo 0)

    # Get total N from first DTI file
    n_total_pheno=$(($(wc -l < "${DTI_FILES[0]}") - 1))
    pct_flagged=0
    if [[ ${n_total_pheno} -gt 0 ]]; then
        pct_flagged=$(echo "scale=1; ${n_flagged_all}*100/${n_total_pheno}" | bc 2>/dev/null || echo "0")
    fi

    echo "  Total subjects in phenotype file: ${n_total_pheno}"
    echo "  Default flagged (SUMMARY):        ${n_flagged_all} (${pct_flagged}%)"
    echo ""

    # Decide: Pass 1 or Pass 2?
    use_pass2=0
    if [[ ${n_total_pheno} -gt 0 ]] && (( n_flagged_all * 100 >= n_total_pheno * 5 )); then
        use_pass2=1
    fi

    if [[ ${use_pass2} -eq 0 ]]; then
        # --- Pass 1 result: <5% flagged, remove all ---
        cp "${FLAGGED_ALL_FILE}" "${OUTLIER_IDS_FILE}"
        n_flagged=${n_flagged_all}
        touch "${FLAGGED_REVIEW_FILE}"

        echo "  Pass 1 result: <5% flagged. Removing all ${n_flagged} flagged subjects."
        echo ""
        echo "  ┌─────────────────────────────────────────────────────────────┐"
        echo "  │ STATISTICAL QC: ${n_flagged} subjects will be removed"
        echo "  └─────────────────────────────────────────────────────────────┘"
        if [[ ${n_flagged} -gt 0 ]]; then
            echo "  Subjects to be REMOVED:"
            cat "${OUTLIER_IDS_FILE}" | sed 's/^/    /'
        else
            echo "  No outliers detected. All subjects pass QC."
        fi

    else
        # --- Pass 2: >=5% flagged, tighten criteria ---
        echo "  ╔═══════════════════════════════════════════════════════════════╗"
        echo "  ║  >=5% of sample flagged by default criteria (${pct_flagged}%)."
        echo "  ║  Applying stricter removal rule:                             ║"
        echo "  ║    Remove = UNION of:                                        ║"
        echo "  ║      • Summary (Across Column Total) = 2                     ║"
        echo "  ║      • Any value in '5.0% of measures >0.001 or <0.999'      ║"
        echo "  ║  Other flagged subjects will be kept for review only.         ║"
        echo "  ╚═══════════════════════════════════════════════════════════════╝"
        echo ""

        SUPP_CSV=$(find "${qc_stat_dir}" -name "*subjects_flagged_byCriteria*" -type f 2>/dev/null | head -1)

        if [[ -n "${SUPP_CSV}" ]] && [[ -f "${SUPP_CSV}" ]]; then
            echo "  Parsing byCriteria CSV: $(basename ${SUPP_CSV})"
            echo ""

            python3 -c "
import csv, sys

csv_file = '${SUPP_CSV}'
out_remove = '${OUTLIER_IDS_FILE}'
out_review = '${FLAGGED_REVIEW_FILE}'
flagged_file = '${FLAGGED_ALL_FILE}'

# Read flagged IDs from Pass 1
with open(flagged_file) as f:
    flagged_ids = set(line.strip() for line in f if line.strip())

# This CSV has a 3-row header:
#   Row 1: subjectID | Summary (Across Column Total) | Across Diagnosis | ...
#   Row 2:           |                               | Across Sites     | ...
#   Row 3:           |                               | 10.0% of ...     | 5.0% of ...
with open(csv_file, 'r') as f:
    reader = csv.reader(f)
    rows = list(reader)

if len(rows) < 4:
    print('  [ERROR] CSV has fewer than 4 rows. Cannot parse.')
    sys.exit(1)

row1, row2, row3 = rows[0], rows[1], rows[2]
n_cols = len(row1)

# Merge 3 header rows
combined_headers = []
for i in range(n_cols):
    parts = []
    for r in [row1, row2, row3]:
        val = r[i].strip() if i < len(r) else ''
        if val:
            parts.append(val)
    combined_headers.append(' | '.join(parts) if parts else f'col_{i}')

print('  Combined column headers (from 3-row header):')
for i, h in enumerate(combined_headers):
    print(f'    [{i}] {h}')
print()

# Find columns
summary_idx = None
strict_idx = None
for i, h in enumerate(combined_headers):
    if 'summary' in h.lower() and 'across' in h.lower():
        summary_idx = i
    if '0.001' in h and '0.999' in h:
        strict_idx = i

if summary_idx is not None:
    print(f'  Summary column:  [{summary_idx}] {combined_headers[summary_idx]}')
else:
    print('  [ERROR] Cannot find Summary column.')
    sys.exit(1)
if strict_idx is not None:
    print(f'  Strict column:   [{strict_idx}] {combined_headers[strict_idx]}')
else:
    print('  [WARN] Cannot find 5.0% >0.001/<0.999 column. Using Summary only.')
print()

# Parse data rows (row index 3+)
remove_ids = set()
remove_summary2 = set()
remove_strict = set()
all_csv_ids = set()

for row in rows[3:]:
    if not row or not row[0].strip():
        continue
    subj = row[0].strip().strip('\"')
    all_csv_ids.add(subj)

    # Condition A: Summary == 2
    if summary_idx is not None and summary_idx < len(row):
        if row[summary_idx].strip() == '2':
            remove_summary2.add(subj)
            remove_ids.add(subj)

    # Condition B: strict column has any non-empty, non-zero value
    if strict_idx is not None and strict_idx < len(row):
        val = row[strict_idx].strip()
        if val and val != '0' and val.lower() not in ('nan', 'na'):
            remove_strict.add(subj)
            remove_ids.add(subj)

# Write removal list
with open(out_remove, 'w') as f:
    for sid in sorted(remove_ids):
        f.write(sid + '\n')

# Review = flagged by SUMMARY but not meeting strict removal
review = flagged_ids - remove_ids
with open(out_review, 'w') as f:
    for sid in sorted(review):
        f.write(sid + '\n')

print(f'  Subjects in byCriteria CSV:        {len(all_csv_ids)}')
print(f'  Summary (Across Column Total) = 2: {len(remove_summary2)}')
if strict_idx is not None:
    print(f'  5.0% measures >0.001/<0.999:       {len(remove_strict)}')
print(f'  UNION (to be REMOVED):             {len(remove_ids)}')
print(f'  Default-flagged kept for review:   {len(review)}')
" 2>&1 | tee "${LOGDIR}/${COHORT}_qc_criteria_parse.log"

            sed -i '/^$/d' "${OUTLIER_IDS_FILE}" 2>/dev/null || true
            n_flagged=$(wc -l < "${OUTLIER_IDS_FILE}" 2>/dev/null || echo 0)
            n_review=$(wc -l < "${FLAGGED_REVIEW_FILE}" 2>/dev/null || echo 0)

        else
            # No byCriteria CSV — fall back to removing all from SUMMARY
            echo "  [WARN] No byCriteria CSV found. Removing all ${n_flagged_all} flagged subjects."
            cp "${FLAGGED_ALL_FILE}" "${OUTLIER_IDS_FILE}"
            n_flagged=${n_flagged_all}
            n_review=0
            touch "${FLAGGED_REVIEW_FILE}"
        fi

        echo ""
        if [[ ${n_flagged} -gt 0 ]]; then
            echo "  ╔═══════════════════════════════════════════════════════════════╗"
            echo "  ║  WARNING: ${n_flagged} SUBJECTS WILL BE REMOVED"
            echo "  ║                                                               ║"
            echo "  ║  These subjects met at least one of:                          ║"
            echo "  ║    • Summary (Across Column Total) = 2                        ║"
            echo "  ║    • Any value in '5.0% of measures >0.001 or <0.999'         ║"
            echo "  ║                                                               ║"
            echo "  ║  ${n_review} additional subjects flagged by default criteria"
            echo "  ║  are KEPT but noted for visual review.                        ║"
            echo "  ║                                                               ║"
            echo "  ║  PLEASE VISUALLY INSPECT flagged subjects before proceeding.  ║"
            echo "  ╚═══════════════════════════════════════════════════════════════╝"
            echo ""
            echo "  Subjects to be REMOVED:"
            cat "${OUTLIER_IDS_FILE}" | sed 's/^/    /'
            if [[ ${n_review} -gt 0 ]]; then
                echo ""
                echo "  Subjects for REVIEW only (NOT removed):"
                cat "${FLAGGED_REVIEW_FILE}" | sed 's/^/    /'
            fi
        else
            echo "  ┌─────────────────────────────────────────────────────────────┐"
            echo "  │ No subjects meet the strict removal criteria.               │"
            echo "  │ ${n_flagged_all} subjects flagged for review only.          │"
            echo "  └─────────────────────────────────────────────────────────────┘"
            touch "${OUTLIER_IDS_FILE}"
        fi
    fi

    # Deactivate conda if it was activated
    conda deactivate 2>/dev/null || true

else
    echo "  Statistical QC skipped (run_statistical_qc=0)."
    n_flagged=0
    touch "${OUTLIER_IDS_FILE}"
fi
echo ""

#=============================================================================
# 4b. Visual QC on flagged subjects (optional)
#=============================================================================
if [[ ${run_visual_qc} -eq 1 ]] && [[ ${n_flagged} -gt 0 ]]; then
    echo "=================================================================="
    echo " [4b/8] Visual QC on flagged subjects (eDTI_visualQC)"
    echo "=================================================================="

    if [[ -z "${TBSS_dir}" ]] || [[ ! -d "${TBSS_dir}" ]]; then
        echo "  [SKIP] TBSS_dir not set or not found: '${TBSS_dir}'"
        echo "         Set TBSS_dir in user settings to enable visual QC."
        echo "         Visual QC requires the TBSS directory from the ENIGMA-DTI pipeline."
    else
        # Check if the visual_qc script exists in the cloned repo
        edti_qc_dir="${outFolder}/_eDTI_additional_QC_install"
        VIS_QC_SCRIPT=""

        # Look in the cloned repo
        if [[ -d "${edti_qc_dir}/eDTI_additional_QC/visual_qc" ]]; then
            VIS_QC_SCRIPT=$(find "${edti_qc_dir}/eDTI_additional_QC/visual_qc" -name "*.sh" -o -name "*.py" 2>/dev/null | head -1)
        fi

        # If not cloned yet, clone it
        if [[ -z "${VIS_QC_SCRIPT}" ]]; then
            mkdir -p "${edti_qc_dir}"
            if [[ ! -d "${edti_qc_dir}/eDTI_additional_QC" ]]; then
                echo "  Cloning eDTI_additional_QC for visual QC..."
                git clone https://github.com/lizhaddad/eDTI_additional_QC.git \
                    "${edti_qc_dir}/eDTI_additional_QC" 2>&1 | tail -3
            fi
            VIS_QC_SCRIPT=$(find "${edti_qc_dir}/eDTI_additional_QC/visual_qc" -name "*.sh" -o -name "*.py" 2>/dev/null | head -1)
        fi

        if [[ -n "${VIS_QC_SCRIPT}" ]]; then
            vis_qc_outdir="${qc_stat_dir}/visual_qc_output"
            mkdir -p "${vis_qc_outdir}"

            echo "  Visual QC script: ${VIS_QC_SCRIPT}"
            echo "  TBSS directory: ${TBSS_dir}"
            echo "  Flagged subjects: ${n_flagged}"
            echo "  Output: ${vis_qc_outdir}"
            echo ""

            # Run visual QC for each flagged subject
            if [[ "${VIS_QC_SCRIPT}" == *.sh ]]; then
                chmod +x "${VIS_QC_SCRIPT}"
                while IFS= read -r subj_id; do
                    [[ -z "${subj_id}" ]] && continue
                    echo "  Running visual QC for: ${subj_id}"
                    bash "${VIS_QC_SCRIPT}" \
                        "${TBSS_dir}" \
                        "${subj_id}" \
                        "${vis_qc_outdir}" \
                        2>&1 | sed 's/^/    /' || \
                        echo "    [WARN] Visual QC failed for ${subj_id}"
                done < "${OUTLIER_IDS_FILE}"
            elif [[ "${VIS_QC_SCRIPT}" == *.py ]]; then
                while IFS= read -r subj_id; do
                    [[ -z "${subj_id}" ]] && continue
                    echo "  Running visual QC for: ${subj_id}"
                    python3 "${VIS_QC_SCRIPT}" \
                        --tbss_dir "${TBSS_dir}" \
                        --subject "${subj_id}" \
                        --output "${vis_qc_outdir}" \
                        2>&1 | sed 's/^/    /' || \
                        echo "    [WARN] Visual QC failed for ${subj_id}"
                done < "${OUTLIER_IDS_FILE}"
            fi

            n_pngs=$(find "${vis_qc_outdir}" -name "*.png" 2>/dev/null | wc -l)
            echo ""
            echo "  Visual QC complete: ${n_pngs} PNG images generated."
            echo "  Please inspect these images in: ${vis_qc_outdir}/"
            echo "  Include them with your QC files to send to the ENIGMA-DTI team."
        else
            echo "  [WARN] No visual QC script found in eDTI_additional_QC/visual_qc/"
            echo "         Check: https://github.com/lizhaddad/eDTI_additional_QC"
        fi
    fi
    echo ""
elif [[ ${run_visual_qc} -eq 1 ]] && [[ ${n_flagged} -eq 0 ]]; then
    echo "  Visual QC: skipped (no subjects flagged by statistical QC)."
    echo ""
fi

#=============================================================================
# 5. Relatedness check & filter
#=============================================================================
echo "=================================================================="
echo " [5/8] Relatedness check"
echo "=================================================================="

if [[ ${cohort_unrelated} -eq 0 ]]; then
    echo "  [STOP] You indicated this is a RELATED/FAMILY cohort."
    echo "         Please contact enigma.DTIgenetics@gmail.com for the"
    echo "         family-based association protocol."
    echo "         This protocol is for UNRELATED samples only."
    exit 1
fi

echo "  Cohort type: unrelated"
echo "  KING kinship cutoff: 0.0884 (first cousins and closer)"

datafile_QC4="${datafile%_QC3}_QC4"

plink2 --bfile "${datafile}" \
       --king-cutoff 0.0884 \
       --make-bed \
       --out "${datafile_QC4}" \
       2>&1 | tee "${LOGDIR}/${COHORT}_plink2_king.log"

awk '{print $1, $2}' "${datafile_QC4}.fam" > "${datafile_QC4}_keeplist.txt"

n_after_king=$(wc -l < "${datafile_QC4}.fam")
n_removed_king=$((n_geno - n_after_king))
pct_king=$(echo "scale=1; ${n_removed_king}*100/${n_geno}" | bc 2>/dev/null || echo "?")

echo ""
echo "  Before: ${n_geno}  |  After: ${n_after_king}  |  Removed: ${n_removed_king} (${pct_king}%)"

if (( n_removed_king * 10 > n_geno )); then
    echo "  [WARNING] >10% removed. Contact enigma.DTIgenetics@gmail.com if unexpected."
fi

cat > "${qc_stat_dir}/${COHORT}_relatedness_stats.txt" <<EOF
ENIGMA-DTI GWAS Relatedness Check
Date: $(date)
Genotype file: ${datafile}
KING cutoff: 0.0884
Before: ${n_geno}
After: ${n_after_king}
Removed: ${n_removed_king} (${pct_king}%)
EOF
echo ""

#=============================================================================
# 6. Convert imputed VCF to PGEN
#=============================================================================
echo "=================================================================="
echo " [6/8] Converting VCF to PGEN"
echo "=================================================================="

# Set pgen_dir from user setting (default to ${imputed_dir}/pgen if empty)
if [[ -n "${pgen_outdir}" ]]; then
    pgen_dir="${pgen_outdir}"
else
    pgen_dir="${imputed_dir}/pgen"
fi
mkdir -p "${pgen_dir}"
echo "  PGEN output directory: ${pgen_dir}"

chmod +x "${run_directory}/convert_dosage_vcf_to_pgen.sh"
"${run_directory}/convert_dosage_vcf_to_pgen.sh" \
    "${imputed_dir}" "${pgen_dir}" "${datafile_QC4}_keeplist.txt" "${datafile_QC4}.fam" "${genome_build}" \
    2>&1 | tee "${LOGDIR}/${COHORT}_vcf_to_pgen.log"
echo ""

#=============================================================================
# 7. Merge DTI phenotype files (excluding stat QC outliers)
#=============================================================================
echo "=================================================================="
echo " [7/8] Merging DTI phenotypes (with QC filtering)"
echo "=================================================================="

localfamFILE="${datafile_QC4}.fam"

DTI_FILES_ARG=$(IFS=';'; echo "${DTI_FILES[*]}")
DTI_METRICS_ARG=$(IFS=';'; echo "${DTI_METRICS[*]}")

(  "${Rbin}" --no-save --slave --args \
    "${DTI_FILES_ARG}" \
    "${DTI_METRICS_ARG}" \
    "${localfamFILE}" \
    "${outFolder}" \
    "${OUTLIER_IDS_FILE}" \
    "${subjID_column}" \
    < "${run_directory}/ENIGMA_DTI_pheno_merge.R"  ) 2>&1 | tee "${LOGDIR}/${COHORT}_pheno_merge.log"
echo ""

#=============================================================================
# 8. Site-level statistics & final summary
#=============================================================================
echo "=================================================================="
echo " [8/8] Site-level statistics & final summary"
echo "=================================================================="

n_final=0
if [[ -f "${outFolder}/Phenotypes.csv" ]]; then
    n_final=$(($(wc -l < "${outFolder}/Phenotypes.csv") - 1))
fi

# Site-level stats
if [[ -n "${site_dummies}" ]] || [[ -n "${siteColumn}" ]]; then
    echo "  Generating site-level statistics..."
    export PHENO_FILE="${outFolder}/Phenotypes.csv"
    export COVAR_CSV="${covariatesCSV}"
    export SITE_DUMMIES="${site_dummies}"
    export SITE_COLUMN="${siteColumn}"
    export SUBJID_COL="${subjID_column}"
    export COHORT
    export OUT_DIR="${qc_stat_dir}"

    "${Rbin}" --no-save --slave 2>&1 <<'RSITE' | tee "${LOGDIR}/${COHORT}_site_stats.log"
pheno_file   <- Sys.getenv("PHENO_FILE")
covar_csv    <- Sys.getenv("COVAR_CSV")
site_dummies <- Sys.getenv("SITE_DUMMIES")
site_col     <- Sys.getenv("SITE_COLUMN")
id_col       <- Sys.getenv("SUBJID_COL")
cohort       <- Sys.getenv("COHORT")
out_dir      <- Sys.getenv("OUT_DIR")

pheno <- read.csv(pheno_file, header=TRUE, stringsAsFactors=FALSE)
covar <- read.csv(covar_csv, header=TRUE, stringsAsFactors=FALSE,
                   na.strings=c("NA","na","","N/A","-999"))

site_vec <- NULL
if (nchar(site_col)>0 && site_col %in% names(covar)) {
  if ("IID" %in% names(covar)) m <- match(pheno$IID, covar$IID)
  else m <- match(pheno$IID, covar[[id_col]])
  site_vec <- covar[[site_col]][m]
} else if (nchar(site_dummies)>0) {
  dummy_cols <- unlist(strsplit(trimws(site_dummies), "\\s+"))
  dummy_cols <- dummy_cols[dummy_cols %in% names(covar)]
  if (length(dummy_cols)>0) {
    if ("IID" %in% names(covar)) m <- match(pheno$IID, covar$IID)
    else m <- match(pheno$IID, covar[[id_col]])
    site_mat <- covar[m, dummy_cols, drop=FALSE]
    site_vec <- apply(site_mat, 1, function(row) {
      idx <- which(row==1)
      if (length(idx)>0) return(dummy_cols[idx[1]])
      return("Reference")
    })
  }
}

if (!is.null(site_vec) && length(site_vec)==nrow(pheno)) {
  pheno$SITE <- site_vec
  sites <- sort(unique(pheno$SITE))
  cat(sprintf("  Sites: %d\n", length(sites)))
  site_counts <- table(pheno$SITE)
  cat("\n  Site counts:\n")
  for (s in names(site_counts)) {
    flag <- ifelse(site_counts[s]<10, " *** WARNING: N<10 ***", "")
    cat(sprintf("    %s: %d%s\n", s, site_counts[s], flag))
  }
  pheno_cols <- setdiff(names(pheno), c("FID","IID","SITE"))
  results <- list()
  for (s in sites) {
    sub <- pheno[pheno$SITE==s, ]
    for (pc in pheno_cols) {
      vals <- as.numeric(sub[[pc]]); vals <- vals[!is.na(vals)]
      if (length(vals)>0) {
        q <- quantile(vals, c(0.25,0.75))
        results[[length(results)+1]] <- data.frame(
          Site=s, N=nrow(sub), Phenotype=pc,
          Mean=mean(vals), SD=sd(vals), Min=min(vals), Max=max(vals),
          Q1=q[1], Q3=q[2], IQR=q[2]-q[1], stringsAsFactors=FALSE)
      }
    }
  }
  if (length(results)>0) {
    res_df <- do.call(rbind, results)
    out_file <- file.path(out_dir, paste0(cohort,"_site_summary_statistics.csv"))
    write.csv(res_df, file=out_file, row.names=FALSE)
    cat(sprintf("\n  Written: %s\n", out_file))
    small <- names(site_counts[site_counts<10])
    if (length(small)>0) cat(sprintf("\n  [WARNING] Sites with N<10: %s\n", paste(small,collapse=", ")))
  }
} else {
  cat("  No site info available.\n")
}
RSITE
    echo ""
fi

# Final summary
echo "=================================================================="
echo " Step 0 — COMPLETE"
echo " Finished: $(date)"
echo "=================================================================="
echo ""
echo " SUBJECT FLOW"
echo " ────────────────────────────────────────────────────────"
echo "  Genotype subjects (QC3 .fam):          ${n_geno}"
echo "  Statistical QC outliers removed:      - ${n_flagged}"
echo "  Related individuals removed (KING):   - ${n_removed_king}"
echo "  Final in Phenotypes.csv:                ${n_final}"
echo " ────────────────────────────────────────────────────────"

if [[ ${n_final} -lt 75 ]]; then
    echo ""
    echo " ╔═══════════════════════════════════════════════════════════════╗"
    echo " ║  WARNING: Only ${n_final} subjects available for analysis.         ║"
    echo " ║  Below recommended minimum of 75. Contact ENIGMA-DTI team.  ║"
    echo " ╚═══════════════════════════════════════════════════════════════╝"
fi

echo ""
echo " Logs: ${LOGDIR}/"
ls -1 "${LOGDIR}/" | sed 's/^/   /'
echo ""
echo " QC files to send: ${qc_stat_dir}/"
ls -1 "${qc_stat_dir}/" 2>/dev/null | sed 's/^/   /'
echo ""
echo " Next: Edit and run run1_ENIGMA-DTI_GWAS_step1.sh"
echo "=================================================================="
