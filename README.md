# ENIGMA-DTI GWAS Pipeline

**PLINK2-based genome-wide association analysis for DTI white matter phenotypes**

This pipeline runs GWAS on diffusion tensor imaging (DTI) phenotypes — fractional anisotropy (FA), mean diffusivity (MD), axial diffusivity (AD), and radial diffusivity (RD) — across 25 white matter ROIs from the JHU Atlas. It is distributed by the ENIGMA-DTI Genetics Working Group and designed to be run independently at each cohort site before meta-analysis.

The pipeline handles phenotype QC, covariate formatting, VCF-to-PGEN conversion, and GWAS execution, producing up to 600 association models per cohort.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Installation](#installation)
3. [Directory Setup](#directory-setup)
4. [Preparing Your Input Files](#preparing-your-input-files)
5. [Step 0 — Phenotype QC & Genotype Prep (run0)](#step-0--phenotype-qc--genotype-prep-run0)
6. [Step 1 — Covariate Formatting & QC (run1)](#step-1--covariate-formatting--qc-run1)
7. [QC Checkpoint — Send Results for Review](#qc-checkpoint--send-results-for-review)
8. [Step 2 — Run GWAS (run2)](#step-2--run-gwas-run2)
9. [Packaging Results](#packaging-results)
10. [GWAS Models](#gwas-models)
11. [Troubleshooting](#troubleshooting)
12. [Repository Structure](#repository-structure)
13. [Contact & Citation](#contact--citation)

---

## Prerequisites

Before starting, you should have completed:

1. **ENIGMA-DTI TBSS pipeline** — producing `combinedROItable_FA.csv`, `combinedROItable_MD.csv`, `combinedROItable_AD.csv`, `combinedROItable_RD.csv`
   - [ENIGMA-DTI TBSS Protocol](https://github.com/ENIGMA-git/ENIGMA-DTI-TBSS-Protocol)

2. **Pre-imputation genotype QC** — producing QC'd PLINK files (`${datafile}_QC3.bed/.bim/.fam`) and ancestry PCA (`${datafile}_PCACovariates.eigenvec`)
   - [ENIGMA Genetics QC Protocol](https://github.com/ENIGMA-git/Genetics)

3. **Genotype imputation** — producing imputed VCF files (`chr*.dose.vcf.gz`) from the Michigan Imputation Server or TOPMed

### Required Software

| Software | Purpose | Link |
|----------|---------|------|
| PLINK2 | GWAS execution | https://www.cog-genomics.org/plink/2.0/ |
| R ≥ 4.0 | Statistical processing | https://cran.r-project.org/ |
| bcftools + tabix | VCF manipulation | https://www.htslib.org/ |
| vcftools | VCF utilities | https://vcftools.github.io/ |
| Conda/Miniconda | Environment management | https://docs.conda.io/ |
| git | Cloning this repo | Standard on most systems |

> **Tip:** The included `setup_environment.sh` script installs everything into a single conda environment. See [Installation](#installation).

---

## Installation

### Option A: Automated setup (recommended)

The `setup_environment.sh` script creates a conda environment with all dependencies:

```bash
# Clone the repository
git clone https://github.com/ravibot93/ENIGMA-DTI-GWAS.git SCRIPTS
cd SCRIPTS

# Run the setup script
bash setup_environment.sh
```

This creates a conda environment called `enigma_dti_gwas` containing R, Python, PLINK2, bcftools, vcftools, tabix, and the `eDTI_outliers` statistical QC tool. The setup takes about 5–10 minutes.

To activate the environment before running the pipeline:

```bash
conda activate enigma_dti_gwas
```

**Custom install location:** If you need the environment on a scratch filesystem (common on HPC), edit `ENV_PREFIX` at the top of `setup_environment.sh`:

```bash
ENV_PREFIX="/scratch/mylab/envs/enigma_dti_gwas"
```

### Option B: Manual installation

If you prefer to manage dependencies yourself, ensure that `plink2`, `R`, `bcftools`, `tabix`, `vcftools`, and `git` are all available on your `PATH`. The pipeline's `run0` script checks for these at startup and will exit with an error message if any are missing.

Install the `eDTI_outliers` Python package manually:

```bash
pip install git+https://github.com/lizhaddad/eDTI_additional_QC.git
```

---

## Directory Setup

Create a working directory and organize your input files:

```bash
mkdir -p $HOME/enigma/DTI_GWAS
cd $HOME/enigma/DTI_GWAS
```

Your directory should look like this before running the pipeline:

```
$HOME/enigma/DTI_GWAS/
├── combinedROItable_FA.csv              # DTI phenotype (from TBSS)
├── combinedROItable_MD.csv              # DTI phenotype (from TBSS)
├── combinedROItable_AD.csv              # DTI phenotype (from TBSS)
├── combinedROItable_RD.csv              # DTI phenotype (from TBSS)
├── Covariates.csv                       # You create this (see below)
├── MYCOHORT_EUR_AB_20260101_QC3.bed     # Genotype files (from QC)
├── MYCOHORT_EUR_AB_20260101_QC3.bim
├── MYCOHORT_EUR_AB_20260101_QC3.fam
├── MYCOHORT_EUR_AB_20260101_PCACovariates.eigenvec  # Ancestry PCA
├── IMPUTED/                             # Imputed VCF files
│   ├── chr1.dose.vcf.gz
│   ├── chr2.dose.vcf.gz
│   ├── ...
│   └── chr22.dose.vcf.gz
└── SCRIPTS/                             # This repository
    ├── ENIGMA_DTI/
    │   ├── run0_ENIGMA-DTI_GWAS_format.sh
    │   ├── run1_ENIGMA-DTI_GWAS_step1.sh
    │   └── run2_ENIGMA-DTI_GWAS_step2.sh
    ├── enigma_backend/
    ├── setup_environment.sh
    └── ...
```

### File naming convention

All output files use the pattern: `COHORT_ANCESTRY_ANALYST_DATE`

- **COHORT** — short cohort abbreviation, no spaces (e.g., `HCP`, `TOP`, `QTIM`)
- **ANCESTRY** — ancestry group matching your genotype QC (e.g., `EUR`, `AFR`, `EAS`)
- **ANALYST** — your initials (e.g., `AB`)
- **DATE** — auto-generated as `YYYYMMDD`

Use the same COHORT, ancestry, and analyst values across all three scripts.

---

## Preparing Your Input Files

### DTI phenotype files (combinedROItable_*.csv)

These are the output of the ENIGMA-DTI TBSS pipeline's ROI extraction step. Each is a CSV with one row per subject and columns for each bilateral-averaged ROI:

```
SubjID, FA_ACR, FA_ALIC, FA_BCC, FA_CC, ...
sub-100206, 0.4521, 0.5632, 0.6123, 0.5234, ...
sub-100307, 0.4401, 0.5521, 0.6001, 0.5102, ...
```

You need at minimum the FA file. Comment out any metric you don't have in the run0 settings.

> **Note on subject IDs:** The pipeline automatically handles prefix mismatches between files. For example, if your DTI files use `sub-100206` but your genotype files use `100206`, the pipeline will match them correctly by stripping common prefixes (`sub-`, `SUB-`, `Sub-`).

### Covariates.csv

Create this file yourself in Excel or a text editor. Required columns:

| Column | Description | Values |
|--------|-------------|--------|
| SubjID | Subject identifier | Must match DTI files (prefix differences are handled automatically) |
| Age | Age in years | Numeric |
| Sex | Biological sex | `M` = male, `F` = female |
| AffectionStatus | Disease status | `A` = affected/case, `C` = control/healthy |

**Cohort type coding:**

| Cohort type | AffectionStatus coding |
|------------|----------------------|
| Population / healthy-only | All subjects = `C` |
| Case-control | Cases = `A`, controls = `C` |
| Case-only | All subjects = `A` |

**Additional columns you may need:**

- **Site dummy variables** — if your cohort spans multiple MRI scanners or sites, add dummy-coded columns (e.g., `Dummy1_MR`, `Dummy2_MR`). With *k* sites, add *k−1* dummy columns.
- **FID / IID** — optional. If your Covariates.csv already has FID and IID columns matching your `.fam` file, the pipeline will use them directly.

**Important:** No missing values are allowed in the covariates. Remove any subjects with incomplete data before running the pipeline.

Example `Covariates.csv`:

```csv
SubjID,Age,Sex,AffectionStatus
100206,29,M,C
100307,26,F,C
100408,31,M,A
100610,28,F,C
```

---

## Step 0 — Phenotype QC & Genotype Prep (run0)

**Script:** `ENIGMA_DTI/run0_ENIGMA-DTI_GWAS_format.sh`

### What it does

1. **Checks software** — verifies plink2, R, vcftools, bcftools, tabix, conda are available
2. **Validates input files** — checks all files exist, detects available DTI metrics, validates column names
3. **Harmonizes subject IDs** — matches subjects across phenotype, covariate, and genotype files
4. **Runs statistical QC** — uses `eDTI_outliers` (Nir, Haddad et al.) to flag and remove outlier subjects based on extreme ROI values
5. **Filters related individuals** — removes related pairs (KING kinship > 0.0885) if your cohort is unrelated
6. **Converts VCF to PGEN** — converts imputed VCF files to PLINK2 PGEN format
7. **Merges DTI phenotypes** — creates a single `Phenotypes.csv` excluding outliers and subjects without genotype data

### Settings to edit

Open `ENIGMA_DTI/run0_ENIGMA-DTI_GWAS_format.sh` and edit the `USER SETTINGS` section at the top:

```bash
# Cohort identifiers (use the same values in all three scripts)
COHORT="MYCOHORT"
ancestry="EUR"
analyst="AB"

# Path to the enigma_backend/ folder
run_directory=$HOME/enigma/DTI_GWAS/SCRIPTS/enigma_backend

# DTI phenotype CSV files (comment out any you don't have)
csvFILE_FA=$HOME/enigma/DTI_GWAS/combinedROItable_FA.csv
csvFILE_MD=$HOME/enigma/DTI_GWAS/combinedROItable_MD.csv
csvFILE_AD=$HOME/enigma/DTI_GWAS/combinedROItable_AD.csv
csvFILE_RD=$HOME/enigma/DTI_GWAS/combinedROItable_RD.csv

# QCed genotype file prefix (no .bed/.bim/.fam extension)
datafile=$HOME/enigma/DTI_GWAS/MYCOHORT_EUR_AB_20260101_QC3

# Path to Covariates.csv
covariatesCSV=$HOME/enigma/DTI_GWAS/Covariates.csv

# Directory with imputed VCF files
imputed_dir=$HOME/enigma/DTI_GWAS/IMPUTED

# Genome build of your imputed data
genome_build="hg19"    # or "hg38"

# Output directory
outFolder=$HOME/enigma/DTI_GWAS/formatted

# Subject ID column name in your DTI CSV files
subjID_column="SubjID"

# Statistical QC settings
run_statistical_qc=1      # 1 = run QC, 0 = skip
```

### Running

```bash
conda activate enigma_dti_gwas
cd $HOME/enigma/DTI_GWAS/SCRIPTS/ENIGMA_DTI
chmod 755 run0_ENIGMA-DTI_GWAS_format.sh
./run0_ENIGMA-DTI_GWAS_format.sh
```

All output is logged to `formatted/logs/${COHORT}_run0.log`.

### What it produces

| Output | Used by |
|--------|---------|
| `formatted/Phenotypes.csv` | run1 |
| `IMPUTED/pgen/chr*.{pgen,pvar,psam}` | run2 |
| `formatted/QC_statistical/` | Send to ENIGMA team |
| `formatted/logs/${COHORT}_run0.log` | Send to ENIGMA team |

---

## Step 1 — Covariate Formatting & QC (run1)

**Script:** `ENIGMA_DTI/run1_ENIGMA-DTI_GWAS_step1.sh`

### What it does

1. **Reformats covariates** — recodes Sex (`M`→`1`, `F`→`2`) and AffectionStatus (`C`→`1`, `A`→`2`) for PLINK2
2. **Generates derived covariates** — Age², Age×Sex, Age²×Sex interaction terms
3. **Merges ancestry PCs** — adds the first 4 (if N ≤ 1000) or 10 (if N > 1000) principal components
4. **Merges global DTI metrics** — adds FA/MD/AD/RD GlobalAverage from Phenotypes.csv
5. **Validates site correction** — checks that site dummies or scanner PCA are included
6. **Generates QC output** — histograms, summary statistics, and run notes

### Settings to edit

Open `ENIGMA_DTI/run1_ENIGMA-DTI_GWAS_step1.sh` and edit the `USER SETTINGS`:

```bash
# Path to backend scripts (same as run0)
run_directory=$HOME/enigma/DTI_GWAS/SCRIPTS/enigma_backend

# Input files
csvFILE=$HOME/enigma/DTI_GWAS/Covariates.csv                          # Your covariates
combinedROItableFILE=$HOME/enigma/DTI_GWAS/formatted/Phenotypes.csv    # From run0
pcaFILE=$HOME/enigma/DTI_GWAS/MYCOHORT_EUR_AB_20260101_PCACovariates.eigenvec
localfamFILE=$HOME/enigma/DTI_GWAS/MYCOHORT_EUR_AB_20260101_QC4.fam

# Column name mappings (must match your Covariates.csv exactly)
ageColumnHeader="Age"
sexColumnHeader="Sex"
maleIndicator="M"               # What value represents males in your data?
CaseControlCohort=0              # 1 = case-control, 0 = healthy-only
affectionStatusColumnHeader="AffectionStatus"
affectedIndicator="A"

# Cohort name and output directory
COHORT="MYCOHORT"
outDir=$HOME/enigma/DTI_GWAS/formatted

# Site/scanner correction (choose one or neither)
site_dummies=""                  # e.g., "Dummy1_MR Dummy2_MR"
scannerPCA_file=""               # Alternative: path to scanner PCA file
scannerPCA_n=6                   # Number of scanner PCs to include
```

### Running

```bash
./run1_ENIGMA-DTI_GWAS_step1.sh
```

### What it produces

| Output | Used by |
|--------|---------|
| `formatted/${COHORT}_enigma_dti_gwas.pheno` | run2 |
| `formatted/${COHORT}_enigma_dti_gwas.covar` | run2 |
| `formatted/QC_ENIGMA/*.pdf` | Send to ENIGMA team |
| `formatted/logs/${COHORT}_run1.log` | Send to ENIGMA team |

---

## QC Checkpoint — Send Results for Review

**Before running the GWAS**, email the following files to **enigma.DTIgenetics@gmail.com**:

- All files from `formatted/logs/`
- All files from `formatted/QC_statistical/`
- All files from `formatted/QC_ENIGMA/` (PDF plots, summary statistics)
- `formatted/QC_statistical/phenotype_merge_QC_log.txt`

Include a brief note describing:
- Your cohort (name, N, ancestry, case-control status)
- Any additional covariates you plan to include (site dummies, scanner PCA)
- Any issues or questions

**Wait for confirmation from the ENIGMA team before proceeding to Step 2.**

---

## Step 2 — Run GWAS (run2)

**Script:** `ENIGMA_DTI/run2_ENIGMA-DTI_GWAS_step2.sh`

### What it does

Runs all GWAS models using PLINK2 `--glm` on imputed dosages. For each of 100 DTI phenotypes (25 ROIs × 4 metrics), up to 6 models are tested, producing up to 600 GWAS runs for case-control cohorts.

### Settings to edit

Open `ENIGMA_DTI/run2_ENIGMA-DTI_GWAS_step2.sh` and edit the `USER SETTINGS`:

```bash
# Directories (same as run0/run1)
run_directory=$HOME/enigma/DTI_GWAS/SCRIPTS/enigma_backend
ENIGMA_DTI_GWAS_dir=$HOME/enigma/DTI_GWAS

# Cohort identifiers (must match run0 and run1 exactly)
COHORT="MYCOHORT"
ancestry="EUR"
analyst="AB"

# Input files (all created by previous steps)
pgen_dir=$HOME/enigma/DTI_GWAS/IMPUTED/pgen
phenoFILE=$HOME/enigma/DTI_GWAS/formatted/${COHORT}_enigma_dti_gwas.pheno
covarFILE=$HOME/enigma/DTI_GWAS/formatted/${COHORT}_enigma_dti_gwas.covar

# Execution mode
mode="qsub"    # "qsub" = submit to cluster, "run" = run locally, "manual" = generate commands only

# Cluster settings (only if mode="qsub")
qsub_mem="8G"
qsub_queue=""

# Chromosomes
chrStart=1; chrEnd=22; doChrX=1
```

### Execution modes

| Mode | Description | Best for |
|------|-------------|----------|
| `qsub` | Submits jobs to SGE/PBS cluster | HPC environments (recommended) |
| `run` | Runs all GWAS locally in sequence | Small cohorts or single-node machines |
| `manual` | Generates command files without running | Review commands first, then execute manually |

### Running

```bash
./run2_ENIGMA-DTI_GWAS_step2.sh
```

---

## Packaging Results

After all GWAS runs complete:

```bash
cd $HOME/enigma/DTI_GWAS
mkdir ENIGMA-DTI-GWAS_output
cp formatted/*.pdf formatted/*.txt formatted/*.csv \
   GWAS_output/combined/*.txt.gz GWAS_output/meta/*.meta.gz \
   *.log ENIGMA-DTI-GWAS_output/

zip -r ${COHORT}_${ancestry}_ENIGMA-DTI-GWAS_${date}_${analyst}.zip \
   ENIGMA-DTI-GWAS_output/
```

Contact **enigma.DTIgenetics@gmail.com** for upload credentials.

---

## GWAS Models

For each of 100 DTI phenotypes (25 ROIs × 4 metrics), the following models are run:

| # | Model | Covariates | Run when |
|---|-------|------------|----------|
| 1 | Full cohort + global average | GlobalAvg + Sex + Age + Age² + Age×Sex + Age²×Sex + PCs + AffectionStatus + site dummies | All cohorts |
| 2 | Controls only + global average | Same as Model 1 minus AffectionStatus | Case-control cohorts |
| 3 | Cases only + global average | Same as Model 1 minus AffectionStatus | Case-control and case-only |
| 4 | Males only + global average | GlobalAvg + Age + Age² + PCs + AffectionStatus + site dummies | All (if N_males ≥ 100) |
| 5 | Females only + global average | Same as Model 4 | All (if N_females ≥ 100) |
| 6 | Full cohort without global average | Same as Model 1 minus GlobalAvg | All cohorts |

This results in up to **600 GWAS runs** for case-control cohorts (500 for population/healthy-only cohorts).

---

## Troubleshooting

### Conda activation crashes the script

**Symptom:** The script exits immediately after trying to activate conda, with no clear error message.

**Cause:** Conda's shell hook references internal variables (`CONDA_SHLVL`, `_CE_CONDA`) that may not exist yet. The pipeline uses `set -u` (exit on unset variables), which kills the script.

**Fix:** This is already handled in the current version of the pipeline. If you're using an older version, wrap the conda activation block:

```bash
set +eu
eval "$(conda shell.bash hook)" 2>/dev/null || \
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate enigma_dti_gwas 2>/dev/null || true
set -eu
```

### Subject IDs don't match between files

**Symptom:** The phenotype merge reports 0 matched subjects, or many subjects are listed as unmatched.

**Cause:** Different tools may use different ID formats. For example, the TBSS pipeline might output `sub-100206` while the genotype files use `100206`.

**Fix:** The pipeline automatically handles common prefix mismatches (`sub-`, `SUB-`, `Sub-`, `sub_`, `SUB_`, `Sub_`). It tries exact matching first, then falls back to normalized matching with prefixes stripped. The QC log reports when normalization was used.

If your IDs differ in a way the pipeline doesn't handle (e.g., completely different ID schemes), you will need to create a mapping file and harmonize the IDs in your input files before running the pipeline.

### PLINK2 GLIBC errors on HPC

**Symptom:** `plink2: /lib64/libc.so.6: version 'GLIBC_2.29' not found`

**Fix:** Use the conda environment (`conda activate enigma_dti_gwas`), which bundles a compatible PLINK2 binary. Alternatively, download a statically-linked PLINK2 binary from https://www.cog-genomics.org/plink/2.0/.

### VCF conversion fails or is slow

**Symptom:** `run0` hangs or crashes during VCF-to-PGEN conversion.

**Tips:**
- Ensure your VCF files are bgzipped (`.vcf.gz`) and tabix-indexed (`.vcf.gz.tbi`)
- Check available disk space — PGEN files can be large
- On HPC, ensure you're running on a compute node, not a login node

### "No space left on device"

**Fix:** The PGEN conversion and GWAS output can require significant disk space. Use a scratch filesystem if available. Check with `df -h`.

---

## Repository Structure

```
├── ENIGMA_DTI/                              # User-facing scripts (EDIT THESE)
│   ├── run0_ENIGMA-DTI_GWAS_format.sh       #   Step 0: Phenotype QC & genotype prep
│   ├── run1_ENIGMA-DTI_GWAS_step1.sh        #   Step 1: Covariate formatting & QC
│   └── run2_ENIGMA-DTI_GWAS_step2.sh        #   Step 2: GWAS execution
│
├── enigma_backend/                          # Backend scripts (DO NOT EDIT)
│   ├── ENIGMA_DTI_pheno_merge.R             #   Merge DTI phenotype CSVs
│   ├── ENIGMA_DTI_create_formatted_files.R  #   Format covariates for PLINK2
│   ├── ENIGMA_DTI_GWAS_QC_plots_stats.R     #   QC plots & summary statistics
│   ├── ENIGMA_unrelatedGWAS_plink_DTI.sh    #   GWAS command generation
│   ├── convert_dosage_vcf_to_pgen.sh        #   VCF → PGEN conversion
│   ├── ENIGMA_DTI_combine_one.sh            #   Combine per-chr results
│   └── ENIGMA_DTI_plot_one.sh               #   Manhattan/QQ plots
│
├── setup_environment.sh                     # Conda environment setup
├── environment.yaml                         # Conda env specification
├── ENIGMA-DTI-GWAS_SOP.md                   # Full protocol document
└── README.md                                # This file
```

---

## Contact & Citation

**ENIGMA-DTI Genetics Support Team:** enigma.DTIgenetics@gmail.com

- [ENIGMA website](http://enigma.ini.usc.edu/)
- [ENIGMA mailing list](http://enigma.ini.usc.edu/members/join/)

### Citation

If you use this pipeline, please cite:

> ENIGMA Genetics Support Team. ENIGMA-DTI-GWAS Association Analysis Protocol [Online]. The Enhancing Neuroimaging Genetics through Meta-Analysis (ENIGMA) Consortium. http://enigma.ini.usc.edu/genetics-protocols/
