# ENIGMA-DTI GWAS Association Analysis Protocol (PLINK2)

**Authors:** ENIGMA-DTI Genetics Working Group  
**Date:** February 2026  
**Adapted from:** ENIGMA-3 Cortical GWAS Protocol (ENIGMA3_SOP_20151224)

---

## Overview

This document provides a GWAS analysis plan for diffusion tensor imaging (DTI) phenotypes using PLINK2. It is designed for **unrelated samples** that have completed the ENIGMA-DTI TBSS pipeline and pre-imputation QC.

For groups with less experience or those who imputed using the ENIGMA protocol, follow the protocols on GitHub:  
https://github.com/ENIGMA-git/ENIGMA-DTI-GWAS

For questions: enigma.DTIgenetics@gmail.com

**Citation:** ENIGMA Genetics Support Team. ENIGMA-DTI-GWAS Association Analysis Protocol [Online]. The Enhancing Neuroimaging Genetics through Meta-Analysis (ENIGMA) Consortium. http://enigma.ini.usc.edu/genetics-protocols/

---

## Prerequisites

You should have:
1. Completed the ENIGMA-DTI TBSS pipeline (producing `combinedROItable_*.csv` files)
2. Completed the pre-imputation QC protocol (producing `${datafile}_QC3` PLINK files)
3. Imputed your genetic data (producing `chr*.dose.vcf.gz` files)

### Required software

- **PLINK2** (https://www.cog-genomics.org/plink2) — **Not PLINK1.9**
- **R** (http://cran.stat.ucla.edu/)
- **vcftools** (https://vcftools.github.io/)
- **bcftools** and **tabix** (https://www.htslib.org/)
- **git**, **wget** (standard Linux tools)
- **Conda/Miniconda** (https://docs.conda.io/projects/conda/en/latest/user-guide/install/) — for eDTI_outliers statistical QC

The `run0` script checks for these and attempts to install the eDTI_outliers package automatically.

---

## User workflow: 3 scripts

We will run **three scripts** in order to set up and perform the GWAS:

| Script | Purpose |
|--------|---------|
| `run0_ENIGMA-DTI_GWAS_format.sh` | Format DTI phenotypes and covariates |
| `run1_ENIGMA-DTI_GWAS_step1.sh` | Generate summary statistics, create analysis-ready files |
| `run2_ENIGMA-DTI_GWAS_step2.sh` | Run all GWAS models |

For each step, edit the user-configurable parameters at the top of the script, then run it. The backend scripts in `enigma_backend/` should **not** be modified.

---

## Step 1 — Organise your input files

Create a working directory and gather the required input files:

```bash
mkdir -p $HOME/enigma/DTI_GWAS
cd $HOME/enigma/DTI_GWAS
```

You will need the following files in your working directory:

1. **combinedROItable_FA.csv** — Mean FA within each JHU Atlas ROI from TBSS
2. **combinedROItable_MD.csv** — Mean MD within each JHU Atlas ROI from TBSS
3. **combinedROItable_AD.csv** — Mean AD within each JHU Atlas ROI from TBSS
4. **combinedROItable_RD.csv** — Mean RD within each JHU Atlas ROI from TBSS
5. **Covariates.csv** — You must create this file (see below)
6. **`${datafile}_QC3` PLINK files** — From pre-imputation QC (`.bed/.bim/.fam`)
7. **`${datafile}_PCACovariates.eigenvec`** — Ancestry PCA from pre-imputation QC
8. **Imputed VCF files** — `chr*.dose.vcf.gz` from imputation

### Creating Covariates.csv

Using Excel or a text editor, create a CSV file with these columns (at minimum):

| Column | Description | Coding |
|--------|-------------|--------|
| FID | Family ID | Must match genotype files |
| IID | Individual ID | Must match genotype files |
| Sex | Biological sex | M = male, F = female |
| Age | Age in years | Numeric |
| AffectionStatus | Disease status | A = affected/case, C = control/healthy |

**Notes:**
- **Case-control cohorts:** Code cases as `A`, controls as `C`
- **Case-only cohorts:** All individuals coded as `A`
- **Population/healthy-only cohorts:** All individuals coded as `C`
- **No missing values allowed.** Remove subjects with any missing covariate data.
- Add additional dummy covariates (e.g., `Dummy1_MR`, `Dummy2_MR`) if needed for scanner/site correction.

### File naming convention

```
COHORT_ANCESTRY_ANALYST_DATE
```
- `COHORT` = cohort abbreviation (use the abbreviation from Grasby et al. 2020 Supp. Table 2 if applicable)
- `ANCESTRY` = e.g., EUR
- `ANALYST` = initials
- `DATE` = YYYYMMDD

---

## Step 2 — Download the scripts

```bash
cd $HOME/enigma/DTI_GWAS
git clone https://github.com/ENIGMA-git/ENIGMA-DTI-GWAS.git SCRIPTS
chmod -R 755 SCRIPTS/
```

Your working directory should look like:

```
$HOME/enigma/DTI_GWAS/
├── combinedROItable_FA.csv
├── combinedROItable_MD.csv
├── combinedROItable_AD.csv
├── combinedROItable_RD.csv
├── Covariates.csv
├── ${datafile}_QC3.bed/.bim/.fam
├── ${datafile}_PCACovariates.eigenvec
└── SCRIPTS/
    ├── run0_ENIGMA-DTI_GWAS_format.sh
    ├── run1_ENIGMA-DTI_GWAS_step1.sh
    ├── run2_ENIGMA-DTI_GWAS_step2.sh
    └── enigma_backend/
        ├── ENIGMA_DTI_pheno_merge.R
        ├── ENIGMA_DTI_create_formatted_files.R
        ├── ENIGMA_DTI_GWAS_QC_plots_stats.R
        ├── ENIGMA_unrelatedGWAS_plink_DTI.sh
        └── convert_dosage_vcf_to_pgen.sh
```

---

## Step 3 — Edit and run `run0_ENIGMA-DTI_GWAS_format.sh`

This script performs phenotype statistical QC, merges the four DTI phenotype files, reformats them, converts imputed VCF to PGEN format, and removes related individuals. All output is logged to `${COHORT}_run0.log`.

Open `SCRIPTS/run0_ENIGMA-DTI_GWAS_format.sh` and set the following parameters:

```bash
#-- Cohort information
COHORT="MYCOHORT"
ancestry="EUR"
analyst="AB"

#-- Backend scripts directory (do not change if using default layout)
run_directory=$HOME/enigma/DTI_GWAS/SCRIPTS/enigma_backend

#-- Full path to R binary (find with: which R)
Rbin=/usr/local/bin/R

#-- Full paths to the four DTI phenotype CSV files
csvFILE_FA=$HOME/enigma/DTI_GWAS/combinedROItable_FA.csv
csvFILE_MD=$HOME/enigma/DTI_GWAS/combinedROItable_MD.csv
csvFILE_AD=$HOME/enigma/DTI_GWAS/combinedROItable_AD.csv
csvFILE_RD=$HOME/enigma/DTI_GWAS/combinedROItable_RD.csv

#-- Full path to your QCed genotype file prefix (no extension)
datafile=$HOME/enigma/DTI_GWAS/COHORT_EUR_AB_20260101_QC3

#-- Full path to directory with imputed VCF files
imputed_dir=$HOME/enigma/DTI_GWAS/IMPUTED

#-- Output directory for formatted phenotype/covariate files
outFolder=$HOME/enigma/DTI_GWAS/formatted

#-- Subject ID column name in your combinedROItable CSV files
subjID_column="SubjID"

#-- STATISTICAL QC SETTINGS (eDTI_outliers) --
run_statistical_qc=1      # 1=run QC, 0=skip
demogCSV=""                # Path to demographics CSV (optional)
dxColumn=""                # Diagnosis column name in demogCSV (optional)
siteColumn=""              # Site column name in demogCSV (optional)
```

Then run:
```bash
cd $HOME/enigma/DTI_GWAS/SCRIPTS
chmod 755 run0_ENIGMA-DTI_GWAS_format.sh
./run0_ENIGMA-DTI_GWAS_format.sh
```

**What this script does (6 steps, all logged):**

1. **Software check** — Verifies plink2, R, vcftools, bcftools, tabix, git, conda are available
2. **Input file check** — Validates all input files exist, logs line counts
3. **Statistical QC (eDTI_outliers)** — Installs the `eDTI_additional_QC` package (Nir, Haddad et al.) into a conda environment, then runs the `eDTI_outliers` statistical QC tool on all four DTI CSV files. Subjects are flagged as outliers using two thresholds:
   - *Threshold 1:* ≥5% of ROI values outside 0.1st/99.9th percentile
   - *Threshold 2:* ≥10% of ROI values outside 0.5th/99.5th percentile
   
   If a demographics CSV is provided with diagnosis and/or site columns, outlier detection is also performed within diagnostic groups and within sites. **Flagged subjects are automatically removed** from all downstream phenotype files. The full list of removed subjects is logged.
4. **Relatedness filtering** — Removes related individuals (KING kinship cutoff 0.0885 = first cousins and closer) → `${datafile}_QC4`
5. **VCF to PGEN conversion** — Converts imputed VCF to PGEN format (required by PLINK2)
6. **Phenotype merge** — Merges the four combinedROItable CSV files into a single `Phenotypes.csv`, excluding statistical QC outliers and subjects without genotype data. A complete subject flow log is produced.

**Subject flow summary** (printed at end of run0):
```
Genotype subjects (pre-QC3):                  N
Removed by relatedness (KING 0.0885):       - X
Genotype subjects (post-relatedness, QC4):    N-X
Statistical QC outliers removed:             - Y
Final subjects in Phenotypes.csv:              N-X-Y
```

---

## Step 4 — Edit and run `run1_ENIGMA-DTI_GWAS_step1.sh`

This script reformats covariates, merges in PCA components and global DTI metrics, generates derived age covariates, and produces QC plots and summary statistics.

Open `SCRIPTS/run1_ENIGMA-DTI_GWAS_step1.sh` and set the following parameters:

```bash
#-- Backend scripts directory
run_directory=$HOME/enigma/DTI_GWAS/SCRIPTS/enigma_backend

#-- Full path to R binary
Rbin=/usr/local/bin/R

#-- Full path to Covariates.csv
csvFILE=$HOME/enigma/DTI_GWAS/Covariates.csv

#-- Full path to merged phenotypes (output from run0)
combinedROItableFILE=$HOME/enigma/DTI_GWAS/formatted/Phenotypes.csv

#-- Full path to ancestry PCA eigenvec file
pcaFILE=$HOME/enigma/DTI_GWAS/COHORT_EUR_AB_20260101_PCACovariates.eigenvec

#-- Full path to QCed .fam file (from run0)
localfamFILE=$HOME/enigma/DTI_GWAS/COHORT_EUR_AB_20260101_QC4.fam

#-- Covariate coding
ageColumnHeader="Age"          # Column header for age
sexColumnHeader="Sex"          # Column header for sex
maleIndicator="M"              # Code for males (M? 1? check carefully!)
CaseControlCohort=0            # 0=population/healthy-only, 1=case-control or case-only
affectionStatusColumnHeader="AffectionStatus"
affectedIndicator="A"          # Code for affected/cases

#-- Cohort name (used in output file naming)
COHORT="MYCOHORT"

#-- Output directory
outDir=$HOME/enigma/DTI_GWAS/formatted
```

Then run:
```bash
./run1_ENIGMA-DTI_GWAS_step1.sh
```

**What this script does:**
1. Recodes Sex (M→1, F→2) and AffectionStatus (C→1, A→2) for PLINK2
2. Generates derived covariates: AgeCsq, AgexSex, AgeCsqxSex
3. Merges in the first 4 (N≤1000) or 10 (N>1000) ancestry PCs
4. Merges in global DTI metrics (FA/MD/AD/RD GlobalAverage)
5. Produces QC histograms and summary statistics
6. Outputs: `${COHORT}_enigma_dti_gwas.pheno` and `${COHORT}_enigma_dti_gwas.covar`

### ⚠️ Checkpoint: Send QC files for review

**Before proceeding to Step 5**, send the following files to enigma.DTIgenetics@gmail.com:

- All `*.pdf`, `*.txt`, and `*.csv` QC log files from the `formatted/` directory
- Relatedness and ancestry stats files
- A summary of which additional covariates you plan to include

**Wait for confirmation before running the GWAS.**

---

## Step 5 — Edit and run `run2_ENIGMA-DTI_GWAS_step2.sh`

This script runs all GWAS models using PLINK2 linear regression on imputed dosages.

Open `SCRIPTS/run2_ENIGMA-DTI_GWAS_step2.sh` and set the following parameters:

```bash
#-- Backend scripts directory
run_directory=$HOME/enigma/DTI_GWAS/SCRIPTS/enigma_backend

#-- Working directory
ENIGMA_DTI_GWAS_dir=$HOME/enigma/DTI_GWAS

#-- Cohort information
COHORT="MYCOHORT"
ancestry="EUR"
analyst="AB"
date="20260101"

#-- Path to imputed PGEN files (output from run0)
pgen_dir=$HOME/enigma/DTI_GWAS/IMPUTED/pgen

#-- Path to phenotype and covariate files (output from run1)
phenoFILE=$HOME/enigma/DTI_GWAS/formatted/${COHORT}_enigma_dti_gwas.pheno
covarFILE=$HOME/enigma/DTI_GWAS/formatted/${COHORT}_enigma_dti_gwas.covar

#-- Unrelated sample (leave as 0)
related=0

#-- Execution mode
mode="manual"   # "manual" = outputs command list; "run" = submits jobs to cluster
Nnodes=1        # Number of parallel nodes (set 1 for local/manual)
```

Then run:
```bash
./run2_ENIGMA-DTI_GWAS_step2.sh
```

### GWAS models

For each DTI metric (FA, MD, AD, RD) × 25 ROIs = **100 phenotypes**, the following models are run:

**Model 1 — Full cohort with global average correction** (all cohort types):
```
Phenotype ~ SNP + GlobalAverage + Sex + Age + AgeCsq + AgexSex + AgeCsqxSex
            + PC1..PC4/10 + AffectionStatus + dummy covariates
```
Output: `${COHORT}_${ancestry}_ENIGMA-DTI_fullcohort_w${metric}_${pheno}_${CHR}_${DATE}_${analyst}.out.gz`

**Model 2 — Controls-only** (case-control cohorts only):
```
Phenotype ~ SNP + GlobalAverage + Sex + Age + AgeCsq + AgexSex + AgeCsqxSex
            + PC1..PC4/10 + dummy covariates
```
Output: `..._Conly_w${metric}_...`

**Model 3 — Cases-only** (case-control and case-only cohorts):
```
Phenotype ~ SNP + GlobalAverage + Sex + Age + AgeCsq + AgexSex + AgeCsqxSex
            + PC1..PC4/10 + dummy covariates
```
Output: `..._Aonly_w${metric}_...`

**Model 4 — Males-only** (all cohort types, if N_males ≥ 100):
```
Phenotype ~ SNP + GlobalAverage + Age + AgeCsq + PC1..PC4/10
            + AffectionStatus + dummy covariates
```
Output: `..._males_w${metric}_...`

**Model 5 — Females-only** (all cohort types, if N_females ≥ 100):
```
Same as Model 4 (no Sex covariates)
```
Output: `..._females_w${metric}_...`

**Model 6 — Full cohort without global average** (all cohort types):
```
Phenotype ~ SNP + Sex + Age + AgeCsq + AgexSex + AgeCsqxSex
            + PC1..PC4/10 + AffectionStatus + dummy covariates
```
Output: `..._fullcohort_wo_...`

This results in up to **600 GWAS runs** for case-control cohorts (500 for population cohorts).

**Important:** Use imputed dosages, not hard calls. Do **not** rescale phenotypes.

---

## Step 6 — Organise summary statistics and upload

### Summary statistics format

For each GWAS, retain these columns:

| Column | Description |
|--------|-------------|
| CHR | Chromosome number |
| BP | Base pair position |
| Effect_allele | Coded/modelled allele |
| NonEffect_allele | The other allele |
| Beta | Beta estimate (≥5 decimal places; 'NA' if unavailable) |
| SE | Standard error (≥5 decimal places; 'NA' if unavailable) |
| Pval | P-value ('NA' if unavailable) |
| N | Sample size (only if varies across SNPs) |

Do **not** include allele frequency, Rsq/INFO, RS number, or Imputed/Genotyped — these go in the meta file.

### Meta data file

Upload one meta file per chromosome: `${COHORT}_${ancestry}_ENIGMA-DTI_${CHR}_${DATE}_${analyst}.meta.gz`

| Column | Description |
|--------|-------------|
| CHR | Chromosome |
| BP | Base pair position |
| RS_number | rsID (if available) |
| Effect_allele | Coded allele |
| NonEffect_allele | Other allele |
| Freq_EA | Effect allele frequency |
| MAF | Minor allele frequency |
| Rsq | Imputation accuracy (Rsq/INFO) |
| Imputed | 1=imputed, 0=genotyped |

### Phenotypic summary statistics

Rename and include the `RUN_NOTES.txt` file:
```bash
mv formatted/RUN_NOTES.txt ${COHORT}_RUN_NOTES_${DATE}_${analyst}.txt
```

### Package and upload

```bash
cd $HOME/enigma/DTI_GWAS
mkdir ENIGMA-DTI-GWAS_output
cp formatted/*.pdf formatted/*.txt formatted/*.csv *.out.gz *.meta.gz *.log \
   ENIGMA-DTI-GWAS_output/
zip -r ${COHORT}_${ancestry}_ENIGMA-DTI-GWAS_${DATE}_${analyst}.zip \
   ENIGMA-DTI-GWAS_output/
```

Upload to the ENIGMA upload server. Contact enigma.DTIgenetics@gmail.com for upload credentials.

---

## Appendix: List of 100 phenotypes

25 ROIs × 4 metrics (FA, MD, AD, RD) = 100 phenotypes (bilateral averages):

| FA | MD | AD | RD |
|----|----|----|-----|
| FA_GlobalAverage | MD_GlobalAverage | AD_GlobalAverage | RD_GlobalAverage |
| FA_ACR | MD_ACR | AD_ACR | RD_ACR |
| FA_ALIC | MD_ALIC | AD_ALIC | RD_ALIC |
| FA_BCC | MD_BCC | AD_BCC | RD_BCC |
| FA_CC | MD_CC | AD_CC | RD_CC |
| FA_CGC | MD_CGC | AD_CGC | RD_CGC |
| FA_CGH | MD_CGH | AD_CGH | RD_CGH |
| FA_CR | MD_CR | AD_CR | RD_CR |
| FA_CST | MD_CST | AD_CST | RD_CST |
| FA_EC | MD_EC | AD_EC | RD_EC |
| FA_FX | MD_FX | AD_FX | RD_FX |
| FA_FXST | MD_FXST | AD_FXST | RD_FXST |
| FA_GCC | MD_GCC | AD_GCC | RD_GCC |
| FA_IC | MD_IC | AD_IC | RD_IC |
| FA_IFO | MD_IFO | AD_IFO | RD_IFO |
| FA_PCR | MD_PCR | AD_PCR | RD_PCR |
| FA_PLIC | MD_PLIC | AD_PLIC | RD_PLIC |
| FA_PTR | MD_PTR | AD_PTR | RD_PTR |
| FA_RLIC | MD_RLIC | AD_RLIC | RD_RLIC |
| FA_SCC | MD_SCC | AD_SCC | RD_SCC |
| FA_SCR | MD_SCR | AD_SCR | RD_SCR |
| FA_SFO | MD_SFO | AD_SFO | RD_SFO |
| FA_SLF | MD_SLF | AD_SLF | RD_SLF |
| FA_SS | MD_SS | AD_SS | RD_SS |
| FA_UNC | MD_UNC | AD_UNC | RD_UNC |
