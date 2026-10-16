###############################################################################
# ENIGMA_DTI_create_formatted_files.R
#
# Backend script: Reformat covariates for PLINK2 GWAS analysis.
# - Recodes Sex and AffectionStatus to numeric (1/2)
# - Generates derived age covariates (AgeCsq, AgexSex, AgeCsqxSex)
# - Merges in ancestry PCs (4 for N≤1000, 10 for N>1000)
# - Merges in global DTI metrics (FA/MD/AD/RD GlobalAverage)
# - Outputs .pheno and .covar files
#
# Called by run1_ENIGMA-DTI_GWAS_step1.sh — do not run directly.
#
# Arguments (positional via --args):
#   1  = csvFILE (Covariates.csv)
#   2  = combinedROItableFILE (Phenotypes.csv from run0)
#   3  = pcaFILE (eigenvec file)
#   4  = ageColumnHeader
#   5  = sexColumnHeader
#   6  = maleIndicator
#   7  = CaseControlCohort (0 or 1)
#   8  = affectionStatusColumnHeader
#   9  = affectedIndicator
#   10 = outDir
#   11 = eName (project name for file naming)
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 11) {
  stop("Insufficient arguments. This script is called by run1_ENIGMA-DTI_GWAS_step1.sh")
}

csvFILE                      <- args[1]
combinedROItableFILE         <- args[2]
pcaFILE                      <- args[3]
ageColumnHeader              <- args[4]
sexColumnHeader              <- args[5]
maleIndicator                <- args[6]
CaseControlCohort            <- as.numeric(args[7])
affectionStatusColumnHeader  <- args[8]
affectedIndicator            <- args[9]
outDir                       <- args[10]
eName                        <- args[11]

# Optional: scanner PCA file
scannerPCA_file <- if (length(args) >= 12 && args[12] != "") args[12] else ""
scannerPCA_n    <- if (length(args) >= 13) as.numeric(args[13]) else 0
scannerPCA_idcol <- if (length(args) >= 14 && args[14] != "") args[14] else "eid"

#=============================================================================
# ID normalization helper
#   Strips common prefixes (sub-, SUB-, Sub-, sub_, SUB_, Sub_) so that
#   IDs from different sources (DTI CSVs, Covariates, .fam) can be matched
#   even when one source uses "sub-100206" and another uses "100206".
#=============================================================================
normalize_id <- function(ids) {
  gsub("^[Ss][Uu][Bb][-_]", "", trimws(as.character(ids)))
}

cat("=== ENIGMA-DTI Covariate Formatting ===\n")

#=============================================================================
# 1. Read input files
#=============================================================================
cat("Reading input files...\n")

covar <- read.csv(csvFILE, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Covariates: %d subjects, %d columns\n", nrow(covar), ncol(covar)))

pheno <- read.csv(combinedROItableFILE, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Phenotypes: %d subjects, %d columns\n", nrow(pheno), ncol(pheno)))

# Read PCA file (tab or space delimited)
# Eigenvec files often have NO header — detect and handle both cases
pca <- tryCatch(
  read.table(pcaFILE, header = FALSE, stringsAsFactors = FALSE),
  error = function(e) {
    read.csv(pcaFILE, header = FALSE, stringsAsFactors = FALSE)
  }
)

# Check if first row looks like a header (contains non-numeric values like "FID", "IID", "PC1")
first_row <- as.character(pca[1, ])
looks_like_header <- any(grepl("^(#?FID|IID|PC|C[0-9])", first_row, ignore.case = TRUE))

if (looks_like_header) {
  # Re-read with header
  pca <- tryCatch(
    read.table(pcaFILE, header = TRUE, stringsAsFactors = FALSE, comment.char = ""),
    error = function(e) {
      read.csv(pcaFILE, header = TRUE, stringsAsFactors = FALSE)
    }
  )
  cat(sprintf("  PCA: %d subjects, %d columns (with header)\n", nrow(pca), ncol(pca)))
} else {
  # No header — assign column names: FID, IID, PC1, PC2, ..., PCN
  n_pcs_in_file <- ncol(pca) - 2
  colnames(pca) <- c("FID", "IID", paste0("PC", 1:n_pcs_in_file))
  cat(sprintf("  PCA: %d subjects, %d PCs (no header, assigned names)\n", nrow(pca), n_pcs_in_file))
}

# Ensure FID and IID are character
pca$FID <- as.character(pca$FID)
pca$IID <- as.character(pca$IID)

#=============================================================================
# 2. Create FID/IID if not present, validate required columns
#=============================================================================
cat("Validating columns...\n")

# If Covariates.csv doesn't have FID/IID, create them by matching
# against Phenotypes.csv (which has FID/IID from the .fam merge in run0)
if (!"FID" %in% names(covar) || !"IID" %in% names(covar)) {
  cat("  FID/IID not found in Covariates.csv — creating from Phenotypes.csv match.\n")

  # Find the subject ID column in covariates
  id_candidates <- c("subjectID", "SubjID", "subjID", "ID", "Subject", "SID")
  covar_id_col <- NULL
  for (cand in id_candidates) {
    if (cand %in% names(covar)) { covar_id_col <- cand; break }
  }
  if (is.null(covar_id_col)) {
    covar_id_col <- names(covar)[1]
    cat(sprintf("  [WARN] No standard ID column found, using first column: '%s'\n", covar_id_col))
  }
  cat(sprintf("  Using '%s' as subject ID column in Covariates.csv\n", covar_id_col))

  # Clean IDs
  covar[[covar_id_col]] <- trimws(as.character(covar[[covar_id_col]]))

  # Phenotypes.csv already has FID and IID (created by run0 pheno_merge)
  if ("FID" %in% names(pheno) && "IID" %in% names(pheno)) {
    # Try matching covar IDs to pheno IID, then FID_IID, then FID
    pheno$FID <- as.character(pheno$FID)
    pheno$IID <- as.character(pheno$IID)
    pheno$FID_IID <- paste(pheno$FID, pheno$IID, sep = "_")

    covar$FID <- NA_character_
    covar$IID <- NA_character_

    # Strategy 1: covar ID matches pheno IID
    m1 <- match(covar[[covar_id_col]], pheno$IID)
    matched1 <- !is.na(m1)

    # Strategy 2: covar ID matches pheno FID_IID
    m2 <- match(covar[[covar_id_col]], pheno$FID_IID)
    matched2 <- !is.na(m2)

    # Strategy 3: normalize covar IDs (strip sub-/SUB-/Sub- prefix) and match to IID
    covar_ids_norm <- normalize_id(covar[[covar_id_col]])
    m3 <- match(covar_ids_norm, pheno$IID)
    matched3 <- !is.na(m3)

    # Strategy 4: normalize pheno IID and match raw covar ID
    pheno_ids_norm <- normalize_id(pheno$IID)
    m4 <- match(covar[[covar_id_col]], pheno_ids_norm)
    matched4 <- !is.na(m4)

    # Strategy 5: normalize both sides
    m5 <- match(covar_ids_norm, pheno_ids_norm)
    matched5 <- !is.na(m5)

    # Pick the strategy with the most matches
    counts <- c(sum(matched1), sum(matched2), sum(matched3), sum(matched4), sum(matched5))
    names(counts) <- c("IID", "FID_IID", "strip_covar->IID", "covar->strip_pheno", "strip_both")
    cat(sprintf("  Match attempts: %s\n", paste(paste(names(counts), counts, sep="="), collapse=", ")))

    best <- which.max(counts)
    m_best <- list(m1, m2, m3, m4, m5)[[best]]
    cat(sprintf("  Best match: %s (%d matches)\n", names(counts)[best], counts[best]))

    for (i in which(!is.na(m_best))) {
      covar$FID[i] <- pheno$FID[m_best[i]]
      covar$IID[i] <- pheno$IID[m_best[i]]
    }

    # Remove unmatched subjects
    n_unmatched <- sum(is.na(covar$FID))
    if (n_unmatched > 0) {
      cat(sprintf("  Removing %d covariate subjects with no phenotype/genotype match\n", n_unmatched))
      covar <- covar[!is.na(covar$FID), ]
    }
    cat(sprintf("  Covariates after matching: %d subjects\n", nrow(covar)))
  } else {
    # Phenotypes.csv also doesn't have FID/IID — use subjectID as both
    cat("  [WARN] Phenotypes.csv also lacks FID/IID. Using subject ID as both.\n")
    covar$FID <- covar[[covar_id_col]]
    covar$IID <- covar[[covar_id_col]]
  }
}

required_cols <- c("FID", "IID", ageColumnHeader, sexColumnHeader)
if (CaseControlCohort == 1 || affectionStatusColumnHeader %in% names(covar)) {
  required_cols <- c(required_cols, affectionStatusColumnHeader)
}

missing_cols <- setdiff(required_cols, names(covar))
if (length(missing_cols) > 0) {
  # Try case-insensitive match as fallback
  for (mc in missing_cols) {
    ci_match <- names(covar)[tolower(names(covar)) == tolower(mc)]
    if (length(ci_match) > 0) {
      cat(sprintf("  [NOTE] Using '%s' for '%s' (case mismatch)\n", ci_match[1], mc))
      names(covar)[names(covar) == ci_match[1]] <- mc
      missing_cols <- setdiff(missing_cols, mc)
    }
  }
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in Covariates.csv: %s\n  Available columns: %s",
                 paste(missing_cols, collapse = ", "),
                 paste(names(covar), collapse = ", ")))
  }
}

# Check for missing values
na_count <- sum(!complete.cases(covar))
if (na_count > 0) {
  cat(sprintf("  WARNING: Removing %d subjects with missing covariate values\n", na_count))
  covar <- covar[complete.cases(covar), ]
}

#=============================================================================
# 3. Recode Sex: male→1, female→2
#=============================================================================
cat("Recoding Sex...\n")

sex_orig <- covar[[sexColumnHeader]]
sex_new <- rep(NA, length(sex_orig))

if (is.numeric(sex_orig)) {
  sex_new[sex_orig == as.numeric(maleIndicator)] <- 1
  sex_new[sex_orig != as.numeric(maleIndicator) & !is.na(sex_orig)] <- 2
} else {
  sex_new[sex_orig == maleIndicator] <- 1
  sex_new[sex_orig != maleIndicator & sex_orig != "" & !is.na(sex_orig)] <- 2
}

n_males <- sum(sex_new == 1, na.rm = TRUE)
n_females <- sum(sex_new == 2, na.rm = TRUE)
cat(sprintf("  Males: %d, Females: %d\n", n_males, n_females))

if (n_males == 0 || n_females == 0) {
  warning("One sex has 0 subjects! Check your maleIndicator setting.")
}

covar$Sex <- sex_new

#=============================================================================
# 4. Recode AffectionStatus: control→1, affected→2
#=============================================================================
cat("Recoding AffectionStatus...\n")

if (affectionStatusColumnHeader %in% names(covar)) {
  aff_orig <- covar[[affectionStatusColumnHeader]]
  aff_new <- rep(NA, length(aff_orig))

  if (is.numeric(aff_orig)) {
    aff_new[aff_orig == as.numeric(affectedIndicator)] <- 2
    aff_new[aff_orig != as.numeric(affectedIndicator) & !is.na(aff_orig)] <- 1
  } else {
    aff_new[aff_orig == affectedIndicator] <- 2
    aff_new[aff_orig != affectedIndicator & aff_orig != "" & !is.na(aff_orig)] <- 1
  }

  n_cases <- sum(aff_new == 2, na.rm = TRUE)
  n_controls <- sum(aff_new == 1, na.rm = TRUE)
  cat(sprintf("  Controls: %d, Cases: %d\n", n_controls, n_cases))

  covar$AffectionStatus <- aff_new
} else {
  # No affection status — code all as controls
  covar$AffectionStatus <- 1
  cat("  No AffectionStatus column — all coded as controls (1)\n")
}

#=============================================================================
# 5. Generate derived age covariates
#=============================================================================
cat("Generating derived age covariates...\n")

age <- as.numeric(covar[[ageColumnHeader]])
meanAge <- mean(age, na.rm = TRUE)
meanSex <- mean(covar$Sex, na.rm = TRUE)

# Center Age and Sex before computing interaction terms
# This prevents multicollinearity (VIF inflation) between
# Sex, AgexSex, and Age
age_centered <- age - meanAge
sex_centered <- covar$Sex - meanSex

covar$Age <- age
covar$AgeCsq <- age_centered^2
covar$AgexSex <- age_centered * sex_centered
covar$AgeCsqxSex <- (age_centered^2) * sex_centered

cat(sprintf("  Mean age: %.2f\n", meanAge))
cat(sprintf("  Mean sex (1=M,2=F): %.2f\n", meanSex))
cat(sprintf("  Min age: %.1f, Max age: %.1f\n", min(age), max(age)))

#=============================================================================
# 6. Merge ancestry PCs
#=============================================================================
cat("Merging ancestry PCs...\n")

# Standardise PCA column names
pca_names <- names(pca)
# Look for FID/IID columns
fid_col <- grep("^(FID|#FID)$", pca_names, value = TRUE)[1]
iid_col <- grep("^IID$", pca_names, value = TRUE)[1]

if (is.na(fid_col)) fid_col <- pca_names[1]
if (is.na(iid_col)) iid_col <- pca_names[2]

names(pca)[names(pca) == fid_col] <- "FID"
names(pca)[names(pca) == iid_col] <- "IID"

# Select PC columns
pc_cols <- grep("^(PC|C)[0-9]+$", names(pca), value = TRUE)
n_pcs_available <- length(pc_cols)

# Use 4 PCs for N≤1000, 10 PCs for N>1000
N <- nrow(covar)
n_pcs <- ifelse(N > 1000, min(10, n_pcs_available), min(4, n_pcs_available))
pc_cols_use <- pc_cols[1:n_pcs]

cat(sprintf("  N = %d → using %d PCs\n", N, n_pcs))

# Standardise PC names to PC1, PC2, ...
pca_subset <- pca[, c("FID", "IID", pc_cols_use)]
names(pca_subset)[3:ncol(pca_subset)] <- paste0("PC", 1:n_pcs)

covar <- merge(covar, pca_subset, by = c("FID", "IID"), all.x = TRUE)

n_no_pca <- sum(is.na(covar$PC1))
if (n_no_pca > 0) {
  cat(sprintf("  WARNING: %d subjects missing PCA data — removing\n", n_no_pca))
  covar <- covar[!is.na(covar$PC1), ]
}

#=============================================================================
# 6b. Merge scanner PCA components (optional)
#=============================================================================
if (scannerPCA_file != "" && file.exists(scannerPCA_file)) {
  cat("\nMerging scanner PCA components...\n")
  cat(sprintf("  File: %s\n", scannerPCA_file))
  cat(sprintf("  Components to use: %d\n", scannerPCA_n))
  cat(sprintf("  ID column: %s\n", scannerPCA_idcol))

  scanner_pca <- read.csv(scannerPCA_file, header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("  Scanner PCA file: %d subjects, %d columns\n", nrow(scanner_pca), ncol(scanner_pca)))

  # Clean ID column
  scanner_pca[[scannerPCA_idcol]] <- trimws(gsub('"', '', as.character(scanner_pca[[scannerPCA_idcol]])))

  # Select columns: ID + first N scanner PCA columns
  scanner_cols <- grep("scanner_pca_[0-9]+$", names(scanner_pca), value = TRUE)
  scanner_cols <- sort(scanner_cols)
  if (length(scanner_cols) > scannerPCA_n) {
    scanner_cols <- scanner_cols[1:scannerPCA_n]
  }
  cat(sprintf("  Using columns: %s\n", paste(scanner_cols, collapse = ", ")))

  # Rename to scanPC1, scanPC2, ... (shorter names for PLINK2)
  scanner_subset <- scanner_pca[, c(scannerPCA_idcol, scanner_cols), drop = FALSE]
  new_names <- paste0("scanPC", 1:length(scanner_cols))
  names(scanner_subset)[2:ncol(scanner_subset)] <- new_names
  cat(sprintf("  Renamed to: %s\n", paste(new_names, collapse = ", ")))

  # Match scanner PCA IDs to covar IDs
  # Try: IID, FID (exact and normalized)
  scanner_ids <- scanner_subset[[scannerPCA_idcol]]
  scanner_ids_norm <- normalize_id(scanner_ids)

  n_by_iid <- sum(covar$IID %in% scanner_ids)
  n_by_fid <- sum(covar$FID %in% scanner_ids)
  # Also try normalized matching
  n_by_iid_norm <- sum(normalize_id(covar$IID) %in% scanner_ids_norm)
  n_by_fid_norm <- sum(normalize_id(covar$FID) %in% scanner_ids_norm)

  if (n_by_iid >= n_by_fid && n_by_iid > 0) {
    cat(sprintf("  Matching on IID (%d matches)\n", n_by_iid))
    names(scanner_subset)[1] <- "IID"
    covar <- merge(covar, scanner_subset, by = "IID", all.x = TRUE)
  } else if (n_by_fid > 0) {
    cat(sprintf("  Matching on FID (%d matches)\n", n_by_fid))
    names(scanner_subset)[1] <- "FID"
    covar <- merge(covar, scanner_subset, by = "FID", all.x = TRUE)
  } else if (n_by_iid_norm > 0) {
    cat(sprintf("  Matching on normalized IID (%d matches, prefix-stripped)\n", n_by_iid_norm))
    # Build lookup via normalized IDs
    covar$norm_id <- normalize_id(covar$IID)
    scanner_subset$norm_id <- normalize_id(scanner_subset[[scannerPCA_idcol]])
    scanner_subset[[scannerPCA_idcol]] <- NULL
    covar <- merge(covar, scanner_subset, by = "norm_id", all.x = TRUE)
    covar$norm_id <- NULL
  } else if (n_by_fid_norm > 0) {
    cat(sprintf("  Matching on normalized FID (%d matches, prefix-stripped)\n", n_by_fid_norm))
    covar$norm_id <- normalize_id(covar$FID)
    scanner_subset$norm_id <- normalize_id(scanner_subset[[scannerPCA_idcol]])
    scanner_subset[[scannerPCA_idcol]] <- NULL
    covar <- merge(covar, scanner_subset, by = "norm_id", all.x = TRUE)
    covar$norm_id <- NULL
  } else {
    cat("  [WARN] No direct ID match. Trying via original subject IDs...\n")
    # Try matching via the original covar subjectID column
    # The covar might have been created with FID/IID from .fam, but the scanner PCA
    # uses the original subject ID (like from Covariates.csv)
    orig_covar <- read.csv(csvFILE, header = TRUE, stringsAsFactors = FALSE)
    id_candidates <- c("subjectID", "SubjID", "subjID", "ID", "Subject", "eid")
    orig_id_col <- NULL
    for (cand in id_candidates) {
      if (cand %in% names(orig_covar)) { orig_id_col <- cand; break }
    }
    if (!is.null(orig_id_col)) {
      orig_ids <- trimws(as.character(orig_covar[[orig_id_col]]))
      n_via_orig <- sum(orig_ids %in% scanner_ids)
      cat(sprintf("  Trying via '%s': %d matches\n", orig_id_col, n_via_orig))
      if (n_via_orig > 0) {
        # Build FID/IID -> original_id lookup
        # covar already has FID/IID; orig_covar has subjectID
        # pheno has FID/IID and was matched from subjectID
        lookup <- data.frame(
          orig_id = orig_ids,
          stringsAsFactors = FALSE
        )
        # Match orig_covar rows to covar rows via phenotype matching
        # Use IID from covar
        if (orig_id_col %in% names(orig_covar)) {
          names(scanner_subset)[1] <- "orig_id"
          # Add orig_id to covar for merging
          # Try finding orig_id via the pheno file which has both FID/IID and was built from subjectID
          pheno_ids <- pheno[, c("FID", "IID")]
          # The pheno file's row order preserves the original subjectID order
          # But we need a direct mapping... let's try numeric IID matching
          covar$orig_id <- covar$IID  # default
          n_match <- sum(covar$orig_id %in% scanner_ids)
          if (n_match == 0) {
            covar$orig_id <- as.character(as.numeric(covar$IID))
            n_match <- sum(covar$orig_id %in% scanner_ids)
          }
          cat(sprintf("  Matching orig_id: %d matches\n", n_match))
          if (n_match > 0) {
            covar <- merge(covar, scanner_subset, by = "orig_id", all.x = TRUE)
            covar$orig_id <- NULL
          }
        }
      }
    }
  }

  n_with_scanner <- sum(!is.na(covar[[new_names[1]]]))
  n_missing <- sum(is.na(covar[[new_names[1]]]))
  cat(sprintf("  Scanner PCA matched: %d subjects\n", n_with_scanner))
  if (n_missing > 0) {
    cat(sprintf("  [WARN] %d subjects missing scanner PCA — removing\n", n_missing))
    covar <- covar[!is.na(covar[[new_names[1]]]), ]
  }
  cat(sprintf("  Covariates after scanner PCA merge: %d subjects\n", nrow(covar)))

} else if (scannerPCA_file != "") {
  cat(sprintf("\n[WARN] Scanner PCA file not found: %s\n", scannerPCA_file))
}

#=============================================================================
# 7. Merge global DTI metrics from Phenotypes.csv
#=============================================================================
cat("Merging global DTI metrics...\n")

global_cols <- c("FA_GlobalAverage", "MD_GlobalAverage", "AD_GlobalAverage", "RD_GlobalAverage")
available_globals <- intersect(global_cols, names(pheno))

# If standard names not found, try to find and rename Average/Global columns
# Standard ENIGMA output: AverageFA, AverageMD, AverageAD, AverageRD
# After prefix: FA_AverageFA, MD_AverageMD, etc.
if (length(available_globals) < 4) {
  for (metric in c("FA", "MD", "AD", "RD")) {
    std_name <- paste0(metric, "_GlobalAverage")
    if (!std_name %in% names(pheno)) {
      # Priority 1: AverageFA (standard ENIGMA output)
      p1_name <- paste0(metric, "_Average", metric)
      p1 <- if (p1_name %in% names(pheno)) p1_name else character(0)
      # Priority 2: Average_${measure} pattern
      p2 <- grep(paste0("^", metric, "_Average_"), names(pheno), value = TRUE)
      # Priority 3: any Average/Global
      p3 <- grep(paste0("^", metric, "_(Average|Global)"), names(pheno), value = TRUE)

      chosen <- NULL
      if (length(p1) > 0) { chosen <- p1[1] }
      else if (length(p2) > 0) { chosen <- p2[1] }
      else if (length(p3) > 0) { chosen <- p3[1] }

      if (!is.null(chosen)) {
        cat(sprintf("  Renaming '%s' -> '%s'\n", chosen, std_name))
        names(pheno)[names(pheno) == chosen] <- std_name
      }
    }
  }
  available_globals <- intersect(global_cols, names(pheno))
}

if (length(available_globals) == 0) {
  cat("  [WARN] No global DTI metrics found. GWAS will run without global correction.\n")
  cat(sprintf("  Available columns matching 'Average' or 'Global': %s\n",
              paste(grep("Average|Global", names(pheno), value = TRUE), collapse = ", ")))
} else if (length(available_globals) < 4) {
  warning(sprintf("Only %d of 4 expected global DTI metrics found in Phenotypes.csv", length(available_globals)))
}

pheno_globals <- pheno[, c("FID", "IID", available_globals)]
covar <- merge(covar, pheno_globals, by = c("FID", "IID"), all.x = TRUE)

cat(sprintf("  Added: %s\n", paste(available_globals, collapse = ", ")))

#=============================================================================
# 8. Clean up and reorder columns
#=============================================================================

# Remove original age/sex column names if different from standardised
# Also remove any leftover ID columns that shouldn't be covariates
cols_to_drop <- c()
if (ageColumnHeader != "Age" && ageColumnHeader %in% names(covar)) {
  cols_to_drop <- c(cols_to_drop, ageColumnHeader)
}
if (sexColumnHeader != "Sex" && sexColumnHeader %in% names(covar)) {
  cols_to_drop <- c(cols_to_drop, sexColumnHeader)
}
if (affectionStatusColumnHeader != "AffectionStatus" && affectionStatusColumnHeader %in% names(covar)) {
  cols_to_drop <- c(cols_to_drop, affectionStatusColumnHeader)
}
# Drop any leftover subject ID columns (these should NOT be covariates)
id_drop_patterns <- c("subjectID", "SubjID", "subjID", "SID", "sid",
                       "Subject", "subject", "ID")
for (idp in id_drop_patterns) {
  if (idp %in% names(covar) && !idp %in% c("FID", "IID")) {
    cols_to_drop <- c(cols_to_drop, idp)
  }
}
if (length(cols_to_drop) > 0) {
  cat(sprintf("  Dropping from covar output: %s\n", paste(cols_to_drop, collapse=", ")))
}
covar <- covar[, !names(covar) %in% cols_to_drop]

#=============================================================================
# 8b. Age range check
#=============================================================================
cat("\n--- Age range check ---\n")
age_vals <- covar$Age
cat(sprintf("  Min age: %.1f\n  Max age: %.1f\n  Mean age: %.1f\n",
            min(age_vals), max(age_vals), mean(age_vals)))

if (min(age_vals) < 0) {
  cat("  ╔═══════════════════════════════════════════════════════════╗\n")
  cat("  ║  WARNING: Negative age detected! Check your data.       ║\n")
  cat("  ╚═══════════════════════════════════════════════════════════╝\n")
}
if (max(age_vals) > 120) {
  cat("  ╔═══════════════════════════════════════════════════════════╗\n")
  cat("  ║  WARNING: Age > 120 detected! Check your data.          ║\n")
  cat("  ╚═══════════════════════════════════════════════════════════╝\n")
}
if (min(age_vals) < 5) {
  cat("  [NOTE] Very young subjects detected (age < 5). Verify this is expected.\n")
}
if (max(age_vals) > 100) {
  cat("  [NOTE] Very old subjects detected (age > 100). Verify this is expected.\n")
}

# Reorder: FID, IID, Sex, Age, AgeCsq, AgexSex, AgeCsqxSex, AffectionStatus, PCs, Globals, dummies
standard_cols <- c("FID", "IID", "Sex", "Age", "AgeCsq", "AgexSex", "AgeCsqxSex", "AffectionStatus")
pc_names <- paste0("PC", 1:n_pcs)
all_standard <- c(standard_cols, pc_names, available_globals)
dummy_cols <- setdiff(names(covar), all_standard)

covar <- covar[, c(all_standard, dummy_cols)]

#=============================================================================
# 9. Write output files
#=============================================================================

# Extract cohort name from eName
cohort_name <- sub("_enigma_dti_gwas", "", eName)

# Write covariate file (tab-delimited for PLINK2)
covar_outfile <- file.path(outDir, paste0(cohort_name, "_enigma_dti_gwas.covar"))
write.table(covar, file = covar_outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Covariate file: %s (%d subjects, %d columns)\n",
            covar_outfile, nrow(covar), ncol(covar)))

# Write phenotype file (subset of pheno matching covar subjects)
pheno_matched <- pheno[pheno$FID %in% covar$FID & pheno$IID %in% covar$IID, ]
pheno_outfile <- file.path(outDir, paste0(cohort_name, "_enigma_dti_gwas.pheno"))
write.table(pheno_matched, file = pheno_outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Phenotype file: %s (%d subjects, %d columns)\n",
            pheno_outfile, nrow(pheno_matched), ncol(pheno_matched)))

# Write subject ID lists for subgroup analyses
# Controls only
if (any(covar$AffectionStatus == 1)) {
  controls <- covar[covar$AffectionStatus == 1, c("FID", "IID")]
  write.table(controls, file = file.path(outDir, paste0(cohort_name, "_Conly.ids")),
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  cat(sprintf("  Controls ID list: %d subjects\n", nrow(controls)))
}

# Cases only
if (any(covar$AffectionStatus == 2)) {
  cases <- covar[covar$AffectionStatus == 2, c("FID", "IID")]
  write.table(cases, file = file.path(outDir, paste0(cohort_name, "_Aonly.ids")),
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  cat(sprintf("  Cases ID list: %d subjects\n", nrow(cases)))
}

# Males only
males <- covar[covar$Sex == 1, c("FID", "IID")]
write.table(males, file = file.path(outDir, paste0(cohort_name, "_males.ids")),
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
cat(sprintf("  Males ID list: %d subjects\n", nrow(males)))

# Females only
females <- covar[covar$Sex == 2, c("FID", "IID")]
write.table(females, file = file.path(outDir, paste0(cohort_name, "_females.ids")),
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
cat(sprintf("  Females ID list: %d subjects\n", nrow(females)))

#=============================================================================
# 10. Write RUN_NOTES.txt with summary statistics
#=============================================================================
cat("Writing RUN_NOTES.txt...\n")

notes_file <- file.path(outDir, "RUN_NOTES.txt")
sink(notes_file)

cat("ENIGMA-DTI GWAS — RUN NOTES\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Project: %s\n\n", eName))

cat("=== Sample Sizes ===\n")
cat(sprintf("Total N: %d\n", nrow(covar)))
cat(sprintf("Males: %d\n", n_males))
cat(sprintf("Females: %d\n", n_females))
cat(sprintf("Controls: %d\n", sum(covar$AffectionStatus == 1)))
cat(sprintf("Cases: %d\n", sum(covar$AffectionStatus == 2)))
cat(sprintf("Ancestry PCs used: %d\n\n", n_pcs))

cat("=== Age Statistics ===\n")
cat(sprintf("Mean: %.2f\n", mean(covar$Age)))
cat(sprintf("SD: %.2f\n", sd(covar$Age)))
cat(sprintf("Range: %.1f - %.1f\n\n", min(covar$Age), max(covar$Age)))

cat("=== Phenotype Summary Statistics ===\n")
pheno_cols <- setdiff(names(pheno_matched), c("FID", "IID"))
for (col in pheno_cols) {
  vals <- as.numeric(pheno_matched[[col]])
  vals <- vals[!is.na(vals)]
  if (length(vals) > 0) {
    cat(sprintf("%s: N=%d, Mean=%.6f, SD=%.6f, Range=[%.6f, %.6f]\n",
                col, length(vals), mean(vals), sd(vals), min(vals), max(vals)))
  }
}

cat("\n=== Covariate Columns ===\n")
cat(paste(names(covar), collapse = "\n"))
cat("\n")

sink()
cat(sprintf("  RUN_NOTES.txt: %s\n", notes_file))

cat("\n=== Covariate formatting complete ===\n")
