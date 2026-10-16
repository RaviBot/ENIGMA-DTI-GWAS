###############################################################################
# ENIGMA_DTI_pheno_merge.R
#
# Backend: Merge available DTI phenotype CSVs into Phenotypes.csv.
# Handles 1-4 metric files. Excludes statistical QC outliers and
# subjects not in the .fam file. Logs every step.
#
# Arguments (positional via --args):
#   1 = semicolon-separated list of CSV file paths
#   2 = semicolon-separated list of metric names (FA;MD;AD;RD)
#   3 = localfamFILE
#   4 = outFolder
#   5 = outlierIDsFile
#   6 = subjID_column
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("Insufficient arguments.")

csv_files_str  <- args[1]
metrics_str    <- args[2]
localfamFILE   <- args[3]
outFolder      <- args[4]
outlierIDsFile <- args[5]
subjID_column  <- args[6]

csv_files <- unlist(strsplit(csv_files_str, ";"))
metrics   <- unlist(strsplit(metrics_str, ";"))

#=============================================================================
# ID normalization helper
#   Strips common prefixes (sub-, SUB-, Sub-, sub_, SUB_, Sub_) so that
#   IDs from different sources (DTI CSVs, Covariates, .fam) can be matched
#   even when one source uses "sub-100206" and another uses "100206".
#=============================================================================
normalize_id <- function(ids) {
  gsub("^[Ss][Uu][Bb][-_]", "", trimws(as.character(ids)))
}

cat("=== ENIGMA-DTI Phenotype Merge ===\n")
cat(sprintf("  Date: %s\n", Sys.time()))
cat(sprintf("  Metrics: %s\n", paste(metrics, collapse=", ")))
cat(sprintf("  Subject ID column: %s\n\n", subjID_column))

#=============================================================================
# 1. Read phenotype files
#=============================================================================
cat("--- Reading phenotype files ---\n")
dti_list <- list()
for (i in seq_along(csv_files)) {
  df <- read.csv(csv_files[i], header=TRUE, stringsAsFactors=FALSE,
                  na.strings=c("NA","na","","N/A","-999","NaN"),
                  strip.white=TRUE, check.names=FALSE)
  # Find ID column
  id_col <- subjID_column
  if (!id_col %in% names(df)) {
    possible <- c("SubjID","subjectID","subjID","IID","ID","Subject")
    found <- possible[possible %in% names(df)]
    if (length(found)>0) { id_col <- found[1] } else { id_col <- names(df)[1] }
  }
  # Clean IDs: remove quotes, whitespace
  df[[id_col]] <- trimws(gsub('"', '', as.character(df[[id_col]])))
  df <- df[df[[id_col]] != "" & !is.na(df[[id_col]]), ]

  cat(sprintf("  %s: %d subjects, %d columns (ID col: %s)\n",
              metrics[i], nrow(df), ncol(df), id_col))

  # Drop known non-ROI columns BEFORE prefixing
  # These are demographic/metadata columns that should NOT be in the phenotype file
  non_roi_patterns <- c("^age$", "^sex$", "^Sex$", "^AGE$", "^gender$", "^Gender$",
                        "^diagnosis$", "^Diagnosis$", "^dx$", "^DX$",
                        "^site$", "^Site$", "^SITE$", "^scanner$", "^Scanner$",
                        "^SID$", "^sid$", "^subjectID$", "^SubjID$", "^subjID$",
                        "^IID$", "^FID$", "^ID$", "^Subject$",
                        "^AverageROI")  # AverageROI is NOT the global avg we want
  cols_to_drop <- c()
  for (pat in non_roi_patterns) {
    matches <- grep(pat, names(df), value = TRUE)
    cols_to_drop <- c(cols_to_drop, matches)
  }
  # Don't drop the ID column (we need it for merging)
  cols_to_drop <- setdiff(cols_to_drop, id_col)
  if (length(cols_to_drop) > 0) {
    cat(sprintf("    Dropping non-ROI columns: %s\n", paste(cols_to_drop, collapse=", ")))
    df <- df[, !names(df) %in% cols_to_drop, drop=FALSE]
  }

  # Rename remaining columns with metric prefix (skip ID col)
  col_names <- names(df)
  for (j in which(col_names != id_col)) {
    old <- col_names[j]
    new_name <- paste0(metrics[i], "_", old)
    col_names[j] <- new_name
  }
  # Standardise global average name
  # ENIGMA-DTI pipeline outputs: AverageFA, AverageMD, AverageAD, AverageRD
  # Some sites may have: Average_${measure}, AverageROI_${measure}
  # After metric prefix: FA_AverageFA, FA_Average_${measure}, etc.
  # All should become: FA_GlobalAverage
  target_name <- paste0(metrics[i], "_GlobalAverage")
  if (!target_name %in% col_names) {
    # Priority:
    #   1. AverageFA/AverageMD etc. (standard ENIGMA output)
    #   2. Average_${measure} or Average_ pattern
    #   3. Any remaining Average/Global column

    avg_metric_pat <- paste0(metrics[i], "_Average", metrics[i])  # FA_AverageFA
    avg_measure_pat <- paste0(metrics[i], "_Average_")             # FA_Average_*

    avg_metric_idx <- which(col_names == avg_metric_pat)
    avg_measure_idx <- grep(gsub("\\$", "\\\\$", paste0("^", avg_measure_pat)), col_names)
    avg_any_idx <- grep(paste0("^", metrics[i], "_(Average|Global)"), col_names)

    if (length(avg_metric_idx) > 0) {
      old_name <- col_names[avg_metric_idx[1]]
      col_names[avg_metric_idx[1]] <- target_name
      cat(sprintf("    Renamed '%s' -> '%s'\n", old_name, target_name))
    } else if (length(avg_measure_idx) > 0) {
      old_name <- col_names[avg_measure_idx[1]]
      col_names[avg_measure_idx[1]] <- target_name
      cat(sprintf("    Renamed '%s' -> '%s'\n", old_name, target_name))
    } else if (length(avg_any_idx) > 0) {
      old_name <- col_names[avg_any_idx[1]]
      col_names[avg_any_idx[1]] <- target_name
      cat(sprintf("    Renamed '%s' -> '%s' (fallback)\n", old_name, target_name))
    } else {
      cat(sprintf("    [WARN] No global average column found for %s\n", metrics[i]))
    }
  }
  names(df) <- col_names

  dti_list[[metrics[i]]] <- df
}
n_start <- nrow(dti_list[[1]])

#=============================================================================
# 2. Statistical QC exclusion
#=============================================================================
cat("\n--- Statistical QC exclusion ---\n")
outlier_ids <- character(0)
if (file.exists(outlierIDsFile) && file.info(outlierIDsFile)$size > 0) {
  outlier_ids <- trimws(readLines(outlierIDsFile))
  outlier_ids <- outlier_ids[outlier_ids != ""]
}
cat(sprintf("  Outlier IDs to exclude: %d\n", length(outlier_ids)))

if (length(outlier_ids) > 0) {
  # Normalize outlier IDs so "sub-100206" matches "100206" and vice versa
  outlier_ids_norm <- normalize_id(outlier_ids)
  for (m in names(dti_list)) {
    n_before <- nrow(dti_list[[m]])
    dti_norm <- normalize_id(dti_list[[m]][[id_col]])
    dti_list[[m]] <- dti_list[[m]][!dti_norm %in% outlier_ids_norm, ]
    n_after <- nrow(dti_list[[m]])
    cat(sprintf("  %s: %d -> %d (removed %d)\n", m, n_before, n_after, n_before-n_after))
  }
  cat("  Removed IDs:\n")
  for (oid in outlier_ids) cat(sprintf("    - %s\n", oid))
}

#=============================================================================
# 3. Merge all available metrics
#=============================================================================
cat("\n--- Merging ---\n")
merged <- dti_list[[1]]
if (length(dti_list) > 1) {
  for (i in 2:length(dti_list)) {
    merged <- merge(merged, dti_list[[i]], by=id_col, all=FALSE)
    cat(sprintf("  After merging %s: %d subjects\n", names(dti_list)[i], nrow(merged)))
  }
}
cat(sprintf("  Merged total: %d subjects, %d columns\n", nrow(merged), ncol(merged)))

#=============================================================================
# 4. Match to .fam
#=============================================================================
cat("\n--- Matching to genotype .fam ---\n")
fam <- read.table(localfamFILE, header=FALSE, stringsAsFactors=FALSE)
names(fam) <- c("FID","IID","PID","MID","Sex","Pheno")
cat(sprintf("  .fam subjects: %d\n", nrow(fam)))

# Pre-compute normalized versions of .fam IDs for prefix-tolerant matching
fam_IID_norm  <- normalize_id(fam$IID)
fam_FID_norm  <- normalize_id(fam$FID)
fam_FIDIID    <- paste(fam$FID, fam$IID, sep="_")
fam_FIDIID_norm <- normalize_id(fam_FIDIID)

# Check if normalization changes any IDs (for logging)
dti_ids_raw  <- merged[[id_col]]
dti_ids_norm <- normalize_id(dti_ids_raw)
norm_active  <- any(dti_ids_raw != dti_ids_norm) || any(fam$IID != fam_IID_norm)
if (norm_active) {
  cat("  [NOTE] ID normalization active (stripping sub-/SUB-/Sub- prefixes)\n")
}

merged$FID <- NA_character_
merged$IID <- NA_character_

n_matched <- 0; n_unmatched <- 0
for (i in 1:nrow(merged)) {
  subj      <- merged[i, id_col]
  subj_norm <- dti_ids_norm[i]

  # 1. Exact match on IID
  idx <- which(fam$IID == subj)
  if (length(idx)==1) { merged$FID[i]<-fam$FID[idx]; merged$IID[i]<-fam$IID[idx]; n_matched<-n_matched+1; next }
  # 2. Exact match on FID_IID
  idx <- which(fam_FIDIID == subj)
  if (length(idx)==1) { merged$FID[i]<-fam$FID[idx]; merged$IID[i]<-fam$IID[idx]; n_matched<-n_matched+1; next }
  # 3. Exact match on FID
  idx <- which(fam$FID == subj)
  if (length(idx)==1) { merged$FID[i]<-fam$FID[idx]; merged$IID[i]<-fam$IID[idx]; n_matched<-n_matched+1; next }
  # 4. Normalized match on IID (e.g., sub-100206 vs 100206)
  idx <- which(fam_IID_norm == subj_norm)
  if (length(idx)==1) { merged$FID[i]<-fam$FID[idx]; merged$IID[i]<-fam$IID[idx]; n_matched<-n_matched+1; next }
  # 5. Normalized match on FID_IID
  idx <- which(fam_FIDIID_norm == subj_norm)
  if (length(idx)==1) { merged$FID[i]<-fam$FID[idx]; merged$IID[i]<-fam$IID[idx]; n_matched<-n_matched+1; next }
  # 6. Normalized match on FID
  idx <- which(fam_FID_norm == subj_norm)
  if (length(idx)==1) { merged$FID[i]<-fam$FID[idx]; merged$IID[i]<-fam$IID[idx]; n_matched<-n_matched+1; next }

  n_unmatched <- n_unmatched + 1
}
cat(sprintf("  Matched: %d | Unmatched: %d\n", n_matched, n_unmatched))

if (n_unmatched > 0) {
  unmatched <- merged[[id_col]][is.na(merged$FID)]
  cat("  Unmatched IDs (removed):\n")
  for (u in unmatched) cat(sprintf("    - %s\n", u))
  merged <- merged[!is.na(merged$FID), ]
}

n_final <- nrow(merged)
cat(sprintf("  Final: %d subjects\n", n_final))

#=============================================================================
# 5. Reorder and write
#=============================================================================
pheno_cols <- setdiff(names(merged), c("FID","IID",id_col))
merged <- merged[, c("FID","IID", pheno_cols)]

outFile <- file.path(outFolder, "Phenotypes.csv")
write.csv(merged, file=outFile, row.names=FALSE, quote=FALSE)
cat(sprintf("\n  Written: %s (%d x %d)\n", outFile, nrow(merged), ncol(merged)))

#=============================================================================
# 6. Write QC log
#=============================================================================
logFile <- file.path(outFolder, "QC_statistical", "phenotype_merge_QC_log.txt")
sink(logFile)
cat(sprintf("ENIGMA-DTI Phenotype Merge QC Log\nDate: %s\n\n", Sys.time()))
cat("=== SUBJECT FLOW ===\n")
cat(sprintf("Starting subjects (first metric): %d\n", n_start))
cat(sprintf("Statistical QC outliers removed:  %d\n", length(outlier_ids)))
cat(sprintf("After merge across metrics:       %d\n", nrow(merged)+n_unmatched))
cat(sprintf("No genotype match removed:        %d\n", n_unmatched))
cat(sprintf("FINAL:                            %d\n", n_final))
cat(sprintf("\nMetrics included: %s\n", paste(metrics, collapse=", ")))
cat(sprintf("Phenotype columns: %d\n", length(pheno_cols)))
if (norm_active) {
  cat("\nNOTE: ID normalization was active (sub-/SUB-/Sub- prefix stripping).\n")
  cat("      Some source files used prefixed IDs that were matched after normalization.\n")
}
if (length(outlier_ids)>0) {
  cat(sprintf("\n=== OUTLIERS REMOVED (%d) ===\n", length(outlier_ids)))
  for (o in outlier_ids) cat(sprintf("  %s\n", o))
}
sink()
cat(sprintf("  QC log: %s\n", logFile))
cat("\n=== Phenotype merge complete ===\n")
