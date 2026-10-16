###############################################################################
# ENIGMA_DTI_GWAS_QC_plots_stats.R
#
# Backend: Generate QC histograms and summary statistics for phenotype
# and covariate data.
#
# Called by run1_ENIGMA-DTI_GWAS_step1.sh — do not run directly.
#
# Arguments (named via --argname value):
#   --cohort       Cohort name
#   --covarFILE    Path to .covar file
#   --phenoFILE    Path to .pheno file
#   --outDir       Output directory for QC files
#   --icvFILE      (optional) Path to ICV CSV file
#   --icvColumnHeader (optional) Column name for ICV
###############################################################################

# Parse named arguments
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list()
  i <- 1
  while (i <= length(args)) {
    if (grepl("^--", args[i])) {
      key <- sub("^--", "", args[i])
      if (i + 1 <= length(args) && !grepl("^--", args[i + 1])) {
        params[[key]] <- args[i + 1]
        i <- i + 2
      } else {
        params[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  return(params)
}

params <- parse_args(args)

cohort    <- params$cohort
covarFILE <- params$covarFILE
phenoFILE <- params$phenoFILE
outDir    <- params$outDir

cat("=== ENIGMA-DTI GWAS QC Plots & Statistics ===\n")
cat(sprintf("  Cohort: %s\n", cohort))

dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

#=============================================================================
# Read data
#=============================================================================
covar <- read.table(covarFILE, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
pheno <- read.table(phenoFILE, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

cat(sprintf("  Covariates: %d subjects\n", nrow(covar)))
cat(sprintf("  Phenotypes: %d subjects\n", nrow(pheno)))

prefix <- paste0(cohort, "_ENIGMA_DTI_GWAS")

#=============================================================================
# Age histograms
#=============================================================================
cat("Generating age histograms...\n")

pdf(file.path(outDir, paste0(prefix, "_Age_histograms.pdf")), width = 10, height = 6)
par(mfrow = c(1, 2))

hist(covar$Age, breaks = 30, main = paste(cohort, "- Age Distribution"),
     xlab = "Age (years)", col = "steelblue", border = "white")

if ("AffectionStatus" %in% names(covar) && length(unique(covar$AffectionStatus)) > 1) {
  age_ctrl <- covar$Age[covar$AffectionStatus == 1]
  age_case <- covar$Age[covar$AffectionStatus == 2]
  hist(age_ctrl, breaks = 30, main = "Age by Affection Status",
       xlab = "Age (years)", col = rgb(0, 0.5, 1, 0.5), border = "white",
       xlim = range(covar$Age))
  hist(age_case, breaks = 30, col = rgb(1, 0.3, 0.3, 0.5), border = "white", add = TRUE)
  legend("topright", c("Controls", "Cases"),
         fill = c(rgb(0, 0.5, 1, 0.5), rgb(1, 0.3, 0.3, 0.5)))
} else {
  plot.new()
  text(0.5, 0.5, "Single affection group\n(no stratification)", cex = 1.2)
}
dev.off()

# Age stats
sink(file.path(outDir, paste0(prefix, "_Age_stats.txt")))
cat(sprintf("Age statistics for %s\n", cohort))
cat(sprintf("N: %d\n", length(covar$Age)))
cat(sprintf("Mean: %.2f\n", mean(covar$Age)))
cat(sprintf("SD: %.2f\n", sd(covar$Age)))
cat(sprintf("Median: %.2f\n", median(covar$Age)))
cat(sprintf("Min: %.1f\n", min(covar$Age)))
cat(sprintf("Max: %.1f\n", max(covar$Age)))
if ("Sex" %in% names(covar)) {
  males <- covar$Age[covar$Sex == 1]
  females <- covar$Age[covar$Sex == 2]
  if (length(males) > 0)
    cat(sprintf("\nMales (N=%d): Mean=%.2f, SD=%.2f\n", length(males), mean(males), sd(males)))
  if (length(females) > 0)
    cat(sprintf("Females (N=%d): Mean=%.2f, SD=%.2f\n", length(females), mean(females), sd(females)))
}
sink()

#=============================================================================
# PC histograms
#=============================================================================
cat("Generating PC histograms...\n")

pc_cols <- grep("^PC[0-9]+$", names(covar), value = TRUE)

if (length(pc_cols) > 0) {
  pdf(file.path(outDir, paste0(prefix, "_PC_histograms.pdf")), width = 12, height = 8)
  nr <- ceiling(sqrt(length(pc_cols)))
  nc <- ceiling(length(pc_cols) / nr)
  par(mfrow = c(nr, nc))

  for (pc in pc_cols) {
    hist(covar[[pc]], breaks = 30, main = pc, xlab = pc,
         col = "steelblue", border = "white")
  }
  dev.off()

  sink(file.path(outDir, paste0(prefix, "_PC_stats.txt")))
  cat(sprintf("PC statistics for %s\n\n", cohort))
  for (pc in pc_cols) {
    vals <- covar[[pc]]
    cat(sprintf("%s: Mean=%.6f, SD=%.6f, Range=[%.6f, %.6f]\n",
                pc, mean(vals), sd(vals), min(vals), max(vals)))
  }
  sink()
}

#=============================================================================
# ROI histograms by metric
#=============================================================================
cat("Generating ROI histograms...\n")

pheno_cols <- setdiff(names(pheno), c("FID", "IID"))

for (metric in c("FA", "MD", "AD", "RD")) {
  metric_cols <- grep(paste0("^", metric, "_"), pheno_cols, value = TRUE)
  if (length(metric_cols) == 0) next

  # All subjects
  pdf(file.path(outDir, paste0(prefix, "_ROIs_histograms_", metric, ".pdf")),
      width = 16, height = 12)
  nr <- ceiling(sqrt(length(metric_cols)))
  nc <- ceiling(length(metric_cols) / nr)
  par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1))

  for (col in metric_cols) {
    vals <- as.numeric(pheno[[col]])
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      hist(vals, breaks = 30, main = gsub(paste0(metric, "_"), "", col),
           xlab = "", col = "steelblue", border = "white", cex.main = 0.8)
    }
  }
  dev.off()

  # By affection status
  if ("AffectionStatus" %in% names(covar) && length(unique(covar$AffectionStatus)) > 1) {
    merged_af <- merge(pheno, covar[, c("FID", "IID", "AffectionStatus")], by = c("FID", "IID"))

    pdf(file.path(outDir, paste0(prefix, "_ROIs_histograms_", metric, "_by_AffectionStatus.pdf")),
        width = 16, height = 12)
    par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1))

    for (col in metric_cols) {
      vals_c <- as.numeric(merged_af[[col]][merged_af$AffectionStatus == 1])
      vals_a <- as.numeric(merged_af[[col]][merged_af$AffectionStatus == 2])
      vals_c <- vals_c[!is.na(vals_c)]
      vals_a <- vals_a[!is.na(vals_a)]

      if (length(vals_c) > 0 || length(vals_a) > 0) {
        all_vals <- c(vals_c, vals_a)
        brk <- seq(min(all_vals), max(all_vals), length.out = 31)
        hist(vals_c, breaks = brk, main = gsub(paste0(metric, "_"), "", col),
             xlab = "", col = rgb(0, 0.5, 1, 0.5), border = "white", cex.main = 0.8)
        if (length(vals_a) > 0)
          hist(vals_a, breaks = brk, col = rgb(1, 0.3, 0.3, 0.5), border = "white", add = TRUE)
      }
    }
    dev.off()
  }
}

#=============================================================================
# Summary statistics CSV
#=============================================================================
cat("Generating summary statistics...\n")

stats_list <- list()
for (col in pheno_cols) {
  vals <- as.numeric(pheno[[col]])
  vals <- vals[!is.na(vals)]
  if (length(vals) > 0) {
    stats_list[[col]] <- data.frame(
      Phenotype = col, N = length(vals),
      Mean = mean(vals), SD = sd(vals),
      Median = median(vals), Min = min(vals), Max = max(vals),
      N_NA = sum(is.na(pheno[[col]])),
      stringsAsFactors = FALSE
    )
  }
}

stats_df <- do.call(rbind, stats_list)
write.csv(stats_df,
          file = file.path(outDir, paste0(prefix, "_ROIs_SummaryStats.csv")),
          row.names = FALSE)

#=============================================================================
# Cohen's d (if case-control)
#=============================================================================
if ("AffectionStatus" %in% names(covar) && length(unique(covar$AffectionStatus)) > 1) {
  cat("Computing Cohen's d...\n")
  merged_cd <- merge(pheno, covar[, c("FID", "IID", "AffectionStatus")], by = c("FID", "IID"))

  d_list <- list()
  for (col in pheno_cols) {
    vals_c <- as.numeric(merged_cd[[col]][merged_cd$AffectionStatus == 1])
    vals_a <- as.numeric(merged_cd[[col]][merged_cd$AffectionStatus == 2])
    vals_c <- vals_c[!is.na(vals_c)]
    vals_a <- vals_a[!is.na(vals_a)]

    if (length(vals_c) > 1 && length(vals_a) > 1) {
      pooled_sd <- sqrt(((length(vals_c) - 1) * var(vals_c) +
                          (length(vals_a) - 1) * var(vals_a)) /
                         (length(vals_c) + length(vals_a) - 2))
      if (pooled_sd > 0) {
        d <- (mean(vals_a) - mean(vals_c)) / pooled_sd
        d_list[[col]] <- data.frame(Phenotype = col, Cohens_d = d, stringsAsFactors = FALSE)
      }
    }
  }

  if (length(d_list) > 0) {
    d_df <- do.call(rbind, d_list)
    write.csv(d_df,
              file = file.path(outDir, paste0(prefix, "_Cohensd_results.csv")),
              row.names = FALSE)

    png(file.path(outDir, paste0(prefix, "_Cohensd_bar_graph_AFF-CON.png")),
        width = 1200, height = 600)
    par(mar = c(10, 4, 3, 1))
    barplot(d_df$Cohens_d, names.arg = d_df$Phenotype, las = 2,
            col = ifelse(d_df$Cohens_d > 0, "coral", "steelblue"),
            main = paste(cohort, "- Cohen's d (Cases vs Controls)"),
            ylab = "Cohen's d", cex.names = 0.5)
    abline(h = 0)
    dev.off()
  }
}

#=============================================================================
# ICV plots (optional)
#=============================================================================
if (!is.null(params$icvFILE) && file.exists(params$icvFILE)) {
  cat("Generating ICV plots...\n")
  icv <- read.csv(params$icvFILE, header = TRUE, stringsAsFactors = FALSE)
  icv_col <- params$icvColumnHeader

  if (icv_col %in% names(icv)) {
    # Try matching on IID or SubjID
    if ("SubjID" %in% names(icv)) {
      icv_merged <- merge(covar[, c("FID", "IID", "Sex")], icv, by.x = "IID", by.y = "SubjID")
    } else if ("IID" %in% names(icv)) {
      icv_merged <- merge(covar[, c("FID", "IID", "Sex")], icv, by = "IID")
    } else {
      icv_merged <- merge(covar[, c("FID", "IID", "Sex")], icv,
                          by.x = "IID", by.y = names(icv)[1])
    }

    if (nrow(icv_merged) > 0) {
      pdf(file.path(outDir, paste0(prefix, "_ICV_hist_by_Sex.pdf")), width = 10, height = 6)
      par(mfrow = c(1, 2))
      icv_m <- icv_merged[[icv_col]][icv_merged$Sex == 1]
      icv_f <- icv_merged[[icv_col]][icv_merged$Sex == 2]
      hist(icv_m, breaks = 30, main = "ICV - Males", xlab = "ICV", col = "steelblue", border = "white")
      hist(icv_f, breaks = 30, main = "ICV - Females", xlab = "ICV", col = "coral", border = "white")
      dev.off()

      pdf(file.path(outDir, paste0(prefix, "_ICV_box_by_Sex.pdf")), width = 8, height = 6)
      boxplot(icv_merged[[icv_col]] ~ icv_merged$Sex, names = c("Males", "Females"),
              main = paste(cohort, "- ICV by Sex"), ylab = "ICV",
              col = c("steelblue", "coral"))
      dev.off()
    }
  }
}

cat("\n=== QC plots and statistics complete ===\n")
cat(sprintf("  Output directory: %s\n", outDir))
