#!/bin/bash
###############################################################################
# ENIGMA_DTI_plot_one.sh
#
# Backend: Generate Manhattan + QQ plots for ONE combined GWAS file.
# Called as a qsub array task.
#
# Arguments:
#   $1 = file list (one combined file path per line)
#   $2 = plots output directory
#   $3 = cohort name
#
# Uses SGE_TASK_ID to pick which file to process.
###############################################################################

set -uo pipefail

file_list="$1"
plots_dir="$2"
cohort="$3"

mkdir -p "${plots_dir}"

# Get this task's file
combined_file=$(sed -n "${SGE_TASK_ID}p" "${file_list}")

if [[ -z "${combined_file}" ]] || [[ ! -f "${combined_file}" ]]; then
    echo "ERROR: No file for task ${SGE_TASK_ID}"
    exit 1
fi

grp=$(basename "${combined_file}" _allchr.txt.gz)

# Skip if plots already exist
if [[ -f "${plots_dir}/${grp}_Manhattan.png" ]] && [[ -f "${plots_dir}/${grp}_QQ.png" ]]; then
    echo "[DONE] ${grp} plots already exist"
    exit 0
fi

echo "Plotting: ${grp}"

R --no-save --slave << RSCRIPT
f       <- "${combined_file}"
outdir  <- "${plots_dir}"
cohort  <- "${cohort}"
grp     <- "${grp}"

# Extract phenotype name
pheno <- sub(".*ENIGMA-DTI_[^_]+_", "", grp)
pheno <- sub("_[0-9]{8}_.*", "", pheno)

cat(sprintf("  Reading %s\n", basename(f)))

# Read header to find column positions
con <- gzfile(f, "r")
header_line <- readLines(con, n = 1)
close(con)

header <- strsplit(header_line, "\t")[[1]]
header <- gsub("^#", "", header)

chr_col <- which(header %in% c("CHROM", "CHR"))[1]
bp_col  <- which(header %in% c("POS", "BP"))[1]
p_col   <- which(header == "P")[1]

if (is.na(chr_col) || is.na(bp_col) || is.na(p_col)) {
  cat(sprintf("  [SKIP] Missing CHR/BP/P columns. Have: %s\n", paste(header, collapse=", ")))
  quit(status = 0)
}

# Read only the 3 columns we need
col_classes <- rep("NULL", length(header))
col_classes[chr_col] <- "character"
col_classes[bp_col]  <- "integer"
col_classes[p_col]   <- "numeric"

con <- gzfile(f, "r")
df <- read.table(con, header = TRUE, colClasses = col_classes,
                  comment.char = "", sep = "\t", fill = TRUE)
close(con)

names(df) <- c("CHR", "BP", "P")[1:ncol(df)]
cat(sprintf("  %s variants\n", format(nrow(df), big.mark = ",")))

# --- QQ plot ---
pvals <- df\$P[!is.na(df\$P) & df\$P > 0 & df\$P <= 1]
if (length(pvals) >= 10) {
  observed <- -log10(sort(pvals))
  expected <- -log10(ppoints(length(pvals)))
  chisq <- qchisq(1 - pvals, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)

  qq_file <- file.path(outdir, paste0(grp, "_QQ.png"))
  png(qq_file, width = 600, height = 600, res = 100)
  plot(expected, observed,
       xlab = expression(Expected ~ -log[10](p)),
       ylab = expression(Observed ~ -log[10](p)),
       main = paste(cohort, "-", pheno), pch = 20, cex = 0.4, col = "steelblue")
  abline(0, 1, col = "red", lwd = 2)
  legend("topleft", bty = "n",
         legend = sprintf("lambda = %.3f\nN = %s", lambda, format(length(pvals), big.mark = ",")))
  dev.off()
  cat(sprintf("  QQ: lambda=%.3f\n", lambda))
} else {
  cat("  [SKIP] Too few valid P values for QQ plot\n")
}

# --- Manhattan plot ---
mdf <- df[!is.na(df\$P) & df\$P > 0 & !is.na(df\$CHR) & !is.na(df\$BP), ]
if (nrow(mdf) >= 100) {
  mdf\$CHR_num <- suppressWarnings(as.numeric(ifelse(mdf\$CHR == "X", 23, mdf\$CHR)))
  mdf <- mdf[!is.na(mdf\$CHR_num), ]
  mdf <- mdf[order(mdf\$CHR_num, mdf\$BP), ]

  chr_max <- aggregate(BP ~ CHR_num, mdf, max)
  chr_max <- chr_max[order(chr_max\$CHR_num), ]
  chr_max\$offset <- cumsum(c(0, head(chr_max\$BP, -1)))
  mdf <- merge(mdf, chr_max[, c("CHR_num", "offset")], by = "CHR_num")
  mdf\$BPcum <- mdf\$BP + mdf\$offset

  axis_df <- aggregate(BPcum ~ CHR_num, mdf, function(x) mean(range(x)))
  axis_labels <- ifelse(axis_df\$CHR_num == 23, "X", as.character(axis_df\$CHR_num))
  col_vec <- ifelse(mdf\$CHR_num %% 2 == 0, "steelblue", "coral")

  man_file <- file.path(outdir, paste0(grp, "_Manhattan.png"))
  png(man_file, width = 1400, height = 500, res = 100)
  par(mar = c(4, 4, 3, 1))
  plot(mdf\$BPcum, -log10(mdf\$P), col = col_vec, pch = 20, cex = 0.3,
       xlab = "Chromosome", ylab = expression(-log[10](p)),
       main = paste(cohort, "-", pheno), xaxt = "n",
       ylim = c(0, max(-log10(mdf\$P)) * 1.1))
  axis(1, at = axis_df\$BPcum, labels = axis_labels, las = 2, cex.axis = 0.8)
  abline(h = -log10(5e-8), col = "red", lty = 2, lwd = 1.5)
  abline(h = -log10(1e-5), col = "blue", lty = 3, lwd = 1)
  dev.off()

  cat(sprintf("  Manhattan: minP=%.2e\n", min(mdf\$P, na.rm = TRUE)))
}

# Write summary line
summary_file <- file.path(outdir, paste0(grp, "_summary.txt"))
writeLines(sprintf("%s\t%s\t%d\t%.4f\t%.2e",
    grp, pheno, nrow(df),
    ifelse(exists("lambda"), lambda, NA),
    min(df\$P, na.rm = TRUE)), summary_file)

cat("  Done\n")
RSCRIPT
