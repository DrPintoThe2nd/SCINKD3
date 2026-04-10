#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript mf_depth_ratio_multi.R <male_files> <female_files> <N_chromosomes>")
}

male_files   <- strsplit(args[1], ",")[[1]]
female_files <- strsplit(args[2], ",")[[1]]
n_chr        <- as.integer(args[3])

# Output prefix
out_prefix <- sub("\\.bed\\.gz$|\\.gz$|\\.bed$", "", basename("SCINKD3"))
png_file <- paste0(out_prefix, ".mf_ratio.png")
pdf_file <- paste0(out_prefix, ".mf_ratio.pdf")

# --- FUNCTION TO READ FILE ---
read_depth <- function(f) {
  read.table(gzfile(f),
             header = FALSE,
             sep = "\t",
             col.names = c("chr", "start", "end", "value"),
             stringsAsFactors = FALSE)
}

# --- READ FILES ---
male_list   <- lapply(male_files, read_depth)
female_list <- lapply(female_files, read_depth)

# --- sanity checks ---
n_rows <- nrow(male_list[[1]])

for (df in c(male_list, female_list)) {
  if (nrow(df) != n_rows) {
    stop("All files must have identical number of rows")
  }
}

ref <- male_list[[1]]

check_coords <- function(df) {
  all(df$chr == ref$chr &
      df$start == ref$start &
      df$end == ref$end)
}

if (!all(sapply(c(male_list, female_list), check_coords))) {
  stop("All files must have identical coordinates/order")
}

# --- BUILD MATRICES ---
male_mat   <- do.call(cbind, lapply(male_list, function(df) df$value))
female_mat <- do.call(cbind, lapply(female_list, function(df) df$value))

# --- LIMIT TO FIRST N CHROMOSOMES ---
chr_order <- unique(ref$chr)
keep_chr <- chr_order[seq_len(min(n_chr, length(chr_order)))]

keep_idx <- ref$chr %in% keep_chr

# Subset consistently
ref        <- ref[keep_idx, ]
male_mat   <- male_mat[keep_idx, , drop = FALSE]
female_mat <- female_mat[keep_idx, , drop = FALSE]

# --- ROW MEANS ---
male_mean   <- rowMeans(male_mat, na.rm = TRUE)
female_mean <- rowMeans(female_mat, na.rm = TRUE)

# --- MEDIAN NORMALIZATION ---
male_norm   <- male_mean / median(male_mean, na.rm = TRUE)
female_norm <- female_mean / median(female_mean, na.rm = TRUE)

# --- LOG2 M/F ---
log2_mf <- log2(male_norm / female_norm)

# Replace Inf / NA with 0 (for plotting stability if needed)
log2_mf[!is.finite(log2_mf)] <- 0

# --- X positions (preserve order) ---
x <- seq_along(log2_mf)

# --- ALTERNATING COLORS BY CHROMOSOME (NO REORDERING) ---
chr <- ref$chr
chr_index <- as.integer(factor(chr, levels = unique(chr)))
cols <- ifelse(chr_index %% 2 == 1, "black", "grey")

plot_ratio <- function() {
  plot(x, log2_mf,
       col = cols,
       pch = 16,
       cex = 1,
       xlab = "Genomic Window (input order)",
       ylab = "log2(M/F Read Depth)",
       main = paste("log2(M/F) Depth (",
                    length(male_files), "M vs",
                    length(female_files), "F )"))
  
  abline(h = 0, col = "red", lty = 2)
  box()
}

# PNG
png(png_file, width = 1600, height = 600)
plot_ratio()
dev.off()

# PDF
pdf(pdf_file, width = 12, height = 5)
plot_ratio()
dev.off()

cat("Plots written to:\n", png_file, "\n", pdf_file, "\n")