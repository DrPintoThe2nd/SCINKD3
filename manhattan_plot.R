#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript manhattan_plot.R input.tsv.gz NUM_CHROMOSOMES")
}

infile <- args[1]
max_chr <- as.integer(args[2])

# Output prefix from input filename
out_prefix <- sub("\\.bed\\.gz$|\\.gz$|\\.bed$", "", basename(infile))
png_file <- paste0(out_prefix, ".manhattan.png")
pdf_file <- paste0(out_prefix, ".manhattan.pdf")

# Read gzipped TSV
df <- read.table(gzfile(infile),
                 header = FALSE,
                 sep = "\t",
                 col.names = c("chr", "start", "end", "value"),
                 stringsAsFactors = FALSE)

# Keep chromosome order as in file
chr_levels <- unique(df$chr)

# Restrict to first X number of chromosomes
chr_levels <- chr_levels[1:min(max_chr, length(chr_levels))]
df <- df[df$chr %in% chr_levels, ]

# Midpoint
df$pos <- (df$start + df$end) / 2

# Preserve order
df$chr <- factor(df$chr, levels = chr_levels)
df <- df[order(df$chr, df$pos), ]

# Chromosome index
df$chr_index <- as.integer(df$chr)

# Chromosome lengths
chr_lengths <- tapply(df$pos, df$chr, max)

# Cumulative offsets
chr_offsets <- c(0, cumsum(chr_lengths[-length(chr_lengths)]))
names(chr_offsets) <- names(chr_lengths)

df$cum_pos <- df$pos + chr_offsets[as.character(df$chr)]

# Alternate colors
df$col <- ifelse(df$chr_index %% 2 == 1, "black", "grey")

# X-axis centers
chr_centers <- chr_offsets + chr_lengths / 2

plot_manhattan <- function() {
  plot(df$cum_pos, log2(df$value),
       col = df$col,
       pch = 16,
#       cex = 0.6,
       xaxt = "n",
       xlab = "Chromosome",
       ylab = "log2(Read Depth)",
       main = paste("Manhattan Plot (first", length(chr_levels), "chromosomes)"))
  
  axis(1, at = chr_centers,
       labels = names(chr_centers),
       las = 2,
       cex.axis = 0.7)
  
  box()
}

# PNG output
png(png_file, width = 1600, height = 600)
plot_manhattan()
dev.off()

# PDF output
pdf(pdf_file, width = 12, height = 5)
plot_manhattan()
dev.off()

cat("Plots written to:\n", png_file, "\n", pdf_file, "\n")
