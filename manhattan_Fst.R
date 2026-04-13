#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3){
  stop("Usage: Rscript manhattan_plot.R <input_file> <output_prefix> <N_chromosomes>")
}

infile <- args[1]
out_prefix <- args[2]
Nchrom <- as.numeric(args[3])

#-----------------------------
# Load data
#-----------------------------
dt <- fread(infile)

# Check required columns
required_cols <- c("CHROM", "MEAN_FST")
missing_cols <- setdiff(required_cols, colnames(dt))
if(length(missing_cols) > 0){
  stop(paste("Missing required columns:", paste(missing_cols, collapse=", ")))
}

#-----------------------------
# Keep only first N chromosomes (in order of appearance)
#-----------------------------
chrom_levels <- unique(dt$CHROM)

if(Nchrom > length(chrom_levels)){
  warning("Requested more chromosomes than exist. Using all available.")
  Nchrom <- length(chrom_levels)
}

keep_chroms <- chrom_levels[1:Nchrom]
dt <- dt[CHROM %in% keep_chroms]

# Preserve original order
dt$CHROM <- factor(dt$CHROM, levels = keep_chroms)

#-----------------------------
# Create cumulative positions for Manhattan plotting
#-----------------------------
dt[, pos_index := .I]

#-----------------------------
# Assign alternating colors
#-----------------------------
chrom_index <- as.numeric(dt$CHROM)
colors <- ifelse(chrom_index %% 2 == 1, "black", "grey")

#-----------------------------
# Plot
#-----------------------------
png(paste0(out_prefix, ".png"), width=1600, height=600)
par(mar=c(5,5,2,2))

plot(
  dt$pos_index,
  dt$MEAN_FST,
  pch=19,
  col=colors,
  xlab="Chromosome",
  ylab="Fst",
  main="Manhattan Plot",
  ylim=range(0, 1)
)

#-----------------------------
# Add chromosome boundaries
#-----------------------------
chrom_breaks <- dt[, .(mid = mean(pos_index)), by=CHROM]

axis(1, at=chrom_breaks$mid, labels=chrom_breaks$CHROM, las=2, cex.axis=0.7)

#-----------------------------
# Mean + significance lines
#-----------------------------
mean_val <- mean(dt$MEAN_FST, na.rm=TRUE)
sig_val <- mean_val + 4*sd(dt$MEAN_FST, na.rm=TRUE)

abline(h = mean_val, col = "grey")
abline(h = sig_val, col = "red")

legend(
  "topright",
  legend=c("Mean", "99.7% CI"),
  col=c("grey", "red"),
  lty=1,
  cex=0.9,
  text.font=4
)

dev.off()

#-----------------------------
# Optional PDF output
#-----------------------------
pdf(paste0(out_prefix, ".pdf"), width=12, height=5)

plot(
  dt$pos_index,
  dt$MEAN_FST,
  pch=19,
  col=colors,
  xlab="Chromosome",
  ylab="Fst",
  main="Manhattan Plot",
  ylim=range(0, 1)
)

axis(1, at=chrom_breaks$mid, labels=chrom_breaks$CHROM, las=2, cex.axis=0.7)
abline(h = mean_val, col = "grey")
abline(h = sig_val, col = "red")

legend(
  "topright",
  legend=c("Mean", "99.7% CI"),
  col=c("grey", "red"),
  lty=1,
  cex=0.9,
  text.font=4
)

dev.off()