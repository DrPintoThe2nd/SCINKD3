#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<4) stop("Usage: Rscript chrom_dotplot.R <prefix> <file1> <file2> <N>")

prefix <- args[1]
f1 <- args[2]
f2 <- args[3]
N  <- as.integer(args[4])

d1 <- read.delim(f1, stringsAsFactors=FALSE)
d2 <- read.delim(f2, stringsAsFactors=FALSE)

# Get total mean values
total1 <- d1$mean[d1$chrom == "total"]
total2 <- d2$mean[d2$chrom == "total"]

# Remove *_region and keep chrom, length, mean
d1 <- d1[!grepl("_region$", d1$chrom) & d1$chrom != "total", c("chrom","length","mean")]
d2 <- d2[!grepl("_region$", d2$chrom) & d2$chrom != "total", c("chrom","length","mean")]

# Normalize mean by total
d1$mean_norm <- d1$mean / total1
d2$mean_norm <- d2$mean / total2

# Limit to first N chromosomes
chrs <- unique(d1$chrom)[1:min(N, length(unique(d1$chrom)))]

d1 <- d1[d1$chrom %in% chrs, ]
d2 <- d2[d2$chrom %in% chrs, ]

d1$Dataset <- "file1"
d2$Dataset <- "file2"

d <- rbind(d1, d2)

p <- ggplot(d, aes(x=length, y=mean_norm, color=Dataset)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE) +
  labs(x="Chromosome length",
       y="Mean / total(mean)",
       title=prefix) +
  theme_minimal()

ggsave(paste0(prefix,".dotplot.png"), p, width=12, height=5, dpi=150)
ggsave(paste0(prefix,".dotplot.pdf"), p, width=12, height=5)

