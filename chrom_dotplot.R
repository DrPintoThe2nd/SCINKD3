#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<4) stop("Usage: Rscript chrom_dotplot.R <prefix> <file1> <file2> <N>")

prefix <- args[1]; f1 <- args[2]; f2 <- args[3]; N <- as.integer(args[4])

d1 <- read.delim(f1, stringsAsFactors=FALSE)
d2 <- read.delim(f2, stringsAsFactors=FALSE)

d1 <- d1[!grepl("_region$", d1$chrom), c("chrom","mean")]
d2 <- d2[!grepl("_region$", d2$chrom), c("chrom","mean")]

chrs <- unique(d1$chrom)[1:min(N, length(unique(d1$chrom)))]

d1 <- d1[d1$chrom %in% chrs, ]; d2 <- d2[d2$chrom %in% chrs, ]
d1$Dataset <- "file1"; d2$Dataset <- "file2"

d <- rbind(d1, d2)
d$chrom <- factor(d$chrom, levels=chrs)
d$x <- as.numeric(d$chrom)

p <- ggplot(d, aes(x=x, y=mean, color=Dataset)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE) +
  scale_x_continuous(breaks=seq_along(chrs), labels=chrs) +
  labs(x="Chromosome", y="Mean", title=prefix) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90,hjust=1))

ggsave(paste0(prefix,".dotplot.png"), p, width=12, height=5, dpi=150)
ggsave(paste0(prefix,".dotplot.pdf"), p, width=12, height=5)
