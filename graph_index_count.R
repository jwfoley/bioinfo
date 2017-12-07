#! /usr/bin/Rscript

# graph the output (STDERR) of index_count.py


library(ggplot2)
library(scales)

args <- commandArgs(trailingOnly = T)
if (length(args) != 1) stop("usage: collision_analysis.R graph.pdf < counts.tsv")
out.file = args[1]

index.counts <- read.table(file("stdin"), sep = "\t")
colnames(index.counts) <- c("sample", "index", "sequence", "frequency")

index.counts$sequence <- factor(index.counts$sequence, levels = with(index.counts, sequence[order(frequency, decreasing = T)]))
index.counts$used <- index.counts$sample != ""

result.plot <- ggplot(index.counts, aes(sequence, frequency, fill = used)) +
	geom_col() +
	scale_y_continuous(label = comma) +
	xlab("index sequence") +
	ylab("reads assigned naively") +
	theme(
		legend.justification = c(1, 1),
		legend.position = c(1, 1),
		panel.grid.major.x = element_blank(),
		axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, family = "mono"),
	)

ggsave(out.file, result.plot, device = "pdf", width = 8, height = 6)

