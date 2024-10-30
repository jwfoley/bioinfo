#! /usr/bin/Rscript

# known limitation: this requires that all lanes are combined together (as if "--no-lane-splitting true") and therefore have identical SampleID lists
# known limitation: this requires a specific order of column names even though several aren't used and could be absent (e.g. Sample_Project)

library(ggplot2)
library(scales)
library(optparse)

DEFAULT.OUTPUT <- "Demultiplex_Stats.pdf"
DEFAULT.DIMS <- "20,36"
SEQUENCE.COLORS <- c(
	"perfect" = "blue3",
	"1 mismatch" = "yellow3",
	"2 mismatches" = "red3"
)

opt <- parse_args(OptionParser(option_list = list(
	make_option(c("-d", "--dim"), action = "store", default = DEFAULT.DIMS, help = paste0("dimensions in cm of output graphs in form W,H (default: ", DEFAULT.DIMS, ")")),
	make_option(c("-k", "--keep_order"), action = "store_true", default = FALSE, dest = "keep.order", help = "keep samples in listed order (default: sort by total reads)"),
	make_option(c("-o", "--output"), action = "store", default = DEFAULT.OUTPUT, help = paste0("path of output PDF file (default: ", DEFAULT.OUTPUT, ")"))
)), positional_arguments = 1)
plot.dims <- as.numeric(strsplit(opt$options$dim, ",")[[1]])

stats.raw <- read.csv(opt$arg, check.names = F)
stopifnot(colnames(stats.raw)[1:8] == c("Lane", "SampleID", "Sample_Project", "Index", "# Reads", "# Perfect Index Reads", "# One Mismatch Index Reads", "# Two Mismatch Index Reads"))
metadata <- subset(stats.raw, Lane == 1)[,2:3]
counts.pooled <- as.matrix(subset(stats.raw, Lane == 1)[,5:8])
for (lane in unique(stats.raw$Lane)[1]) {
	stopifnot(identical(metadata, subset(stats.raw, Lane == lane)[,2:3]))
	counts.pooled <- counts.pooled + as.matrix(subset(stats.raw, Lane == lane)[,5:8])
}
stopifnot(counts.pooled[,1] == rowSums(counts.pooled[,2:4]))
stats.pooled <- data.frame(
	SampleID = factor(rep(metadata$SampleID, each = 3), levels = (if (opt$options$keep.order) metadata$SampleID else metadata$SampleID[order(counts.pooled[,1], decreasing = T)])),
	index.sequence = factor(c("perfect", "1 mismatch", "2 mismatches"), levels = c("perfect", "1 mismatch", "2 mismatches")),
	reads = as.vector(t(counts.pooled[,2:4]))
)

plot.stats <- ggplot(stats.pooled, aes(SampleID, reads, fill = index.sequence)) +
	geom_col(position = "stack", width = 1) +
	scale_x_discrete(limits = rev) +
	scale_y_continuous(label = comma, expand = c(0, 0), position = "right", sec.axis = sec_axis(~ . / sum(stats.pooled$reads), labels = percent)) +
	scale_fill_manual(values = SEQUENCE.COLORS) +
	coord_flip() +
	theme( # note x and y are reversed in here
		panel.grid.major.y = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.text.y = element_text(hjust = 0),
		legend.position = "inside",
		legend.justification = c("right", "bottom")
	)

ggsave(opt$options$output, plot.stats, width = plot.dims[1], height = plot.dims[2], units = "cm")

