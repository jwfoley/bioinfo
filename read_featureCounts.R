# easily read counts from Subread featureCounts output
# omit gene position info, read only a matrix of integer counts
# convert BAM/CRAM filenames to basenames

read.featureCounts <- function(file) {
	raw <- read.table(file,
		header = T,
		row.names = 1,
		check.names = F,
		stringsAsFactors = F,
		sep = "\t"
	)
	stopifnot("unrecognized column headers" = all(colnames(raw)[1:5] == c("Chr", "Start", "End", "Strand", "Length")))
	result <- as.matrix(raw[,6:ncol(raw)])
	colnames(result) <- sub("\\.bam$|\\.cram$", "", basename(colnames(result)))
	result
} 

