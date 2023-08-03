# easily read counts from Subread featureCounts output
# omit gene position info, read only a matrix of integer counts
# optionally convert BAM/CRAM file paths to basenames, but this can create duplicates

read.featureCounts <- function(file, simplify = TRUE) {
	raw <- read.table(file,
		header = T,
		row.names = 1,
		check.names = F,
		stringsAsFactors = F,
		sep = "\t"
	)
	stopifnot("unrecognized column headers" = all(colnames(raw)[1:5] == c("Chr", "Start", "End", "Strand", "Length")))
	result <- as.matrix(raw[,6:ncol(raw)])
	if (simplify) colnames(result) <- sub("\\.bam$|\\.cram$", "", basename(colnames(result)))
	if (any(duplicated(colnames(result)))) warning("duplicate colnames")
	result
} 

