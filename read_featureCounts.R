# easily read counts from Subread featureCounts output
# omit gene position info, read only a matrix of integer counts
# optionally remove directory paths, but this can create duplicates

read.featureCounts <- function(file, remove.dir = TRUE) {
	raw <- read.table(file,
		header = T,
		check.names = F,
		stringsAsFactors = F,
		sep = "\t"
	)
	stopifnot("unrecognized featureCounts format" = identical(c("Geneid", "Chr", "Start", "End", "Strand", "Length"), colnames(raw)[1:6]))
	counts <- as.matrix(raw[,7:ncol(raw)])
	
	duplicate.genes <- unique(raw$Geneid[duplicated(raw$Geneid)])
	if (length(duplicate.genes) > 0) warning(paste(c("duplicate gene IDs:", duplicate.genes), collapse = " "))
	rownames(counts) <- raw$Geneid
	
	colnames(counts) <- sub("\\.bam$|\\.cram$", "", colnames(counts))
	if (remove.dir) colnames(counts) <- basename(colnames(counts))
	duplicate.samples <- unique(colnames(counts)[duplicated(colnames(counts))])
if (length(duplicate.samples) > 0) warning(paste(c("duplicate sample names:", duplicate.samples), collapse = " "))
	
	return(counts)
} 

