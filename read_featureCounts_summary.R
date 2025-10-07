# easily read category counts from Subread featureCounts summary file
# optionally remove directory paths from library names, but this can create duplicates

read.featureCounts.summary <- function(file, remove.dir = TRUE) {
	counts <- as.matrix(read.table(file,
		header = T,
		check.names = F,
		stringsAsFactors = F,
		row.names = 1,
		sep = "\t"
	))
	
	colnames(counts) <- sub("\\.bam$|\\.cram$", "", colnames(counts))
	if (remove.dir) colnames(counts) <- basename(colnames(counts))
	duplicate.samples <- unique(colnames(counts)[duplicated(colnames(counts))])
	if (length(duplicate.samples) > 0) warning(paste(c("duplicate sample names:", duplicate.samples), collapse = " "))
	
	return(counts)
} 

# easily calculate proportion of reads assigned from featureCounts category counts file
get.prop.assigned <- function(file, ...) {
	counts <- read.featureCounts.summary(file, ...)
	return(counts["Assigned",] / colSums(counts))
}

