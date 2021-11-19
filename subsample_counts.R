library(abind)
library(parallel)
options(mc.cores = detectCores())

# subsample the counts in an integer vector to a certain proportion
subsample.count.vec <- function(count.vec, p) sapply(count.vec, function(count) rbinom(1, count, p))

# subsample the counts in an integer matrix to certain proportion(s)
subsample.count.mat <- function(count.mat, p.vec) {
	if (length(p.vec) == 1) p.vec <-rep(p.vec, ncol(count.mat)) else stopifnot(length(p.vec) == ncol(count.mat))
	result <- do.call(cbind, lapply(seq(length(p.vec)), function(i) subsample.count.vec(count.mat[,i], p.vec[i])))	
	dimnames(result) <- dimnames(count.mat)
	result
}

# create alternative count matrices subsampling all the columns to different target totals
# subsample proportion may be calculated as a fraction of column sum or of a user-supplied total (e.g. total sequenced reads including uncounted, unaligned)
subsample.increments <- function(count.mat, totals = colSums(count.mat), increment = 5E5) {
	stopifnot(length(totals) == ncol(count.mat))
	targets <- increment * seq(floor(max(totals) / increment))
	do.call(abind, c(mclapply(targets, function(target) subsample.count.mat(count.mat, target / totals)), along = 3, list(new.names = format(targets, scientific = F))))
}

