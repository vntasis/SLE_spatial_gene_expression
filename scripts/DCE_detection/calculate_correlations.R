# Calclulates the correlation matrix with p-values
cor.prob <- function (X, dfr = nrow(X) - 2, method = "spearman") {
  R <- cor(X, method = method)
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  return(R)
}

# Transforms the correlation matrix calculated from the cor.prob function
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = as.numeric(rownames(m)[row(m)[ut]]),
             j = as.numeric(rownames(m)[col(m)[ut]]),
             cor=t(m)[ut],
             p=m[ut])
}


# This function requires the bins data.frame produced by the calc_bins ('split_in_bins.R'), 
# the bin_counts matrix computed by the find_bin_counts function ('split_in_bins.R'), 
# the chromosome whose correlation matrix you want to calculate (corresponds to the chromosome field of bins data.frame), 
# the method to be used for the correlation coefficient calculation (default is Spearman) 
# and the function to be used for this calculation (cor or cor.prob).
# It outputs a correlation matrix for the bins (with non zero expression) in the corresponding chromosome. 
# Each bin is named after its start coordinate.

calc_cor_matrix <- function(bins, bin_counts, chr, cor_method = "spearman", correlation){
  z <- rowSums(bin_counts)
  cor_matrix <- correlation(t(bin_counts[which(bins$chromosome == chr & z > 0),]), method = cor_method)
  rownames(cor_matrix) <- as.character(bins$start[which(bins$chromosome == chr & z > 0)])
  colnames(cor_matrix) <- as.character(bins$start[which(bins$chromosome == chr & z > 0)])
  return(cor_matrix)
}


# The next function has the same requirements as the preceding one.
# Also, the chrs vector is a character vector containing the names of chromosomes under study. 
# The cor.prob function, the flattenSquareMatrix function and the calc_cor_matrix function need to be in the global environment.
# It calculates the correlation matrices for all chromosomes, flattens them and rbinds them in a single data.frame. It outputs this data.frame.

calc_cor_final <- function(bin_counts){
  final <- data.frame(chrom = character(), i = numeric(), j = numeric(), cor = numeric(), p = numeric())
  for (chr in chrs){
    x <- calc_cor_matrix(bins = bins, bin_counts = bin_counts, chr = chr, cor_method = "spearman", correlation = cor.prob)
    x <- flattenSquareMatrix(x)
    chrom <- rep(chr, nrow(x))
    x <- cbind(chrom, x)
    x$chrom <- as.character(x$chrom)
    final <- rbind(final, x)
    rm(x)
  }
  return(final)
}

