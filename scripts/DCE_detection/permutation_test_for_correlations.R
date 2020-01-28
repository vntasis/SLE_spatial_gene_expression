# This function requires the bin_counts matrix computed by the find_bin_counts function (split_in_bins.R), 
# the cor_final data.frame computed by the calc_cor_final function (calculate_correlations.R),
# the correlation function to be used, the number of permutations to be implemented and the chromosome to be permuted.
# It returns a vector with the permutation p-values of the correlation coefficients of the corresponding chromosome.

permutations <- function(chr, n_of_permutations, correlation, bin_counts, cor_final){
  score <- rep(0,length(cor_final$chrom[which(cor_final$chrom == chr)]))
  
  shuffle <- function(k, chr, correlation, score){
    print(paste0("current number of permutation : ", k))
    x <- cor_final[which(cor_final$chrom == chr),]
    z <- rowSums(bin_counts)
    bin_counts <- bin_counts[which(bins$chromosome == chr & z > 0),]
    bin_counts <- apply(X = bin_counts, MARGIN = 2, FUN = sample)
    
    cor_matrix <- correlation(t(bin_counts), method = "spearman")
    rownames(cor_matrix) <- as.character(bins$start[which(bins$chromosome == chr & z > 0)])
    colnames(cor_matrix) <- as.character(bins$start[which(bins$chromosome == chr & z > 0)])
    cor_matrix <- flattenSquareMatrix(cor_matrix)
    
    pos_indices <- which(x[,"cor"] > 0)
    indices <- which(cor_matrix[pos_indices,"cor"] >= x[pos_indices,"cor"])
    score[pos_indices[indices]] <- score[pos_indices[indices]] + 1
    
    neg_indices <- which(x[,"cor"] < 0)
    indices <- which(cor_matrix[neg_indices,"cor"] <= x[neg_indices,"cor"])
    score[neg_indices[indices]] <- score[neg_indices[indices]] + 1
    return(score)
  }
  
  for(i in 1:n_of_permutations){
    score <- shuffle(k = i, chr = chr, correlation = correlation, score = score)
  }
  score <- score/n_of_permutations
  return(score)
}


# Use parallel computing to calculate the permutation scores for every chromosome (requires the parallel package). It returns a list.
# Here is an example for the healthy group.
# system.time(
#   permutation_scores_healthy <- mclapply(X = chrs, FUN = permutations, n_of_permutations = 1000, correlation = cor, 
#                                          bin_counts = bin_counts_healthy, cor_final = cor_final_healthy, mc.cores = 4)
# )

