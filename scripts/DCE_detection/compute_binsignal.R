binsignal <- function(chr, thres = 0.05, group, window.size){
  
  bin_counts <- get(x = paste0("bin_counts_", group))
  z <- rowSums(bin_counts)
  permutation_scores <- get(x = paste0("permutation_scores_", group))
  x <- get(x = paste0("cor_final_", group))
  x <- x[which(x$chrom == chr),]
  x$cor[which(permutation_scores[[chr]] > thres)] <- 0
  map <- diag(1, nrow = nrow(bins[which(bins$chromosome == chr & z > 0),]), ncol = nrow(bins[which(bins$chromosome == chr & z > 0),]))
  map[upper.tri(map)] <- x$cor
  map <- t(map)
  map[upper.tri(map)] <- x$cor
  
  rank_test <- function(x,y, ...){
    tryCatch(
      wilcox.test(x,y, ...),
      error = function(e){
        return(list(p.value = NaN))
      }
    )
  }
  
  n_bins <- nrow(map)
  mean.cf <- rep(0, n_bins)
  pvalue_left <- rep(1, n_bins)
  pvalue_right <- rep(1, n_bins)
 
  ##Diamond function from TOPDOM for creating binsignal
  for(i in 1:(n_bins-1))
  {
    lowerbound = max(1, i-window.size+1)
    upperbound = min(i+window.size, n_bins)
    diamond = map[lowerbound:i, (i+1):upperbound]
    mean.cf[i] = mean(diamond)

    pvalue_right[i] <- rank_test(diamond, map[lowerbound:i, lowerbound:i], 'less')$p.value
    pvalue_left[i] <- rank_test(diamond, map[(i+1):upperbound, (i+1):upperbound], 'less')$p.value
  }
  
  for (i in 1:(n_bins-1))
  {
    if(!is.nan(pvalue_right[i])) pvalue_right[i]=pvalue_right[i]
    else pvalue_right[i]=1

    if(!is.nan(pvalue_left[i])) pvalue_left[i]=pvalue_left[i]
    else pvalue_left[i]=1
  }
  
  ##Binsignal table
  bin.signal <- matrix(data = c(mean.cf, pvalue_left, pvalue_right),nrow=n_bins,ncol=3)
  rownames(bin.signal) <- as.character(bins$start[which(bins$chromosome == chr & z > 0)])
  colnames(bin.signal) <- c("binsignal", "p-value_L", "p-value_R")
  return(bin.signal)
}

