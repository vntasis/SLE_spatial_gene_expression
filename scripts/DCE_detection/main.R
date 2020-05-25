########Main script of the COD detection pipeline########

# load required libraries
library(parallel)
library(features)

# source all the required functions and prepare the data
filenames <- c('Initial_filtering.R',
               'data_preparation.R',
               'calculate_correlations.R',
               'permutation_test_for_correlations.R',
               'compute_binsignal.R',
               'detect_DCEs.R',
               'merge_DCEs.R')
sapply(filenames, source, echo = T)

# Calculate correlation matrix
chrs <- unique(bins$chromosome)

cor_final_healthy <- calc_cor_final(bin_counts = bin_counts_healthy)
cor_final_da1 <- calc_cor_final(bin_counts = bin_counts_da1)
cor_final_da2 <- calc_cor_final(bin_counts = bin_counts_da2)
cor_final_da3 <- calc_cor_final(bin_counts = bin_counts_da3)

# Calculate permutation scores
system.time(
  permutation_scores_healthy <- mclapply(X = chrs, FUN = permutations, n_of_permutations = 1000, correlation = cor,
                                         bin_counts = bin_counts_healthy, cor_final = cor_final_healthy, mc.cores = 4)
)


system.time(
  permutation_scores_da1 <- mclapply(X = chrs, FUN = permutations, n_of_permutations = 1000, correlation = cor,
                                     bin_counts = bin_counts_da1, cor_final = cor_final_da1, mc.cores = 4)
)

system.time(
  permutation_scores_da2 <- mclapply(X = chrs, FUN = permutations, n_of_permutations = 1000, correlation = cor,
                                     bin_counts = bin_counts_da2, cor_final = cor_final_da2, mc.cores = 4)
)

system.time(
  permutation_scores_da3 <- mclapply(X = chrs, FUN = permutations, n_of_permutations = 1000, correlation = cor,
                                     bin_counts = bin_counts_da3, cor_final = cor_final_da3, mc.cores = 4)
)




# compute intra cod coexpression
compute_intra_cod_coexp <- function(cods, chr, thres = 0.05, group){
  bin_counts <- get(x = paste0("bin_counts_", group))
  map <- calc_cor_matrix(bins = bins, bin_counts = bin_counts, chr = chr, cor_method = "spearman", correlation = cor)
  rnames <- rownames(map)
  cnames <- colnames(map)
  if (thres < 1){
    permutation_scores <- get(x = paste0("permutation_scores_", group))
    x <- get(x = paste0("cor_final_", group))
    x <- x[which(x$chrom == chr),]
    x$cor[which(permutation_scores[[chr]] > thres)] <- 0
    map <- diag(1, nrow = nrow(map), ncol = ncol(map))
    map[upper.tri(map)] <- x$cor
    map <- t(map)
    map[upper.tri(map)] <- x$cor
    rownames(map) <- rnames
    colnames(map) <- cnames
  }
  average_intra_cod <- numeric()
  for (i in 1:nrow(cods)){
    cod <- map[cods[i,"START"]:cods[i,"END"], cods[i,"START"]:cods[i,"END"]]
    average_intra_cod[i] = mean(cod[upper.tri(cod)])
  }
  return(average_intra_cod)
}

# check intra cod coexpression using different window sizes to define cods
test_for_window <- function(windows, chr, group){
  intra_cod <- list()
  for (w in windows){
    print(w)
    cods <- codfinder(chr = chr, group = group, window.size = w, bin_signal_thres = 0.25,
                      perm_pvalue_thres = 0.05)
    cods <- cods$CODs
    cods <- cod_merge(cods, interfere = 2, group = group, chr = chr, permut_thres = 0.05)
    intra_cod[[as.character(w)]] <- compute_intra_cod_coexp(cods, chr, thres = 0.05, group)
  }
  return(intra_cod)
}

check_w <- test_for_window(3:20, "chr1", "healthy")
str(check_w)
sapply(check_w, function(x) length(x))
sapply(check_w, function(x) summary(x))


# Find CODs for all chrs and all groups

bin_signals <- list()
CODs <- list()
make_COD_list <- function(group, window.size, chrs){
  print(group)
  cods <- list()
  bin_signals <- get('bin_signals', envir = .GlobalEnv)
  bin_signals[[group]] <- list()
  for (chr in chrs){
    print(chr)
    c <- codfinder(chr, group, window.size, bin_signal_thres = 0.25,
                   perm_pvalue_thres = 0.05)
    cods[[chr]] <- c$CODs
    bin_signals[[group]][[chr]] <- c$bin_signal
  }
  assign('bin_signals', value = bin_signals, envir = .GlobalEnv)
  cods
}
CODs[['healthy']] <- make_COD_list('healthy', 3, chrs)

## Merge cods with a gap of at most two bins between them
merge_cods <- function(cod_list, group){
  for (i in 1:length(cod_list)){
    print(names(cod_list)[[i]])
    cod_list[[i]] <- cod_merge(cod_list[[i]], 2, group, names(cod_list)[[i]], 0.05)
  }
  cod_list
}
CODs$healthy <- merge_cods(CODs$healthy, 'healthy')

CODs[['da1']] <- make_COD_list('da1', 3, chrs)
CODs$da1 <- merge_cods(CODs$da1, 'da1')

CODs[['da2']] <- make_COD_list('da2', 3, chrs)
CODs$da2 <- merge_cods(CODs$da2, 'da2')

CODs[['da3']] <- make_COD_list('da3', 3, chrs)
CODs$da3 <- merge_cods(CODs$da3, 'da3')

## Find the cod coordinates
coCODs <- function(cods){
  c_CODs <- list()
  for (group in names(cods)){
    c_CODs[[group]] <- list()
    for (chr in names(cods[[1]])){
      c_CODs[[group]][[chr]] <- matrix(0, nrow(cods[[group]][[chr]]), ncol(cods[[group]][[chr]]), dimnames = list(NULL, c("START", "END")))
      bin_starts <- bin_starts <- as.numeric(rownames(bin_signals[[group]][[chr]]))

      c_CODs[[group]][[chr]][,"START"] <- bin_starts[cods[[group]][[chr]][,"START"]]
      c_CODs[[group]][[chr]][,"END"] <- bin_starts[cods[[group]][[chr]][,"END"]]+9999
    }
  }
  return(c_CODs)
} #find the cod coordinates
c_CODs <- coCODs(CODs)
