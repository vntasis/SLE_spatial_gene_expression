########Number of CODs per chromosome and per group########
noCpch <- function(cods){
  n_cods_per_chr <- list()
  for (group in names(cods)){
    n_cods_per_chr[[group]] <- sapply(cods[[group]], nrow)
  }
  return(n_cods_per_chr)
}

n_cods_per_chr <- noCpch(CODs)
n_cods_per_chr
sapply(n_cods_per_chr, sum)

########Number of bins found inside a COD########
nob <- function(cods){
  n_bins <- list()
  for (group in names(cods)){
    for (chr in names(cods[[names(cods)[1]]])){
      x <- cods[[group]][[chr]]
      n <- x[,"END"]-x[,"START"]+1
      n_bins[[group]][[chr]] <- unname(n)
    }
  }
  return(n_bins)
}
n_bins <- nob(CODs)
sapply(n_bins, function(x) summary(unlist(x , use.names = F)))


########Percentage of bins inside CODs########
no_non_zero_bins <- function(groups, chrs){
  n_bins <- list()
  for (group in groups){
    bin_counts <- get(x = paste0("bin_counts_", group))
    for (chr in chrs){
      z <- rowSums(bin_counts[which(bins$chromosome == chr),])
      n_bins[[group]][[chr]] <- sum(z > 0)
    }
  }
  return(n_bins)
}
no_nz_bins <- no_non_zero_bins(c("healthy", "da1", "da2", "da3"), chrs = chrs)
percentages <- function(groups, chrs, n_bins, no_nz_bins){
  perc <- matrix(0,1,length(chrs))
  for (group in groups){
    perc <- rbind(perc, sapply(chrs, function(x) (sum(n_bins[[group]][[x]])/no_nz_bins[[group]][x])*100))
  }
  perc <- perc[-1,]
  dimnames(perc)<- list(groups, chrs)
  return(perc)
}
perc <- percentages(groups = c("healthy", "da1", "da2", "da3"), chrs, n_bins, no_nz_bins)


########Size of CODs########
soCODs <- function(cods){
  s_CODs <- list()
  for (group in names(cods)){
    for (chr in chrs){
      bin_starts <- as.numeric(rownames(bin_signals[[group]][[chr]]))
      sizes <- bin_starts[cods[[group]][[chr]][,"END"]]-bin_starts[cods[[group]][[chr]][,"START"]]+10000
      s_CODs[[group]][[chr]] <- sizes
    }
  }
  return(s_CODs)
}
s_CODs <- soCODs(CODs)
sapply(s_CODs, function(x)summary(unlist(x, use.names = F)))
s_CODs <- lapply(s_CODs, unlist, use.names = F)
s_CODs <- lapply(s_CODs, function(x)x/1000000)
s_CODs1 <- data.frame(size =c(s_CODs$healthy, s_CODs$da1, s_CODs$da2, s_CODs$da3), 
                      group = c(rep("Healthy", length(s_CODs$healthy)), rep("DA1", length(s_CODs$da1)), 
                                rep("DA2", length(s_CODs$da2)), rep("DA3", length(s_CODs$da3))))

########Find intra cod coexpression########

compute_intra_cod_coexp <- function(cods, thres = 0.05){
  intra_cod <- list()
  for (group in names(cods)){
    intra_cod[[group]] <- list()
    print(group)
    
    xfeth <- feather::feather(paste0("RData/correlations_", group,".feather"))
    for (chr in names(cods[[1]])){
      print(chr)
      bin_counts <- get(x = paste0("bin_counts_", group))
      z <- rowSums(bin_counts)
      n <- sum(bins$chromosome == chr & z > 0)
      nam <- as.character(bins$start[which(bins$chromosome == chr & z > 0)])
      rm(bin_counts, z)
      
      
      x <- xfeth[xfeth$chrom == chr,]
      
      if (thres < 1){
        permutation_scores <- get(x = paste0("permutation_scores_", group))
        x$cor[permutation_scores[[chr]] > thres] <- 0
        rm(permutation_scores)
        gc()
      } 
      map <- diag(1, nrow = n, ncol = n)
      map[upper.tri(map)] <- x$cor
      map <- t(map)
      map[upper.tri(map)] <- x$cor
      rownames(map) <- nam
      colnames(map) <- nam
      rm(n, x)
      gc()
      codss <- cods[[group]][[chr]]
      
      intra_cod[[group]][[chr]] <- numeric(length(nrow(codss)))
      if (nrow(codss) >= 1){
        for (i in 1:nrow(codss)){
          cod <- map[codss[i,"START"]:codss[i,"END"], codss[i,"START"]:codss[i,"END"]]
          intra_cod[[group]][[chr]][i] <- mean(cod[upper.tri(cod)])
        }
      }
    }
  }
  return(intra_cod)
}

intra_cod <- compute_intra_cod_coexp(CODs)
sapply(intra_cod, function(x)summary(unlist(x, use.names = F)))

########Find inter cod coexpression########

compute_inter_cod_coexp <- function(cods, thres = 0.05){
  inter_cod <- list()
  for (group in names(cods)){
    inter_cod[[group]] <- list()
    print(group)
    
    xfeth <- feather::feather(paste0("RData/correlations_", group,".feather"))
    for (chr in names(cods[[1]])){
      print(chr)
      bin_counts <- get(x = paste0("bin_counts_", group))
      z <- rowSums(bin_counts)
      n <- sum(bins$chromosome == chr & z > 0)
      nam <- as.character(bins$start[which(bins$chromosome == chr & z > 0)])
      rm(bin_counts, z)
      
      
      x <- xfeth[xfeth$chrom == chr,]
      
      if (thres < 1){
        permutation_scores <- get(x = paste0("permutation_scores_", group))
        x$cor[permutation_scores[[chr]] > thres] <- 0
        rm(permutation_scores)
        gc()
      } 
      map <- diag(1, nrow = n, ncol = n)
      map[upper.tri(map)] <- x$cor
      map <- t(map)
      map[upper.tri(map)] <- x$cor
      rownames(map) <- nam
      colnames(map) <- nam
      rm(n, x)
      gc()
      codss <- cods[[group]][[chr]]
      
      inter_cod[[group]][[chr]] <- numeric(length(nrow(codss)))
      if (nrow(codss) > 1){
        combinations <- combn(1:nrow(codss), 2)
        inter_cod[[group]][[chr]] <- numeric(length(ncol(combinations)))
        for (i in 1:ncol(combinations)){
          cod <- map[codss[combinations[1,i],"START"]:codss[combinations[1,i],"END"], 
                     codss[combinations[2,i],"START"]:codss[combinations[2,i],"END"]]
          inter_cod[[group]][[chr]][i] <- mean(cod)
        }
      }
    }
  }
  return(inter_cod)
}

inter_cod <- compute_inter_cod_coexp(CODs)
sapply(inter_cod, function(x)summary(unlist(x, use.names = F)))
