cod_merge <- function(cods, interfere, group, chr, permut_thres ){
  
  # calculate cod gaps
  if (nrow(cods) > 1){
    dists <- sapply(1:(nrow(cods)-1), function(x) cods[x+1,1]-cods[x,2]-1)
  }else return(cods)
  
  # find regions with continuous gaps
  find_regions <- function(inds){
    regions <- list()
    if (length(inds) > 0) {
      regions[[1]] <- inds[1] 
      if (length(inds) > 1){
        p <- 1
        for (i in 2:length(inds)){
          if (inds[i] == inds[i-1]+1) regions[[p]] <- c(regions[[p]], inds[i]) 
          else{
            p <- p + 1
            regions[[p]] <- inds[i]
          }
        }
      } 
    }
    regions
  } 
  regions <- find_regions(which(dists <= interfere))
  
  # find map of correlations for the corresponding sample and group from the environment
  find_map <- function(group, chr, thres){
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
    map
  }
  map <- find_map(group, chr, permut_thres)
  
  if (length(regions) > 0){
    opt_merges <- list()
    
    # for each region find the optimum way to merge cods
    for (l in 1:length(regions)){
      region <- regions[[l]]
      interfere_index <- numeric()
      combs <- list()
      combs_score <- numeric()
      
      # calculate initial merging and score the different scenarios
      ds <- dists[region]
      region_cods <- cods[append(region, (region[length(region)] + 1)),]
      for (i in 2:nrow(region_cods)){
        x <- unlist(apply(combn(1:nrow(region_cods),i), 2, find_regions), recursive = F)
        if (length(x) > 0) x <- x[sapply(x, function(n) length(n) == i)]
        if (length(x) > 0) x <- x[sapply(x, function(k) sum(ds[k[-length(k)]])) <= interfere]
        if (length(x) > 0){
          combs <- c(combs, lapply(x, function(k) paste(k, collapse = '-')))
          interfere_index <- c(interfere_index, sapply(x, function(k) sum(ds[k[-length(k)]])))
          for (c in x){
            cors <- map[region_cods[min(c),'START']:region_cods[max(c), 'END'], 
                        region_cods[min(c),'START']:region_cods[max(c), 'END']]
            cors <- mean(cors[upper.tri(cors)])
            rest_cods <- region_cods[-c, , drop=F]
            if (nrow(rest_cods) > 0) for (j in 1:nrow(rest_cods)){
              cod_cor <- map[rest_cods[j,'START']:region_cods[j, 'END'], 
                             rest_cods[j,'START']:region_cods[j, 'END']]
              
              cors <- c(cors, mean(cod_cor[upper.tri(cod_cor)]))
            }
            combs_score <- c(combs_score, (mean(cors) + (length(c)-1)/(nrow(region_cods) -1)))
          }
        }
      }
      
      # combine the different merge events (in case the size of the region allows that) and score them
      repeat{
        new_combs <- list(0)
        for (i in 1:length(combs)){
          x <- unlist(strsplit(combs[[i]], split = '-'))
          ints <- sapply(combs, function(k) length(intersect(x, unlist(strsplit(k, split = '-')))))
          if (sum(ints == 0) > 0) to_pair <- combs[ints == 0] else next
          for  (j in 1:length(to_pair)){
            pair <- c(combs[[i]], to_pair[[j]])
            if (!(any(sapply(combs, function(k) identical(pair[order(pair)],k[order(k)])))) &&
                !(any(sapply(new_combs, function(k) identical(pair[order(pair)],k[order(k)]))))){
              new_combs[[length(new_combs) + 1]] <- pair
            }
          }
        }
        new_combs <- new_combs[-1]
        if (length(new_combs) > 0){
          combs <- c(combs, new_combs)
          for (new_comb in new_combs){
            cod_cors <- numeric()
            for (pr in new_comb){
              ids <- as.integer(unlist(strsplit(pr,'-')))
              rest_cods <- region_cods[-ids, , drop=F]
              cod_cor <- map[region_cods[min(ids),'START']:region_cods[max(ids), 'END'], 
                             region_cods[min(ids),'START']:region_cods[max(ids), 'END']]
              cod_cors <- c(cod_cors, mean(cod_cor[upper.tri(cod_cor)]))
            }
            
            if (nrow(rest_cods) > 0) for (j in 1:nrow(rest_cods)){
              cod_cor <- map[rest_cods[j,'START']:region_cods[j, 'END'], 
                             rest_cods[j,'START']:region_cods[j, 'END']]
              cod_cors <- c(cod_cors, mean(cod_cor[upper.tri(cod_cor)]))
            }
            combs_score <- c(combs_score, (mean(cod_cors) + (length(unlist(strsplit(new_comb, '-'))) - length(new_comb)) / (nrow(region_cods) -1)))
          }
        } else break
      }
      
      # choose the solution with the maximum score
      opt_comb <- which(combs_score == max(combs_score))
      opt_merges[[l]] <- combs[[opt_comb]]
    }
    
    # Built new cod matrix
    new_cods <- matrix(0, 1, 2, dimnames = list(NULL, c("START", "END")))
    cods_in_regions <- unlist(sapply(regions, function(k) append(k, (k[length(k)] + 1))))
    # sapply(cods_in_regions, 
    #        function(k) which(sapply(regions, function(region) (k %in% region || (k-1) %in% region))))
    cods_not_in_regions <- setdiff(1:nrow(cods), cods_in_regions)
    
    for (idx in cods_not_in_regions) new_cods <- rbind(new_cods, cods[idx, ])
    for (i in 1:length(regions)){
      region <- (function(k) append(k, (k[length(k)] + 1)))(regions[[i]])
      for (j in 1:length(opt_merges[[i]])){
        idxs_to_merge <- as.integer(unlist(strsplit(opt_merges[[i]][j], '-')))
        new_cods <- rbind(new_cods, c(cods[region[min(idxs_to_merge)], "START"], 
                                      cods[region[max(idxs_to_merge)], "END"]))
      }
      not_merged <- setdiff(1:length(region), 
                            as.integer(unlist(strsplit(opt_merges[[i]], '-'))))
      if (length(not_merged) > 0) for (m in not_merged){
        new_cods <- rbind(new_cods, cods[region[m], ])
      }
    }
    new_cods <- new_cods[-1, , drop=F]
    new_cods[,1] <- new_cods[order(new_cods[,1], decreasing = F),1]
    new_cods[,2] <- new_cods[order(new_cods[,2], decreasing = F),2]
  } else return(cods)
  new_cods
}
