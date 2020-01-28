#######Detection of Disruptors - Genes located in the gap between two DCEs that were united in the Healthy group########
find_genes <- function(chr, start, end, coordinates, DA_group){
  
  data <- data_normalized[,which(DA == DA_group)]
  z <- rowSums(data)
  
  genes <- intersect(rownames(coordinates)[coordinates[, "start"] >= start &  coordinates[, "start"] <= end & coordinates[, "chromosome"] == chr],
                     rownames(data)[z > 0])
  coordinates[genes,]
}


find_genes_of_split_cods <- function(overlap_cutoff=0){
  
  genes <- list()
  for (group in names(rearrangements)){
    
    print(group)
    gene_name <- character()
    cod_control <- integer()
    cod_test <- integer()
    start <- integer()
    end <- integer()
    chromosome <- character()
    
    for (chr in names(rearrangements[[group]])){
      
      cod_idxs_control <- unique(rearrangements[[group]][[chr]]$split[, "control"])
      if (length(cod_idxs_control) > 0){
        
        print(chr)
        
        for (j in cod_idxs_control){
          
          cod_idxs_test <- rearrangements[[group]][[chr]]$split[
            rearrangements[[group]][[chr]]$split[,"control"] == j ,]
          if (all(cod_idxs_test[, "overlap"] < overlap_cutoff)) next
          cod_idxs_test <- cod_idxs_test[, "test"]
          
          split_cods_control <- c_CODs[["healthy"]][[chr]][j, , drop = F]
          split_cods_test <- c_CODs[[group]][[chr]][cod_idxs_test, , drop = F]
          
          genes_control <- find_genes(chr,
                                      split_cods_control[,1],
                                      split_cods_control[,2],
                                      coordinates$coordinates,
                                      "Healthy")
          
          
          test_cod <- rep(0, nrow(genes_control))
          
          for (i in 1:nrow(split_cods_test)){
            
            gene_test <- rownames(find_genes(chr, 
                                     split_cods_test[i,1],
                                     split_cods_test[i,2], 
                                     coordinates$coordinates,
                                     toupper(group)))
            common <- intersect(rownames(genes_control), gene_test)
            test_cod[match(common, rownames(genes_control))] <- cod_idxs_test[i]
          }
          
          gene_name <- c(gene_name, rownames(genes_control))
          cod_control <- c(cod_control, rep(j, nrow(genes_control)))
          cod_test <- c(cod_test, test_cod)
          chromosome <- c(chromosome, genes_control[,"chromosome"])
          start <- c(start, genes_control[,"start"])
          end <- c(end, genes_control[,"end"])
        }
      }
    }
    if (!(length(gene_name) > 0)) next
    genes[[group]] <- data.frame("gene_name" = gene_name, 
                                 "cod_control" = cod_control,
                                 "cod_test" = cod_test,
                                 "chromosome" = chromosome,
                                 "start" = start,
                                 "end" = end, stringsAsFactors = F)
  }
  genes
}
genes_of_split_cods <- find_genes_of_split_cods(0.1)



trim <- function(l){
  repeat if (l[[1]][1] == 0) l <- l[-1,] else break
  l
}

disruptors <- 
  genes_of_split_cods %>% 
  lapply(FUN = function(x){
  unique(x$chromosome) %>% 
    lapply(FUN = function(k) {x %>% filter(chromosome == k) %$% 
        unique(cod_control) %>%
        lapply(FUN = function(t){
          x %>% filter(chromosome == k) %>%
            filter(.$cod_control == t) %$%
            data.frame(cod_test, gene_name, stringsAsFactors = F) %>% trim %>% 
            .[rev(rownames(.)),] %>% trim %>% 
            .[rev(rownames(.)),] %$% gene_name[cod_test == 0]
        })}) %>% unlist
})
