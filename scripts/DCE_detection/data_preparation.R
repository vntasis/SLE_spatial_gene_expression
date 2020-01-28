#Split each chromosome in bins and calculate their coordinates
#This function requires the coordinates data.frame produced by the previous function and a bin_zise,
#splits the genome in equally sized bins and outputs a data.frame containing the coordinates of the bins 
calc_bins <- function(coordinates, bin_size = 10000){
  bins <- data.frame(chromosome = character(), start = numeric(), end = numeric(), stringsAsFactors = F)
  for (chr in unique(coordinates$chromosome)){
    x <- coordinates[which(coordinates$chromosome == chr),]
    start <- min(x$start)
    end <- max(x$end)
    if (((end-start+1) %% bin_size) == 0){
      bin_start <- seq(from = start, to = end, by = bin_size)
      bin_end <- bin_start -1
      bin_start <- bin_start[-length(bin_start)]
      bin_end <- bin_end[-1]
      bin_end[length(bin_end)] <- bin_end[length(bin_end)]+1
    }else{
      bin_start <- c(seq(from = start, to = end, by = bin_size), end)
      bin_end <- bin_start -1
      bin_start <- bin_start[-length(bin_start)]
      bin_end <- bin_end[-1]
      bin_end[length(bin_end)] <- bin_end[length(bin_end)]+1
    }
    t <- data.frame(chromosome = rep(chr, length(bin_start)), start = bin_start, end = bin_end, stringsAsFactors = F)
    bins <- rbind(bins, t)
  }
  return(bins)
}

bins <- calc_bins(coordinates = coordinates$coordinates)

#make DA groups
sledai <- assist_info$Clinical.SLEDAI
names(sledai) <- assist_info$Sample.ID
DA <- character()
DA[which(is.na(sledai))] <- "Healthy"
DA[which(sledai < 3)] <- "DA1"
DA[which(sledai < 9 & sledai >= 3)] <- "DA2"
DA[which(sledai >= 9)] <- "DA3"
names(DA) <- names(sledai)
DA <- as.factor(DA)
summary(DA)


# Assign counts to the bins. DA is factor that includes the patient disease activity (DA) class
data_healthy <- data_normalized[,which(DA == "Healthy")]
data_da1 <- data_normalized[,which(DA == "DA1")]
data_da2 <- data_normalized[,which(DA == "DA2")]
data_da3 <- data_normalized[,which(DA == "DA3")]

#reverse coordinates for genes of '-' strand
reverse_coordinates <- function(coordinates){
  for (i in 1:nrow(coordinates)){
    start <- coordinates[i, 'start']
    end <- coordinates[i, 'end']
    if (coordinates[i, "strand"] == '-'){
      coordinates[i, 'start'] <- end
      coordinates[i, 'end'] <- start
    }
  }
  coordinates
}

coordinates$coordinates <- reverse_coordinates(coordinates$coordinates)

library(parallel)
#This function requires the bins data.frame produced by the calc_bins, the coordinates data.frame produced by the find_coordinates,
#the data matrix, containing the normalized counts and the gene_names as its rownames,
#the index of the row (refering the bins data.frame) of the last bin needed to calculate its counts and the number of cores to be utilized (it implements parallel computing - requires the parallel package).
#It outputs a matrix with the bin_counts that correspond to the bins in the bins data.frame
find_bin_counts <- function(bins, data, gene_coordinates, lastbin = nrow(bins), cores = 3){
  
  
  assign_counts <- function(bin, bins, data, gene_coordinates){
    counts <- rep(0, ncol(data))
    x <- gene_coordinates[which(gene_coordinates$chromosome == bins[bin,"chromosome"] 
                                & gene_coordinates$start >= bins[bin, "start"] 
                                & gene_coordinates$start <= bins[bin, "end"]),]
    if (nrow(x) > 0){
      t <- data[rownames(x), , drop = F]
      counts <- colMeans(t)
    }
    return(counts)
  }
  bin_counts <- mclapply(X = 1:lastbin, FUN = assign_counts, mc.cores = cores, bins = bins, data = data, gene_coordinates = gene_coordinates)
  bin_counts <- matrix(data = unlist(bin_counts, use.names = F), byrow = T, nrow = length(bin_counts), ncol = ncol(data))
  colnames(bin_counts) <- colnames(data)
  return(bin_counts)
}

bin_counts_healthy <- find_bin_counts(bins = bins, data = data_healthy, gene_coordinates = coordinates$coordinates, lastbin = nrow(bins))
bin_counts_da1 <- find_bin_counts(bins = bins, data = data_da1, gene_coordinates = coordinates$coordinates, lastbin = nrow(bins))
bin_counts_da2 <- find_bin_counts(bins = bins, data = data_da2, gene_coordinates = coordinates$coordinates, lastbin = nrow(bins))
bin_counts_da3 <- find_bin_counts(bins = bins, data = data_da3, gene_coordinates = coordinates$coordinates, lastbin = nrow(bins))
