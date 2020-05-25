library(tidyverse)
library(magrittr)

#import data
raw_counts <- read.table(file.choose(), header = T, sep = "\t", stringsAsFactors = F)
rownames(raw_counts) <- raw_counts[,1]
raw_counts <- raw_counts[,-1]

#import annotation file
annotation_data  <- read.table(file = file.choose(), header = F, sep = "\t", stringsAsFactors = F)
head(annotation_data)
colnames(annotation_data) <- c("Chromosome", "Start", "End", "gene_id", "transcript_id", "gene_type", "gene_name", "strand")

table(annotation_data$gene_type)
#gene_names_occurencies <- table(annotation_data$gene_name)

#set excluded gene type
gene_types <- as.character(unique(annotation_data$gene_type))
gene_types
excluded_types <- gene_types[c(1,5,7,10,12,13,18,22,23,25,26,27,29,30)]

#discard genes of excluded type
excluded_genes <- character()

for (excluded_type in excluded_types){
  excluded_genes <-
    annotation_data %>%
    filter(gene_type == excluded_type) %$%
    gene_name %>%
    append(x = excluded_genes)
}
excluded_genes %<>% unique %>% intersect(rownames(raw_counts))
length(excluded_genes)

raw_counts_filtered <- raw_counts[!(rownames(raw_counts) %in% excluded_genes),]



#discard genes with entries of multiple chromosomes or multiple strands
annotation_data_filtered <- annotation_data %>% filter(!gene_type %in% excluded_types)


#find the chromosomal coordinates of the genes from the annotation file
#This function reuires a matrix with the gene_names of interest as its rownames (data)
#and an annotation data.frame containing gene_names and gene_coordinates
#It gives a list as an output, containing the coordinates of the coresponding genes (in coordinates)
#and the names of the genes (included in data) having multiple entries in different chromosomes  (in error)
find_coordinates <- function(data, annotation){
  g <-character()
  chromosome <- character()
  start <- numeric()
  end <- numeric()
  strand <- character()
  error <- character()
  for (gene in rownames(data)){
    x <- annotation[annotation$gene_name == gene,]
    chr <- as.character(unique(x$Chromosome))
    stran <- unique(x$strand)
    if (length(chr) == 1 && length(stran) == 1){
      chromosome <- append(x = chromosome, values = chr)
      g <- append(x = g, values = gene)
      start <- append(x = start, values = min(x$Start))
      end <- append(x = end, values = max(x$End))
      strand <- append(x = strand, values = stran)
    }else{
      error <- append(x = error, values = gene)
      #print(paste0("Error: Gene ", gene, " is assigned to more than one chromosomes"))
    }
  }
  return(list(error = error, coordinates = data.frame(chromosome = chromosome, start = start, end = end, strand = strand, row.names  = g, stringsAsFactors = F)))
}


coordinates <- find_coordinates(raw_counts_filtered, annotation_data_filtered)
length(coordinates$error)
raw_counts_filtered <- raw_counts_filtered[!(rownames(raw_counts_filtered) %in% coordinates$error),]
dim(raw_counts_filtered)

#discard mitochondrial genes
Mt_genes <-
  coordinates$coordinates %>%
  rownames_to_column('gene') %>%
  filter(chromosome == "chrM") %$%
  gene

raw_counts_filtered <- raw_counts_filtered[!(rownames(raw_counts_filtered) %in% Mt_genes),]


# Normalization for library size
library(MDSeq)
data_normalized <- normalize.counts(counts = raw_counts_filtered, method = "RLE", verbose = T)

coordinates <- find_coordinates(raw_counts_filtered, annotation_data_filtered)
# Normalization for gene length
lengths <- coordinates$coordinates$end - coordinates$coordinates$start + 1
lengths %<>% `/`(1000)
for (i in 1:nrow(data_normalized)){
  print(i)
  data_normalized[i,] <- data_normalized[i,]/lengths[i]
}
