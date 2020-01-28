########Enhancer Analysis########

#load data - An example with CD4 data
eg_association <- read.table(file = "RData/external_data/enhancer_atlas/cd4", 
                             header = F, sep = "\t", quote = "", stringsAsFactors = F)


#tidy data
eg_association %<>% 
  rename(chrom_Enh = V1, 
         chromStart = V2, 
         chromEnd = V3,
         Gene_ID = V4,
         chrom_Gene = V5,
         TSS = V6,
         Transcript = V7,
         signalValue = V8,
         EP_Score = V9) %>%
  arrange(chrom_Enh)


#change transcript id to gene name
mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                path="/biomart/martservice",dataset="hsapiens_gene_ensembl")

id2name <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcription_start_site", 
                                "chromosome_name", "start_position", "end_position",
                                "external_gene_name", "transcript_gencode_basic"), mart = mart)

rownames(id2name) <- id2name$ensembl_transcript_id

eg_association %<>%
  mutate(gene_name = id2name[eg_association$Transcript,"external_gene_name"])


#check enhancers and tss for involvment in DCEs
f_cods <- function(chr, coor, group, type){
  result <- 0
  cods <- c_CODs[[group]][[chr]]
  if (type == "enh"){
    res <- which(cods[,1] <= coor[1] & cods[,2] >= coor[2])
    if (length(res) > 0) result <- res
  }else if (type == "tss"){
    res <- which(cods[,1] <= coor & cods[,2] >= coor)
    if (length(res) > 0) result <- res
  }
  result
}


#--------------------------------
## Repeat this for patient groups
#--------------------------------
eg_association <- 
  eg_association %$%
  sapply(1:length(chrom_Enh), function(x) f_cods(chr = chrom_Enh[x],
                                       coor = c(chromStart[x], chromEnd[x]), 
                                       group = "healthy", 
                                       type = "enh")) %>%
  mutate(.data = eg_association, cod_healthy_enh = .)


eg_association <- 
  eg_association %$%
  sapply(1:length(chrom_Gene), function(x) f_cods(chr = chrom_Gene[x],
                                                 coor = TSS[x], 
                                                 group = "healthy", 
                                                 type = "tss")) %>%
  mutate(.data = eg_association, cod_healthy_tss = .)
#--------------------------------


#discard genes not expressed in the dataset
groups <- c("Healthy", "DA1", "DA2", "DA3")
for (i in 1:length(groups)){
  data <- data_normalized[,which(DA == groups[i])]
  z <- rowSums(data)
  
  k <- c(11,13,15,17)
  for (j in 1:nrow(eg_association)){
    if (!(is.na(z[eg_association$gene_name[j]])) && z[eg_association$gene_name[j]] > 0) next()
    eg_association[j,k[i]:(k[i]+1)] <- 0
  }
}


#--------------------------------
## Repeat this for patient groups
#--------------------------------

#check enhancer tss interactions - they are accumulated inside the same DCE

eg_association <-
  eg_association %$%
  sapply(1:length(chrom_Enh), function(x) identical(cod_healthy_enh[x], cod_healthy_tss[x])) %>%
  replace(eg_association$cod_healthy_enh == 0, FALSE) %>%
  mutate(.data = eg_association, cod_healthy_interaction = .)

apply(eg_association[,19:22], 2, summary)

#--------------------------------


#------------------------------------------------
## Repeat this for the rest of the patient groups 
#------------------------------------------------

#Isolate depleted and emerged interactions for each patient group

depleted_interactions <- list()
depleted_interactions[["da1"]] <- 
  eg_association %>%
  filter(cod_healthy_interaction == T & cod_da1_interaction == F) %$%
  gene_name

sapply(depleted_interactions, length)

emerged_interactions <- list()
emerged_interactions[["da1"]] <- 
  eg_association %>%
  filter(cod_healthy_interaction == F & cod_da1_interaction == T) %$%
  gene_name

sapply(emerged_interactions, length)

#------------------------------------------------



# Functional analysis

#library(clusterProfiler)

depleted_functions <-
  gProfileR::gprofiler(depleted_interactions)


emerged_functions <-
  gProfileR::gprofiler(emerged_interactions)



