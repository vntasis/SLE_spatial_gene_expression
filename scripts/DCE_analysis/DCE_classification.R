########Find cods with at least 80% overlap########
find_overlap <- function(group, c_CODS = c_CODs, control_group){
  print(group)
  
  overlap <- list()
  for (chr in names(c_CODS[[group]])){
    print(chr)
    overlap[[chr]] <- matrix(0,1,3, dimnames = list(NULL, c("control", "test", "overlap")))
    x <- c_CODS[[group]][[chr]]
    co <- control_group[[chr]]
    
    if(nrow(x) > 0){
      for (i in 1:nrow(x)){
        
        k <- 1
        repeat {
          # print(chr)
          # print(paste0("i = ", i))
          # print(paste0("k = ", k))
          if (k > nrow(co)) break
          else if (x[i, "START"] >= co[k, "END"]) k <- k+1
          else if (x[i, "END"] <= co[k, "START"]) break
          else if (((x[i, "START"] >= co[k, "START"] && x[i, "START"] < co[k, "END"]) ||  (x[i, "END"] > co[k, "START"] && x[i, "END"] <= co[k, "END"])) ||
                   ((co[k, "START"] >= x[i, "START"] && co[k, "START"] < x[i, "END"]) ||  (co[k, "END"] > x[i, "START"] && co[k, "END"] <= x[i, "END"]))) {
            
            o <- length(intersect(x[i, "START"]:x[i, "END"], co[k, "START"]:co[k, "END"]))/length(union(x[i, "START"]:x[i, "END"], co[k, "START"]:co[k, "END"]))
            overlap[[chr]] <- rbind(overlap[[chr]], c(k,i,o))
            k <- k+1
          }
        }
      }
    }
    overlap[[chr]] <- overlap[[chr]][-1, , drop = F]
  }
  return(overlap)
}

overlap <- lapply(X = c("da1", "da2", "da3"), FUN = find_overlap, control_group = c_CODs$healthy)
names(overlap) <- c("da1", "da2", "da3")
lapply(overlap, function(y) sapply(y, function(x)summary(x[, "overlap"])))

N_cods <- function(overlap, threshold){
  output <- matrix(0,1, length(names(overlap[[1]])))
  for (group in names(overlap)){
    #print(group)
    n_cods <- numeric()
    for (chr in names(overlap[[group]])){
      #print(chr)
      # x <- overlap[[group]][[chr]][,"overlap"]
      # n_cods <- c(n_cods, sum(x >= threshold))
      x <- unique(overlap[[group]][[chr]][which(overlap[[group]][[chr]][,"overlap"] >= threshold),"test"])
      n_cods <- c(n_cods, length(x))
    }
    output <- rbind(output, n_cods)
  }
  output <- output[-1,]
  dimnames(output) <- list(names(overlap), names(overlap[[1]]))
  return(output)
}
N_cods(overlap, 0)
N_cods(overlap, 0.8)
n_cods_per_chr <- noCpch(CODs)
P_cods <- function(number_of_cods, overlap, threshold){
  x <- N_cods(overlap, threshold)
  for (i in 1:nrow(x)){
    x[i,] <- x[i,]/number_of_cods[[rownames(x)[i]]]
  }
  return(x)
}
P_cods(n_cods_per_chr, overlap, 0.8)
P_cods(n_cods_per_chr, overlap, 0)


########Find and categorize COD rearrangements (split, merge, emerged, depleted, shrinkage, enlargment, displacement, intact, one border altered, both borders altered)########
find_rearrangements <- function(group, overlap, control_group, c_CODS){
  
  output <- list()
  for (chr in names(overlap[[group]])){
    #overlap[[group]][[chr]] <- overlap[[group]][[chr]][overlap[[group]][[chr]][, "overlap"] >= 0.1, , drop = F]
    output[[chr]] <- list()
    
    cod_row_indices_control <- seq_along(control_group[[chr]][,1])
    cod_row_indices_test <- seq_along(c_CODS[[group]][[chr]][,1])
    output[[chr]][["depleted"]] <- setdiff(cod_row_indices_control, overlap[[group]][[chr]][,"control"])
    output[[chr]][["emerged"]] <- setdiff(cod_row_indices_test, overlap[[group]][[chr]][,"test"])
    
    
    freq_control <- table(overlap[[group]][[chr]][,"control"])
    freq_test <- table(overlap[[group]][[chr]][,"test"])
    
    categorized <- numeric()
    output[[chr]][["split"]] <- matrix(0,1,3, dimnames = list(NULL, c("control", "test", "overlap")))
    if(sum(freq_control > 1) > 0){
      for (spl in names(freq_control)[freq_control > 1]){
        output[[chr]][["split"]] <- rbind(output[[chr]][["split"]], overlap[[group]][[chr]][which(overlap[[group]][[chr]][,"control"] == spl), ])
        categorized <- c(categorized, which(overlap[[group]][[chr]][,"control"] == spl))
      }
    }
    output[[chr]][["split"]] <- output[[chr]][["split"]][-1, , drop = F]
    
    output[[chr]][["merged"]] <- matrix(0,1,3, dimnames = list(NULL, c("control", "test", "overlap")))
    if(sum(freq_test > 1) > 0){
      for (spl in names(freq_test)[freq_test > 1]){
        output[[chr]][["merged"]] <- rbind(output[[chr]][["merged"]], overlap[[group]][[chr]][which(overlap[[group]][[chr]][,"test"] == spl), ])
        categorized <- c(categorized, which(overlap[[group]][[chr]][,"test"] == spl))
      }
    }
    output[[chr]][["merged"]] <- output[[chr]][["merged"]][-1, , drop = F]
    
    output[[chr]][["intact"]] <- overlap[[group]][[chr]][overlap[[group]][[chr]][, "overlap"] == 1 , "control"]
    categorized <- c(categorized, which(overlap[[group]][[chr]][, "overlap"] == 1))
    
    #print(categorized)
    ifelse(length(categorized) > 0, non_categorized <- overlap[[group]][[chr]][-categorized, , drop = F],
           non_categorized <- overlap[[group]][[chr]])
    border_alteration <- data.frame("control" = non_categorized[,"control"], "test" = non_categorized[, "test"],
                                    "overlap" = non_categorized[, "overlap"], "border_altered" = rep("", nrow(non_categorized)),
                                    stringsAsFactors = F)
    if (nrow(non_categorized) > 0){
      for (i in 1:nrow(non_categorized)){
        control <- non_categorized[i, "control"]
        test <- non_categorized[i, "test"]
        
        if (control_group[[chr]][control, "START"] != c_CODS[[group]][[chr]][test, "START"]) border_alteration[["border_altered"]][i] <- "left"
        if (control_group[[chr]][control, "END"] != c_CODS[[group]][[chr]][test, "END"]){
          x <- "right"
          if (border_alteration[["border_altered"]][i] == "left") x <- "both"
          border_alteration[["border_altered"]][i] <- x
        }
      }
    }
    output[[chr]][["border_alteration"]] <- border_alteration
  }
  return(output)
}
rearrangements <- lapply(c("da1", "da2", "da3"), find_rearrangements, overlap = overlap, control_group = c_CODs$healthy, c_CODS = c_CODs)
names(rearrangements) <- c("da1", "da2", "da3")


# recategorize border_shift
for (group in names(rearrangements)){
  for (chr in names(rearrangements[[group]])){
    
    rearrangements[[group]][[chr]][["expanded"]] <- 
      matrix(0, 1, 3, dimnames = list(NULL, c("control", "test", "overlap")))
    rearrangements[[group]][[chr]][["contracted"]] <- 
      matrix(0, 1, 3, dimnames = list(NULL, c("control", "test", "overlap")))
    rearrangements[[group]][[chr]][["shift"]] <- 
      matrix(0, 1, 3, dimnames = list(NULL, c("control", "test", "overlap")))
    
    if (nrow(rearrangements[[group]][[chr]][["border_alteration"]]) > 0){
      for (i in 1:nrow(rearrangements[[group]][[chr]][["border_alteration"]])){
        #print(i)
        size_control <- s_CODs[["healthy"]][[chr]][
          rearrangements[[group]][[chr]][["border_alteration"]][i, "control"]
          ]
        size_test <- s_CODs[[group]][[chr]][
          rearrangements[[group]][[chr]][["border_alteration"]][i, "test"]
          ]
        # print(size_control)
        # print(size_test)
        if (size_test > size_control){
          
          rearrangements[[group]][[chr]][["expanded"]] <-
            rbind(rearrangements[[group]][[chr]][["expanded"]], 
                  rearrangements[[group]][[chr]][["border_alteration"]][i, 1:3])
          
        } else if (size_test < size_control){
          
          rearrangements[[group]][[chr]][["contracted"]] <-
            rbind(rearrangements[[group]][[chr]][["contracted"]],
                  rearrangements[[group]][[chr]][["border_alteration"]][i, 1:3])
          
        } else if (size_test == size_control){
          
          rearrangements[[group]][[chr]][["shift"]] <-
            rbind(rearrangements[[group]][[chr]][["shift"]],
                  rearrangements[[group]][[chr]][["border_alteration"]][i, 1:3])
        }
      }
    }
    rearrangements[[group]][[chr]][["expanded"]] <- 
      rearrangements[[group]][[chr]][["expanded"]][-1, , drop = F]
    rearrangements[[group]][[chr]][["contracted"]] <- 
      rearrangements[[group]][[chr]][["contracted"]][-1, , drop = F]
    rearrangements[[group]][[chr]][["shift"]] <- 
      rearrangements[[group]][[chr]][["shift"]][-1, , drop = F]
  }
}
