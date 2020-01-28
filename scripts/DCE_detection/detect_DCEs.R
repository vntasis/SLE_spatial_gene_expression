codfinder <- function(chr, group, window.size, bin_signal_thres = 0.25, perm_pvalue_thres){
  
  # calculate bin signal
  bin_signal <- binsignal(chr = chr, thres = perm_pvalue_thres, group = group, window.size = window.size)
  
  # apply spline interpolation
  binsignal_spline <- features((1:nrow(bin_signal)), bin_signal[,1], 
                             smoother="smooth.spline",  
                             all.knots = T)

  # refine local minima and local maxima
  local_min <- round(binsignal_spline$cpts[binsignal_spline$curvature > 0])
  signal <- bin_signal[local_min, 1]
  local_min <- local_min[signal < bin_signal_thres]

  local_max <- round(binsignal_spline$cpts[binsignal_spline$curvature < 0])
  signal <- bin_signal[local_max, 1]
  local_max <- local_max[signal >= bin_signal_thres]

  rm(signal)



  # detect initial cod regions
  cod_matrix <- matrix(0,1,2)
  cod <- c(0,0)
  #interfere_index <- 0
  for (lmax in local_max){

      if (!(lmax %in% cod_matrix[nrow(cod_matrix),1]:cod_matrix[nrow(cod_matrix),2])){
          
          if (length(local_min[local_min < lmax]) > 0) cod[1] <- max(local_min[local_min < lmax]) else{
              cod[1] <- 1
          }
          if (length(local_min[local_min > lmax]) > 0) cod[2] <- min(local_min[local_min > lmax]) else{
              cod[2] <- nrow(bin_signal)
          }
          cod_matrix <- rbind(cod_matrix, cod)

      }
      cod <- c(0,0)
  }

  rm(cod)

  # Refine cod borders
  discard <- integer()
  for (i in 2:nrow(cod_matrix)){
      start <- 0
      end <- 0

      for (j in cod_matrix[i,1]:cod_matrix[i,2]){
          
          if (j == 1) start <- j
          else if (start == 0 && bin_signal[j,1] >= bin_signal_thres && bin_signal[j-1,2] <= 0.05) start <- j

      }
      for (j in cod_matrix[i,2]:cod_matrix[i,1]){

          if (start != 0 && end == 0 && bin_signal[j,1] >= bin_signal_thres && bin_signal[j,3] <= 0.05) end <- j

      }

      if (start == 0 || end == 0) discard <- c(discard, i)
      else{

          cod_matrix[i,1] <- start
          cod_matrix[i,2] <- end

      }

  }

  # discard cods with insignificant bordrers 
  if (length(discard) > 0) cod_matrix <- cod_matrix[-discard, , drop = F]
  rm(discard, start, end)
  
  
  colnames(cod_matrix) <- c("START", "END")
  cod_matrix <- cod_matrix[-1, , drop = F]
  if(nrow(cod_matrix) >= 1){
   
    # filter out cods smaller than window.size bins
    size <- cod_matrix[,2]-cod_matrix[,1]+1
    cod_matrix <- cod_matrix[size > (window.size - 1), , drop = F]
    
    
  }
    
  list(bin_signal = bin_signal, CODs = cod_matrix)
}

