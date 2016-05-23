checkScaleFree <- function (k, nBreaks = 10, removeFirst = FALSE) 
{
  discretized.k = cut(k, nBreaks)
  dk = tapply(k, discretized.k, mean)
  p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
  breaks1 = seq(from = min(k), to = max(k), length = nBreaks + 1)
  hist1 = hist(k, breaks = breaks1, plot = FALSE, right = TRUE)
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  p.dk = ifelse(is.na(p.dk), 0, p.dk)
  log.dk = as.vector(log10(dk))
  if (removeFirst) {
    p.dk = p.dk[-1]
    log.dk = log.dk[-1]
  }
  log.p.dk = as.numeric(log10(p.dk + 1e-09))
  lm1 = lm(log.p.dk ~ log.dk)
  lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
  datout = data.frame(Rsquared.SFT = summary(lm1)$r.squared, slope.SFT = summary(lm1)$coefficients[2, 1], truncatedExponentialAdjRsquared = summary(lm2)$adj.r.squared)
  datout
}

R2_connectivity <- function(adj_mat) {
  #--------------------------For scale Free-----------------------------
  k <- rowSums(adj_mat)
  r2 <- checkScaleFree(k)$Rsquared.SFT
  r2_conn <- list(r2=r2, mean_conn=mean(k), median_conn=median(k), max_conn=max(k), min_conn=min(k))
  
  #------------------For connectons-----------------------
  tissue_indexed <- c(0, as.numeric(table(unlist(lapply(strsplit(colnames(adj_mat), split = '_'), function(x) {return(x[1])})))))
  for(i in 2:length(tissue_indexed)) tissue_indexed[i] <- tissue_indexed[i] <- sum(tissue_indexed[i:(i-1)])
  
  TS_conn_mat <- matrix(NA, nrow=(length(tissue_indexed)-1), ncol=4)
  colnames(TS_conn_mat) <- c('mean', 'median', 'max', 'min')
  
  CT_conn_mat <- matrix(NA, nrow=((((length(tissue_indexed)-1)*(length(tissue_indexed)-1))-(length(tissue_indexed)-1))/2), ncol=4)
  colnames(CT_conn_mat) <- c('mean', 'median', 'max', 'min')
  CT_counter <- 1
  
  for(i in 1:(length(tissue_indexed)-1)) {
    temp_adj_mat <- adj_mat[(tissue_indexed[i]+1):tissue_indexed[i+1], (tissue_indexed[i]+1):tissue_indexed[i+1]]
    k <- rowSums(temp_adj_mat)
    TS_conn_mat[i, ] <- c(mean(k), median(k), max(k), min(k))
    
    if(i < length(tissue_indexed)-1) {
      for(j in 1:(length(tissue_indexed)-1-i)) {
        temp_adj_mat <- adj_mat[(tissue_indexed[i]+1):tissue_indexed[i+1], (tissue_indexed[i+j]+1):tissue_indexed[i+j+1]]
        k <- rowSums(temp_adj_mat)
        CT_conn_mat[CT_counter, ] <- c(mean(k), median(k), max(k), min(k))
        CT_counter <- CT_counter + 1
      }
    }
    
  }
  
  
  scale_free_conn_list <- list(r2_conn=r2_conn, TS_conn_mat=TS_conn_mat, CT_conn_mat=CT_conn_mat)
  # save(scale_free_conn_list, file='Scale_free_conn_list.RData')
  
  
  R2_Conn <- list(r2=r2_conn$r2, mean_conn=r2_conn$mean_conn, TS_mean_conn=paste(TS_conn_mat [ ,'mean'], collapse = ', '), mean_TS_mean_conn=mean(TS_conn_mat [ ,'mean']),
                  max_TS_mean_conn=max(TS_conn_mat [ ,'mean']), min_TS_mean_conn=min(TS_conn_mat [ ,'mean']), CT_mean_conn=paste(CT_conn_mat[ ,'mean'], collapse = ', '),
                  mean_CT_mean_conn=mean(CT_conn_mat[ ,'mean']), max_CT_mean_conn=max(CT_conn_mat[ ,'mean']), min_CT_mean_conn=min(CT_conn_mat[ ,'mean']))
  
  return(R2_Conn)
}


