
# -------------------Parameter Checking -----------------------
ParameterChecking <- function(tissue_names = NULL, tissue_expr_file_names = NULL, MV_sd_thr = 0.5, TS_power = 6, CT_power = 3) {
  
  parameter_checking_success <- FALSE
  
  if(!is.null(tissue_names) & !is.null(tissue_expr_file_names) & is.double(MV_sd_thr) & is.double(TS_power) & is.double(CT_power)) {
    
    if(is.vector(tissue_names) & is.character(tissue_names)
       & is.vector(tissue_expr_file_names) & is.character(tissue_expr_file_names)) {
      
      if(length(tissue_names) != length(tissue_expr_file_names)) {stop("Missmatch Tissues name and its file names ....")}
      
      parameter_checking_success <- TRUE
      
    } 
  } 
  
  return(parameter_checking_success)
}

# -------------------Load expression data and filtered most variant genes --------------
LoadExprData<-function(tissue_name, tissue_file_name, MV_sd_thr) {
  
  
  if(substr(tissue_file_name, start = nchar(tissue_file_name)-5, stop = nchar(tissue_file_name)) == ".RData") {
    
    datExpr <- get(load(tissue_file_name))
  } else if(substr(tissue_file_name, start = nchar(tissue_file_name)-3, stop = nchar(tissue_file_name)) == ".txt") {
    
    datExpr <- read.delim(tissue_file_name, header = TRUE, sep = '\t', quote = "")
    datExpr <- as.matrix(datExpr)
    rownames(datExpr) <- datExpr[ ,1]
    datExpr <- datExpr[ ,-1]
    datExpr <- apply(datExpr, c(1,2), as.numeric)
  } else stop("Unsupported input file !!!")
  
  
  colnames(datExpr) <- paste(tissue_name,"_",colnames(datExpr),sep="")
  
  MV_datExpr <- apply(datExpr,2,sd)
  
  MV_sd_thr <- abs(MV_sd_thr)
  if(MV_sd_thr >= 0.0 & MV_sd_thr < 1.0) {
    MV_datExpr <- MV_datExpr[MV_datExpr>MV_sd_thr]
  } else {
    MV_sd_thr <- round(MV_sd_thr)
    MV_datExpr <- sort(MV_datExpr, decreasing=TRUE)
    MV_datExpr <- MV_datExpr[1:MV_sd_thr]
  }
  
  datExpr <- datExpr[,names(MV_datExpr)]
  
  return(datExpr)
  
}

# -------------------Similarity matrix from all tissue expression data ----------------------------------
SimilarityFromExpr <- function(tissue_names = NULL, tissue_expr_file_names = NULL, MV_sd_thr = 0.5, cor_method = "pearson") {
  
  cor_mat <- AdjacencyFromExpr(tissue_names, tissue_expr_file_names, MV_sd_thr, TS_power = 1, CT_power = 1, cor_method)
  
  return(cor_mat)
}

# -------------------Adjacency matrix from all tissue expression data ----------------------------------------- 
AdjacencyFromExpr <- function(tissue_names = NULL, tissue_expr_file_names = NULL, MV_sd_thr = 0.5, TS_power = 6, CT_power = 3, cor_method = "pearson") {
  
  if(!ParameterChecking(tissue_names, tissue_expr_file_names, MV_sd_thr, TS_power = 6, CT_power = 3)) stop("Missing or Invalid Parameters ! ! !")
  
  total_tissue <- length(tissue_names)
  
  vector_expr<-vector(mode="list", length=total_tissue)
  
  ##----This is to track tissue in whole adjacency matrix--------------------
  tissue_index_adj<-vector(mode="integer", length=(total_tissue+1))
  rc_names<-c()
  
  for(i in 1:total_tissue) {
    
    TS_expr_data<-LoadExprData(tissue_names[i], tissue_expr_file_names[i], MV_sd_thr)
    vector_expr[[i]]<-TS_expr_data
    
    tissue_index_adj[i+1]<- tissue_index_adj[i]+ncol(vector_expr[[i]])
    rc_names<-c(rc_names,colnames(vector_expr[[i]]))
    
  }
  
  adj_mat<-matrix(0,nrow=tissue_index_adj[total_tissue+1], ncol=tissue_index_adj[total_tissue+1])
  rownames(adj_mat)<-rc_names
  colnames(adj_mat)<-rc_names
  
  for(i in 1:(total_tissue-1)) {
    adj_mat[(tissue_index_adj[i]+1):tissue_index_adj[i+1], (tissue_index_adj[i]+1):tissue_index_adj[i+1]] <- abs(cor(vector_expr[[i]], method = cor_method))^TS_power
    
    for(j in (i+1):total_tissue) {
      common_Samples <- intersect(rownames(vector_expr[[i]]),rownames(vector_expr[[j]]))
      adj_mat[(tissue_index_adj[i]+1):tissue_index_adj[i+1],(tissue_index_adj[j]+1):tissue_index_adj[j+1]] <- abs(cor(vector_expr[[i]][common_Samples,],vector_expr[[j]][common_Samples,], method = cor_method))^CT_power
      adj_mat[(tissue_index_adj[j]+1):tissue_index_adj[j+1],(tissue_index_adj[i]+1):tissue_index_adj[i+1]] <- t(adj_mat[(tissue_index_adj[i]+1):tissue_index_adj[i+1],(tissue_index_adj[j]+1):tissue_index_adj[j+1]])
    }
  }
  
  ##----------------This is for final(last) Tissue specific Adj matrix---------------------------------------------------
  adj_mat[(tissue_index_adj[total_tissue]+1):tissue_index_adj[total_tissue+1],(tissue_index_adj[total_tissue]+1):tissue_index_adj[total_tissue+1]] <- abs(cor(vector_expr[[total_tissue]], method = cor_method))^TS_power
  
  return(adj_mat)
}

# ------------------TOM from Adjacency matrix ----------------------- 

Cross_Tissue_TOM <- function(adj_mat) {
  
  l <- adj_mat %*% t(adj_mat)
  print("Matrix Multiplication done...")
  rsums<-rowSums(adj_mat)
  D_val<-matrix(rsums,nrow=nrow(adj_mat),ncol=ncol(adj_mat))
  
  total_col <- ncol(D_val)
  if(total_col <= 5000) col_limit <- total_col else col_limit <- 5000
  start_col <- 1
  end_col <- col_limit
  
  tmp_D <- matrix(0, nrow = total_col, ncol = total_col)
  
  repeat {
    tmp_D[ ,start_col:end_col] <- pmin(D_val[ ,start_col:end_col], t(D_val)[ ,start_col:end_col])
    if(end_col==total_col) break
    start_col <- end_col + 1
    end_col <- end_col + col_limit
    if(end_col > total_col) end_col <- total_col
  }  
  
  rm(D_val)
  rm(rsums)
  gc()
  print("Matrix pmin done...")
  
  TOM_mat <- matrix(NA, nrow=nrow(adj_mat), ncol=ncol(adj_mat))
  rownames(TOM_mat) <- rownames(adj_mat)
  colnames(TOM_mat) <- colnames(adj_mat)
  
  TOM_mat <- (adj_mat + l) / (tmp_D + 1 - adj_mat)
  
  print("Tom done...")
  return(TOM_mat) 
  
}

# ------------------Cross-Tissue Clusters Table from TOM matrix --------------

Clusters_Table <- function(TOM_mat, minClusterSize = 20) {
  
  library('dynamicTreeCut')
  dist_TOM_mat <- 1 - TOM_mat
  rm(TOM_mat)
  gc()
  h_TOM <- hclust(as.dist(dist_TOM_mat), method="average")
  dynamicMods = cutreeDynamic(h_TOM, method = "tree", deepSplit = TRUE, minClusterSize = minClusterSize);
  
  all_gene_names <- h_TOM$labels
  total_clusters <- length(table(dynamicMods)) - 1
  vector_clusters <- vector(mode="list", length=(total_clusters))
  for(i in 1:total_clusters) {
    
    index <- which(dynamicMods == i)
    vector_clusters[[i]] <- all_gene_names[index]
    names(vector_clusters[[i]]) <- index
    
  }
  
  clusters_table <- do.call(rbind, lapply(seq_along(vector_clusters), function(i) {data.frame(Cluster_ID=i, Gene_Symbol=vector_clusters[[i]])}))
  
  clusters_table <- as.matrix(clusters_table)
  clusters_table <- cbind(clusters_table, '')
  colnames(clusters_table) <- c('Cluster ID', 'Tissue', 'Gene Symbol')
  
  clusters_table[ ,c(2,3)] <- cbind(unlist(lapply(strsplit(clusters_table[ ,2], split = '_'), function(x) {return(x[1])})), unlist(lapply(strsplit(clusters_table[ ,2], split = '_'), function(x) {return(x[2])})))
  
  return(clusters_table)
}

# -----------------Clusters Deatils Information from clusters table ---------------

Clusters_Details <- function(clusters_table, cluster_type_thr = 0.95) {
  
  tissue_names <- names(table(clusters_table[ ,"Tissue"]))
  total_clusters <- max(as.numeric(clusters_table[ ,"Cluster ID"]))
  
  clusters_table_details <- matrix(0, nrow = total_clusters, ncol = (5+length(tissue_names)))
  colnames(clusters_table_details) <- c('Cluster ID', 'Cluster Size', 'Cluster Type', 'Cluster Tissues', tissue_names, 'Dominant Tissue')
  
  for(i in 1:total_clusters) {
    
    temp_cluster <- clusters_table[as.numeric(clusters_table[ ,"Cluster ID"]) == i, ]
    
    clusters_table_details[i, "Cluster ID"] <- i
    
    clusters_table_details[i, "Cluster Size"] <- nrow(temp_cluster)
    
    if(max(round(table(temp_cluster[ ,"Tissue"])/nrow(temp_cluster), 2)) >= cluster_type_thr) clusters_table_details[i, "Cluster Type"] <- 'TS' else clusters_table_details[i, "Cluster Type"] <- 'CT'
    
    clusters_table_details[i, "Cluster Tissues"] <- paste(names(table(temp_cluster[ ,"Tissue"])), collapse = ',')
    
    clusters_table_details[i, names(table(temp_cluster[ ,"Tissue"]))] <- table(temp_cluster[ ,"Tissue"])
    
    clusters_table_details[i, "Dominant Tissue"] <- names(which.max(round(table(temp_cluster[ ,"Tissue"]))))
    
  }
  
  return(clusters_table_details)
}

# ------------------Cross-Tissue Clusters from gene expression data ------------------------

XWGCNA_Clusters <- function(tissue_names = NULL, tissue_expr_file_names = NULL, MV_sd_thr = 0.5, cluster_type_thr = 0.95, minClusterSize = 20) {
  
  adj_mat <- AdjacencyFromExpr(tissue_names, tissue_expr_file_names, MV_sd_thr)
  
  TOM_mat <- Cross_Tissue_TOM(adj_mat)
  
  clusters_table <- Clusters_Table(TOM_mat, minClusterSize = minClusterSize)
  
  clusters_details <- Clusters_Details(clusters_table, cluster_type_thr = cluster_type_thr)
  
  write.table(clusters_table, file = 'Clusters_table.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  write.table(clusters_details, file = 'Clusters_details.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  XWGCNA_clusters = list(clusters_table=clusters_table, clusters_details=clusters_details)
  
  return(XWGCNA_clusters)
}