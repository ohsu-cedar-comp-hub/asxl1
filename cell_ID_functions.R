
#Automate ID of clusters from a Seurat object and the haemosphere bone marrow data
#haemoshere bone marrow data: download from https://www.haemosphere.org/datasets/show  
#Haemopedia-Mouse-RNASeq version 1.5

#Author: Trevor Enright

####example on use:

###either:
# download data from haemosphere and run: 
# Haemopedia <- read.delim("/PATH/TO/Haemopedia-Mouse-RNASeq_raw.txt")
# Haem_prep <- preprocess_downloaded_raw_counts(Haemopedia_matrix = Haemopedia) 
# Haem_data <- get_Haem_data(Continue_with = Haem_prep, assay = "RNA", scale_max = 3) 
#to prepare to matrix of normalized and scaled counts of haemopoesis cells
#
### OR (recommended)
# use the file from the github repo: https://github.com/ohsu-cedar-comp-hub/asxl1
# Haem_data <- readRDS(/PATH/TO/FILE/FILE.RDS)
# preprocess a seurat object through to obtaining cluster information ie, After FindNeighbors() and FindClusters(), 
# seuratObj <- add_cell_haem_data(seuratObj = seuratObj, markers = NULL, assay = "RNA", res = "seurat_clusters", genes_data = Haem_data, top_unique = 1, num_genes = 20, lineage = NULL, quick = TRUE, max_pval=0.05)


#Required Libraries
library(plyr)
library(dplyr)
library(Seurat) #version 3.1.0
library(biomaRt)

#Helper Function
substrRight <- function(x, n){
  #returns n number of characters from character x
  substr(x, nchar(x)-n+1, nchar(x))
}


#Helper function
get_intersection <- function(seuratObj,external, seurat.markers, num_genes = 20, assay = "RNA", max_pval = 1) {
  # returns genes from external (haemopedia data, with common names) that of DE markers "seurat.markers" from FindAllMarkers that intersect
  # These genes will be used for averaging cell expression of those specific genes in the get_cell_haem_data
  
  maxCluster <- max(as.numeric(levels(seurat.markers$cluster)))
  #remove genes of p value adjusted statsistically significant
  seurat.markers <- seurat.markers[which(seurat.markers$p_val_adj <= max_pval),]
  
  if (is.null(num_genes)){num_genes <- nrow(seurat.markers)}
  
  Top <- NULL
  Bottom <- NULL
  #Create data frames of Top genes per cluster and bottom genes per cluster
  for (i in 0:maxCluster){
    Top <- rbind(Top,seurat.markers[which(seurat.markers$cluster == i),] %>% top_n(n = num_genes, wt = avg_logFC))
    Bottom <- rbind(Bottom,seurat.markers[which(seurat.markers$cluster == i),] %>% top_n(n = num_genes, wt = -avg_logFC))
  }
  #ensuring genes are in haemopedia
  genes_of_interest <- intersect(external[["external_gene_name"]],Top$gene)
  #ensuring genes are in scaled.data
  genes_of_interest <- genes_of_interest[which(genes_of_interest %in% rownames(seuratObj@assays[[assay]]@scale.data))]  
  
  return(Top[which(genes_of_interest %in% Top$gene),])
  
}

#cluster ID function
get_cell_type <- function (seuratObj, cluster, n = NULL, res = "integrated_snn_res.0.1", tops, gene_data, top_unique = 1, show_ranking = FALSE, lineage = NULL, assay = "integrated", quick = TRUE) {

  # seuratObj is the seurat object to use for the clustering
  # res is the resolution / Ident used for which meta data that clustered the cells
  # n is the number of cells to sample for per cluster default all cells per cluster
  # tops is the a sorted or otherwise filtered output from FindAllMarkers() on the seuratObj.
  # gene_data is the output from get_Haem_data() or a data_matrix with cell names as columns and entries as normalized and scaled count. and the last column is gene names.
  #ensure number of cells is acceptable
  # lineage is a lookup table (2 column dataframe with "celltype" as one column and "Lineage" in another). default is NULL which returns cell type 
  # quick accepts a boolean value to either run without attempting to find more markers if the number of genes is insufficient to calculate an average or to simply return the cluster number (default)
  #
  # averages gene expression across cells and finds highest correlation among cell expression of cell types in the haemopedia data (gene_data)
  #returns a list of names associated with 
  max_n <- length(which(seuratObj@meta.data[[res]] == cluster))
  if (is.null(n)){
    n <- max_n
  } 
  else if (n > max_n){
    warning(sprintf("number of cells is too large! Continuing with all %s cells in cluster %s",max_n,cluster))
    n <- max_n
  }
  else if (n < 2) {
    #ensures rowMeans works
    n <- 2
  }
  #vector of numbers
  calcuated_mean <- NULL
  correlate_calc <- NULL
  calcuated_med <- NULL
  #sample cells from 
  cells_to_check <- sample(which(seuratObj@meta.data[[res]] == cluster),size = n, replace = F)
  
  #get genes that are in both datasets (ie: in bulk RNA hematopedia and DE genes for the given cluster)
  gene_names <- intersect(gene_data[["external_gene_name"]],tops$gene[which(tops$cluster == cluster)])
  
  
  if ((length(gene_names) <= 5) & (quick == FALSE)){
    old_length <- length(gene_names)
    Idents(seuratObj) <- res
    DefaultAssay(seuratObj) <- assay
    tops <- FindMarkers(seuratObj,logfc.threshold = 0.1,min.pct = 0.01, ident.1 = cluster)
    gene_names <- intersect(gene_data[["external_gene_name"]],rownames(tops))
    warning(sprintf("cluster: %s experienced a low length of genes to test at only %s. Searching for DE genes in this cluster and relaxing parameters.  If you do not desire this feature, input quick = TRUE.  Now currently ran at %s number of genes", cluster,old_length, length(gene_names)))
  } 
  if (length(gene_names) <= 1){
    warning(sprintf("cluster: %s experienced a low length of genes to test at only %s.  Returning the cluster number to prevent errors", cluster, length(gene_names)))
    if (quick == TRUE){warning(sprintf("FindMarkers with relaxed parameters did not find sufficient DE genes for the cluster %s", cluster))}
    else {warning("Consider relaxing parameters for FindAllMarkers or rerun with 'quick = TRUE'. Or perhaps this is not an important cluster and the resolution is too high") }
    return(cluster)
  }
  
  average_expression_cells <- rowMeans(seuratObj@assays[[assay]]@scale.data[gene_names,colnames(seuratObj[,cells_to_check])])
  Haem_genes <- which(gene_data$external_gene_name %in% gene_names)
  
  df <- data.frame(genes = gene_names, average_expression = average_expression_cells)
  
  #loop through all (cell types) but the last column (common gene names) in the normalized and scaled haemosphere data
  for (i in 1:(ncol(gene_data)-1)){
    correlate_calc[i] <- cor(x = gene_data[Haem_genes,i],y = average_expression_cells, method = "pearson")
    
    calcuated_med[i] <- median((gene_data[Haem_genes,i] - average_expression_cells)**2)
    calcuated_mean[i] <- mean((gene_data[Haem_genes,i] - average_expression_cells)**2)
    
  }
  
  #get return data top unique names
  data_return <- data.frame(cor = correlate_calc, median = calcuated_med, average = calcuated_mean, cell_name = rev(rev(colnames(gene_data))[-1]))
  
  #filter out 90% of high average difference
  data_return <- filter(data_return, data_return$average < quantile(data_return$average,0.10))
  ##sort by correlation
  data_return <- data_return[order(-data_return$cor),]
  
  #get the number of unique top matches and keep only the amount requested (typically 1 to 3)
  data_return$cell_type <- gsub("\\..*","",data_return$cell_name)
  if (is.null(lineage)){
    name_unique <- as.data.frame(table(data_return$cell_type))
    name_unique <- name_unique[which(name_unique$Var1 %in% head(unique(data_return$cell_type, fromLast = T), n = top_unique)),]
    
    if (show_ranking) {
      data_return <- paste(name_unique$Var1,name_unique$Freq,sep="_",collapse = "/")
    }
    else {
      data_return <- paste(name_unique$Var1,collapse = "/")
    }
    
    #return the column name that best fits the profile without the extra stuff after the period
    #return(gsub("\\..*","",colnames(gene_data)[which(calcuated_mean == min(calcuated_mean))]))
    return(data_return)
  }
  else {
    data_return <- merge(x = data_return,y = lineage,by.x = "cell_type", by.y = "celltype", all.x = T)
    data_return <- data_return[order(-data_return$cor),]
    #data_return$cell_lineage <- gsub("\\ .*","",data_return$cell_lineage)
    data_return <- head(unique(data_return$Lineage, fromLast = T), n = top_unique)
    return(paste(data_return,collapse = "/"))
  }
}



add_cell_haem_data <- function (seuratObj, markers = NULL, assay = NULL, genes_data, res = NULL, num_genes = 20, k = NULL, top_unique = 1, lineage = NULL, sub_name = NULL, max_pval = 0.05, quick = TRUE) {
  ##
  # Adds meta data to a seuratObj.  takes "markers" if already computed from FindAllMarkers(seuratObj,only.pos = TRUE) otherwise will compute
  # "res" (resolution) which is a cluster column from meta.data such as seurat_clusters (default), recommended for mouse hematopoetic integrated wt and mut objects is "integrated_snn_res.0.1"
  # or a number to find a meta.data column, or a string
  # "lineage" is a lookup table (2 column dataframe with "celltype" as one column and "Lineage" in another). default is NULL which returns cell type 
  # "genes_data" which is for the external cell count reference (haemopedia - output from get_Haem_data)
  # "num_genes" is the number of genes to use for averaging gene expression from FindAllMarkers(),  if num_genes = NULL then it takes the max number of genes per cluster
  # "k" is the number of cells per cluster to average across to test (default all cells per cluster)
  # "top_unique" is the number of cell types to return for each clusters. 
  # "sub_name" takes a name to use as the column name for the new meta data containing the cell ID.  this can be used to help understand the same clusters at different parameters within the same object
  # "k" is the amount of cells within the cluster to sample and perform average expression (default all the cells within a cluster)
  # "max_pval" is the maximum p_val_adj (from markers) allowed to perform average expression. increase to 1 to get all genes
  # "quick" accepts a boolean value to either run without attempting to find more markers if the number of genes is insufficient to calculate an average gene expression for a given cluster OR to simply return the cluster number (default : TRUE)
  if (is.null(assay)) {
    assay <- DefaultAssay(seuratObj)  
  }
  seurat.markers <- markers
  current_res <- res
  nameIt <- "new_cell_type"
  if (!(is.null(sub_name))){
    nameIt <- paste(sub_name,nameIt,sep = "_")
  }
  
  if (is.character(current_res)){
    nameIt <- gsub(current_res,sprintf("cell_type_%s_%s",assay,substrRight(current_res,4)),current_res)
    if (!(is.null(sub_name))){
      nameIt <- paste(sub_name,nameIt,sep = "_")
    }
  }
  else if(is.numeric(current_res)){
    current_res <- grep(sprintf("^%s_snn_res.%s$",assay,current_res),colnames(seuratObj@meta.data), value = T)
    
    nameIt <- sprintf("cell_type_%s_%s",assay,current_res)
    if (!(is.null(sub_name))){
      nameIt <- paste(sub_name,nameIt,sep = "_")
    }
  }
  else if(is.null(res)){
    current_res <- "seurat_clusters"
    nameIt <- sprintf("cell_type_%s_%s",assay,current_res)
    if (!(is.null(sub_name))){
      nameIt <- paste(sub_name,nameIt,sep = "_")
    }
    if (!(current_res %in% colnames(seuratObj@meta.data))) {stop("Seurat object needs to have clusters, as in snn or otherwise.  update 'res' variable")}
    
  }
  seuratObj <- SetIdent(object = seuratObj, value = current_res)
  
  if (is.null(markers)){
    print("Finding DE genes between clusters")
    seurat.markers <- FindAllMarkers(object = seuratObj, 
                                     only.pos = TRUE, 
                                     verbose = FALSE)
  } else {
    seurat.markers <- markers
  }
  
  ## get number of clusters as per current res
  maxCluster <- max(c(unique(seuratObj@meta.data[[current_res]])))-1
  
  
  #get genes to compare against the data and the cell information
  Top <- get_intersection(seuratObj = seuratObj,external = genes_data, seurat.markers = seurat.markers,num_genes = num_genes, assay = assay, max_pval = max_pval)
  
  
  
  cluster_id <- NULL
  for (cluster in 0:maxCluster) {
    print(paste("Calculating cluster", cluster, sep = ": "))
    cluster_id[cluster+1] <- get_cell_type(seuratObj = seuratObj,cluster = cluster, gene_data = genes_data , n = k,res = current_res, tops = Top, top_unique = top_unique, lineage = lineage, assay = assay, quick = quick)
  }
  #set up for adding cell type to seurat object
  
  cell_cluster <- data.frame(cell = seuratObj@meta.data[[current_res]],
                             row.names = rownames(seuratObj@meta.data))
  
  ###################### MAP VALUES ##################################
  ########## to add meta data column of "cell" which is cell type ####
  cell_cluster$cell <- mapvalues(cell_cluster$cell,
                                 from = as.character(c(0:maxCluster)),
                                 to = cluster_id)
  
  
  
  # Finally, add meta data
  seuratObj <- AddMetaData(seuratObj, metadata = cell_cluster, col.name = nameIt)
  #update Idents in object to the new column name
  Idents(seuratObj) <- nameIt
  return(seuratObj)
}




preprocess_downloaded_raw_counts <- function(Haemopedia_matrix){
  # Adds common names to "raw counts summarized as genes" (transcipts per million) from Hemosphere website
  # used function read.delim()
  # Load common names to use instead of Ensembl gene ID
  # This may take a while.  Recommend saving rds of "gene_tranfser" lookup table
  
  ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
  filters = listFilters(ensembl)
  gene_transfer <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                     filters =  "ensembl_gene_id", 
                     values = Haemopedia_matrix$geneId, #this assumes the ensembl gene IDs are in a column named "geneId"; possible failure point for updated versions
                     mart = ensembl)
  
  
  #incorporate common names into the Haemopedia dataset
  Haemopedia <- merge(gene_transfer,Haemopedia_matrix,by.x = "ensembl_gene_id", by.y = "geneId", all.y = T)
  return(Haemopedia)
}

get_Haem_data <- function(Continue_with, assay = "RNA", scale_max = 3){
  #Continue_with is a data.frame / matrix from output of function 'preprocess_downloaded_raw_counts'
  # returns a scaled matrix of the raw transcripts per million, with cell names as column names and entries are the normalized and scaled count. and the last column is gene names.
  # The assay can be "SCT" but is recommended to use "RNA" (default)
  # scale_max is the max value to scale the data to. 3 (default). Please ensure this is identical to the scale.max found
  # in ScaleData() of your Seurat object.
  
  Haem_it_up <- Continue_with[,-(1:2)]
  
  #get column names which are for a cell_id
  col_names_type <- gsub("\\..*","",colnames(Haem_it_up))
  
  #Creating additional cells by averaging cell types with only two samples
  #to ensure Seurat scaling/normalizing will work
  for (i in unique(col_names_type)){
    mean_calc <- NULL
    columns_to_compute <- Haem_it_up[,which(col_names_type == i)]
    if (is.data.frame(columns_to_compute)) {
      mean_calc <- rowMeans(columns_to_compute)
      Haem_it_up[[i]] <- as.integer(mean_calc)
    }
  }
  
  
  #get matrix of Haemopedia
  Haem_mat <- as.matrix(Haem_it_up)
  
  rownames(Haem_mat) <- Continue_with$ensembl_gene_id
  #create seurat object of the data / normalize and scale
  Haem_BM <- CreateSeuratObject(Haem_mat,project = "BM_Haem")
  Haem_BM[["percent.mt"]] <- PercentageFeatureSet(Haem_BM, pattern = "^mt-")
  if(assay == "SCT"){
    Haem_BM <- SCTransform(Haem_BM, vars.to.regress = "percent.mt", return.only.var.genes = FALSE)
  }
  Haem_BM <- NormalizeData(object = Haem_BM, assay = "RNA")
  Haem_BM <- ScaleData(object = Haem_BM, assay = "RNA", features = rownames(Haem_BM),scale.max = scale_max)
  Haem_data <- Haem_BM@assays[[assay]]@scale.data
  
  # merge by ensembl gene ID to get common names
  Haem_data <- merge(gene_transfer,Haem_data,by.x = "ensembl_gene_id", by.y = 0, all.y = T)
  
  #save common name
  Haem_names <- Haem_data$external_gene_name
  #obtain only scaled data
  Haem_data <- Haem_data[,-(1:2)]
  #Create row names
  rownames(Haem_data) <- make.names(Haem_names,unique = T)
  # adds to the end of the matrix
  Haem_data$external_gene_name <- Haem_names
  
  #Haem_data only consists of data from Haemopedia not added data
  Haem_data <- Haem_data[,which(grepl("\\.",colnames(Haem_data)) | colnames(Haem_data) == "external_gene_name")]
  return(Haem_data)
}