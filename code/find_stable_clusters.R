
library(cluster)
library(parallelDist)
library(doMC)
registerDoMC(20)

### Seurat ###

FindStableClustersSeurat <- function(object, resolution, dims) {

  # Compute graph
  object <- FindNeighbors(object, reduction = "pca", verbose = F)


  # Compute distance matrix
  pc <- object@reductions$pca@cell.embeddings[, dims]
  distance <- parDist(pc, method = "euclidean")
  
  stability <- foreach(i = 1:length(resolution), .combine = cbind) %dopar% {
    
    object <- FindClusters(object = object, resolution = resolution[i],
                           verbose = F)
    cluster.id <- as.numeric(as.character(object@meta.data[, paste0("integrated_snn_res.", resolution[i])]))
    
    silhouette <- silhouette(x = cluster.id, dist = distance)[,3]
    
    data.frame(silhouette)

  }
  
  colnames(stability) <- resolution
  rownames(stability) <- colnames(object)
  
  return(stability)
}



### SnapATAC ###

FindStableClustersSnapATAC <- function(object, resolution, dims) {

  # Compute distance matrix
  pc <- object@smat@dmat[, dims]
  distance <- parDist(pc, method = "euclidean")

  # Compute stability of different clusterings
  stability <- foreach(i = 1:length(resolution), .combine = rbind) %dopar% {
    
    object = runClusterRep(obj = object, tmp.folder=tempdir(),
                           louvain.lib="leiden", resolution = resolution[i], seed.use = 10)
    

    cluster.id <- as.numeric(as.character(object@cluster))
    silhouette <- silhouette(x = cluster.id, dist = distance)[,3]

    data.frame(resolution = i, silhouette = mean(silhouette),
               clusters = length(unique(object@cluster)))
  }

  return(stability)
}
