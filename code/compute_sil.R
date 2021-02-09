
library(doMC)
registerDoMC(20)


compute.sil <- function(x, dist) {
  dist <- as.matrix(dist)
  
  data.b <- foreach(idx.point = 1:length(x), .combine = rbind) %dopar% {
    
    
    b <- Inf
    nearest.cluster <- "NA"
    
    for (cluster in unique(x)) {
      
      if (cluster != x[idx.point]) {
        
        idx.cluster <- which(x == cluster)
        b.try <- mean(dist[idx.point, idx.cluster])
        nearest.cluster <- c(nearest.cluster, cluster)[which.min(c(b, b.try))]
        b <- min(b, b.try)
      }
    }
    
    data.frame(b = b, nearest.cluster = nearest.cluster)
  }
  
  
  data.a <- foreach (idx.point = 1:length(x), .combine = rbind) %dopar% {
    
    cluster <- x[idx.point]
    idx.cluster <- which(x == cluster)[-idx.point]
    a <- mean(dist[idx.point, idx.cluster])
    data.frame(a = a)  
  }
  
  sil <- (data.b$b - data.a$a) / mapply(max, data.a$a, data.b$b)
  return(sil)
  
}