runClusterRep <- function(obj, tmp.folder, louvain.lib, resolution, seed.use, ...) {
  UseMethod("runClusterRep", obj);
}

runClusterRep.default <- function(
  obj,
  tmp.folder,
  louvain.lib=c("R-igraph", "leiden"),
  resolution=1.0,
  seed.use=10,
  ...
){
  cat("Epoch: checking input parameters\n", file = stderr())

  if(missing(obj)){
    stop("obj is missing");
  }else{
    if(!is.snap(obj)){
      stop("obj is not a snap object");
    }
  }

  if(missing(tmp.folder)){
    stop("tmp.folder is missing")
  }else{
    if(!dir.exists(tmp.folder)){
      stop("tmp.folder does not exist");
    }
  }

  louvain.lib = match.arg(louvain.lib);

  if(louvain.lib == "leiden"){
    if (!requireNamespace("leiden", quietly = TRUE)) {
      stop("Please install leiden - learn more at https://github.com/TomKellyGenetics/leiden")
    }
  }

  if(isKgraphEmpty(obj@graph)){
    stop("knn graph is empty, run 'runKNN' first")
  }

  if(is.na(as.numeric(resolution))){
    stop("resolution must be numeric class!")
  }

  if(louvain.lib == "R-igraph"){
    data.use = getGraph(obj@graph);
    data.use = data.use + t(data.use);
    g = graph_from_adjacency_matrix(data.use, weighted=TRUE, mode="undirected");
    cat("Epoch: finding clusters using R-igraph\n", file = stderr())
    set.seed(seed.use);
    cl = cluster_louvain(g);
    obj@cluster = factor(cl$membership);
  }else if(louvain.lib == "leiden"){
    cat("Epoch: finding clusters using leiden\n", file = stderr())
    data.use = getGraph(obj@graph);
    set.seed(seed.use);
    obj@cluster <- factor(leiden(data.use, resolution=resolution, seed = 10, ...));
  }else{
    stop("unrecognized louvain.lib option")
  }
  gc();
  return(obj);
}

isKgraphEmpty <- function(obj){
  if(is.null(obj)){
    stop("obj is empty")
  }else{
    if(!is(obj, "kgraph")){
      stop("obj is not a kgraph object");
    }else{
      if((x = nrow(obj@mat)) > 0L){
        return(FALSE)
      }
      if((x=length(obj@file) > 0L)){
        if(file.exists(obj@file)){
          return(FALSE)
        }
      }
    }
  }
  return(TRUE);
}


getGraph <- function(obj){
  if(is.null(obj)){
    stop("obj is empty")
  }else{
    if(!is(obj, "kgraph")){
      stop("obj is not a kgraph object");
    }
    if(isKgraphEmpty(obj)){
      stop("obj is empty");
    }
  }

  if((x=nrow(obj@mat) != 0L)){
    return(obj@mat);
  }

  if(file.exists(obj@file)){
    edgeList = read.table(obj@file, header=FALSE);
    if((x=ncol(edgeList)) != 3L){
      stop(paste(obj@file, " does not have 3 columns"));
    }
    num.node = max(edgeList[,c(1,2)]);
    M1 = sparseMatrix(i=edgeList[,1], j=edgeList[,2], x=edgeList[,3], dims=c(num.node,num.node));
    M2 = sparseMatrix(i=edgeList[,2], j=edgeList[,1], x=edgeList[,3], dims=c(num.node,num.node));
    M = M1 + M2;
    rm(M1, M2);
    gc();
    return(M);
  }
}


readBins <- function(file, bin.size=5000){
  if(exists('h5closeAll', where='package:rhdf5', mode='function')){
    rhdf5::h5closeAll();
  }else{
    rhdf5::H5close();
  }
  if(!file.exists(file)){stop(paste("Error @addBmat: ", file, " does not exist!", sep=""))};
  if(!isSnapFile(file)){stop(paste("Error @addBmat: ", file, " is not a snap-format file!", sep=""))};
  if(!(bin.size %in% showBinSizes(file))){stop(paste("Error @addBmat: bin.size ", bin.size, " does not exist in ", file, "\n", sep=""))};
  options(scipen=999);
  binChrom = tryCatch(binChrom <- h5read(file, paste("AM", bin.size, "binChrom", sep="/")), error = function(e) {stop(paste("Warning @readaddBmatSnap: 'AM/bin.size/binChrom' not found in ",file))})
  binStart = tryCatch(binStart <- h5read(file, paste("AM", bin.size, "binStart", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/binStart' not found in ",file))})
  if(bin.size == 0){
    binEnd = tryCatch(binEnd <- h5read(file, paste("AM", bin.size, "binEnd", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/binStart' not found in ",file))})
  }else{
    binEnd = binStart + as.numeric(bin.size) -1;
  }
  if((length(binChrom) == 0) || (length(binStart) == 0)){stop("Error @addBmat: bin is empty! Does not support empty snap file")}
  if(length(binChrom) != length(binStart)){
    stop(paste("Error @addBmat: ", "binChrom and binStart has different length!", sep=""))
  }else{
    nBin = length(binChrom);
  }
  bins = GRanges(binChrom, IRanges(as.numeric(binStart),binEnd), name=paste(paste(binChrom, binStart, sep=":"), binEnd, sep="-"));
  rm(binChrom, binStart);

  # Updated
  bins = sortSeqlevels(bins);
  bins = sort(bins);
  return(bins)
}


addBmatToSnapSingle <- function(obj, file, bin.size=5000){
  # close the previously opened H5 file
  if(exists('h5closeAll', where='package:rhdf5', mode='function')){
    rhdf5::h5closeAll();
  }else{
    rhdf5::H5close();
  }

  if(missing(obj)){
    stop("obj is missing")
  }else{
    if(!is(obj, "snap")){
      stop("obj is not a snap object")
    }
  }

  if(!file.exists(file)){stop(paste("Error @addBmat: ", file, " does not exist!", sep=""))};
  if(!isSnapFile(file)){stop(paste("Error @addBmat: ", file, " is not a snap-format file!", sep=""))};
  if(!(bin.size %in% showBinSizes(file))){stop(paste("Error @addBmat: bin.size ", bin.size, " does not exist in ", file, "\n", sep=""))};
  obj@bmat = Matrix(0,0,0);
  ############################################################################
  barcode = as.character(tryCatch(barcode <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @addBmat: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));

  bin.sizeList = showBinSizes(file);
  if(length(bin.sizeList) == 0){stop("Error @addBmat: bin.sizeList is empty! Does not support reading empty snap file")}
  if(!(bin.size %in% bin.sizeList)){stop(paste("Error @addBmat: ", bin.size, " does not exist in bin.sizeList, valid bin.size includes ", toString(bin.sizeList), "\n", sep=""))}

  bins = readBins(file, bin.size);
  obj@feature = bins;
  idx = as.numeric(tryCatch(idx <- h5read(file, paste("AM", bin.size, "idx", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/idx' not found in ",file))}));
  idy = as.numeric(tryCatch(idy <- h5read(file, paste("AM", bin.size, "idy", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/idy' not found in ",file))}));
  count = as.numeric(tryCatch(count <- h5read(file, paste("AM", bin.size, "count", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/count' not found in ",file))}));

  if(!all(sapply(list(length(idx),length(idy),length(count)), function(x) x == length(count)))){stop("Error: idx, idy and count has different length in the snap file")}
  nBarcode = length(barcode);
  nBin = length(obj@feature);
  M = sparseMatrix(i=idx, j =idy, x=count, dims=c(nBarcode, nBin));
  rownames(M) = barcode;
  obj@bmat = M[match(obj@barcode, rownames(M)),]
  rm(idx, idy, count, M);
  if(exists('h5closeAll', where='package:rhdf5', mode='function')){
    rhdf5::h5closeAll();
  }else{
    rhdf5::H5close();
  }
  gc();
  return(obj);
}

addBmatToSnap <- function(obj, bin.size, do.par, num.cores){
  UseMethod("addBmatToSnap", obj);
}

#' @export
addBmatToSnap.default <- function(obj, bin.size=5000, do.par=FALSE, num.cores=1){
  # close the previously opened H5 file
  if(exists('h5closeAll', where='package:rhdf5', mode='function')){
    rhdf5::h5closeAll();
  }else{
    rhdf5::H5close();
  }
  if(missing(obj)){
    stop("obj is missing")
  }else{
    if(!is(obj, "snap")){
      stop("obj is not a snap object")
    }
  }

  if(!is.numeric(num.cores)){
    stop("num.cores is not an integer")
  }
  num.cores = round(num.cores);

  fileList = as.list(unique(obj@file));

  # check if snap files exist
  if(any(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)){
    idx = which(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)
    print("error: these files does not exist")
    print(fileList[idx])
    stop()
  }

  # check if files are all snap files
  if(any(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)){
    idx = which(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)
    print("error: these files are not snap file")
    print(fileList[idx])
    stop()
  }

  # check if BM session exist
  if(any(do.call(c, lapply(fileList, function(x){ "AM" %in% h5ls(x, recursive=1)$name  })) == FALSE)){
    idx = which(do.call(c, lapply(fileList, function(x){ "AM" %in% h5ls(x, recursive=1)$name  })) == FALSE)
    print("error: the following nsap files do not contain AM session")
    print(fileList[idx])
    stop()
  }

  if(any(do.call(c, lapply(fileList, function(x){(bin.size %in% showBinSizes(x))})) == FALSE)){
    idx = which(do.call(c, lapply(fileList, function(x){(bin.size %in% showBinSizes(x))})) == FALSE)
    print("error: chosen bin size does not exist in the following snap files")
    print(fileList[idx])
    stop()
  }

  # check if bins match
  bin.list = lapply(fileList, function(x){
    readBins(x, bin.size=bin.size)
  })

  if(!all(sapply(bin.list, FUN = identical, bin.list[[1]]))){
    stop("bins does not match between snap files, please regenerate the cell-by-bin matrix by snaptools")
  }

  # read the snap object
  message("Epoch: reading cell-bin count matrix session ...");
  if(do.par){
    obj.ls = mclapply(fileList, function(file){
      idx = which(obj@file == file)
      addBmatToSnapSingle(obj[idx,], file, bin.size=bin.size);
    }, mc.cores=num.cores);
  }else{
    obj.ls = lapply(fileList, function(file){
      idx = which(obj@file == file)
      addBmatToSnapSingle(obj[idx,], file, bin.size=bin.size);
    });
  }

  # combine
  if((x=length(obj.ls)) == 1L){
    res = obj.ls[[1]]
  }else{
    res = Reduce(snapRbind, obj.ls);
  }
  obj@feature = res@feature;
  obj@bmat = res@bmat;
  rm(res, obj.ls);
  gc()
  return(obj);
}
