library(lmtest)

compute.gene.expr.rank <- function(path, name, TFs) {

  ## Load ES metrices
  # Load gene enrichment score
  ges <- read.table(gzfile(paste0(path, "/", name, "/", "ges.esw.csv.gz")),
                    header = T, sep = ",")

  rownames(ges) <- ges$gene
  ges <- ges[,-1]

  # Load normalized specificity index
  si <- read.table(gzfile(paste0(path, "/", name, "/", "nsi.esw.csv.gz")),
                   header = T, sep = ",")
  rownames(si) <- si$gene
  si <- si[,-1]

  # Load expression proportion
  ep <- read.table(gzfile(paste0(path, "/", name, "/", "ep.esw.csv.gz")),
                   header = T, sep = ",")
  rownames(ep) <- ep$gene
  ep <- ep[,-1]

  # Load t-statistics
  tstat <- read.table(gzfile(paste0(path, "/", name, "/", "det.esw.csv.gz")),
                      header = T, sep = ",")
  rownames(tstat) <- tstat$gene
  tstat <- tstat[,-1]

  ges <- ges[na.omit(match(TFs, rownames(ges))),]
  ges.rank <- apply(ges, 2, rank)
  si <- si[na.omit(match(TFs, rownames(si))),]
  si.rank <- apply(si, 2, rank)
  ep <- ep[na.omit(match(TFs, rownames(ep))),]
  ep.rank <- apply(ep, 2, rank)
  tstat <- tstat[na.omit(match(TFs, rownames(tstat))),]
  tstat.rank <- apply(tstat, 2, rank)

  rank <- data.frame(matrix(NA, nrow = nrow(ges), ncol = ncol(ges)))
  dimnames(rank) <- dimnames(ges)

  for (i in 1:nrow(rank)) {
    for (j in 1:ncol(rank)) {
      mean <- mean(ges.rank[i,j], si.rank[i,j],
                   ep.rank[i,j], tstat.rank[i,j])
      rank[i,j] <- mean
    }
  }
  return(rank)
}

compute.average.access <- function(object, ident, TFs) {

  average.access <- data.frame(matrix(NA, nrow = length(TFs),
                                      ncol = length(unique(object@meta.data[,ident]))))
  rownames(average.access) <- TFs
  colnames(average.access) <- unique(object@meta.data[,ident])

  for (i in colnames(average.access)) {
    average <- apply(object@assays$chromvar@data[, which(object@meta.data[,ident] == i)], 1, mean)
    average.access[, i] <- average
  }
  return(average.access)
}


compute.modality.cor<- function(gene.expr.rank, average.access) {

  gene.expr.rank <- gene.expr.rank[, colnames(average.access)]

  cor <- data.frame(matrix(NA, nrow = nrow(gene.expr.rank), ncol = 2))
  colnames(cor) <- c("Rho", "p.value")

  for (i in 1:nrow(gene.expr.rank)) {

    cor.test <- cor.test(as.numeric(gene.expr.rank[i, ]),
                         as.numeric(average.access[rownames(average.access) == rownames(gene.expr.rank)[i], ]),
                         method = "spearman", exact = F, alternative="greater")

    cor$Rho[i] <- cor.test$estimate
    cor$p.value[i] <- cor.test$p.value
  }
  return(cor)

}

find.markers <- function(object, features, group, ident.1, ident.2 = NULL) {

  if(is.null(ident.2)) {

    ident.2 <- unique(object@meta.data[, group][!(object@meta.data[, group] %in% ident.1)])
  }

  data.use <- object@assays$chromvar@data[features, ]
  cells.1 <- colnames(object)[which(object@meta.data[, group] %in% ident.1)]
  cells.2 <- colnames(object)[which(object@meta.data[, group] %in% ident.2)]

  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])

  data.use <- data.use[, rownames(group.info), drop = FALSE]

  da.motifs <- data.frame(matrix(NA, nrow = nrow(data.use),
                                 ncol = 4))
  rownames(da.motifs) <- rownames(data.use)
  colnames(da.motifs) <- c("average.zscore1", "average.zscore2",
                           "zscore.diff", "p")

  for (i in 1:nrow(da.motifs)) {

    average.zscore1 <- mean(data.use[i, which(group.info$group == "Group1")])
    average.zscore2 <- mean(data.use[i, which(group.info$group == "Group2")])
    zscore.diff <- average.zscore1 - average.zscore2

    # LR test
    model.data <- cbind(GENE = data.use[i, ], group.info)
    fmla <- as.formula(object = "group ~ GENE")
    fmla2 <- as.formula(object = "group ~ 1")

    model1 <- glm(formula = fmla, data = model.data, family = "binomial",
                  control = list(maxit = 50))
    model2 <- glm(formula = fmla2, data = model.data, family = "binomial",
                  control = list(maxit = 50))
    lrtest <- lrtest(model1, model2)

    da.motifs$average.zscore1[i] <- average.zscore1
    da.motifs$average.zscore2[i] <- average.zscore2
    da.motifs$zscore.diff[i] <- zscore.diff

    da.motifs$p[i] <- lrtest$Pr[2]
  }
  return(da.motifs)

}