.testClusters <- function(locCor, FDR.cluster){
                                        # Procedure 1 at level q
                                        # Weights:
                                        # I)  sum(w.cluster) = m
                                        # II) w.cluster = fac * length.cluster
  Z.cluster <- locCor$Z.cluster
  length.cluster <- locCor$length.cluster
  sigma.cluster <- locCor$sigma.cluster
  pValsList <- locCor$pValsList
  
  q <- FDR.cluster
  m <- length(Z.cluster)
  fac <- m/sum(length.cluster)
  w.cluster <- fac * length.cluster

  p.cluster <- 1 - pnorm(Z.cluster / sigma.cluster)
  ind.sort <- order(p.cluster)
  p.cluster.sorted <- p.cluster[ind.sort]
  w.cluster.sorted <- w.cluster[ind.sort]
  k.ind <- which(p.cluster.sorted <= cumsum(w.cluster.sorted)/m * q)
  if(length(k.ind) == 0){
    stop("No CpG clusters rejected.\n")
  }
  k <- max(k.ind)
  u.1 <- p.cluster.sorted[k]

  ind.reject <- which(p.cluster <= u.1)
  n.reject <- length(ind.reject)
  if(n.reject == 0){
    cat("No CpG clusters rejected.\n")
  } else{
    cat(n.reject, "CpG clusters rejected.\n")
    CpGs.clust.reject <- pValsList[ind.reject]
    CpGs.clust.not.reject <- pValsList[-ind.reject]
    cluster.chr <- sapply(pValsList, function(x) x$chr[1])
    cluster.id <- names(cluster.chr)
    cluster.start <- sapply(pValsList, function(x) x$pos[1])
    cluster.end <- sapply(pValsList, function(x) x$pos[nrow(x)])
    clusters.GR <- GRanges(seqnames = cluster.chr,
                           ranges = IRanges(start=cluster.start,
                             end=cluster.end))
    elementMetadata(clusters.GR) <- cluster.id
    clusters.reject <- clusters.GR[ind.reject]
    clusters.not.reject <- clusters.GR[-ind.reject]
    sigma.clusters.reject <- sigma.cluster[ind.reject]
    return(list(FDR.cluster = FDR.cluster,
                CpGs.clust.reject = CpGs.clust.reject,
                CpGs.clust.not.reject = CpGs.clust.not.reject,
                clusters.reject = clusters.reject,
                clusters.not.reject = clusters.not.reject,
                sigma.clusters.reject = sigma.clusters.reject, 
                variogram=locCor$variogram,
                m = m,
                k = k,
                u.1 = u.1))
                
  }
}


setMethod("testClusters",
          signature=c(locCor = "list", FDR.cluster = "numeric"),
          .testClusters)


setMethod("testClusters",
          signature=c(locCor = "list", FDR.cluster = "missing"),
          function(locCor) {
            .testClusters(locCor, FDR.cluster = 0.05)
          })
