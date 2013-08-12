.trimClusters <- function(clusters.rej, FDR.loc){
  
  # 2. a) Definitions
  q <- clusters.rej$FDR.cluster
  m.0 <- (clusters.rej$m-clusters.rej$k)/(1-q)
  all.clusters <- c(clusters.rej$CpGs.clust.reject,
                    clusters.rej$CpGs.clust.not.reject)
  z.li <- unlist(sapply(all.clusters, function(x) x$z.score))
  c.i <- sapply(all.clusters, nrow)
  mu.i <- (sum(z.li)/sum(c.i)) / clusters.rej$sigma.clusters.reject # standardized expectation of the cluster test statistic under H_1i

 # 2. b) Estimate p values

  integrand <- function(u){
    out <- ( (m.0/m) * (pnorm( (qnorm(u.1, lower.tail=FALSE) - rho.li*u) / sqrt(1-rho.li^2) , lower.tail=FALSE)) 
    + (1-m.0/m) * (pnorm( (qnorm(u.1, lower.tail=FALSE) - rho.li*u - mu.i[i]) / sqrt(1-rho.li^2), lower.tail=FALSE)) ) * dnorm(u)
  return(out)
  }

  CpGs.clust.reject <- clusters.rej$CpGs.clust.reject
  sigma.clusters.reject <- clusters.rej$sigma.clusters.reject
  m <- clusters.rej$m
  u.1 <- clusters.rej$u.1
  variogram <- clusters.rej$variogram
  rownames(variogram) <- as.character(variogram[,"h"])
  for(i in seq(along=CpGs.clust.reject)){
    part <- clusters.rej$CpGs.clust.reject[[i]]
    c.i <- nrow(part)
    z.li <- part$z.score
    for(loc in 1:c.i){
      diffs <- abs(part$pos[loc] - part$pos[-loc])
      ind.diffs <- sapply(diffs, function(x) which(x == variogram[, "h"]))
      rho.lk <- 1 - variogram[as.character(diffs), "v.sm"]
      rho.li <- (1 + sum(rho.lk)) / (c.i*sigma.clusters.reject[i])
      rho.li <- min(rho.li, 1) # Sometimes, smoothing of the variogram results in 0 for very small distances. To avoid correlations > 1.
      integral <- integrate(integrand, lower=z.li[loc], upper=Inf)
      p <- integral$value * 1 / (m.0/m * u.1 + (1-m.0/m) * pnorm((qnorm(u.1, lower.tail=FALSE)-mu.i[i]), lower.tail=FALSE))
      CpGs.clust.reject[[i]]$p.li[loc] <- p
    }
  }


  # 2. c) two-stage procedure at level FDR.loc
  q.2 <- FDR.loc
  clusters.reject.unlist <- do.call("rbind", CpGs.clust.reject)
  p.loc <- clusters.reject.unlist$p.li
  ind.sort <- order(p.loc)
  p.loc.sorted <- p.loc[ind.sort]
  m.2 <- length(p.loc)
  q.2. <- q.2/(1+q.2)
  p.threshold.1 <- (seq(along=p.loc.sorted)/m.2)*q.2.
  k.1 <- max(which(p.loc.sorted <= p.threshold.1))
  m.02 <- m.2 - k.1
  p.threshold.2 <- (seq(along=p.loc.sorted)/m.02)*q.2.
  k.2 <- max(which(p.loc.sorted <= p.threshold.2))
  p.threshold <- (k.2/m.02)*q.2.
  ind.reject <- which(p.loc <= p.threshold)

  clusters.trimmed <- clusters.reject.unlist[ind.reject,]
  # clusters.trimmed$p.val <- clusters.trimmed$p.li
  
  return(clusters.trimmed)

}

setMethod("trimClusters",
          signature=c(clusters.rej = "list", FDR.loc = "numeric"),
          .trimClusters)


setMethod("trimClusters",
          signature=c(clusters.rej = "list", FDR.loc = "missing"),
          function(clusters.rej) {
            .trimClusters(clusters.rej, FDR.loc = 0.2)
          })
