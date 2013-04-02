.estLocCor <- function(vario.sm){
  cl.p.list <- vario.sm$pValsList
  variogram <- vario.sm$variogram[,c("h", "v.sm")]
  sigma.cluster <- numeric(length=length(cl.p.list))
  Z.cluster <- numeric(length=length(cl.p.list))
  length.cluster <- numeric(length=length(cl.p.list))
  for(i in seq(along=cl.p.list)){
    n.na <- 0
    part <- cl.p.list[[i]]
    c.i <- nrow(part)
    if( c.i == 1 ){
      sigma.cluster[i] <- 1
    } else {
      comb <- combn(1:c.i, 2)
      dist.comb <- part$pos[comb[2,]] - part$pos[comb[1,]]
      dist.comb <- sort(dist.comb)
      ind <- which(is.element(variogram[,"h"], dist.comb))
      dist.comb.tab <- table(dist.comb)
      # if variogram is not available for a distance, the missing correlation is subsetted by its median
      if(length(ind) < length(dist.comb.tab)){
        ind.na <- which(!is.element(dist.comb, variogram[,"h"]))
        n.na <- length(ind.na)
        dist.comb <- dist.comb[-ind.na]
        ind <- which(is.element(variogram[,"h"], dist.comb))
        dist.comb.tab <- table(dist.comb)
      }
      ind.vario <- rep(ind, times=as.vector(dist.comb.tab))
      cor <- 1 - variogram[ind.vario,"v.sm"]
      if(n.na > 0){
        median.cor <- median(cor)
        cor <- c(cor, rep(median.cor, n.na))
      }
      sigma.cluster[i] <- 1/c.i * sqrt( c.i + 2*sum(cor, na.rm=TRUE) )
    }
    Z.cluster[i] <- mean(part$z.score, na.rm=TRUE)
    length.cluster[i] <- max(part$pos) - min(part$pos)
  }
  return(list(variogram=vario.sm$variogram,
              pValsList=vario.sm$pValsList,
              sigma.cluster=sigma.cluster,
              Z.cluster=Z.cluster,
              length.cluster=length.cluster))
}


setMethod("estLocCor",
          signature=c(vario.sm = "list"),
          .estLocCor)
