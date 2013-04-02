.variogram <- function(geo.data, positions){
  ind.ord <- order(positions)
  geo.data.n <- geo.data[ind.ord,,drop=FALSE]
  positions.n <- positions[ind.ord]
  all.na <- apply(geo.data.n, 1, function(x) all(is.na(x)))
  geo.data.n <- geo.data.n[!all.na,,drop=FALSE]
  positions.n <- positions.n[!all.na]
  d <- as.matrix(dist(cbind(positions.n)))
  d[upper.tri(d, diag=TRUE)] <- NA
  h <- sort(unique(as.vector(d)))

  v <- matrix(numeric(), ncol=2, nrow=length(h))
  colnames(v) <- c("h", "v")
  i.h <- 0
  pb <- txtProgressBar(min=0, max=length(h), style=3)
  for(i in seq(along=h)){
    h.i <- h[i]
    ind <- d==h.i
    ind.1 <- apply(ind, 1, any, na.rm=TRUE)
    ind.2 <- apply(ind, 2, any, na.rm=TRUE)
    v[i,2] <- median( (geo.data.n[ind.1,,drop=FALSE] - geo.data.n[ind.2,,drop=FALSE] )^2, na.rm=TRUE ) / 0.455 / 2
    v[i,1] <- h.i
    i.h <- i.h + 1
    setTxtProgressBar(pb, value=i.h)
  }
  close(pb)

  return(v)
}

