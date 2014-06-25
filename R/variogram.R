.variogram <- function(geo.data, positions){
  positions.sample <- as.integer(rownames(geo.data))
  ind.ord <- order(positions.sample)
  geo.data.n <- geo.data[ind.ord,,drop=FALSE]
  positions.sample.n <- positions.sample[ind.ord]
  all.na <- apply(geo.data.n, 1, function(x) all(is.na(x)))
  geo.data.n <- geo.data.n[!all.na,,drop=FALSE]
  positions.sample.n <- positions.sample.n[!all.na]
  d.sample <- as.matrix(dist(cbind(positions.sample.n)))
  d.sample[upper.tri(d.sample, diag=TRUE)] <- NA
  h <- sort(unique(as.vector(d.sample)))

  if(identical(positions, positions.sample)){
    h.est <- h
  } else {
    ind.ord <- order(positions)
    positions.n <- positions[ind.ord]
    d <- as.matrix(dist(cbind(positions.n)))
    d[upper.tri(d, diag=TRUE)] <- NA
    h.est <- sort(unique(as.vector(d)))
  }

  v <- matrix(numeric(), ncol=2, nrow=length(h))
  colnames(v) <- c("h", "v")
  i.h <- 0
  pb <- txtProgressBar(min=0, max=length(h), style=3)
  for(i in seq(along=h)){
    h.i <- h[i]
    ind <- d.sample==h.i
    ind.1 <- apply(ind, 1, any, na.rm=TRUE)
    ind.2 <- apply(ind, 2, any, na.rm=TRUE)
    v[i,2] <- median( (geo.data.n[ind.1,,drop=FALSE] - geo.data.n[ind.2,,drop=FALSE] )^2, na.rm=TRUE ) / 0.455 / 2
    v[i,1] <- h.i
    i.h <- i.h + 1
    setTxtProgressBar(pb, value=i.h)
  }
  close(pb)

  return(list(v=v, h.est=h.est))
}

