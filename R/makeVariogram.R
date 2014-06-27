.makeVariogram <- function(test.out, make.variogram, sample.clusters, max.dist){
  test.out$z.score <- qnorm(test.out$p.val, lower.tail=FALSE)
  ind <- which(abs(test.out$z.score) == Inf)
  test.out$z.score[ind] <- NA
  cl.p <- test.out[!is.na(test.out$z.score),]
  cl.p.list <- split(cl.p, cl.p$cluster.id, drop=TRUE)
  rm(cl.p)
  if(is.numeric(sample.clusters)){
      ind.sample <- sample(seq(along=cl.p.list), size = sample.clusters, replace = FALSE)
  }
  # positions within clusters:
  pos.new <- lapply(cl.p.list,
                    function(x){
                      start <- min(x$pos)
                      x$pos <- x$pos - min(x$pos) + 1
                    }
                    )
  cl.p.list <- mapply(function(x,y){
      x$pos.new <- y
      return(x)},
                      cl.p.list, pos.new,
                      SIMPLIFY=FALSE)
  if(is.null(sample.clusters)){
    cl.p.list.sample <- cl.p.list
  }
  if(is.numeric(sample.clusters)){
    cl.p.list.sample <- cl.p.list[ind.sample]
  }
  if(make.variogram == TRUE){
    data.list <- lapply(cl.p.list.sample,
                        function(x){
                          y <- x$z.score
                          names(y) <- x$pos.new
                          return(y)
                        }
                        )
    positions.sample <- sort(unique(do.call("c", lapply(cl.p.list.sample, function(x) x$pos.new))))
    positions <- sort(unique(do.call("c", lapply(cl.p.list, function(x) x$pos.new))))
    geo.data <- matrix(numeric(), ncol=length(cl.p.list.sample), nrow=length(positions.sample))
    rownames(geo.data) <- positions.sample
    for(i in seq(along=data.list)){
      x <- data.list[[i]]
      geo.data[names(x), i] <- x
    }
      # geo.data: z-score for each position relative to sample clusters
      # positions: all positions relative to all clusters
    vario <- .variogram(geo.data, positions, max.dist)
    vario$v <- vario$v[!is.na(vario$v[,"v"]),]
    return(list(variogram=vario, pValsList=cl.p.list))
  } else{
    return(list(variogram=NULL, pValsList=cl.p.list))
  }
}


setMethod("makeVariogram",
          signature=c(test.out = "data.frame", make.variogram = "logical", sample.clusters = "numeric", max.dist = "numeric"),
          .makeVariogram)


setMethod("makeVariogram",
          signature=c(test.out = "data.frame", make.variogram = "missing", sample.clusters = "missing", max.dist = "missing"),
          function(test.out) {
            .makeVariogram(test.out, make.variogram=TRUE, sample.clusters = NULL, max.dist = 500)
          })

setMethod("makeVariogram",
          signature=c(test.out = "data.frame", make.variogram = "missing", sample.clusters = "numeric", max.dist = "numeric"),
          function(test.out, sample.clusters, max.dist) {
            .makeVariogram(test.out, make.variogram=TRUE, sample.clusters = sample.clusters, max.dist = max.dist)
          })

setMethod("makeVariogram",
          signature=c(test.out = "data.frame", make.variogram = "logical", sample.clusters = "missing", max.dist = "missing"),
          function(test.out, make.variogram) {
            .makeVariogram(test.out, make.variogram=make.variogram, sample.clusters = NULL, max.dist = 500)
          })

setMethod("makeVariogram",
          signature=c(test.out = "data.frame", make.variogram = "logical", sample.clusters = "missing", max.dist = "numeric"),
          function(test.out, make.variogram, max.dist) {
            .makeVariogram(test.out, make.variogram=make.variogram, sample.clusters = NULL, max.dist = max.dist)
          })

setMethod("makeVariogram",
          signature=c(test.out = "data.frame", make.variogram = "logical", sample.clusters = "numeric", max.dist = "missing"),
          function(test.out, make.variogram, sample.clusters) {
            .makeVariogram(test.out, make.variogram=make.variogram, sample.clusters = sample.clusters, max.dist = 500)
          })
