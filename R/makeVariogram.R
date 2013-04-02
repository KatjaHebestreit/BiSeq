.makeVariogram <- function(test.out, make.variogram){
  test.out$z.score <- qnorm(test.out$p.val, lower.tail=FALSE)
  cl.p <- test.out[!is.na(test.out$z.score),]
  cl.p.list <- split(cl.p, cl.p$cluster.id, drop=TRUE)
  rm(cl.p)
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
  data.list <- lapply(cl.p.list,
                      function(x){
                        y <- x$z.score
                        names(y) <- x$pos.new
                        return(y)
                      }
                      )
  positions <- sort(unique(do.call("c", sapply(cl.p.list, function(x) x$pos.new))))
  geo.data <- matrix(numeric(), ncol=length(cl.p.list), nrow=length(positions))
  rownames(geo.data) <- positions
  for(i in seq(along=data.list)){
    x <- data.list[[i]]
    geo.data[names(x), i] <- x
  }
  if(make.variogram == TRUE){
    vario <- .variogram(geo.data, positions)
    vario <- vario[!is.na(vario[,"v"]),]
    return(list(variogram=vario, pValsList=cl.p.list))
  } else{
    return(list(variogram=NULL, pValsList=cl.p.list))
  }
}


setMethod("makeVariogram",
          signature=c(test.out = "data.frame", make.variogram = "logical"),
          .makeVariogram)


setMethod("makeVariogram",
          signature=c(test.out = "data.frame", make.variogram = "missing"),
          function(test.out) {
            .makeVariogram(test.out, make.variogram=TRUE)
          })
