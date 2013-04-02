.smoothVariogram <- function(variogram, sill, bandwidth){
  if(is(variogram, "list")){
    vario <- variogram$variogram
  } else {
    vario <- variogram
  }
  pred <- lokerns(x = vario[,"h"],
                  y = vario[,"v"],
                  x.out=vario[,"h"],
                  bandwidth = bandwidth)
  pred$est[ pred$est < 0 ] <- 0
  vario.sm <- data.frame(vario, v.sm = pred$est)
  if(!missing(sill)){
    ind <- which(vario.sm[,"v.sm"] > sill)
    if(length(ind) > 0){
      ind.first <- min(ind)
      vario.sm[ind.first:nrow(vario.sm), "v.sm"] <- sill
    }
  }
  if(is(variogram, "list")){
    return(list(variogram=vario.sm, pValsList=variogram$pValsList))
  } else {
    return(vario.sm)
  }
}

setMethod("smoothVariogram",
          signature=c(variogram = "matrix", sill = "numeric", bandwidth = "numeric"),
          .smoothVariogram)

setMethod("smoothVariogram",
          signature=c(variogram = "matrix", sill = "numeric", bandwidth = "missing"),
          function(variogram, sill){
            bandwidth <- seq(10,1000, length.out=nrow(variogram))
            .smoothVariogram(variogram, sill, bandwidth)
          }
          )

setMethod("smoothVariogram",
          signature=c(variogram = "list", sill = "numeric", bandwidth = "numeric"),
          .smoothVariogram)


setMethod("smoothVariogram",
          signature=c(variogram = "list", sill = "numeric", bandwidth = "missing"),
          function(variogram, sill){
            bandwidth <- seq(10,1000, length.out=nrow(variogram$variogram))
            .smoothVariogram(variogram, sill, bandwidth)
          }
          )
