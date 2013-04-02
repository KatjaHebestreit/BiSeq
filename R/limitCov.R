.limitCov <- function(object, maxCov){
  indCov <- totalReads(object) > maxCov
  fraction <- methReads(object)[indCov] / totalReads(object)[indCov]
  totalReads(object)[indCov] <- as.integer(maxCov)
  methReads(object)[indCov] <- as.integer(round(fraction * maxCov))
  return(object)
}

setMethod("limitCov",
          signature=c(object = "BSraw", maxCov="numeric"),
          .limitCov)
