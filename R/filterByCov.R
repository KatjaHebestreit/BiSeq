.filterByCov <- function(object, minCov=10, global=FALSE){

  if(global){
    freq.cov <- apply(totalReads(object), 1, function(x) sum(x >= minCov))
    ind <- freq.cov == ncol(object)
    rowData.new <- rowData(object)[ind]
    totalReads.new <- totalReads(object)[ind,]
    methReads.new <- methReads(object)[ind,]
  }else{
    ind <- totalReads(object) >= minCov
    ind.out <- apply(ind, 1, function(x) sum(x) == 0)
    ind.dim <- which(!ind, arr.ind=TRUE)
    rowData.new <- rowData(object)[!ind.out]
    totalReads.new <- totalReads(object)
    totalReads.new[ind.dim] <- 0L
    totalReads.new <- totalReads.new[!ind.out,]
    methReads.new <- methReads(object)
    methReads.new[ind.dim] <- 0L
    methReads.new <- methReads.new[!ind.out,]
  }

  return(BSraw(colData=colData(object), 
               rowData=rowData.new, 
               totalReads=totalReads.new, 
               methReads=methReads.new))
}


setMethod("filterByCov",
          signature=c(object = "BSraw", minCov = "numeric", global = "logical"),
          .filterByCov)

setMethod("filterByCov",
          signature=c(object = "BSraw", minCov = "missing", global = "missing"),
          function(object) {
            .filterByCov(object, minCov = 10, global = FALSE)
          })

setMethod("filterByCov",
          signature=c(object = "BSraw", minCov = "numeric", global = "missing"),
          function(object, minCov) {
            .filterByCov(object, minCov, global = FALSE)
          })

setMethod("filterByCov",
          signature=c(object = "BSraw", minCov = "missing", global = "logical"),
          function(object, global) {
            .filterByCov(object, minCov = 10, global)
          })
