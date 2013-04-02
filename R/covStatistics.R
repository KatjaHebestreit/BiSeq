.covStatistics.rel <- function(object){
  sites <- apply(methLevel(object), 2, function(x) sum(!is.na(x)))
  return(list(Covered_CpG_sites = sites))
}

.covStatistics.raw <- function(object){
  t <- totalReads(object)
  ind.sites <-  !is.na(t) & t != 0 
  sites <- colSums(ind.sites)
  cov.median <- integer(length=ncol(object))
  names(cov.median) <- colnames(object)
  for(i in 1:ncol(object)){
    cov.median[i] <- median(t[ind.sites[,i],i])
  }
  return(list(Covered_CpG_sites = sites,
              Median_coverage = cov.median))
}

setMethod("covStatistics",
    signature=c(object = "BSraw"),
    .covStatistics.raw)

setMethod("covStatistics",
    signature=c(object = "BSrel"),
    .covStatistics.rel)
