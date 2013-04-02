.covBoxplots <- function(object, ...){
  
  t <- totalReads(object)
  ind.sites <-  !is.na(t) & t != 0 
  sites <- colSums(ind.sites)
  cov.median <- list()
  for(i in 1:ncol(object)){
    cov.median[[i]] <- t[ind.sites[,i],i]
  }
  names(cov.median) <- colnames(object)

  args = list(...)

  args$x <- cov.median

  if (!is.element("ylab", names(args))) {
    args$ylab="Coverage"
  }

  if (!is.element("names", names(args))) {
    args$names=colnames(object)
  }
  
  do.call(boxplot, args)
  
}


setMethod("covBoxplots",
          signature=c(object = "BSraw"),
          .covBoxplots)
