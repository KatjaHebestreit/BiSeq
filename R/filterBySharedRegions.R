.filterBySharedRegions_BSraw <- function(object, groups, perc.samples, no.samples, minCov){

  if(missing(groups)){
  
    if(!missing(perc.samples)){
      if(perc.samples <= 1 & perc.samples > 0){
        no.samples <- round(perc.samples*ncol(object))
      } else{
        stop("perc.samples has to be a numeric between 0 and 1.")
      }
    }

    if(no.samples > ncol(object)){
      stop("no.samples has to be smaller or equal to the number of samples of the object.")
    }

    if(any(no.samples < 1)){
      stop("no.samples has to be an integer greater or equal to 1.")
    }
    
    freq.pos <- rowSums(totalReads(object) >= minCov, na.rm=TRUE)
    ind <- freq.pos >= no.samples
    
  } else{
    
    g <- groups
    n.g <- as.vector(table(g))

    if(!missing(perc.samples)){
      if(all(perc.samples <= 1 & perc.samples > 0)){
        l.p <- length(perc.samples)
        if( (l.p == length(levels(g))) | l.p == 1 ){
          no.samples <- round( perc.samples * n.g )
        }
        else{
          stop("Length of perc.samples should be 1 or the same as the number of group levels.") 
        }
      } else{
        stop("perc.samples has to be a numeric between 0 and 1.")
      }
    }

    if(any(no.samples > n.g)){
      stop("no.samples has to be smaller or equal to the number of samples of the group.")
    }

    if(length(no.samples) > length(levels(g))){
      stop("Length of no.samples should be 1 or the same as the number of group levels.")
    }

    freq.pos <- matrix(nrow=nrow(object), ncol=length(n.g))
    for(j in seq(along=n.g)){
      ind.part <- g == levels(g)[j]
      part <- totalReads(object)[,ind.part]
      freq.pos[,j] <- rowSums(part >= minCov, na.rm=TRUE)
    }

    ind.m <- freq.pos >= no.samples
    ind <- apply(ind.m, 1, all)

  }

  rowData.new <- rowData(object)[ind]
  totalReads.new <- totalReads(object)[ind,]
  methReads.new <- methReads(object)[ind,]
  object.new <- BSraw(
                    colData = colData(object),
                    rowData = rowData.new,
                    totalReads = totalReads.new,
                    methReads = methReads.new
                    )
 
  return(object.new)
}

.filterBySharedRegions_BSrel <- function(object, groups, perc.samples, no.samples){

  if(missing(groups)){
  
    if(!missing(perc.samples)){
      if(perc.samples <= 1 & perc.samples > 0){
        no.samples <- round(perc.samples*ncol(object))
      } else{
        stop("perc.samples has to be a numeric between 0 and 1.")
      }
    }

    if(no.samples > ncol(object)){
      stop("no.samples has to be smaller or equal to the number of samples of the object.")
    }

    if(any(no.samples < 1)){
      stop("no.samples has to be an integer greater or equal to 1.")
    }
    
    freq.pos <- rowSums(!is.nan(methLevel(object)))
    ind <- freq.pos >= no.samples
    
  } else{
    
    g <- groups
    n.g <- as.vector(table(g))

    if(!missing(perc.samples)){
      if(all(perc.samples <= 1 & perc.samples > 0)){
        l.p <- length(perc.samples)
        if( (l.p == length(levels(g))) | l.p == 1 ){
          no.samples <- round( perc.samples * n.g )
        }
        else{
          stop("Length of perc.samples should be 1 or the same as the number of group levels.") 
        }
      } else{
        stop("perc.samples has to be a numeric between 0 and 1.")
      }
    }
    
    if(any(no.samples > n.g)){
      stop("no.samples has to be smaller or equal to the number of samples of the group.")
    }

    if(length(no.samples) > length(levels(g))){
      stop("Length of no.samples should be 1 or the same as the number of group levels.")
    }

    freq.pos <- matrix(nrow=nrow(object), ncol=length(n.g))
    for(j in seq(along=n.g)){
      ind.part <- g == levels(g)[j]
      part <- methLevel(object)[,ind.part]
      freq.pos[,j] <- rowSums(!is.nan(part))
    }

    ind.m <- freq.pos >= no.samples
    ind <- apply(ind.m, 1, all)

  }

  rowData.new <- rowData(object)[ind]
  methLevel.new <- methLevel(object)[ind,]
  object.new <- BSrel(
                    colData = colData(object),
                    rowData = rowData.new,
                    methLevel = methLevel.new
                    )
 
  return(object.new)
}


setMethod("filterBySharedRegions",
          signature=c(object = "BSraw", perc.samples = "missing", no.samples = "numeric", minCov = "numeric"),
          function(object, groups, no.samples, minCov) {
            .filterBySharedRegions_BSraw(object, groups, no.samples=no.samples, minCov=minCov)
          })

setMethod("filterBySharedRegions",
          signature=c(object = "BSraw", perc.samples = "numeric", no.samples = "missing", minCov = "numeric"),
          function(object, groups, perc.samples, minCov) {
            .filterBySharedRegions_BSraw(object, groups, perc.samples, minCov=minCov)
          })

setMethod("filterBySharedRegions",
          signature=c(object="BSraw", perc.samples = "numeric", no.samples = "missing", minCov = "missing"),
          function(object, groups, perc.samples) {
            .filterBySharedRegions_BSraw(object, groups, perc.samples, minCov=1)
          })

setMethod("filterBySharedRegions",
          signature=c(object="BSraw", perc.samples = "missing", no.samples = "numeric", minCov = "missing"),
          function(object, groups, no.samples) {
            .filterBySharedRegions_BSraw(object, groups, no.samples=no.samples, minCov=1)
          }
          )

setMethod("filterBySharedRegions",
          signature=c(object="BSraw", groups = "missing", perc.samples = "missing", no.samples = "numeric", minCov = "missing"),
          function(object, no.samples) {
            .filterBySharedRegions_BSraw(object, no.samples=no.samples, minCov=1)
          }
          )

setMethod("filterBySharedRegions",
          signature=c(object="BSraw", groups = "missing", perc.samples = "numeric", no.samples = "missing", minCov = "missing"),
          function(object, perc.samples) {
            .filterBySharedRegions_BSraw(object, perc.samples=perc.samples, minCov=1)
          }
          )

setMethod("filterBySharedRegions",
          signature=c(object="BSraw", groups = "missing", perc.samples = "missing", no.samples = "numeric", minCov = "numeric"),
          function(object, no.samples, minCov) {
            .filterBySharedRegions_BSraw(object, no.samples=no.samples, minCov=minCov)
          }
          )

setMethod("filterBySharedRegions",
          signature=c(object="BSraw", group = "missing", perc.samples = "numeric", no.samples = "missing", minCov = "numeric"),
          function(object, perc.samples, minCov) {
            .filterBySharedRegions_BSraw(object, perc.samples=perc.samples, minCov=minCov)
          }
          )

setMethod("filterBySharedRegions",
          signature=c(object = "BSraw", groups = "missing", perc.samples = "missing", no.samples = "missing", minCov = "numeric"),
          function(object, minCov) {
            .filterBySharedRegions_BSraw(object, perc.samples=1, minCov=minCov)
          })

setMethod("filterBySharedRegions",
          signature=c(object = "BSraw", groups = "missing", perc.samples = "missing", no.samples = "missing", minCov = "missing"),
          function(object) {
            .filterBySharedRegions_BSraw(object, perc.samples=1, minCov=1)
          })




setMethod("filterBySharedRegions",
          signature=c(object = "BSrel", perc.samples = "missing", no.samples = "numeric"),
          function(object, groups, no.samples) {
            .filterBySharedRegions_BSrel(object, groups, no.samples=no.samples)
          })

setMethod("filterBySharedRegions",
          signature=c(object = "BSrel", perc.samples = "numeric", no.samples = "missing"),
          function(object, groups, perc.samples) {
            .filterBySharedRegions_BSrel(object, groups, perc.samples)
          })

setMethod("filterBySharedRegions",
          signature=c(object="BSrel",  perc.samples = "numeric", no.samples = "missing"),
          function(object, groups, perc.samples) {
            .filterBySharedRegions_BSrel(object, groups, perc.samples)
          })

setMethod("filterBySharedRegions",
          signature=c(object="BSrel", perc.samples = "missing", no.samples = "numeric"),
          function(object, groups, no.samples) {
            .filterBySharedRegions_BSrel(object, groups, no.samples=no.samples)
          }
          )

setMethod("filterBySharedRegions",
          signature=c(object="BSrel", groups = "missing", perc.samples = "missing", no.samples = "numeric"),
          function(object, no.samples) {
            .filterBySharedRegions_BSrel(object, no.samples=no.samples)
          }
          )

setMethod("filterBySharedRegions",
          signature=c(object="BSrel", groups = "missing", perc.samples = "numeric", no.samples = "missing"),
          function(object, perc.samples) {
            .filterBySharedRegions_BSrel(object, perc.samples=perc.samples)
          }
          )

setMethod("filterBySharedRegions",
          signature=c(object = "BSrel", groups = "missing", perc.samples = "missing", no.samples = "missing"),
          function(object) {
            .filterBySharedRegions_BSrel(object, perc.samples=1)
          })
