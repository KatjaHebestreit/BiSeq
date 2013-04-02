.summarizeRegions_BSraw <- function(object, regions, outputAll=FALSE) {

  if (is.null(names(regions))) {
    names(regions) = as.character(1:length(regions))
  }
  
  totalReads = matrix(0L, nrow=length(regions), ncol=ncol(object),
    dimnames=list(names(regions), colnames(object)))
  methReads = matrix(0L, nrow=length(regions), ncol=ncol(object),
    dimnames=list(names(regions), colnames(object)))

  for (n in levels(seqnames(object))) {
    ind = which(seqnames(regions) == n)
    objectChr = object[as.character(seqnames(object)) == n,  ]
    cpgs = ranges(rowData(objectChr))

    ov = as.data.frame(as.matrix(findOverlaps(ranges(regions[ind, ]), cpgs)))
    m = split.data.frame(ov, ov$query)
    
    if (length(m) > 0) {
      for (s in colnames(objectChr)) {
        totalReads[ind[as.integer(names(m))], s] = as.integer(sapply(m, function(x) {
          return(sum(totalReads(objectChr)[x$subject, s]))
        }))
        methReads[ind[as.integer(names(m))], s] = as.integer(sapply(m, function(x) {
          return(sum(methReads(objectChr)[x$subject, s]))
        }))
      }
    }
  }

  rrbs = BSraw(colData=colData(object), 
               rowData=regions,
               totalReads=totalReads, 
               methReads=methReads)

  if (outputAll == FALSE) {
    ind = apply(totalReads(rrbs), 1, sum) != 0
    rrbs = rrbs[ind, ]
  }

  return(rrbs)
}



.summarizeRegions_BSrel <- function(object, regions, outputAll=FALSE) {

  if (is.null(names(regions))) {
    names(regions) = as.character(1:length(regions))
  }
  
  methLevel = matrix(0L, nrow=length(regions), ncol=ncol(object),
    dimnames=list(names(regions), colnames(object)))

  for (n in levels(seqnames(object))) {
    ind = which(seqnames(regions) == n)
    objectChr = object[as.character(seqnames(object)) == n,  ]
    cpgs = ranges(rowData(objectChr))

    ov = as.data.frame(as.matrix(findOverlaps(ranges(regions[ind, ]), cpgs)))
    m = split.data.frame(ov, ov$query)
    
    if (length(m) > 0) {
      for (s in colnames(object)) {
        methLevel[ind[as.integer(names(m))], s] = as.integer(sapply(m, function(x) {
          return(mean(methLevel(objectChr)[x$subject, s], na.rm=TRUE))
        }))
      }
    }
  }

  rrbs = BSrel(colData=colData(object), 
               rowData=regions,
               methLevel=methLevel)

  if (outputAll == FALSE) {
    ind = !apply(is.na(methLevel(rrbs)), 1, all)
    rrbs = rrbs[ind, ]
  }

  return(rrbs)
}



setMethod("summarizeRegions",
  signature=c(object="BSraw", regions="GRanges", outputAll="logical"),
  .summarizeRegions_BSraw)

setMethod("summarizeRegions",
  signature=c(object="BSraw", regions="GRanges", outputAll="missing"),
  function(object, regions) {
    .summarizeRegions_BSraw(object, regions, outputAll=FALSE)
  })

setMethod("summarizeRegions",
  signature=c(object="BSrel", regions="GRanges", outputAll="logical"),
  .summarizeRegions_BSrel)

setMethod("summarizeRegions",
  signature=c(object="BSrel", regions="GRanges", outputAll="missing"),
  function(object, regions) {
    .summarizeRegions_BSrel(object, regions, outputAll=FALSE)
  })
