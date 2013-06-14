.readBismark <- function(files, colData) {

  if (nrow(colData) != length(files)) {
    stop("Row number of colData must equal length of files.")
  }

  methData = list()

  for (i in 1:length(files)) {

    cat(paste("Processing sample ", rownames(colData)[i], " ... \n", sep=""))

    bismark <- scan(files[i], skip=0, sep="\t",
                    what=list("character", integer(), NULL, NULL, integer(), integer()))
    
    methData[[i]] = GRanges(
              seqnames=bismark[[1]],
              ranges=IRanges(start=bismark[[2]], width=1),
              methylated=bismark[[5]],
              reads=bismark[[5]]+bismark[[6]])
    
    rm(bismark)
  }

  cat("Building BSraw object.\n")

  fData <- methData[[1]]
  
  if(length(methData) > 1){
    for(i in seq(along=methData)[-1]){
      fData <- unique(c(fData, methData[[i]]))
    }
  }
  elementMetadata(fData) <- NULL
  names(fData) <- as.character(1:length(fData))

  tReads <- matrix(integer(length = length(fData) * length(methData)), nrow=length(fData))
  mReads <- matrix(integer(length = length(fData) * length(methData)), nrow=length(fData))

  for(i in seq(along=methData)){
    mtch <- findOverlaps(fData, methData[[i]])
    mtch.m <-  as.matrix(mtch)
    ind <- mtch.m[, 1]
    tReads[ind, i] <- elementMetadata(methData[[i]])$reads
    mReads[ind, i] <- elementMetadata(methData[[i]])$methylated
  }

  colnames(tReads) <- rownames(colData)
  colnames(mReads) <- rownames(colData)
  rownames(tReads) <- names(fData)
  rownames(mReads) <- names(fData)

  rrbs = BSraw(
    colData = colData,
    rowData = fData,
    totalReads = tReads,
    methReads = mReads)
  
  return(rrbs)
}

setMethod("readBismark",
    signature=c(files="character", colData="DataFrame"),
    .readBismark)
          
setMethod("readBismark",
    signature=c(files="character", colData="data.frame"),
    function(files, colData) {
      colData = as(colData, "DataFrame")
      .readBismark(files, colData)
    })

setMethod("readBismark",
    signature=c(files="character", colData="character"),
    function(files, colData) {
      colData = DataFrame(row.names=colData)
      .readBismark(files, colData)
    })
