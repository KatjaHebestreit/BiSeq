.readBismark <- function(otFiles, obFiles, colData) {

  if (length(otFiles) != length(obFiles)) {
    stop("Arguments otFiles and obFiles must be vectors of equal length.")
  }
  if (nrow(colData) != length(otFiles)) {
    stop("Row number of colData must equal length of otFiles and obFiles.")
  }

  methData = list()

  for (i in 1:length(otFiles)) {

    cat(paste("Processing sample ", rownames(colData)[i], " ... ", sep=""))
    
#    dfT <- read.table(otFiles[i], skip=1, header=FALSE, sep="\t",
#                      stringsAsFactors=FALSE, comment.char="")
#    dfB <- read.table(obFiles[i], skip=1, header=FALSE, sep="\t",
#                      stringsAsFactors=FALSE, comment.char="")
#    dfT$Strand <- "+"
#    dfB$Strand <- "-"
    
    listT <- scan(otFiles[i], skip=1, sep="\t",
                  what=list(NULL, "character", "character", "numeric", NULL))
    listB <- scan(obFiles[i], skip=1, sep="\t",
                  what=list(NULL, "character", "character", "numeric", NULL))

    dfT <- data.frame(Meth=listT[[2]], Chr=listT[[3]], Pos=as.integer(listT[[4]]), Strand="+",
                      stringsAsFactors=FALSE)
    dfB <- data.frame(Meth=listB[[2]], Chr=listB[[3]], Pos=as.integer(listB[[4]]), Strand="-",
                      stringsAsFactors=FALSE)
    rm(listT, listB)
    
    df <- rbind(dfT, dfB)
    rm(dfT, dfB)

    df <- split.data.frame(df, df[,2])
    df <- lapply(df, function(x) {
      return(x[order(x[, 3]),])
    })
    counts <- list()
  
    for (n in names(df)) {
      t <- table(df[[n]][,3])
      ind <- match(as.integer(names(t)), df[[n]][,3])
      counts[[n]] <- data.frame(
          position = as.integer(names(t)),
          methylated = 0,
          reads = as.integer(t),
          chrom = n,
          strand = df[[n]][ind, "Strand"])
      t <- table(df[[n]][df[[n]][, 1] == "+", 3])
      if (length(t) > 0) {
        ind <- match(as.integer(names(t)), counts[[n]]$position)
        counts[[n]][ind, "methylated"] <- as.integer(t)
      }
    }
    rm(df)

    methData[[i]] = GRanges(
              seqnames=unlist(lapply(counts, function(x) {x$chrom})),
              ranges=IRanges(start=unlist(lapply(counts, function(x) {x$position})), width=1),
              strand=unlist(lapply(counts, function(x) {x$strand})),
              methylated=unlist(lapply(counts, function(x) {as.integer(x$methylated)})),
              reads=unlist(lapply(counts, function(x) {as.integer(x$reads)}))
              )
    methData[[i]] <- methData[[i]][order(methData[[i]])]
    rm(counts)

    cat("done.\n")
  }
  
  fData <- methData[[1]]
  
  if(length(methData) > 1){
    for(i in seq(along=methData)[-1]){
      fData <- union(fData, methData[[i]])
    }
  }

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
    signature=c(otFiles="character", obFiles="character", colData="DataFrame"),
    .readBismark)
          
setMethod("readBismark",
    signature=c(otFiles="character", obFiles="character", colData="data.frame"),
    function(otFiles, obFiles, colData) {
      colData = as(colData, "DataFrame")
      .readBismark(otFiles, obFiles, colData)
    })

setMethod("readBismark",
    signature=c(otFiles="character", obFiles="character", colData="character"),
    function(otFiles, obFiles, colData) {
      colData = DataFrame(row.names=colData)
      .readBismark(otFiles, obFiles, colData)
    })
