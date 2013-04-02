.writeBED_BSraw <- function(object, name, file) {

  if (ncol(object) != length(name) | ncol(object) != length(file)) {
    stop("Character vectors name and file must have the same length as the object has samples.")
  }

  colFunc <- colorRamp(colors=c("green", "black", "red"))
  strand(object) <- "*"
  object = sort(object)
  
  for (i in 1:ncol(object)) {
    object.i <- object[, i]
    ind.cov <- totalReads(object.i) != 0
    if (sum(ind.cov) > 1) {
      object.i <- object.i[ind.cov, ]
      bed <- rowData(object.i)
      elementMetadata(bed)$score <- methReads(object.i) / totalReads(object.i)
      elementMetadata(bed)$name <- totalReads(object.i)
      m <- colFunc(elementMetadata(bed)$score) / 255
      elementMetadata(bed)$itemRgb <- rgb(m[,1], m[,2], m[,3])

      bed <- as(as(bed, "RangedData"), "UCSCData")
      bed@trackLine@name <- paste("\"", name[i], "\"", sep="")
      export.ucsc(bed, con=file[i], subformat="bed")

    } else { # no CpGs are covered
      warning(paste(name[i], ": No CpG is covered - no bed file has been written", sep=""))
    }
  }
}

.writeBED_BSrel <- function(object, name, file) {

  if (ncol(object) != length(name) | ncol(object) != length(file)) {
    stop("Character vectors name and file must have the same length as the object has samples.")
  }

  colFunc <- colorRamp(colors=c("green", "black", "red"))
  object = sort(object)
  
  for (i in 1:ncol(object)) {
    object.i <- object[, i]
    ind.cov <- !is.na(methLevel(object.i))
    if (sum(ind.cov) > 1) {
      object.i <- object.i[ind.cov, ]
      bed <- rowData(object.i)
      elementMetadata(bed)$score <- methLevel(object.i)
      elementMetadata(bed)$name <- "\"\""
      m <- colFunc(elementMetadata(bed)$score) / 255
      elementMetadata(bed)$itemRgb <- rgb(m[,1], m[,2], m[,3])

      bed <- as(as(bed, "RangedData"), "UCSCData")
      bed@trackLine@name <- paste("\"", name[i], "\"", sep="")
      export.ucsc(bed, con=file[i], subformat="bed")

    } else { # no CpGs are covered
      warning(paste(name[i], ": No CpG is covered - no bed file has been written", sep=""))
    }  
  }
}

setMethod("writeBED",
    signature=c(object="BSraw", name="character", file="character"),
    .writeBED_BSraw)

setMethod("writeBED",
    signature=c(object="BSraw", name="character", file="missing"),
    function(object, name) {
      .writeBED_BSraw(object, name=name, file=paste(colnames(object), ".bed", sep=""))
    })

setMethod("writeBED",
    signature=c(object="BSraw", name="missing", file="character"),
    function(object, file) {
      .writeBED_BSraw(object, name=colnames(object), file=file)
    })

setMethod("writeBED",
    signature=c(object="BSraw", name="missing", file="missing"),
    function(object) {
      .writeBED_BSraw(object, name=colnames(object), file=paste(colnames(object), ".bed", sep=""))
    })



setMethod("writeBED",
    signature=c(object="BSrel", name="character", file="character"),
    .writeBED_BSrel)

setMethod("writeBED",
    signature=c(object="BSrel", name="character", file="missing"),
    function(object, name) {
      .writeBED_BSrel(object, name=name, file=paste(colnames(object), ".bed", sep=""))
    })

setMethod("writeBED",
    signature=c(object="BSrel", name="missing", file="character"),
    function(object, file) {
      .writeBED_BSrel(object, name=colnames(object), file=file)
    })

setMethod("writeBED",
    signature=c(object="BSrel", name="missing", file="missing"),
    function(object) {
      .writeBED_BSrel(object, name=colnames(object), file=paste(colnames(object), ".bed", sep=""))
    })
