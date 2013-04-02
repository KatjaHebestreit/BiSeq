.compareTwoSamples <- function(object, sample1, sample2, minDiff, max.dist){
  
  if(minDiff <= 0){
    stop("minDiff should be greater than 0.")
  }
  if(minDiff > 1){
    stop("minDiff should be smaller than 1.")
  }
  meth.diff <- data.frame(chr = as.character(seqnames(rowData(object))),
                          pos = start(ranges(rowData(object))),
                          meth.group1 = methLevel(object)[,sample1],
                          meth.group2 = methLevel(object)[,sample2]
                          )
  meth.diff$meth.diff <- meth.diff$meth.group1 - meth.diff$meth.group2
  meth.diff$direction <- sign(meth.diff$meth.diff)
  meth.diff$p.val <- -abs(meth.diff$meth.diff)
  dmrs <- findDMRs(test.out=meth.diff, alpha = -minDiff, max.dist = max.dist)
  return(dmrs[,c("median.meth.group1", "median.meth.group2", "median.meth.diff")])
}


setMethod("compareTwoSamples",
          signature=c(object = "BSrel", sample1 = "numeric", sample2 = "numeric", minDiff = "numeric", max.dist = "numeric"),
          .compareTwoSamples)

setMethod("compareTwoSamples",
          signature=c(object = "BSrel", sample1 = "character", sample2 = "character", minDiff = "numeric", max.dist = "numeric"),
          .compareTwoSamples)

setMethod("compareTwoSamples",
          signature=c(object = "BSraw", sample1 = "numeric", sample2 = "numeric", minDiff = "numeric", max.dist = "numeric"),
          function(object, sample1, sample2, minDiff, max.dist) {
            .compareTwoSamples(object=rawToRel(object), sample1, sample2, minDiff, max.dist)
          }
          )

setMethod("compareTwoSamples",
          signature=c(object = "BSraw", sample1 = "character", sample2 = "character", minDiff = "numeric", max.dist = "numeric"),
          function(object, sample1, sample2, minDiff, max.dist) {
            .compareTwoSamples(object=rawToRel(object), sample1, sample2, minDiff, max.dist)
          }
          )

