.annotateGRanges <- function(object, regions, name, regionInfo){

  regions <- unique(regions)
  MetaCol <- ncol(elementMetadata(object)) + 1
  if(missing(regionInfo)){
    ind <- overlapsAny(object, regions)
    elementMetadata(object)[, MetaCol] <- FALSE
    elementMetadata(object)[, MetaCol][ind] <- TRUE
    colnames(elementMetadata(object))[MetaCol] <- name
  }else{
    if(is(elementMetadata(regions)[,regionInfo], "factor")){
      elementMetadata(regions)[, regionInfo] <- as.character(elementMetadata(regions)[, regionInfo])
    }
    overl <- findOverlaps(query=object, subject=regions)
    matches <- as.data.frame(as.matrix(overl))

    helper <- rep(NA, length(object))
    if(nrow(matches) > 0){
      tab <- as.data.frame(table(matches$query)) # how many regions map to region in object
      ind <- tab$Freq > 1
      ind.many <- as.numeric(as.character(tab$Var1[ind]))
      matches.one <- matches[!is.element(matches$query, ind.many), ]

      helper[matches.one$query] <- elementMetadata(regions)[, regionInfo][matches.one$subject]
      
      for(i in ind.many){
        ind.reg <- matches$subject[matches$query == i]
        ids <- paste(elementMetadata(regions)[,regionInfo][ind.reg], collapse=",")
        helper[i] <- ids
      }
    }
    
    elementMetadata(object)[, MetaCol] <- helper
    colnames(elementMetadata(object))[MetaCol] <- name
  }
  return(object)
}

setMethod("annotateGRanges",
    signature=c(object="GRanges", region="GRanges", name="character", regionInfo="character"),
    .annotateGRanges)

setMethod("annotateGRanges",
    signature=c(object="GRanges", region="GRanges", name="character", regionInfo="integer"),
    .annotateGRanges)

setMethod("annotateGRanges",
    signature=c(object="GRanges", region="GRanges", name="character", regionInfo="missing"),
    .annotateGRanges)
