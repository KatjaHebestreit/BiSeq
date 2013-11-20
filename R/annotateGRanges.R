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

      ids <- elementMetadata(regions)[, regionInfo]
      
      if(is(ids, "CompressedCharacterList")){
        ids.l <- sapply(ids, length)
        ids.char <- character(length=length(ids))
        ids.l.1 <- which(ids.l == 1)
        ids.l.n <- which(ids.l > 1)
        if(length(ids.l.1) > 1){
          ids.char[ids.l.1] <- unlist(ids[ids.l.1])
        }
        if(length(ids.l.n) > 1){
          ids.char[ids.l.n] <- sapply(ids[ids.l.n], function(x) paste(x, collapse = ","))
        }
        ids <- ids.char
      } 
      helper[matches.one$query] <- ids[matches.one$subject]
      
      for(i in ind.many){
        ind.reg <- matches$subject[matches$query == i]
        names.reg <- ids[ind.reg]
        names.reg <- sort(unique(names.reg))
        ids.i <- paste(names.reg, collapse=",")
        helper[i] <- ids.i
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
