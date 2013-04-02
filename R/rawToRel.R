.rawToRel <- function(object){
  object.rel <- BSrel(colData = colData(object),
                      rowData = rowData(object),
                      methLevel = methReads(object) / totalReads(object)
                      )
  return(object.rel)
}

setMethod("rawToRel",
          signature=c(object = "BSraw"),
          .rawToRel)
