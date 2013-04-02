setClass("BSrel",
         contains = "SummarizedExperiment"
         )

setClass("BSraw",
         contains = "SummarizedExperiment"
         )


setValidity("BSrel", function(object){
  if(length(assays(object)) != 1)
    return("The assays slot in BSrel object must be of length one.")
  if(!is(assay(object), "matrix"))
    return("The methLevel slot of an BSrel object must be a matrix.")
  if(mode(assay(object)) != "numeric")
    return("The methLevel matrix of an BSrel object must contain numeric data.")
}
            )

setValidity("BSraw", function(object){
  if(length(assays(object)) != 2)
    return("The assays slot in BSraw object must be of length two.")
  if(!(all( is.element(names(assays(object)), c("totalReads", "methReads")) )))
    return("The assays slot in BSraw object must contain totalReads and methReads.")
  if(!all( sapply(assays(object), class) == "matrix" ))
    return("The totalReads and methReads slots of an BSrel object must be matrices.")
  if(!all(sapply(assays(object), typeof) == "integer" ))
    return("The totalReads and methReads matrices of an BSraw object must contain integer data.")
}
            )

