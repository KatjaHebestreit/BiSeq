.globalTest <- function(response, alternative, ...) {

  args <- list(...)

  if (is.element("subsets", names(args)) && is(args$subsets, "GRanges")) {
    subsets <- as.list(findOverlaps(args$subsets, alternative))
    if (!is.null(names(args$subsets))) {
      names(subsets) <- names(args$subsets)
    }
    args$subsets <- subsets
  }

  eSet <- ExpressionSet(methLevel(alternative),
                        phenoData=AnnotatedDataFrame(as.data.frame(colData(alternative))))

  args$response <- response
  args$alternative <- eSet
  return(do.call(gt, args))
}


setMethod("globalTest",
    signature=c(response="ANY", alternative="BSrel"),
    .globalTest)
