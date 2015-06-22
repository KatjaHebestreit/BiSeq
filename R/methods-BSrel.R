setMethod("BSrel", signature(methLevel="matrix",
                             rowRanges="GRanges"),
          function(methLevel,
                   rowRanges,
                   colData = DataFrame(row.names=colnames(methLevel)),
                   metadata = list(),
                   ...)
          {
            new("BSrel", SummarizedExperiment(
                assays = SimpleList(methLevel = methLevel),
                rowRanges = rowRanges,
                colData = colData,
                metadata = metadata))
          })

setMethod("methLevel", signature(object ="BSrel"),
          function(object){
            return(assays(object)$methLevel)
          }
          )

setReplaceMethod("methLevel", signature(object="BSrel", value="matrix"),
                 function(object, value){
                   assays(object)$methLevel <- value
                   return(object)
                 }
                 )


setMethod("combine", signature(x ="BSrel", y = "BSrel"),
          function(x,y){
            if(any(is.element(colnames(x), colnames(y)))){
              stop("The BSrel objects to combine should not have samples in common!")
            }
            colData.new <- rbind(colData(x), colData(y))
            rowRanges.new <- sort(unique(c(rowRanges(x), rowRanges(y))))
            ind.match.x <- match(rowRanges.new, rowRanges(x))
            ind.match.y <- match(rowRanges.new, rowRanges(y))
            nr <- length(rowRanges.new)
            nc <- nrow(colData.new)
            methLevel.new <- matrix(numeric(length = nr*nc),
                                    ncol=nc,
                                    nrow=nr,
                                    dimnames=list(names(rowRanges.new), rownames(colData.new)))
            methLevel.new[,] <- cbind(methLevel(x)[ind.match.x,], methLevel(y)[ind.match.y,])
            methLevel.new[is.na(methLevel.new)] <- NaN
            z <- BSrel(colData = colData.new,
                       rowRanges = rowRanges.new,
                       methLevel = methLevel.new
                       )
            return(z)
          }
          )

