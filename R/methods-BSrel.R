setMethod("BSrel", signature(methLevel="matrix",
                             rowData="GRanges"),
          function(methLevel,
                   rowData,
                   colData = DataFrame(row.names=colnames(methLevel)),
                   exptData = SimpleList(),
                   ...)
          {
            ssla <- new("ShallowSimpleListAssays",
                        data=SimpleList(
                          methLevel = methLevel))
            new("BSrel",
                assays = ssla,
                rowData = rowData,
                colData = colData,
                exptData = exptData)
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
            rowData.new <- sort(unique(c(rowData(x), rowData(y))))
            ind.match.x <- match(rowData.new, rowData(x))
            ind.match.y <- match(rowData.new, rowData(y))
            nr <- length(rowData.new)
            nc <- nrow(colData.new)
            methLevel.new <- matrix(numeric(length = nr*nc),
                                    ncol=nc,
                                    nrow=nr,
                                    dimnames=list(names(rowData.new), rownames(colData.new)))
            methLevel.new[,] <- cbind(methLevel(x)[ind.match.x,], methLevel(y)[ind.match.y,])
            methLevel.new[is.na(methLevel.new)] <- NaN
            z <- BSrel(colData = colData.new,
                       rowData = rowData.new,
                       methLevel = methLevel.new
                       )
            return(z)
          }
          )

