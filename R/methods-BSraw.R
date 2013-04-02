setMethod("BSraw", signature(methReads="matrix",
                             totalReads="matrix",
                             rowData="GRanges"),
          function(methReads,
                   totalReads,
                   rowData,
                   colData = DataFrame(row.names=colnames(methReads)),
                   exptData = SimpleList(),
                   ...)
          {
            ssla <- new("ShallowSimpleListAssays",
                        data=SimpleList(
                          totalReads = totalReads,
                          methReads = methReads))
            new("BSraw",
                assays = ssla,
                rowData = rowData,
                colData = colData,
                exptData = exptData)
          })

setMethod("totalReads", signature(object ="BSraw"),
          function(object){
            return(assays(object)$totalReads)
          }
          )

setReplaceMethod("totalReads", signature(object="BSraw", value="matrix"),
                 function(object, value){
                   assays(object)$totalReads <- value
                   return(object)
                 }
                 )

setMethod("methReads", signature(object ="BSraw"),
          function(object){
            return(assays(object)$methReads)
          }
          )

setReplaceMethod("methReads", signature(object="BSraw", value="matrix"),
                 function(object, value){
                   assays(object)$methReads <- value
                   return(object)
                 }
                 )

setMethod("combine", signature(x ="BSraw", y="BSraw"),
          function(x,y){
            if(any(is.element(colnames(x), colnames(y)))){
              stop("The BSraw objects to combine should not have samples in common!")
            }
            colData.new <- rbind(colData(x), colData(y))
            rowData.new <- sort(unique(c(rowData(x), rowData(y))))
            ind.match.x <- findOverlaps(rowData.new, rowData(x), select="first")
            ind.match.y <- findOverlaps(rowData.new, rowData(y), select="first")
            nr <- length(rowData.new)
            nc <- nrow(colData.new)
            methReads.new <- matrix(integer(length = nr*nc),
                                    ncol=nc,
                                    nrow=nr,
                                    dimnames=list(names(rowData.new), rownames(colData.new)))
            methReads.new[,] <- cbind(methReads(x)[ind.match.x,], methReads(y)[ind.match.y,])
            methReads.new[is.na(methReads.new)] <- 0L
            totalReads.new <- matrix(integer(length = nr*nc),
                                    ncol=nc,
                                    nrow=nr,
                                    dimnames=list(names(rowData.new), rownames(colData.new)))
            totalReads.new[,] <- cbind(totalReads(x)[ind.match.x,], totalReads(y)[ind.match.y,])
            totalReads.new[is.na(totalReads.new)] <- 0L
            z <- BSraw(colData = colData.new,
                       rowData = rowData.new,
                       methReads = methReads.new,
                       totalReads = totalReads.new
                       )
            return(z)
          }
          )
