setMethod("BSraw", signature(methReads="matrix",
                             totalReads="matrix",
                             rowRanges="GRanges"),
          function(methReads,
                   totalReads,
                   rowRanges,
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
                rowRanges = rowRanges,
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
            rowRanges.new <- sort(unique(c(rowRanges(x), rowRanges(y))))
            ind.match.x <- findOverlaps(rowRanges.new, rowRanges(x), select="first")
            ind.match.y <- findOverlaps(rowRanges.new, rowRanges(y), select="first")
            nr <- length(rowRanges.new)
            nc <- nrow(colData.new)
            methReads.new <- matrix(integer(length = nr*nc),
                                    ncol=nc,
                                    nrow=nr,
                                    dimnames=list(names(rowRanges.new), rownames(colData.new)))
            methReads.new[,] <- cbind(methReads(x)[ind.match.x,], methReads(y)[ind.match.y,])
            methReads.new[is.na(methReads.new)] <- 0L
            totalReads.new <- matrix(integer(length = nr*nc),
                                    ncol=nc,
                                    nrow=nr,
                                    dimnames=list(names(rowRanges.new), rownames(colData.new)))
            totalReads.new[,] <- cbind(totalReads(x)[ind.match.x,], totalReads(y)[ind.match.y,])
            totalReads.new[is.na(totalReads.new)] <- 0L
            z <- BSraw(colData = colData.new,
                       rowRanges = rowRanges.new,
                       methReads = methReads.new,
                       totalReads = totalReads.new
                       )
            return(z)
          }
          )
