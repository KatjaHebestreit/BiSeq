.findDMRs <- function(test.out, alpha, max.dist, diff.dir){   

  test.out <- test.out[!is.na(test.out$p.val),]
  test.out$direction[test.out$meth.diff < 0] <- -1
  test.out$direction[test.out$meth.diff > 0] <- 1
  chromosomes <- unique(test.out$chr)
  clustered.pvals.l <- list()
  for(i in seq(along=chromosomes)){  # each chromosome
    chr <- chromosomes[i]
    ind <- which(test.out$chr == chr)
    test.out.chr <- test.out[ind,]

    signif.unordered <- test.out.chr[test.out.chr$p.val <= alpha,]
    signif <- signif.unordered[order(signif.unordered$pos),]

    if(nrow(signif) > 1) {

      pos <- signif$pos

      if(diff.dir){
        dir <- signif$direction
      } else{
        dir <- rep(1, length(pos))
      }
      
      pos.dist <- (pos[-1] - pos[-length(pos)]) * dir[-1] * dir[-length(pos)]  # how close is CpG from the subsequent CpG? multiplied with direction of methylation difference
      pos.close <- pos.dist <= max.dist & pos.dist > 0 # CpGs have to lie cloes to each other and methylation difference must have the same direction
      
      clus.break <- which(!pos.close)
      l <- length(clus.break)

      if(l > 0) {
        
        cluster <- data.frame(chr = rep(chr, l+1),
                              start = numeric(length=l+1),
                              end = numeric(length=l+1),
                              median.p = numeric(length=l+1),
                              median.meth.group1 = numeric(length=l+1),
                              median.meth.group2 = numeric(length=l+1),
                              median.meth.diff = numeric(length=l+1)
                              )
      
        for(b in seq(l)){ 
          break.prev <- clus.break[b-1]
          break.b <- clus.break[b]
          if(b == 1) {
            cluster$start[1] <- pos[1]
            cluster$end[1] <- pos[break.b]
            cluster$median.p[1] <- median(signif$p.val[1:break.b])
            if(is.element("meth.group1", colnames(signif))) {
              cluster$median.meth.group1[1] <- median(signif$meth.group1[1:break.b])
              cluster$median.meth.group2[1] <- median(signif$meth.group2[1:break.b])
            }
            cluster$median.meth.diff[1] <- median(signif$meth.diff[1:break.b])
          } else {
            cluster$start[b] <- pos[break.prev + 1]
            cluster$end[b] <- pos[break.b]
            cluster$median.p[b] <- median(signif$p.val[ (break.prev + 1) : break.b ])
            if(is.element("meth.group1", colnames(signif))) {
              cluster$median.meth.group1[b] <- median(signif$meth.group1[ (break.prev + 1) : break.b ])
              cluster$median.meth.group2[b] <- median(signif$meth.group2[ (break.prev + 1) : break.b ])
            }
            cluster$median.meth.diff[b] <- median(signif$meth.diff[ (break.prev + 1) : break.b ])
          }
        }
        cluster$start[l+1] <- pos[clus.break[l] + 1]
        cluster$end[l+1] <- pos[length(pos)]
        cluster$median.p[l+1] <- median(signif$p.val[ (clus.break[l] + 1) : length(pos) ])
        if(is.element("meth.group1", colnames(signif))) {
          cluster$median.meth.group1[l+1] <- median(signif$meth.group1[ (clus.break[l] + 1) : length(pos) ])
          cluster$median.meth.group2[l+1] <- median(signif$meth.group2[ (clus.break[l] + 1) : length(pos) ])
        }
        cluster$median.meth.diff[l+1] <- median(signif$meth.diff[ (clus.break[l] + 1) : length(pos) ])

      }

      if(l == 0) {  ## all CpGs in one cluster
        cluster <- data.frame(chr = signif$chr[1],
                              start = signif$pos[1],
                              end = rev(signif$pos)[1],
                              median.p = median(signif$p.val),
                              median.meth.group1 = median(signif$meth.group1),
                              median.meth.group2 = median(signif$meth.group2),
                              median.meth.diff = median(signif$meth.diff)
                              )
        
      }
      clustered.pvals.l[[i]] <- cluster
    }
  
    if(nrow(signif) == 1) {
      cluster <- data.frame(chr = signif$chr,
                            start = signif$pos,
                            end = signif$pos,
                            median.p = signif$p.val,
                            median.meth.group1 = signif$meth.group1,
                            median.meth.group2 = signif$meth.group2,
                            median.meth.diff = signif$meth.diff
                            )
      clustered.pvals.l[[i]] <- cluster
    }
    
  }

  ind.null <- sapply(clustered.pvals.l, is.null)
  clustered.pvals.l <- clustered.pvals.l[!ind.null]
  
  cp <- sapply(clustered.pvals.l, function(x) {
    y <- GRanges(seqnames = x$chr,
                 ranges = IRanges(start = x$start, end = x$end)
                 )
    elementMetadata(y)$median.p <- x$median.p
    elementMetadata(y)$median.meth.group1 <- x$median.meth.group1
    elementMetadata(y)$median.meth.group2 <- x$median.meth.group2
    elementMetadata(y)$median.meth.diff <- x$median.meth.diff
    return(y)
  })
  
  clustered.pvals <- do.call(c, cp)
  
  return(sort(clustered.pvals))
}
    

setMethod("findDMRs",
          signature=c(test.out = "data.frame", alpha = "numeric", max.dist="numeric", diff.dir="logical"),
          .findDMRs)

setMethod("findDMRs",
          signature=c(test.out = "data.frame", alpha = "numeric", max.dist="numeric", diff.dir="missing"),
          function(test.out, alpha, max.dist) {
            .findDMRs(test.out, alpha, max.dist, diff.dir = TRUE)
          })

setMethod("findDMRs",
          signature=c(test.out = "data.frame", alpha = "missing", max.dist="numeric", diff.dir="missing"),
          function(test.out, max.dist) {
            .findDMRs(test.out, alpha = Inf, max.dist, diff.dir = TRUE)
          })

setMethod("findDMRs",
          signature=c(test.out = "data.frame", alpha = "missing", max.dist="numeric", diff.dir="logical"),
          function(test.out, max.dist, diff.dir) {
            .findDMRs(test.out, alpha = Inf, max.dist, diff.dir)
          })
