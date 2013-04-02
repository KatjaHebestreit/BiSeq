.predictMeth <- function(object, h, grid.dist, mc.cores){

  strand(object) <- "*"
  object <- sort(object)
  
  rowData.new <- GRanges()
  seqlevels(rowData.new) <- seqlevels(rowData(object))
  methLevels <- matrix(nrow=0, ncol=ncol(object))
  colnames(methLevels) <- colnames(object)

  l.chr <- table(seqnames(object))
  i.chr <- 0
  pb <- txtProgressBar(min=0, max=sum(l.chr), style=3)
  
  for(chr in unique(as.character(seqnames(object)))){ # for each chromosome
    ind.chr <- as.character(seqnames(object)) == chr
    object.chr <- object[ind.chr, ]

    # split object.chr by clusters:
    object.split <- split(object.chr, f = as.factor(elementMetadata(object.chr)$cluster.id), drop=TRUE)
    
    predict.cluster <- function(object.part, h=h, grid.dist=grid.dist){ # for each cluster
      
      clus.chr <- as.character(unique(seqnames(object.part)))
      clus.pos <- start(ranges(object.part))
      if(is.numeric(grid.dist)){
        clus.start <- min(start(ranges(object.part)))
        clus.end <- max(start(ranges(object.part)))
        clus.grid <- seq(clus.start, clus.end, grid.dist)
      } else{
        clus.grid <- clus.pos
      }
      rrbs.predict <- matrix(nrow=length(clus.grid), ncol=ncol(object))
      colnames(rrbs.predict) <- colnames(object)
      for(i in 1:ncol(object.part)){ # for each sample
        part <- object.part[,i]
        rrbs.predict[, i] <- binomLikelihoodSmooth(pred.pos = clus.grid,
                                                   pos = clus.pos,
                                                   m = methReads(part),
                                                   n = totalReads(part),
                                                   h = h)
      }
      f.help <- GRanges(seqnames=clus.chr, ranges=IRanges(start=clus.grid, end=clus.grid))
      elementMetadata(f.help)$cluster.id <- unique(elementMetadata(object.part)$cluster.id)
      seqlevels(f.help) <- seqlevels(rowData(object))

      return(list(f.help, rrbs.predict))
    }

    out.list <- mclapply(object.split, predict.cluster, h=h, grid.dist=grid.dist, mc.cores=mc.cores, mc.preschedule=FALSE)
    
    out1.pre <- lapply(out.list, function(x) x[[1]])
    
    out1 <- out1.pre[[1]]
    if(length(out1.pre) >= 2){
      for(k in 2:length(out1.pre)){
        out1 <- c(out1, out1.pre[[k]])
      }
    }
    out2 <- do.call("rbind", lapply(out.list, function(x) x[[2]]))
    rowData.new <- c(rowData.new, out1)
    methLevels <- rbind(methLevels, out2)

    i.chr <- i.chr + l.chr[chr]
    setTxtProgressBar(pb, value=i.chr)
  }
  close(pb)
  
  rownames(methLevels) <- 1:nrow(methLevels)
  names(rowData.new) <- 1:length(rowData.new)
  
  predictedMeth <- BSrel(colData=colData(object), 
                         rowData=rowData.new, 
                         methLevel=methLevels)
  
  return(predictedMeth)
}

setMethod("predictMeth",
          signature=c(object = "BSraw", h = "numeric", grid.dist = "numeric", mc.cores = "numeric"),
          .predictMeth)

setMethod("predictMeth",
          signature=c(object="BSraw", h = "numeric", grid.dist = "missing", mc.cores = "missing"),
          function(object, h) {
            .predictMeth(object, h, grid.dist=NULL, mc.cores = 1)
          })

setMethod("predictMeth",
          signature=c(object="BSraw", h = "numeric", grid.dist = "missing", mc.cores = "numeric"),
          function(object, h, mc.cores) {
            .predictMeth(object, h, grid.dist=NULL, mc.cores)
          })

setMethod("predictMeth",
          signature=c(object="BSraw", h = "missing", grid.dist = "missing", mc.cores = "numeric"),
          function(object, mc.cores) {
            .predictMeth(object, h=80, grid.dist=NULL, mc.cores=mc.cores)
          })

setMethod("predictMeth",
          signature=c(object="BSraw", h = "missing", grid.dist = "missing", mc.cores = "missing"),
          function(object) {
            .predictMeth(object, h = 80, grid.dist=NULL, mc.cores = 1)
          })
