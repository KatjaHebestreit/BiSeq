.clusterSites <- function(object, groups, perc.samples, min.sites, max.dist, mc.cores, ...){

  if(missing(groups)){
    object <- filterBySharedRegions(object, perc.samples=perc.samples, ...)
  } else {
    object <- filterBySharedRegions(object, groups=groups, perc.samples=perc.samples, ...)
  }
  strand(object) <- "*"
  object <- sort(object)
  
  rowData.new <- new("GRanges")
  totalReads.new <- matrix(nrow=0, ncol=ncol(object))
  colnames(totalReads.new) <- colnames(object)
  methReads.new <- matrix(nrow=0, ncol=ncol(object))
  colnames(methReads.new) <- colnames(object)
  
  # split rrbs object by chromosomes:
  rrbs.split <- split(object, f = as.factor(seqnames(object)), drop = TRUE)
  
  clus.func <- function(x, max.dist, min.sites) {

    pos <- start(x)

    pos.dist <- pos[-1] - pos[-length(pos)] # distance between CpG to the subsequent CpG
    pos.close <- pos.dist <= max.dist

    clus.break <- which(!pos.close)

    cluster <- data.frame(start = numeric(length=length(clus.break)),
                          end = numeric(length=length(clus.break)),
                          size = numeric(length=length(clus.break)))

    if(length(clus.break) > 0){
      for(i in seq(along=clus.break)){ 
        break.prev <- clus.break[i-1]
        break.i <- clus.break[i]
        if(i == 1) {
          cluster$start[1] <- pos[1]
          cluster$end[1] <- pos[break.i]
          cluster$size[1] <- break.i
        } else {
          cluster$start[i] <- pos[break.prev + 1]
          cluster$end[i] <- pos[break.i]
          cluster$size[i] <- break.i - break.prev
        }
      }
      if(break.i != length(pos)-1){
        start <- pos[break.i + 1]
        end <- pos[length(pos)]
        size <- length(pos) - break.i
        cluster <- rbind(cluster, data.frame(start, end, size))
      }

      cluster.keep <- cluster$size >= min.sites
      cluster <- cluster[cluster.keep,]

      if(nrow(cluster) > 0){
        cluster.gr <- GRanges(seqnames = unique(seqnames(x)),
                               ranges = IRanges(start=cluster$start, end=cluster$end))
        elementMetadata(cluster.gr)$cluster.id <- paste(seqnames(cluster.gr), "_", seq(along=cluster$start), sep="")
      
        mtch <- findOverlaps(query = cluster.gr, subject = rowData(x))
        mtch.m <-  as.matrix(mtch)
        rowData.clust <- rowData(x)[mtch.m[,2],]
        elementMetadata(rowData.clust)$cluster.id <- elementMetadata(cluster.gr)$cluster.id[mtch.m[,1]]
        totalReads.clust <- totalReads(x)[mtch.m[,2],]
        methReads.clust <- methReads(x)[mtch.m[,2],]
        return( BSraw(colData = colData(x),
                      rowData = rowData.clust,
                      totalReads = totalReads.clust,
                      methReads = methReads.clust)
               )
      }
    }

  }

  rrbs.clust <- mclapply(rrbs.split, clus.func, max.dist= max.dist, min.sites = min.sites,
                         mc.cores = mc.cores, mc.preschedule = FALSE) # for each chromosome
  
  ind.null <- sapply(rrbs.clust, is.null)
  rrbs.clust <- rrbs.clust[!ind.null]

  names(rrbs.clust) <- NULL

  rowData.clust <- do.call(c, lapply(rrbs.clust, function(x) rowData(x) ))

  totalReads.clust <- do.call(rbind, lapply(rrbs.clust, function(x) totalReads(x) ))

  methReads.clust <- do.call(rbind, lapply(rrbs.clust, function(x) methReads(x) ))

  rownames(totalReads.clust) <- names(rowData.clust)
  rownames(methReads.clust) <- names(rowData.clust)

  object.clust <- BSraw(colData=colData(object),
                        rowData=rowData.clust,
                        totalReads=totalReads.clust,
                        methReads=methReads.clust)
  return(object.clust)
}


setMethod("clusterSites",
          signature=c(object = "BSraw", perc.samples = "numeric", min.sites = "numeric", max.dist = "numeric", mc.cores = "numeric"),
          .clusterSites)

setMethod("clusterSites",
          signature=c(object = "BSraw", groups = "missing", perc.samples = "missing", min.sites = "missing", max.dist = "missing", mc.cores = "missing"),
          function(object, ...){
            .clusterSites(object, perc.samples = 1, min.sites = 20, max.dist = 100, mc.cores = 1, ...)
          })

setMethod("clusterSites",
          signature=c(object = "BSraw", perc.samples = "numeric", min.sites = "missing", max.dist = "missing", mc.cores = "missing"),
          function(object, groups, perc.samples, ...){
            .clusterSites(object, groups, perc.samples, min.sites = 20, max.dist = 100, mc.cores = 1, ...)
          })

setMethod("clusterSites",
          signature=c(object = "BSraw",  perc.samples = "numeric", min.sites = "numeric", max.dist = "missing", mc.cores = "missing"),
          function(object, groups, perc.samples, min.sites, ...){
            .clusterSites(object, groups, perc.samples, min.sites, max.dist = 100, mc.cores = 1, ...)
          })

setMethod("clusterSites",
          signature=c(object = "BSraw",  perc.samples = "numeric", min.sites = "numeric", max.dist = "numeric", mc.cores = "missing"),
          function(object, groups, perc.samples, min.sites, max.dist, ...){
            .clusterSites(object, groups, perc.samples, min.sites, max.dist, mc.cores = 1, ...)
          })

setMethod("clusterSites",
          signature=c(object = "BSraw", groups = "missing", perc.samples = "numeric", min.sites = "numeric", max.dist = "numeric", mc.cores = "numeric"),
          function(object, perc.samples, min.sites, max.dist, mc.cores, ...){
            .clusterSites(object, perc.samples = perc.samples, min.sites = min.sites, max.dist = max.dist, mc.cores = mc.cores, ...)
          })

setMethod("clusterSites",
          signature=c(object = "BSraw", groups = "missing", perc.samples = "numeric", min.sites = "numeric", max.dist = "numeric", mc.cores = "missing"),
          function(object, perc.samples, min.sites, max.dist, mc.cores, ...){
            .clusterSites(object, perc.samples = perc.samples, min.sites = min.sites, max.dist = max.dist, mc.cores = 1, ...)
          })

setMethod("clusterSites",
          signature=c(object = "BSraw", perc.samples = "numeric", min.sites = "missing", max.dist = "numeric", mc.cores = "missing"),
          function(object, groups, perc.samples, max.dist, ...){
            .clusterSites(object, groups, perc.samples, min.sites = 20, max.dist = max.dist, mc.cores = 1, ...)
          })


