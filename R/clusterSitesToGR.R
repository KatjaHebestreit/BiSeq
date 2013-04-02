.clusterSitesToGR <- function(object){

  fData.split <- split(rowData(object), seqnames(rowData(object)), drop=TRUE)
  clusters <- lapply(fData.split, function(x){
    cluster.ids <- unique(elementMetadata(x)$cluster.id)
    df <- data.frame(chr = unique(seqnames(x)),
                     start = numeric(length=length(cluster.ids)),
                     end =  numeric(length=length(cluster.ids)),
                     cluster.id = cluster.ids)
    for(i in seq(along=cluster.ids)){
      id <- cluster.ids[[i]]
      ind <- which(elementMetadata(x)$cluster.id == id)
      start <- min(start(ranges(x))[ind])
      end <- max(start(ranges(x))[ind])
      df$start[i] <- start
      df$end[i] <- end
    }
    return(df)
  }
                     )

  clusters.df <- do.call(rbind, clusters)

  gr <- GRanges(seqnames=clusters.df$chr,
                ranges=IRanges(start=clusters.df$start, end=clusters.df$end))
  elementMetadata(gr)$cluster.id <- clusters.df$cluster.id

  return(gr)
}

setMethod("clusterSitesToGR",
          signature=c(object = "BSraw"),
          .clusterSitesToGR)

setMethod("clusterSitesToGR",
          signature=c(object = "BSrel"),
          .clusterSitesToGR)
