.plotMeth <- function(object.raw, object.rel, region, col.lines, lwd.lines, col.points, ...) {

  object.raw <- subsetByOverlaps( object.raw, region )
  object.rel <- subsetByOverlaps( object.rel, region )

  f.raw <- start(ranges(rowData(object.raw)))
  t <- totalReads(object.raw)

  ind.cov <-  t[,1] > 0 
  pos.cov <- f.raw[ ind.cov ]
  max.cov <- max(t[ind.cov,1])
  min.cov <- min(t[ind.cov,1])
  alpha <- t[ind.cov,1] / max.cov

  args <- list(...)

  if (!is.element("xlim", names(args))) {
    args$xlim <- c(start(ranges(region)), end(ranges(region)))
  }
  if (!is.element("ylim", names(args))) {
    args$ylim <- c(0,1)
  }
  if (!is.element("xlab", names(args))) {
    args$xlab <- "Position"
  }
  if (!is.element("ylab", names(args))) {
    args$ylab <- "Methylation"
  }
  if (!is.element("main", names(args))) {
    args$main <- colnames(object.raw)[1]
  }

  if( missing(col.points) ) {
    args$col <- "blue"
  } else {
    args$col <- col.points
  }
  
  if( missing(col.lines) ) {
    col.lines <- "black"
  }

  if( missing(lwd.lines) ) {
    lwd.lines <- 1
  }

  args$x <- pos.cov
  args$y <- methReads(object.raw)[ind.cov,1] / totalReads(object.raw)[ind.cov,1]
  col.rgb <- col2rgb(args$col)
  args$bg <- rgb(red = col.rgb["red",1]/255,
                 blue = col.rgb["blue",1]/255,
                 green = col.rgb["green",1]/255,
                 alpha = alpha)
  args$pch <- 21

  do.call(plot, args)

  # one curve for each CpG cluster
  object.rel.split <- split(object.rel, elementMetadata(object.rel)$cluster.id)
  pos.rel.split <- lapply(object.rel.split, function(x) start(ranges(rowData(x))))

  for(i in seq(along=pos.rel.split)) {
    lines(pos.rel.split[[i]],
          methLevel(object.rel.split[[i]])[,1],
          col=col.lines,
          lwd=lwd.lines)
  }

  ind.min.cov <- which(t[ind.cov,1] == min.cov)[1]
  ind.max.cov <- which(t[ind.cov,1] == max.cov)[1]

  if (!is.element("cex", names(args))) {
    args$cex <- 1
  }

  legend("topright",
         pch=args$pch,
         col=args$col,
         pt.cex=args$cex,
         pt.bg=args$bg[c(ind.min.cov, ind.max.cov)],
         title="Coverage",
         legend=c(paste(min.cov, "x",sep=""),
           paste(max.cov, "x",sep="")))
  
  if(ncol(object.raw) != 1 | ncol(object.raw) != 1) {
    warning("Only the first sample was plotted.")
  }

}

setMethod("plotMeth",
    signature=c(object.raw="BSraw", object.rel="BSrel", region="GRanges"),
    function(object.raw, object.rel, region, col.lines, lwd.lines, col.points, ...) {
      .plotMeth(object.raw, object.rel, region, col.lines, lwd.lines, col.points, ...)
    })

