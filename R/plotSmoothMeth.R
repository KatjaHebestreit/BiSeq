.plotSmoothMeth <- function(object.rel, region, groups, group.average, ...) {

  object.rel <- sort(object.rel)
  object.rel <- subsetByOverlaps( object.rel, region + 1000 )

  if(nrow(object.rel) == 0) {
    stop("No methylation data in the given region.")
  }
  
  if (missing(groups)) {
    groups <- factor(rep(1, ncol(object.rel)))
  } else {
    groups <- factor(as.character(groups))
  }

  n.g <- length(levels(groups))
  
  if(group.average) {
    object.groups <- list()
    for( i in seq(along=levels(groups)) ) {
      object.groups[[i]] <- object.rel[, groups == levels(groups)[i]]
    }
    meth.ave <- sapply(object.groups, function(x) rowMeans(methLevel(x), na.rm=TRUE))
    colnames(meth.ave) <- levels(groups)
    mLevel <- meth.ave
  } else {
    mLevel <- methLevel(object.rel)
  }
  n <- ncol(mLevel)

  if(is.element("cluster.id", colnames(elementMetadata(object.rel)))){
  # one curve for each CpG cluster
    object.split <- split(object.rel, elementMetadata(object.rel)$cluster.id)
    pos.rel.split <- lapply(object.split, function(x) start(x))
  } else {
    object.split <- list(object.rel)
    pos.rel.split <- list(start(object.rel))
  }

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
  
  # set col  
  if (!is.element("col", names(args))) {
    if (group.average) {
      cols <- rainbow(n.g)
    } else {
      cols <- rainbow(n.g)[ as.numeric(groups) ]
    }
  } else {
    cols <- args$col
    if (group.average) {
      if (length(cols) != n.g) {
        cols <- rep(cols, n.g)
      }
    } else {
      if (length(cols) != n.g) {
        if (length(cols) != n) {
          cols <- rep(cols, n)
        }
      } else {
        cols <- cols[ as.numeric(groups) ]
      }
    }
  }

  if (!is.element("lty", names(args))) {
    if (group.average) {
      lty <- rep(1, n.g)
    } else {
      lty <- rep(1, n)
    }
  } else {
    lty <- args$lty
    if (group.average) {
      if (length(lty) != n.g) {
        lty <- rep(lty, n.g)
      }
    } else {
      if (length(lty) != n.g) {
        if (length(lty) != n) {
          lty <- rep(lty, n)
        }
      } else {
        lty <- lty[ as.numeric(groups) ]
      }
    }
  }

  if (!is.element("lwd", names(args))) {
    if (group.average) {
      lwd <- rep(1, n.g)
    } else {
      lwd <- rep(1, n)
    }
  } else {
    lwd <- args$lwd
    if (group.average) {
      if (length(lwd) != n.g) {
        lwd <- rep(lwd, n.g)
      }
    } else {
      if (length(lwd) != n.g) {
        if (length(lwd) != n) {
          lwd <- rep(lwd, n)
        }
      } else {
        lwd <- lwd[ as.numeric(groups) ]
      }
    }
  }
  

  args$x <- 1

  do.call(plot, args)

  # along samples / groups
  for(i in 1:n) {
    # along CpG clusters
    for(j in seq(along=pos.rel.split)) {
      lines(pos.rel.split[[j]],
            methLevel(object.split[[j]])[,i],
            col = cols[i],
            lty = lty[i],
            lwd = lwd[i])
    }
  }
}

setMethod("plotSmoothMeth",
    signature=c(object.rel="BSrel", region="GRanges", group.average="logical"),
    function(object.rel, region, groups, group.average, ...) {
      .plotSmoothMeth(object.rel, region=region, groups=groups, group.average=group.average, ...)
    })


setMethod("plotSmoothMeth",
    signature=c(object.rel="BSrel", region="GRanges", group.average="missing"),
    function(object.rel, region, groups, ...) {
      .plotSmoothMeth(object.rel, region=region, groups=groups, group.average=FALSE, ...)
    })
