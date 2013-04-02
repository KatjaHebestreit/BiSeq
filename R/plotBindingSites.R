.plotBindingSites <- function(object, regions, width, groups, quantiles, bandwidth, ...) {

  if (missing(width)) {
    width = max(IRanges:::width(regions))
  }
  if (missing(groups)) {
    groups = factor(rep(1, ncol(object)))
  } else {
    groups = factor(as.character(groups))
  }
  if (missing(bandwidth)) {
    bandwidth = width/8
  }
  if (missing(quantiles)) {
    if (ncol(object) > 4) {
      quantiles = c(0.25, 0.75)
    } else {
      quantiles = c()
    }
  }

  smooth.width = width+(2*ceiling(bandwidth))
  
  regions = resize(regions, width=smooth.width, fix="center")

  relPos <- numeric()
  methData <- numeric()
  matchMatrix <- as.matrix(findOverlaps(regions, rowData(object)))
  if (nrow(matchMatrix) == 0) {
    stop("No CpG sites within the given regions.")
  }
  hitInds <- unique(matchMatrix[, 1])

  pb.max <- length(hitInds) + 0.2*length(hitInds)
  pb <- txtProgressBar(min=0, max=pb.max, style=3)
  progBar <- 1

  for (i in hitInds) {
    setTxtProgressBar(pb, value=progBar)
    progBar <- progBar + 1
    ind = matchMatrix[matchMatrix[, 1] == i, 2]
    relPos = c(relPos, start(rowData(object)[ind]) - start(regions[i]) + 1)
    methData = rbind(methData, methLevel(object)[ind, ])
  }
  methData = split(as.data.frame(methData), factor(relPos))
  methData = lapply(methData, function(x) {
    apply(x, 2, mean, na.rm=TRUE) # mean better than median??
  })
  setTxtProgressBar(pb, value=pb.max)
  close(pb)
  
  relPos = as.numeric(names(methData))
  methData = as.matrix(do.call(rbind, args=methData))
  # all(row.names(methData) == as.character(relPos)
  relPos = relPos - (smooth.width / 2)
  x = relPos[1]:relPos[length(relPos)]


  sData = list()
  # for each group
  for (i in 1:length(levels(groups))) {  
    mData = t(apply(methData[, groups == levels(groups)[i]], 1, quantile, probs=c(0.5, quantiles), na.rm=TRUE))
    notNA = !apply(mData, 1, function(r) {any(is.na(r))})    
    sData[[i]] = apply(mData[notNA, ], 2, function(y) {
        ksmooth(x=relPos[notNA], y=y, kernel="normal", bandwidth=bandwidth,
            x.points=x)$y
    })
  }

  
  # prepare plot of first group
  args = list(...)

  args$x = x
  args$y = sData[[1]][, 1]

  # set ylim
  if (!is.element("ylim", names(args))) {
    args$ylim = c(max(0, min(sapply(sData, min, na.rm=TRUE), na.rm=TRUE) - 0.05),
        min(1, max(sapply(sData, max, na.rm=TRUE), na.rm=TRUE) + 0.05))
  }

  # set xlim
  if (!is.element("xlim", names(args))) {
    args$xlim = c(-round(width/2), round(width/2))
  }

  # set xlab and ylab
  if (!is.element("xlab", names(args))) {
    args$xlab="Genomic position relative to BS"
  }
  if (!is.element("ylab", names(args))) {
    args$ylab="Median methylation"
  }

  # set col  
  if (!is.element("col", names(args))) {
    cols = rainbow(length(levels(groups)))
    args$col = cols[1]
  } else {
    cols = args$col
    if (length(cols) != length(levels(groups))) {
      cols = rep(cols, length(levels(groups)))
    }
    args$col = cols[1]
  }
  
  args$type = "l"
  
  do.call(plot, args)
  if (ncol(sData[[1]]) > 1) {
    for (i in 2:ncol(sData[[1]])) {
      lines(x=x, y=sData[[1]][, i], col=rgb(t(col2rgb(cols[1]))/255, alpha=0.6), lty=5)
    }
  }

  # plot other groups
  if (length(sData) > 1) {
    for (g in 2:length(sData)) {
      lines(x=x, y=sData[[g]][, 1], col=cols[g])
      if (ncol(sData[[g]]) > 1) {
        for (i in 2:ncol(sData[[g]])) {
          lines(x=x, y=sData[[g]][, i], col=rgb(t(col2rgb(cols[g]))/255, alpha=0.6), lty=5)
        }
      }
    }
  }
  
#  cat(paste("Generating plot with", nrow(matchMatrix), "CpGs from", length(unique(matchMatrix[, 1])),
#      "distinct regions.\n"))

}

setMethod("plotBindingSites",
    signature=c(object="BSrel", regions="GRanges"),
    function(object, regions, width, groups, quantiles, bandwidth, ...) {
      .plotBindingSites(object, regions, width=width, groups=groups, quantiles=quantiles, bandwidth=bandwidth, ...)
    })

setMethod("plotBindingSites",
    signature=c(object="BSraw", regions="GRanges"),
    function(object, regions, width, groups, quantiles, bandwidth, ...) {
      .plotBindingSites(rawToRel(object), regions, width=width, groups=groups, quantiles=quantiles, bandwidth=bandwidth, ...)
    })
