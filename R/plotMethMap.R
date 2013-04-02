.plotMethMap <- function(object, region, groups, intervals=FALSE, ...) {

  if (length(region) != 1) {
    stop("Argument \"region\" must contain exactly one genomic region.")
  }

  strand(object) <- "*"
  object <- sort(object)
  ind <- overlapsAny(rowData(object), region)
  if (sum(ind) == 0) {
    stop("No data to plot within the given region/samples.")
  }

  md <- methLevel(object)[ind, ]
  pos <- start(rowData(object))[ind]
  rownames(md) <- as.character(pos)
  md <- t(md)
  ind <- !apply(is.na(md), 2, all)
  md <- md[, ind]
  pos <- pos[ind]

  if (intervals) {
    fullMd <- matrix(NA, nrow=nrow(md), ncol=pos[length(pos)]-pos[1]+1)
    rownames(fullMd) <- rownames(md)
    colnames(fullMd) <- as.character(pos[1]:pos[length(pos)])
    fullMd[, colnames(md)] <- md
    md <- fullMd
    rm(fullMd)
  }
  
  ind.na <- apply(md, 1, function(x) all(is.na(x)))
  if (sum(ind.na) > 0) {
    names.na <- names(which(ind.na))
    md <- md[!ind.na, ]
    object <- object[, !ind.na]
    if (is.element("RowSideColors", names(args))) {
      args$RowSideColors <- args$RowSideColors[!ind.na]
    }
    warning("Sample(s) ", names.na, " omitted due to missing data within \"region\".")
  }

  if (intervals) {
    cn <- as.character(pretty(colnames(md), n=5))
    colnames(md)[!is.element(colnames(md), cn)] <- ""
  }
  
  args <- list(...)
  
  # set RowSideColors if groups is given
  if (!missing(groups) & !is.element("RowSideColors", names(args))) {
    pDat <- as.character(groups)
    uPDat <- unique(pDat)
    cols <- .categorialColors(length(uPDat))
    names(cols) <- uPDat
    args <- c(list(RowSideColors=cols[pDat]), args)
  }

  # default scale="none"
  if (!is.element("scale", names(args))) {
    args <- c(args, list(scale="none"))
  }

  # default zlim=c(0,1)
  if (!is.element("zlim", names(args))) {
    args <- c(args, list(zlim=c(0,1)))
  }
  
  # default Colv = NA
  if (!is.element("Colv", names(args))) {
    args <- c(args, list(Colv=NA))
  }

  # default cexCol=0.8 if intervals=TRUE and no labCol given
  if (intervals &
      !is.element("labCol", names(args)) &
      !is.element("cexCol", names(args))) {
    args <- c(args, list(cexCol=0.8))
  }

  # default colors green - black - red
  if (!is.element("col", names(args))) {
    colF <- colorRampPalette(colors=c("green", "black", "red"))
    args <- c(args, list(col=colF(64)))
  }

  args <- c(list(x=md), args)
  do.call(heatmap, args)
}


setMethod("plotMethMap",
    signature=c(object="BSrel", region="GRanges", groups="factor", intervals="logical"),
    .plotMethMap)

setMethod("plotMethMap",
    signature=c(object="BSrel", region="GRanges", groups="missing", intervals="logical"),
    .plotMethMap)

setMethod("plotMethMap",
    signature=c(object="BSrel", region="GRanges", groups="factor", intervals="missing"),
    .plotMethMap)

setMethod("plotMethMap",
    signature=c(object="BSrel", region="GRanges", groups="missing", intervals="missing"),
    .plotMethMap)

setMethod("plotMethMap",
    signature=c(object="BSraw", region="GRanges", groups="factor", intervals="logical"),
    function(object, region, groups, intervals, ...) {
      .plotMethMap(object=rawToRel(object), region=region, groups=groups, intervals=intervals, ...)
    })

setMethod("plotMethMap",
    signature=c(object="BSraw", region="GRanges", groups="missing", intervals="logical"),
    function(object, region, intervals, ...) {
      .plotMethMap(object=rawToRel(object), region=region, intervals=intervals, ...)
    })

setMethod("plotMethMap",
    signature=c(object="BSraw", region="GRanges", groups="factor", intervals="missing"),
    function(object, region, groups, ...) {
      .plotMethMap(object=rawToRel(object), region=region, groups=groups, ...)
    })

setMethod("plotMethMap",
    signature=c(object="BSraw", region="GRanges", groups="missing", intervals="missing"),
    function(object, region, ...) {
      .plotMethMap(object=rawToRel(object), region=region, ...)
    })
