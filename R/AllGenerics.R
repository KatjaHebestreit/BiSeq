setGeneric("annotateGRanges", function(object, regions, name, regionInfo) 
  standardGeneric("annotateGRanges"))

setGeneric("betaRegression", function(formula, link, object, mc.cores, ...) 
    standardGeneric("betaRegression"))

setGeneric("binomLikelihoodSmooth", function(pred.pos, pos, m, n, h) 
  standardGeneric("binomLikelihoodSmooth"))

setGeneric("BSraw", function(methReads,
                             totalReads,
                             rowData,
                             colData = DataFrame(row.names=colnames(methReads)),
                             exptData = SimpleList(),
                             ...
                             )
           standardGeneric("BSraw"),
           signature=c("methReads", "totalReads", "rowData"))

setGeneric("BSrel", function(methLevel,
                             rowData,
                             colData = DataFrame(row.names=colnames(methLevel)),
                             exptData = SimpleList(),
                             ...
                             )
           standardGeneric("BSrel"),
           signature=c("methLevel", "rowData"))

setGeneric("clusterSites", function(object, groups, perc.samples, min.sites, max.dist, mc.cores, ...) 
  standardGeneric("clusterSites"))

setGeneric("clusterSitesToGR", function(object) 
  standardGeneric("clusterSitesToGR"))

setGeneric("compareTwoSamples", function(object, sample1, sample2, minDiff, max.dist) 
  standardGeneric("compareTwoSamples"))

setGeneric("covBoxplots", function(object, ...) 
    standardGeneric("covBoxplots"))

setGeneric("covStatistics", function(object) 
    standardGeneric("covStatistics"))

setGeneric("estLocCor", function(vario.sm) 
    standardGeneric("estLocCor"))

setGeneric("filterByCov", function(object, minCov, global) 
  standardGeneric("filterByCov"))

setGeneric("filterBySharedRegions", function(object, groups, perc.samples, no.samples, minCov) 
  standardGeneric("filterBySharedRegions"))

setGeneric("findDMRs", function(test.out, alpha, max.dist, diff.dir) 
  standardGeneric("findDMRs"))

setGeneric("limitCov", function(object, maxCov) 
    standardGeneric("limitCov"))

setGeneric("logisticRegression", function(formula, link, object, mc.cores) 
    standardGeneric("logisticRegression"))

setGeneric("makeVariogram", function(test.out, make.variogram) 
    standardGeneric("makeVariogram"))

setGeneric("methLevel", function(object) 
    standardGeneric("methLevel"))

setGeneric("methLevel<-", function(object, value) 
    standardGeneric("methLevel<-"))

setGeneric("methReads", function(object) 
    standardGeneric("methReads"))

setGeneric("methReads<-", function(object, value) 
    standardGeneric("methReads<-"))

setGeneric("plotBindingSites", function(object, regions, width, groups, quantiles, bandwidth, ...) 
  standardGeneric("plotBindingSites"))

setGeneric("plotMeth", function(object.raw, object.rel, region, col.lines, lwd.lines, col.points, ...) 
  standardGeneric("plotMeth"))

setGeneric("plotMethMap", function(object, region, groups, intervals, ...) 
  standardGeneric("plotMethMap"))

setGeneric("plotSmoothMeth", function(object.rel, region, groups, group.average, ...) 
  standardGeneric("plotSmoothMeth"))

setGeneric("predictMeth", function(object, h, grid.dist, mc.cores) 
  standardGeneric("predictMeth"))

setGeneric("rawToRel", function(object) 
    standardGeneric("rawToRel"))

setGeneric("readBismark", function(files, colData) 
    standardGeneric("readBismark"))

setGeneric("smoothVariogram", function(variogram, sill, bandwidth) 
  standardGeneric("smoothVariogram"))

setGeneric("summarizeRegions", function(object, regions, outputAll) 
    standardGeneric("summarizeRegions"))

setGeneric("testClusters", function(locCor, FDR.cluster) 
    standardGeneric("testClusters"))

setGeneric("totalReads", function(object) 
    standardGeneric("totalReads"))

setGeneric("totalReads<-", function(object, value) 
    standardGeneric("totalReads<-"))

setGeneric("trimClusters", function(clusters.rej, FDR.loc) 
    standardGeneric("trimClusters"))

setGeneric("writeBED", function(object, name, file) 
    standardGeneric("writeBED"))
