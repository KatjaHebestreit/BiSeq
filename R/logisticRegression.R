.logisticRegression <- function(formula, link, object, mc.cores = 1){
  
  strand(object) <- "*"
  object <- sort(object)

  if(link == "loglog"){
    inv.link <- function(x){
      exp(-exp(-x))
    }
  }
  if(link == "logit"){c
    inv.link <- function(x){
      1 / ( 1 + exp(-x))
    }
  }
  if(link == "probit"){
    inv.link <- function(x){
      pnorm(x)
    }
  }
  if(link == "cloglog"){
    inv.link <- function(x){
      1 - exp(-exp(x))
    }
  }
  if(link == "log"){
    inv.link <- function(x){
      exp(x)
    }
  }

  object.split <- split(object, f = as.factor(seqnames(object)), drop = TRUE)

  logistic.regr <- function(object.chr, formula){
    chr <- as.character(seqnames(rowData(object.chr))[1])
    pred.meth.chr <- methLevel(object.chr)

    p.val <- rep(NA,nrow(object.chr))
    meth.group1 <- rep(NA,nrow(object.chr))
    meth.group2 <- rep(NA,nrow(object.chr))
    meth.diff <- rep(NA,nrow(object.chr))
    direction <- rep(NA,nrow(object.chr))

    for(j in 1:nrow(pred.meth.chr)){
      pred.meth <- pred.meth.chr[j,]
      pred.meth[pred.meth == 0] <- 1*10^-10
      pred.meth[pred.meth == 1] <- 1-1*10^-10
      data <- cbind(pred.meth = pred.meth,
                    meth.odds = log(pred.meth/ (1-pred.meth)),
                    as.data.frame(colData(object)))
      if(length(unique(data$pred.meth))==1){ # all samples have the same methylation level
        p.val[j] <- 1
        meth.group1[j] <- unique(data$pred.meth)
        meth.group2[j] <- unique(data$pred.meth)
        meth.diff[j] <- 0
      } else {
        lmodel <- try(summary(lm(formula = as.formula(paste("meth.odds ~ ", formula[2])), data=data)))
        if(class(lmodel) == "try-error"){
          if( grepl("contrasts can be applied only to factors with 2 or more levels", lmodel[1]) ){  ## Too many samples with no coverage
            p.val[j] <- NA
            meth.diff[j] <- NA              }
        } else{
          p.val[j] <- lmodel$coefficients[2, 4]
          baseline <- lmodel$coefficients[1, 1]
          coef <- lmodel$coefficients[2, 1]
          meth.group1[j] <- inv.link(baseline)
          meth.group2[j] <- inv.link(baseline + coef)
          meth.diff[j] <-  meth.group1[j] - meth.group2[j]
        }
      }
      direction[j] <- sign(meth.diff[j])
    }
    
    out <- data.frame(chr = chr,
                      pos = start(ranges(object.chr)),
                      p.val,
                      meth.group1,
                      meth.group2,
                      meth.diff)
    return(out)
  }
  summary.df.l <- mclapply(object.split, logistic.regr, formula=formula,
                           mc.cores = mc.cores)
  summary.df <- do.call(rbind, summary.df.l)
  summary.df$cluster.id <- elementMetadata(rowData(object))$cluster.id
  return(summary.df)
}



setMethod("logisticRegression",
          signature=c(formula = "formula", link = "character", object="BSrel", mc.cores = "numeric"),
          .logisticRegression)

setMethod("logisticRegression",
          signature=c(formula = "formula", link = "character", object="BSrel", mc.cores = "missing"),
          function(formula, link, object) {
            .logisticRegression(formula, link, object, mc.cores = 1)
          })
