.triangularWeight <- function(pred.pos, pos, h) {
  dist = pos %*% t(rep(1, length(pred.pos)))
  dist = t(t(dist) - pred.pos)
  w = apply(dist, 2, function(u) {
    triW = (1 - abs(u/h))
    triW[triW < 0] = 0
    return(triW)
  })
  return(w)
}


.binomLikelihoodSmooth <- function(pred.pos, pos, m, n, h) {
  w = .triangularWeight(pred.pos, pos, h)
  wm = t(w) %*% m
  wn = t(w) %*% n
  p = wm / wn
  return(p)
}


setMethod("binomLikelihoodSmooth",
    signature=c(pred.pos="ANY", pos="ANY", m="ANY", n="ANY", h="ANY"),
    .binomLikelihoodSmooth)
