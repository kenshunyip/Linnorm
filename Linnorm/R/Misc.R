rowVars <- function(x) {
  return ((rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1) )  )
}
rowSDs <- function(x) {
  return ((rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1) ) ^ 0.5 )
}

gammaloglik <- function(w,x,y,z) {
    .Call('gammaloglik', PACKAGE = 'Linnorm',w,x,y,z)
}
gammaShape <- function(x) {
    .Call('gammaShape', PACKAGE = 'Linnorm', x)
}

LocateLambda <- function(x,y) {
    .Call('LocateLambda', PACKAGE = 'Linnorm', x,y)
}
