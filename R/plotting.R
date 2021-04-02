
betaPlot <- function(pcdiff, cluster) {
  ix <- which(names(out$vals) == cluster)
  betas <- out$vals[[ix]]$beta
  plot(betas,type="h")
}
