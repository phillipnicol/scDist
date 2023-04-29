#' @export
sim_scDist <- function(nn,
                      dist_true,
                      tau,
                      N1,
                      N2,
                      d=20,
                      G,
                      J,
                      augur=FALSE) {


  beta_true <- rep(0,G)
  beta_true[1:nn] <- runif_on_sphere(n=1,d=nn,r=dist_true)[1,]

  y <- matrix(0, nrow=(N1+N2)*J,ncol=G)
  cntr <- 1
  for(i in 1:G) {
    cntr <- 1
    for(j in 1:N1) {
      omega <- rnorm(1,mean=0,sd=tau)
      y[cntr:(cntr+J-1),i] <- omega+rnorm(J)
      cntr <- cntr+J
    }
    for(j in 1:N2) {
      omega <- rnorm(1,mean=0, sd=tau)
      y[cntr:(cntr+J-1),i] <- beta_true[i]+omega+rnorm(J)
      cntr <- cntr+J
    }
  }

  response <- rep(0,(N1+N2)*J)
  response[1:(N1*J)] <- 1

  samples <- c()
  for(i in 1:(N1+N2)) {
    samples <- c(samples, rep(i,J))
  }

  time.augur <- 0
  if(augur) {
    meta.data <- data.frame(cell_type=rep("A",(N1+N2)*J),label=as.character(response))
    expr <- t(y)
    rownames(expr) <- 1:G; colnames(expr) <- 1:(J*(N1+N2))
    start <- Sys.time()
    auc <- calculate_auc(input=expr,meta=meta.data)
    auc.augur <- auc$AUC$auc
    end <- Sys.time()
    time.augur <- end-start
  }
  meta.data <- data.frame(response=response,samples=as.factor(samples),clusters="A")
  start <- Sys.time()
  out <- scDist(t(y),meta.data,fixed.effects=c("response"),
                   random.effects=c("samples"),
                   clusters="clusters",
                   d=d)
  end <- Sys.time()

  results <- list()
  results$dist <- out$results$Dist.
  results$p <- out$results$p.F
  results$beta_true <- beta_true
  results$scDist_obj <- out
  results$scDist_time <- end-start
  if(augur) {
    results$augur <- auc
    results$auc.augur <- auc.augur
    results$augur_time <- time.augur
  }
  results
}
