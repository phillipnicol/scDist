
FullSimulate <- function(nn,
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
  out <- pcDiffPop(y,meta.data,fixed.effects=c("response"),
                   random.effects=c("samples"),
                   clusters="clusters",
                   d=d)
  end <- Sys.time()

  results <- list()
  results$dist <- out$results$Dist.
  results$p <- out$results$p.F
  results$beta_true <- beta_true
  results$pcdp_obj <- out
  results$scDist_time <- end-start
  if(augur) {
    results$augur <- auc
    results$auc.augur <- auc.augur
    results$augur_time <- time.augur
  }
  results
}

J <- c(10,100,10^3,10^4)
reps <- 2
res.scDist <- matrix(0,nrow=length(J),ncol=reps)
res.augur <- res.scDist
for(j in 1:length(J)) {
  for(k in 1:reps) {
    res <- FullSimulate(nn=100,dist_true=10,tau=0.5,N1=5,N2=5,G=1000,J=J[j],augur=FALSE)
    res.scDist[j,k] <- res$scDist_time
    res.augur[j,k] <- 0
  }
}

saveRDS(res.scDist, "scDist_time.RDS")
saveRDS(res.augur,"augur_time.RDS")
