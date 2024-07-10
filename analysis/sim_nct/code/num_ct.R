library(lme4)
library(lmerTest)
library(scDist)
library(Augur)
library(muscat)
library(SingleCellExperiment)

simCellType <- function(D,tau,G=1000,N1=5,N2=5,J=50,label="A",my.pi=0.9) {
  beta_true <- rep(0,G)

  print(my.pi)
  z <- sample(1:2,size=G,replace=TRUE, prob=c(my.pi,1-my.pi))
  beta_true[z == 1] <- rnorm(n=sum(z==1), mean=0, sd=0.1)
  beta_true[z == 2] <- rnorm(n=sum(z==2), mean=0, sd=1)

  beta_true <- D/sqrt(sum(beta_true^2))*beta_true

  print(J)
  y <- matrix(0, nrow=(N1+N2)*J,ncol=G)
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

  meta.data <- data.frame(response=response,patient=as.factor(samples),clusters=label)
  out <- list()
  out$Y <- t(y); out$meta.data <- meta.data
  return(out)
}

count_deg <- function(Y, meta.data, cts) {
  perturb <- rep(0, length(cts))
  for(k in 1:length(cts)) {
    ixs <- which(meta.data$clusters == cts[k])
    Y.sub <- Y[,ixs]
    meta.sub <- meta.data[ixs,]
    response <- meta.sub$response
    patient <- meta.sub$patient
    p.vals <- apply(Y.sub, 1, function(y) {
      fit <- lm(y~response)
      summary(fit)$coefficients[2,4]
    })
    perturb[k] <- sum(p.vals < 0.05)
    print(perturb[k])
  }

  return(perturb)
}

set.seed(1)
nct <- 10; cts <- letters[1:10]
reps <- 10
res <- array(0, dim=c(nct, reps, 3))
D.true <- rexp(n=10,rate=0.05)
pi.true <- rbeta(n=10,shape1=1, shape2=1)
cell.type.size <- round(seq(100, 10^4, length.out=10))
J <- 50
tau <- 1
G <- 1000

for(i in 1:reps) {
  Y <- matrix(0,nrow=G,ncol=0)
  meta.data <- data.frame(response=NULL,
                          patient=NULL,
                          clusters=NULL)
  cell.type.size <- sample(cell.type.size,
                           size=10,
                           replace=FALSE)
  for(k in 1:nct) {
    out <- simCellType(D=D.true[k],J=cell.type.size[k],
                       N1=5,N2=5,
                       label=letters[k],tau=tau,
                       my.pi=pi.true[k])
    Y <- cbind(Y,out$Y)
    meta.data <- rbind(meta.data,out$meta.data)
  }

  out <- scDist(Y,meta.data,fixed.effects = "response",
                random.effects="patient",
                clusters="clusters",d=15)

  res[,i,1] <- out$results$Dist.

  res[,i,2] <- count_deg(Y, meta.data, cts)

  N1 <- 5; N2 <- 5
  meta.data.a <- data.frame(cell_type=meta.data$clusters,label=as.character(meta.data$response))
  rownames(Y) <- 1:G; colnames(Y) <- 1:ncol(Y)
  auc <- calculate_auc(input=Y,meta=meta.data.a)
  res[,i,3] <- as.vector(auc$AUC[order(auc$AUC$cell_type),]$auc)
}

saveRDS(res, "../data/sim_results.RDS")

res <- readRDS("../data/sim_results.RDS")

library(tidyverse)

my.cor <- apply(res, c(2,3), function(x) {
  cor(x, D.true, method="pearson")
})
my.cor[is.na(my.cor)] <- 0
df <- reshape2::melt(my.cor)


