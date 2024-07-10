### Simulated negative control

library(Seurat)
library(scDist)
library(Augur)
library(tidyverse)

### Simulated null data

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

  if(augur) {
    meta.data <- data.frame(cell_type=rep("A",(N1+N2)*J),label=as.character(response))
    expr <- t(y)
    rownames(expr) <- 1:G; colnames(expr) <- 1:(J*(N1+N2))
    auc <- calculate_auc(input=expr,meta=meta.data)
    auc.augur <- auc$AUC$auc
  }
  meta.data <- data.frame(response=response,samples=as.factor(samples),clusters="A")
  out <- scDist(t(y),meta.data,fixed.effects=c("response"),
                random.effects=c("samples"),
                clusters="clusters",
                d=d)

  results <- list()
  results$dist <- out$results$Dist.
  results$p <- out$results$p.F
  results$beta_true <- beta_true
  results$pcdp_obj <- out
  if(augur) {
    results$augur <- auc
    results$auc.augur <- auc.augur
  }
  results
}


## Null example
set.seed(1)
library(ashr)
library(scDist)
library(Augur)
library(lmerTest)
library(uniformly)
reps <- 100
tau <- seq(0,1,by=0.1)
res.augur <- matrix(0,nrow=reps,ncol=length(tau))
res.pc <- matrix(0,nrow=reps,ncol=length(tau))
for(i in 1:length(tau)) {
  for(j in 1:reps) {
    cat(i, " ", j, "\n")
    sim <- FullSimulate(nn=10,
                        dist_true=0,
                        tau=tau[i],
                        N1=5,
                        N2=5,
                        G=1000,
                        J=50,
                        augur=TRUE)
    res.augur[j,i] <- sim$auc.augur
    res.pc[j,i] <- sim$dist
    cat(i, " ", j, " ", res.pc[j,i]," ", res.augur[j,i], "\n")
  }
}

colnames(res.pc) <- seq(0,1,by=0.1)
df.gg <- data.frame(mean=colMeans(res.pc),
                    tau=colnames(res.pc),
                    se=apply(res.pc,2,function(x) sd(x)))
p <- ggplot(data=df.gg,aes(x=tau,y=mean,
                           ymin=mean-se,ymax=mean+se))
p <- p + geom_errorbar()
p <- p + geom_point(color="blue")
p <- p + geom_hline(yintercept=0,colour="red",linetype="dashed")
#p <- p + geom_text(x=0.8, y=0.5,vjust=-1,label="Ground Truth",color="red")
p <- p + theme_bw()
p <- p + xlab("Sample-specific effect")
p <- p + ylab("Distance (scDist)")
p


colnames(res.augur) <- seq(0,1,by=0.1)
df.gg <- data.frame(mean=colMeans(res.augur),
                    tau=colnames(res.augur),
                    se=apply(res.augur,2,function(x) sd(x)))
p <- ggplot(data=df.gg,aes(x=tau,y=mean,
                           ymin=mean-se,ymax=mean+se))
p <- p + geom_errorbar()
p <- p + geom_point(color="blue")
p <- p + geom_hline(yintercept=0.5,colour="red",linetype="dashed")
#p <- p + geom_text(x=0.8, y=0.5,vjust=-1,label="Ground Truth",color="red")
p <- p + theme_bw()
p <- p + xlab("Sample-specific effect")
p <- p + ylab("Distance (scDist)")
p






