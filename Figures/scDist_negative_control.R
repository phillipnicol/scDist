### This is figure 3 in the manuscript

#Load Seurat
library(Seurat)

## COVID_HC_TEST

library(Augur)
library(scDist)
library(ggplot2)

### Load wilk et al COVID-19 data
### Sco <- readRDS("wilk_covid.RDS")

Sco <- Sco[,Sco$Status == "Healthy"]

nc <- length(unique(Sco$cell.type.coarse))
set.seed(4589)
reps <- 20
res.augur <- matrix(0.5,nrow=nc,ncol=reps)
rownames(res.augur) <- unique(Sco$cell.type.coarse)
res.pc <- matrix(0,nrow=nc,ncol=reps)
rownames(res.pc) <- unique(Sco$cell.type.coarse)
expr <- Sco@assays$SCT@scale.data

for(i in 1:reps) {
  print(i)
  group1 <- sample(unique(Sco$Donor),size=3,replace=FALSE)


  meta.data <- data.frame(label=ifelse(Sco$Donor %in% group1, "A", "B"),
                          cell_type=Sco$cell.type.coarse,
                          sample=Sco$Donor)

  expr <- Sco@assays$SCT@scale.data

  out <- scDist(expr,meta.data,fixed.effects="label",
                random.effects="sample",
                clusters="cell_type",d=20)
  res.pc[,i] <- out$results$Dist.
}

covid.dist <- res.pc
covid.dist <- as.data.frame(covid.dist)
df.gg <- reshape2::melt(t(covid.dist))

p <- ggplot(data=df.gg,aes(x=Var2,y=value))
p <- p + geom_boxplot(fill="lightblue",
                      outlier.shape=NA)
p <- p+geom_hline(yintercept=0,color="red",
                  linetype="dashed")
p <- p + theme_bw()
p <- p + ylab("Distance (Posterior Median)")
p <- p + xlab("Cell type") + ylim(0,10)
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



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
    sim <- FullSimulate(nn=10,
                        dist_true=0,
                        tau=tau[i],
                        N1=5,
                        N2=5,
                        G=1000,
                        J=50,
                        augur=FALSE)
    res.augur[j,i] <- 0
    res.pc[j,i] <- sim$dist
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


### Hierarchical clustering



expr <- Sco@assays$SCT[Sco@assays$SCT@var.features,]

cell.types <- unique(Sco$cell.type.fine)

meta.data <- data.frame(cell.type=Sco$cell.type.fine,
                        samples=Sco$Donor,
                        cluster="A")

nc <- length(cell.types)
D <- matrix(0, nrow=nc, ncol=nc)
for(i in 1:nc) {
  for(j in 1:nc) {
    if(i < j) {
      next
    } else if(i == j) {
      next
    }
    ixs <- which(Sco$cell.type.fine %in% c(cell.types[i], cell.types[j]))
    expr.sub <- as.matrix(expr[,ixs])
    meta.sub <- meta.data[ixs,]
    cat(i, " ", j, "\n")
    try({
      out <- scDist(expr.sub,meta.sub,fixed.effects="cell.type",
                    random.effects="samples",
                    clusters="cluster")
      D[i,j] <- out$results$Dist.
    })
  }
}

colnames(D) <- cell.types
rownames(D) <- cell.types



## Hclust
library("ggdendro")

#D <- readRDS("Dist_mat.RDS")
my.dist <- as.dist(D)
hc <- hclust(my.dist)

p_hc <- ggdendrogram(hc, rotate=FALSE, size=2)


p_hc


