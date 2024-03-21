#Load Seurat
library(Seurat)

## COVID_HC_TEST

library(Augur)
library(ggplot2)
library(scDist)

### Load wilk et al COVID-19 data
#Sco <- readRDS("../../data/blish_covid.seu.rds")
#Sco <- UpdateSeuratObject(Sco)
load("../../data/blish_covid_seu.rda")

Sco <- Sco[,Sco$Status == "Healthy"]

nc <- length(unique(Sco$cell.type.coarse))
set.seed(4589)
reps <- 20
res.augur <- matrix(0.5,nrow=nc,ncol=reps)
rownames(res.augur) <- unique(Sco$cell.type.coarse)
res.pc <- matrix(0,nrow=nc,ncol=reps)
res.pv <- matrix(0,nrow=nc,ncol=reps)
rownames(res.pc) <- unique(Sco$cell.type.coarse)
rownames(res.pv) <- unique(Sco$cell.type.coarse)
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
  res.pv[,i] <- out$results$p.val
}

saveRDS(res.pc, "../data/scDist_results.RDS")
saveRDS(res.augur, "../data/augur_results.RDS")
saveRDS(res.pv, "../data/scDist_pval_results.RDS")


res.pc  <- readRDS("../data/scDist_results.RDS")
res.augur <- readRDS("../data/augur_results.RDS")
res.pv <- readRDS("../data/scDist_pval_results.RDS")

library(tidyverse)
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
ggsave(p, filename="../plots/scDist_dist.png")

covid.pv <- res.pv |> as.data.frame()
df.gg <- reshape2::melt(t(covid.pv))
df.gg$value <- -log10(df.gg$value)

p <- ggplot(data=df.gg,aes(x=Var2,y=value))
p <- p + geom_boxplot(fill="lightblue",
                      outlier.shape=NA)
p <- p+geom_hline(yintercept=-log10(0.05),color="red",
                  linetype="dashed")
p <- p + theme_bw()
p <- p + ylab("-log10 p-value (scDist)")
p <- p + xlab("Cell type")
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(p, filename="../plots/scDist_pval.png")



### FMT comparison
res.pv <- readRDS("../data/scDist_pval_results_FMT.RDS")


covid.pv <- res.pv |> as.data.frame()
df.gg <- reshape2::melt(t(covid.pv))
df.gg$value <- -log10(df.gg$value)

p <- ggplot(data=df.gg,aes(x=Var2,y=value))
p <- p + geom_boxplot(fill="lightblue",
                      outlier.shape=NA)
p <- p+geom_hline(yintercept=-log10(0.05),color="red",
                  linetype="dashed")
p <- p + theme_bw()
p <- p + ylab("-log10 p-value (scDist)")
p <- p + xlab("Cell type")
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(p, filename="../plots/scDist_pval_FMT.png")
