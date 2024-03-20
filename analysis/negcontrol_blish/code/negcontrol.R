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