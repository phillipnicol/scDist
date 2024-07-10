### This one does not subsample CD14 Monocytes





library(DuoClustering2018)
library(Seurat)
library(SingleCellExperiment)

sce <- DuoClustering2018::sce_filteredExpr10_Zhengmix8eq()
Y <- sce@assays@data$counts

meta <- sce@colData |> as.data.frame()

Sco <- CreateSeuratObject(counts=Y,meta.data=meta)

Sco <- SCTransform(Sco)

library(scDist)

cell.combos <- expand.grid(unique(meta$phenoid),
                           unique(meta$phenoid))
cell.combos$dist <- 0
cell.combos$auc <- 0
cell.combos$degBF <- 0
cell.combos$degFDR <- 0
cell.combos$degRAW <- 0

Y <- Sco@assays$SCT@scale.data


for(i in 1:nrow(cell.combos)) {
  if(cell.combos[i,1] == cell.combos[i,2]) {
    next
  }

  ixs <- which(Sco@meta.data$phenoid %in% cell.combos[i,1] | Sco@meta.data$phenoid %in% cell.combos[i,2])

  Y.sub <- Y[,ixs]
  meta.sub <- Sco@meta.data[ixs,]
  meta.sub$hold <- rep("A", nrow(meta.sub))

  out <- scDist(normalized_counts = Y.sub,
                meta.data = meta.sub,
                fixed.effects = "phenoid",
                clusters="hold",
                d=10)

  cell.combos$dist[i] <- out$results$Dist.

  #augur <- calculate_auc(input=Y.sub,meta=meta.sub,label_col="phenoid",cell_type_col="hold")
  #cell.combos$augur[i] <- max(augur$AUC$auc - 0.5, 0)

  pvals <- apply(Y.sub, 1, function(x) {
    t.test(x[meta.sub$phenoid == cell.combos[i,1]],
           x[meta.sub$phenoid == cell.combos[i,2]])$p.val
  })
  cell.combos$degBF[i] <- sum(pvals < 0.05/nrow(Y.sub))
  cell.combos$degRAW[i] <- sum(pvals < 0.05)
  cell.combos$degFDR[i] <- sum(p.adjust(pvals,method="fdr") < 0.05)
}


saveRDS(cell.combos, "../data/cell.combos.dist.equal.RDS")
