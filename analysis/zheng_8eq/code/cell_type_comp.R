

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

cell.type.sizes <- c(395,395, 50, 395,395,395,395,395)
ixs <- c()
set.seed(1)
for(i in 1:length(unique(Sco@meta.data$phenoid))) {
  jxs <- which(Sco@meta.data$phenoid == unique(Sco@meta.data$phenoid)[i])
  jxs <- sample(jxs, size=cell.type.sizes[i], replace=FALSE)
  ixs <- c(ixs, jxs)
}

Sco <- Sco[,ixs]

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
                clusters="hold")

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

saveRDS(cell.combos, "../data/cell.combos.dist.RDS")

cell.combos <- readRDS("../data/cell.combos.dist.RDS")
blish_res <- readRDS("../../blish/data/scDist_results.RDS")

blish_dists <- blish_res$Dist.

D <- matrix(cell.combos$dist, nrow=8,ncol=8)
rownames(D) <- unique(cell.combos$Var1)
colnames(D) <- rownames(D)
zheng.dists <- as.vector(D[D>0])

library(ggplot2)

# Create a data frame with the two vectors
data <- data.frame(
  variable = rep(c("zheng.dists", "blish_dists"), times = c(length(zheng.dists), length(blish_dists))),
  value = c(zheng.dists, blish_dists)
)

custom_colors <- c("forestgreen", "skyblue")
# Plot the density
p <- ggplot(data, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  labs(x = "Distance", y = "Density", fill = "Dataset") +
  scale_fill_manual(values = custom_colors,
                    labels = c("Blish Distances", "Zheng Distances"),
                    name = "Dataset") +
  theme_bw()

ggsave(p,
       filename="../plots/blish_v_zheng_dist.png")


## Multidimension scaling
cell.combos$auc[is.na(cell.combos$auc)] <- 0
D.1 <- matrix(cell.combos$dist, nrow=8,ncol=8)
#D.2 <- matrix(cell.combos$auc,nrow=8,ncol=8)
D.2 <- D.1
D.3 <- matrix(cell.combos$degFDR, nrow=8, ncol=8)
rownames(D.1) <- unique(cell.combos$Var1)
colnames(D.1) <- rownames(D.1)
rownames(D.2) <- unique(cell.combos$Var1)
colnames(D.2) <- rownames(D.1)
rownames(D.3) <- unique(cell.combos$Var1)
colnames(D.3) <- rownames(D.3)
my.mds1 <- cmdscale(D.1)
my.mds2 <- cmdscale(D.2)
my.mds3 <- cmdscale(D.3)


#df <- data.frame(x=c(my.mds1[,1],my.mds2[,1],my.mds3[,1]),
#                 y=c(my.mds1[,2],my.mds2[,2],my.mds3[,2]),
#                 text=rep(rownames(my.mds1),3),
#                 method=c(rep("scDist",8),
#                          rep("Augur",8),
#                          rep("nDEG", 8)))

df <- data.frame(x=c(my.mds1[,1],my.mds3[,1]),
                 y=c(my.mds1[,2],my.mds3[,2]),
                 text=rep(rownames(my.mds1),2),
                 method=c(rep("scDist",8),
                          rep("nDEG", 8)))

library(ggrepel)

p <- ggplot(data=df,aes(x=x,y=y,label=text)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "gray50",
                  segment.alpha = 0.5,
                  color = "black",
                  size = 3,
                  max.overlaps = Inf) +
  facet_wrap(~method, nrow=1,scales="free")

ggsave(p, filename="../plots/mds.png")





## Hclust
library("ggdendro")

#D <- readRDS("Dist_mat.RDS")
my.dist <- as.dist(D.1)
hc <- hclust(my.dist)

p_d1 <- ggdendrogram(hc, rotate=FALSE, size=2)


## Hclust
library("ggdendro")

#D <- readRDS("Dist_mat.RDS")
my.dist <- as.dist(D.3)
hc <- hclust(my.dist)

p_d3 <- ggdendrogram(hc, rotate=FALSE, size=2)

library(ggpubr)

p <- ggarrange(p_d3, p_d1, nrow=1, labels=c("nDEG", "scDist"))


