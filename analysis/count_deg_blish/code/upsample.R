
library(Seurat)
library(scDist)
library(SingleCellExperiment)

Sco <- readRDS("../../data/blish_covid.seu.rds")

upsample <- function(Sco,ct.size) {
  ixs <- c()
  for(ct in unique(Sco$cell.type.coarse)) {
    ix <- sample(which(Sco$cell.type.coarse == ct), size=ct.size,
                 replace=TRUE)
    ixs <- c(ixs,ix)
  }
  return(ixs)
}


muscat_deg <- function(Sco.sub) {
  sce.sub <- SingleCellExperiment(counts=Sco@assays$RNA@counts[,ixs],
                                  colData=meta.data[ixs,])
  sce.sub <- prepSCE(sce.sub,
                     kid = "cell.type.coarse", # subpopulation assignments
                     gid = "Status",  # group IDs (ctrl/stim)
                     sid = "Donor",   # sample IDs (ctrl/stim.1234)
                     drop = TRUE)

  pb <- aggregateData(sce.sub,
                      assay = "counts", fun = "sum",
                      by = c("cluster_id", "sample_id"))

  res <- pbDS(pb, verbose = FALSE)

  tbl <- res$table[[1]]
  # filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
  tbl_fil <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
    dplyr::arrange(u, p_adj.loc)
  })

  # nb. of DS genes & % of total by cluster
  n_de <- vapply(tbl_fil, nrow, numeric(1))
  p_de <- format(n_de / nrow(sce.sub) * 100, digits = 3)
  df <- data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)
  df
}


meta.data <- data.frame(response=Sco$Status,
                        sample=Sco$Donor,
                        cluster=Sco$cell.type.coarse)
out.full <- scDist(Sco@assays$SCT@scale.data[Sco@assays$SCT@var.features,],
                   meta.data,
                   fixed.effects="response",
                   random.effects="sample",
                   clusters="cluster")



set.seed(1)
reps <- 20
nc <- length(table(Sco$cell.type.coarse))
ct.size <- 10^{2:5}
res <- array(dim=c(nc,length(ct.size),reps,2))

for(i in 1:length(ct.size)) {
  for(j in 1:reps) {
    ixs <- upsample(Sco,ct.size=ct.size[i])
    meta.sub <- meta.data[ixs,]
    Y.sub <- Sco@assays$SCT@scale.data[,ixs]
    Y.sub <- Y.sub[Sco@assays$SCT@var.features,]

    out.sub <- scDist(Y.sub,
                      meta.sub,
                      fixed.effects="response",
                      random.effects="sample",
                      clusters="cluster")

    res[,i,j,1] <- out.sub$results$Dist.

    sce.sub <- SingleCellExperiment(assays=list(counts=counts.new),
                                    colData=coldata)
    coldata <- DataFrame(meta.sub)

    df <- muscat_deg(Sco, ixs)
    res[,i,j,2] <- df$`#DS`
  }
}
