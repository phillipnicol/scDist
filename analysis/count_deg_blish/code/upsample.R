
library(scDist)
library(SingleCellExperiment)
library(Matrix)
library(muscat)

Y.sct <- readRDS("../../data/blish_SCT_counts.RDS")
Y.counts <- readMM("../../data/blish_raw_counts.txt")
meta.data <- readRDS("../../data/blish_meta.RDS")

upsample <- function(meta.data,ct.size) {
  ixs <- c()
  for(ct in unique(meta.data$cluster)) {
    ix <- sample(which(meta.data$cluster == ct), size=ct.size,
                 replace=TRUE)
    ixs <- c(ixs,ix)
  }
  return(ixs)
}


muscat_deg <- function(sce.sub) {
  sce.sub <- prepSCE(sce.sub,
                     kid = "cluster", # subpopulation assignments
                     gid = "response",  # group IDs (ctrl/stim)
                     sid = "sample",   # sample IDs (ctrl/stim.1234)
                     drop = TRUE)

  pb <- aggregateData(sce.sub,
                      assay = "counts", fun = "sum",
                      by = c("cluster_id", "sample_id"))

  res <- pbDS(pb, verbose = FALSE, min_cells=0)

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


meta.data <- meta.data[,c("Status", "Donor", "cell.type.coarse")]
colnames(meta.data) <- c("response", "sample", "cluster")

out.full <- scDist(Y.sct,
                   meta.data,
                   fixed.effects="response",
                   random.effects="sample",
                   clusters="cluster")

sce <- SingleCellExperiment(assay=list(counts=Y.counts), colData=meta.data)



set.seed(1)
reps <- 20
nc <- length(table(meta.data$cluster))
cts <- unique(meta.data$cluster)
ct.size <- round(10^{seq(log10(200), 4, length.out=10)})
res <- array(dim=c(nc,length(ct.size),reps,2))
time <- array(dim=c(length(ct.size), reps,2))

for(i in 1:length(ct.size)) {
  for(j in 1:reps) {
    cat(i, " ", j, "\n")
    ixs <- upsample(meta.data,ct.size=ct.size[i])
    meta.sub <- meta.data[ixs,]
    Y.sub <- Y.sct[,ixs]

    start <- Sys.time()
    out.sub <- scDist(Y.sub,
                      meta.sub,
                      fixed.effects="response",
                      random.effects="sample",
                      clusters="cluster")
    end <- Sys.time()
    time[i,j,1] <- difftime(end,start,units="secs")

    res[,i,j,1] <- out.sub$results$Dist.

    sce.sub <- sce[,ixs]
    start <- Sys.time()
    df <- muscat_deg(sce.sub)
    end <- Sys.time()
    time[i,j,2] <- difftime(end,start,units="secs")
    res[,i,j,2] <- df$`#DS`
  }
}

saveRDS(res, "../data/upsample_res.RDS")
saveRDS(time, "../data/upsample_time.RDS")

res <- readRDS("../data/upsample_res.RDS")
time <- readRDS("../data/upsample_time.RDS")

library(tidyverse)

#plot of cell types paths

df <- reshape2::melt(res)
df <- df |> group_by(Var1, Var2, Var4) |> summarize(mean = mean(value))

df$Var1 <- cts[df$Var1]
df$Var2 <- ct.size[df$Var2]

p <- df |> ggplot(aes(x=Var2,y=mean,color=Var1)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Var4) +
  theme_bw()

ggsave(plot=p, filename="../plots/upsample_ct_path.png")

## Plot of correlation with ranking on full data
df <- reshape2::melt(res)

scDist.cor <- apply(res[,,,1], c(2,3), function(x) {
  cor(out.full$results$Dist., x)
})

deg.cor <- apply(res[,,,2], c(2,3), function(x) {
  cor(out.full$results$Dist., x)
})

