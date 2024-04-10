
library(scDist)
library(SingleCellExperiment)
library(Matrix)
library(muscat)
library(MAST)
library(tidyverse)

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

mast_deg <- function(Y.norm, meta.sub, cts) {
  res <- data.frame(cts=cts, ndeg=0)

  colnames(Y.norm) <- paste0("cell", 1:ncol(Y.norm))

  for(i in 1:length(cts)) {
    Y.norm.sub <- Y.norm[,meta.sub$cluster == cts[i]]
    meta.sub.sub <- meta.sub[meta.sub$cluster == cts[i],]

    df <- reshape2::melt(Y.norm.sub)
    df <- df |> mutate(response=meta.sub.sub$response[Var2],
                       sample=meta.sub.sub$sample[Var2])

    sca <- FromFlatDF(df, idvars="Var2",
                      primerid="Var1",
                      measurement = "value",
                      cellvars=c("response","sample"))
    try({
    fit <- zlm(~response+(1|sample),
               sca,
               method="glmer",
               ebayes=FALSE,
               strictConvergence=FALSE)

    summaryCond <- suppressMessages(MAST::summary(fit,
                                                  doLRT='responseHealthy'))
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(summaryDt[summaryDt$contrast=='responseHealthy'
                                & summaryDt$component=='logFC', c(1,7,5,6,8)],
                      summaryDt[summaryDt$contrast=='responseHealthy'
                                & summaryDt$component=='H', c(1,4)],
                      by = 'primerid')

    num.deg <- sum(fcHurdle[,6] < 0.05/nrow(Y.norm.sub))
    if(!is.na(num.deg)) {
      res[i,2] <- num.deg
    }   })
  }


  return(res)
}


meta.data <- meta.data[,c("Status", "Donor", "cell.type.coarse")]
colnames(meta.data) <- c("response", "sample", "cluster")

out.full <- scDist(Y.sct,
                   meta.data,
                   fixed.effects="response",
                   random.effects="sample",
                   clusters="cluster")


sce <- SingleCellExperiment(assay=list(counts=Y.counts), colData=meta.data)

deg.full <- muscat_deg(sce)


set.seed(1)
reps <- 20
nc <- length(table(meta.data$cluster))
cts <- rownames(out.full$results)
ct.size <- round(10^{seq(log10(200), 4, length.out=10)})
res <- array(dim=c(nc,length(ct.size),reps,3))
time <- array(dim=c(length(ct.size), reps,3))

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

    #run Muscat
    sce.sub <- sce[,ixs]
    start <- Sys.time()
    df <- muscat_deg(sce.sub)
    end <- Sys.time()
    time[i,j,2] <- difftime(end,start,units="secs")
    res[,i,j,2] <- df$`#DS`

    #Run MAST
    start <- Sys.time()
    mast.res <- mast_deg(Y.sub, meta.sub, cts)
    end <- Sys.time()
    res[,i,j,3] <- mast.res[,2]
    time[i,j,3] <-  difftime(end,start,units="secs")
  }
}

saveRDS(res, "../data/upsample_res.RDS")
saveRDS(time, "../data/upsample_time.RDS")
saveRDS(cts, "../data/cts.RDS")
saveRDS(ct.size, "../data/ct.size.RDS")


res <- readRDS("../data/upsample_res.RDS")
time <- readRDS("../data/upsample_time.RDS")
cts <- readRDS("../data/cts.RDS")
ct.size <- readRDS("../data/ct.size.RDS")

library(tidyverse)

#plot of cell types paths

df <- reshape2::melt(res)
df <- df |> group_by(Var1, Var2, Var4) |> summarize(mean = mean(value))

df$Var1 <- cts[df$Var1]
df$Var2 <- ct.size[df$Var2]

df$Var4 <- c("scDist", "nDEG")[df$Var4]

p <- df |> ggplot(aes(x=Var2,y=mean,color=Var1)) +
  geom_point() +
  geom_line(linetype="dashed") +
  scale_x_log10() +
  facet_wrap(~Var4, nrow=2, scales="free_y") +
  theme_bw() +
  xlab("# of cells per cell type")

ggsave(plot=p, filename="../plots/upsample_ct_path.png")

## Timing

df <- reshape2::melt(time)
df$Var1 <- ct.size[df$Var1]
df <- df |> group_by(Var1, Var3) |>
  summarise(mean=mean(value)/60)
df$Var3 <- c("scDist", "nDEG")[df$Var3]

p <- ggplot(data=df,aes(x=Var1,
                        y=mean,
                        color=Var3)) +
  geom_point() + geom_line() +
  labs(color="Method", x="# of cells per cell type",
       y="Average runtime (minutes)") +
  theme_bw()

ggsave(p, filename="../plots/upsample_runtime.png")
