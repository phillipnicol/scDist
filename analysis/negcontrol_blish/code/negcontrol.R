#Load Seurat
library(Seurat)
library(lme4)

# Trabunzi code
runTrabunzi <- function(expr, meta.data) {
  lmm <- rep(0, length(unique(meta.data$cell_type)))

  for(j in 1:length(unique(meta.data$cell_type))) {
    ct <- unique(meta.data$cell_type)[j]
    print(ct)
    Y.sub <- expr[,meta.data$cell_type==ct]
    G <- nrow(Y.sub)
    meta.sub <- meta.data[meta.data$cell_type == ct,]

    Y.sub <- rbind(Y.sub, meta.sub$label, meta.sub$sample)
    rownames(Y.sub) <- c(paste0("Genes", 1:G), "Response", "Patient")
    Y.sub <- t(Y.sub); Y.sub <- as.data.frame(Y.sub)
    df <- reshape2::melt(Y.sub,id.vars=c("Response", "Patient"))
    df$Patient <- as.factor(df$Patient); df$Response <- as.factor(df$Response)
    df$value <- as.numeric(df$value)
    df$VxR <- factor(df$variable:df$Response)
    print(dim(df))
    fit <- lmer(value ~ (1 | Patient) + (1 | variable) + (1 | VxR),
                data=df)
    random_effects <- VarCorr(fit)

    lmm[j] <- sqrt(random_effects$`VxR`[1])
    print(lmm[j])
  }
  return(lmm)
}


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
res.tb <- res.pc
rownames(res.pc) <- unique(Sco$cell.type.coarse)
rownames(res.pv) <- unique(Sco$cell.type.coarse)
expr <- Sco@assays$SCT@scale.data

runtime <- matrix(0, nrow=reps, ncol=3)

for(i in 1:reps) {
  print(i)
  group1 <- sample(unique(Sco$Donor),size=3,replace=FALSE)


  meta.data <- data.frame(label=ifelse(Sco$Donor %in% group1, "A", "B"),
                          cell_type=Sco$cell.type.coarse,
                          sample=Sco$Donor)

  expr <- Sco@assays$SCT@scale.data

  start <- Sys.time()
  out <- scDist(expr,meta.data,fixed.effects="label",
                random.effects="sample",
                clusters="cell_type",d=20)
  end <- Sys.time()
  runtime[i,1] <- difftime(end,start,units="secs")


  res.pc[,i] <- out$results$Dist.
  res.pv[,i] <- out$results$p.val

  start <- Sys.time()
  res.tb[,i] <- runTrabunzi(expr,meta.data)
  end <- Sys.time()
  runtime[i,3] <- difftime(end,start,units="secs")
}


saveRDS(res.pc, "../data/scDist_results.RDS")
saveRDS(res.augur, "../data/augur_results.RDS")
saveRDS(res.pv, "../data/scDist_pval_results.RDS")
saveRDS(res.tb, "../data/trabunzi_results.RDS")
saveRDS(runtime, "../data/runtime.RDS")

stop("Stopping")


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



## Trabunzi
res.tb <- readRDS("../data/trabunzi_results.RDS")
res.pc <- readRDS("../data/scDist_pval_results.RDS")
runtime <- readRDS("../data/runtime.RDS")
rownames(res.tb) <- rownames(res.pc)

library(tidyverse)

covid.tb <- res.tb |> as.data.frame()
df.gg <- reshape2::melt(t(covid.tb))

p <- ggplot(data=df.gg,aes(x=Var2,y=value,group=Var2))
p <- p + geom_boxplot(fill="lightblue")
p <- p+geom_hline(yintercept=0,color="red",
                  linetype="dashed")
p <- p + theme_bw()
p <- p + ylab("Variance")
p <- p + xlab("Cell type")
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1 <- p


runtime.comp <- runtime[,c(1,3)]
df <- reshape2::melt(runtime.comp)
df$Var2 <- c("scDist", "Trabunzi")[df$Var2]

p <- ggplot(data=df,aes(x=Var2,group=Var2, y=value/60)) +
  geom_boxplot(fill="lightgreen") +
  theme_bw() +
  scale_y_log10() +
  ylab("Runtime in minutes") +
  xlab("")
p2 <- p

library(ggpubr)

p <- ggarrange(p1, p2, nrow=2,labels=c("a","b"),
               heights=c(1.5,1))

ggsave(p, filename="../plots/trabunzi.png",
       width=5, height=4,units="in")
