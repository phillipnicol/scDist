
library(tidyverse)

path <- "../../analysis/count_deg_blish/data"
res <- readRDS(file.path(path, "upsample_res.RDS"))
time <- readRDS(file.path(path, "upsample_time.RDS"))
cts <- readRDS(file.path(path,"cts.RDS"))
ct.size <- readRDS(file.path(path,"ct.size.RDS"))

#plot of cell types paths

df <- reshape2::melt(res)
df <- df |> group_by(Var1, Var2, Var4) |> summarize(mean = mean(value))

df$Var1 <- cts[df$Var1]
df$Var2 <- ct.size[df$Var2]

df <- df |> filter(Var1 == "CD14 Monocyte")

df$Var4 <- c("scDist", "nDEG (Muscat)")[df$Var4]

pa <- df |> filter(Var4 == "scDist") |>
  ggplot(aes(x=Var2,y=mean,color=Var1)) +
  geom_point() +
  geom_line(linetype="dashed") +
  scale_x_log10() +
  theme_bw() +
  ylab("CD14 Monocyte Distance (scDist)")+
  xlab("# of cells") +
  labs(color="Cell type") +
  ylim(0,25) +
  guides(color="none")

pb <-  df |> filter(Var4 == "nDEG (Muscat)") |>
  ggplot(aes(x=Var2,y=mean,color=Var1)) +
  geom_point() +
  geom_line(linetype="dashed") +
  scale_x_log10() +
  theme_bw() +
  ylab("CD14 Monocyte nDEG (Muscat)")+
  xlab("# of cells") +
  labs(color="Cell type") +
  ylim(0,450) +
  guides(color="none")



df <- reshape2::melt(time)
df$Var1 <- ct.size[df$Var1]
df <- df |> group_by(Var1, Var3) |>
  summarise(mean=mean(value)/60)
df$Var3 <- c("scDist", "nDEG (Muscat)")[df$Var3]

p <- ggplot(data=df,aes(x=Var1,
                        y=mean,
                        color=Var3)) +
  geom_point() + geom_line() +
  labs(color="Method", x="# of cells per cell type",
       y="Average runtime (minutes)") +
  theme_bw()
ggsave(p, filename="muscat_scDist_runtime.png")

path <- "../../analysis/sim_samples/data"
res <- readRDS(file.path(path,"sim_results.RDS"))
samples <- c(2,3,5,10,20,50,100)
set.seed(1)
nct <- 10; cts <- letters[1:10]
reps <- 10
D.true <- rexp(n=10,rate=0.05)

library(tidyverse)

#plot of cell types paths

df <- reshape2::melt(res)
df <- df |> group_by(Var1, Var2, Var4) |> summarize(mean = mean(value))

df$Var1 <- LETTERS[1:10][df$Var1]
df$Var2 <- samples[df$Var2]
df$Var4 <- c("scDist", "nDEG (lme4)", "Augur")[df$Var4]

pc <- df |> ggplot(aes(x=Var2,y=mean,color=Var1)) +
  geom_point() +
  geom_line(linetype="dashed") +
  #scale_x_log10() +
  facet_wrap(~Var4, nrow=3, scales="free_y") +
  theme_bw() +
  labs(color="Cell type")+
  xlab("Number of patients") +
  ylab("Perturbation")


sp.cor <- apply(res, 2:4, function(x) {
  cor(x, D.true, method="pearson")
})
sp.cor[is.na(sp.cor)] <- 0
df <- reshape2::melt(sp.cor)

df <- df |> group_by(Var1, Var3) |> summarize(mean=mean(value))

df$Var1 <- samples[df$Var1]
df$Var3 <- c("scDist", "nDEG (lme4)", "Augur")[df$Var3]

pd <- df |> ggplot(aes(x=Var1,y=mean,color=Var3)) +
  geom_point() +
  geom_line(linetype="dashed") +
  #scale_x_log10() +
  theme_bw() +
  labs(color="Method")+
  xlab("Number of patients") +
  ylab("Correlation")

library(ggpubr)
p <- ggarrange(pc, pd, nrow=2, labels=c("a","b"),
               heights=c(1.5,1))
ggsave(p, filename="vary_sample_all_method.png",
       width=9, height=9)

path <- "../../analysis/zheng_8eq/data"
cell.combos <- readRDS(file.path(path,"cell.combos.dist.RDS"))

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


#df <- data.frame(x=c(my.mds1[,1],my.mds2[,1],my.mds3[,1]),
#                 y=c(my.mds1[,2],my.mds2[,2],my.mds3[,2]),
#                 text=rep(rownames(my.mds1),3),
#                 method=c(rep("scDist",8),
#                          rep("Augur",8),
#                          rep("nDEG", 8)))



## Hclust
library("ggdendro")

#D <- readRDS("Dist_mat.RDS")
my.dist <- as.dist(D.1)
hc <- hclust(my.dist)

p_d1 <- ggdendrogram(hc, rotate=FALSE, size=2)
p_d1 <- p_d1 + ggtitle("scDist") + theme(plot.title = element_text(hjust = 0.5))

#D <- readRDS("Dist_mat.RDS")
my.dist <- as.dist(D.3)
hc <- hclust(my.dist)

p_d3 <- ggdendrogram(hc, rotate=FALSE, size=2)
p_d3 <- p_d3 + ggtitle("nDEG") + theme(plot.title = element_text(hjust = 0.5))



##Equal
cell.combos <- readRDS(file.path(path,"cell.combos.dist.equal.RDS"))

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


#df <- data.frame(x=c(my.mds1[,1],my.mds2[,1],my.mds3[,1]),
#                 y=c(my.mds1[,2],my.mds2[,2],my.mds3[,2]),
#                 text=rep(rownames(my.mds1),3),
#                 method=c(rep("scDist",8),
#                          rep("Augur",8),
#                          rep("nDEG", 8)))



## Hclust
library("ggdendro")

#D <- readRDS("Dist_mat.RDS")
my.dist <- as.dist(D.1)
hc <- hclust(my.dist)

p_d1.equal <- ggdendrogram(hc, rotate=FALSE, size=2)
p_d1.equal <- p_d1.equal + ggtitle("scDist") + theme(plot.title = element_text(hjust = 0.5))

#D <- readRDS("Dist_mat.RDS")
my.dist <- as.dist(D.3)
hc <- hclust(my.dist)

p_d3.equal <- ggdendrogram(hc, rotate=FALSE, size=2)
p_d3.equal <- p_d3.equal + ggtitle("nDEG") + theme(plot.title = element_text(hjust = 0.5))

p <- ggarrange(p_d1.equal, p_d3.equal)
ggsave(p, filename="zheng_full.pdf")


library(ggpubr)


p <- ggarrange(pb,pa,p_d3,p_d1,
               labels=c("a","b","c","d"))


ggsave(p, filename="count_deg_updated.pdf",height=8,width=8)


