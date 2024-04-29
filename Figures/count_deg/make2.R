
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
  xlab("# of cells per cell type") +
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
  xlab("# of cells per cell type") +
  labs(color="Cell type") +
  ylim(0,450) +
  guides(color="none")


set.seed(42)
res.1 <- res[,,,1]
res.2 <- res[,,,2]
for(i in 1:20) {
  for(j in 1:13) {
    res.1[j,,i] <- sample(res.1[j,,i], size=10, replace=FALSE)
    res.2[j,,i] <- sample(res.2[j,,i], size=10, replace=FALSE)
  }
}

gt <- rep(0,13); gt[c(2,6)] <- 1
TPR1 <- apply(res.1, 2:3, function(x) {
  pred <- rep(0,13)
  pred[order(x,decreasing=TRUE)[1:3]] <- 1
  #sum(gt[c(2,6)] == pred[c(2,6)])/2
  ifelse(2 %in% order(x,decreasing = TRUE)[1:3], 1, 0)
})

gt <- rep(0,13); gt[1:2] <- 1
TPR2 <- apply(res.2, 2:3, function(x) {
  pred <- rep(0,13)
  pred[order(x,decreasing=TRUE)[1:2]] <- 1
  #sum(gt[c(1,2)] == pred[c(1,2)])/2
  ifelse(2 %in% order(x,decreasing = TRUE)[1:3], 1, 0)
})

df <- data.frame(scDist = as.vector(TPR1),
                 nDEG = as.vector(TPR2))
df <- reshape2::melt(df)

df$value <- as.character(df$value)
# Convert to a new vector with desired values
df$value <- ifelse(df$value == "0", "No matches",
                  ifelse(df$value == "0.5", "1 match", "2 matches"))
df$value <- factor(df$value, levels=c("2 matches", "1 match", "No matches"))

pb <- ggplot(data=df,aes(x=value,fill=variable)) +
  geom_bar(position="dodge") +
  scale_fill_manual(values=c("forestgreen",
                             "firebrick")) +
  #facet_wrap(~variable,nrow=2) +
  theme_bw() +
  #guides(fill="none") +
  #xlab("True positive rate") +
  labs(fill = "Method") +
  xlab("")




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

res <- res[,,,1:2]

#plot of cell types paths

df <- reshape2::melt(res)
df <- df |> group_by(Var1, Var2, Var4) |> summarize(mean = mean(value))

df$Var1 <- LETTERS[1:10][df$Var1]
df$Var2 <- samples[df$Var2]
df$Var4 <- c("scDist", "nDEG (lme4)")[df$Var4]

pc <- df |> ggplot(aes(x=Var2,y=mean,color=Var1)) +
  geom_point() +
  geom_line(linetype="dashed") +
  #scale_x_log10() +
  facet_wrap(~Var4, nrow=3, scales="free_y") +
  theme_bw() +
  labs(color="Cell type")+
  xlab("Number of patients") +
  ylab("Perturbation")


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

p_d1 <- ggdendrogram(hc, rotate=FALSE, size=2)+ theme(axis.text.y = element_blank())
p_d1 <- p_d1 + ggtitle("scDist") + theme(plot.title = element_text(hjust = 0.5))

#D <- readRDS("Dist_mat.RDS")
my.dist <- as.dist(D.3)
hc <- hclust(my.dist)

p_d3 <- ggdendrogram(hc, rotate=FALSE, size=2) + theme(axis.text.y = element_blank())
p_d3 <- p_d3 + ggtitle("nDEG") + theme(plot.title = element_text(hjust = 0.5))

library(ggpubr)

p <- ggarrange(pa,pb,p_d1, p_d3, nrow=2,ncol=2,
               labels=c("a","b","c","d"))




p1 <- ggarrange(pa, pb, pc, pd, nrow=2,ncol=2,labels=c("a","b","c","d"))
p2 <- ggarrange(pe, nrow=1,ncol=1,labels=c("e"))

p <- ggarrange(p1,p2,nrow=2, heights = c(2,1))


ggsave(p, filename="count_deg_updated.png",height=8,width=8)
