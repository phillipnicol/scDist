
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

df$Var4 <- c("scDist", "nDEG")[df$Var4]

pa <- df |> ggplot(aes(x=Var2,y=mean,color=Var1)) +
  geom_point() +
  geom_line(linetype="dashed") +
  scale_x_log10() +
  facet_wrap(~Var4, nrow=2, scales="free_y") +
  theme_bw() +
  xlab("# of cells per cell type")



df <- reshape2::melt(time)
df$Var1 <- ct.size[df$Var1]
df <- df |> group_by(Var1, Var3) |>
  summarise(mean=mean(value)/60)
df$Var3 <- c("scDist", "nDEG")[df$Var3]

pb <- ggplot(data=df,aes(x=Var1,
                        y=mean,
                        color=Var3)) +
  geom_point() + geom_line() +
  labs(color="Method", x="# of cells per cell type",
       y="Average runtime (minutes)") +
  theme_bw()

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
  scale_x_log10() +
  facet_wrap(~Var4, nrow=3, scales="free_y") +
  theme_bw() +
  labs(color="Cell type")+
  xlab("Number of patients") +
  ylab("Perturbation")

## Pearson correlation

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
  scale_x_log10() +
  theme_bw() +
  labs(color="Cell type")+
  xlab("Number of patients") +
  ylab("Pearson correlation")

p <- ggarrange(pa,pb,pc,pd,
               nrow=2, ncol=2)


### MDS


path <- "../../analysis/zheng_8eq/data"
cell.combos <- readRDS(file.path(path,"cell.combos.dist.RDS"))

## Multidimension scaling
cell.combos$auc[is.na(cell.combos$auc)] <- 0
D.1 <- matrix(cell.combos$dist, nrow=8,ncol=8)
D.2 <- matrix(cell.combos$auc,nrow=8,ncol=8)
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


df <- data.frame(x=c(my.mds1[,1],my.mds2[,1],my.mds3[,1]),
                 y=c(my.mds1[,2],my.mds2[,2],my.mds3[,2]),
                 text=rep(rownames(my.mds1),3),
                 method=c(rep("scDist",8),
                          rep("Augur",8),
                          rep("nDEG", 8)))

library(ggrepel)

pe <- ggplot(data=df,aes(x=x,y=y,label=text)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "gray50",
                  segment.alpha = 0.5,
                  color = "black",
                  size = 3,
                  max.overlaps = Inf) +
  facet_wrap(~method, nrow=1,scales="free")

library(ggpubr)
p <- ggarrange(ggarrange(pa,pb,pc,pd,
                         nrow=2,ncol=2,
                         heights=c(3,4),
                         labels=c("a","b","c","d"),
                         widths=c(0.6*16,0.4*16)),
               pe,
               nrow=3,ncol=1, heights=c(7,4),labels=c("","e"))
ggsave(p, filename="count_deg.png",
       units="in",width=16,height=22)

p1 <- plot_grid(pa,pb,pc,pd,align="h",axis="tblr",nrow=2)

p2 <- plot_grid(p1, pe, nrow=2, rel_heights = c(2,1), axis="tblr", align="v")

ggsave(p2, filename="count_deg.png",
       units="in",width=16,height=22)
