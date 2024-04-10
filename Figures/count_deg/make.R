
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


library(ggpubr)
p <- ggarrange(ggarrange(pa,pb,pc,pd,
                         nrow=2,ncol=2,
                         widths=c(0.7,0.3,0.7,0.3)),
               p_deg_sim,
               nrow=2,ncol=1,
               widths=c(1,1), heights=c(1,1))
