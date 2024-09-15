### This is Figure 1 in the manuscript


library(Augur)
library(ggplot2)

### Download the Wilk et al data
###Sco <- readRDS("wilk_covid.RDS")

Sco <- Sco[,Sco$Status == "Healthy"]

nc <- length(unique(Sco$cell.type.coarse))
set.seed(4589)
reps <- 20
res.augur <- matrix(0.5,nrow=nc,ncol=reps)
rownames(res.augur) <- unique(Sco$cell.type.coarse)
res.pc <- matrix(0,nrow=nc,ncol=reps)
rownames(res.pc) <- unique(Sco$cell.type.coarse)

for(i in 1:reps) {
  print(i)
  group1 <- sample(unique(Sco$Donor),size=3,replace=FALSE)


  meta.data <- data.frame(label=ifelse(Sco$Donor %in% group1, "A", "B"),
                          cell_type=Sco$cell.type.coarse,
                          sample=Sco$Donor)

  expr <- Sco@assays$SCT@scale.data

  augur <- calculate_auc(expr,meta=meta.data,n_threads = 4)
  augur$AUC <- as.data.frame(augur$AUC)
  for(j in augur$AUC$cell_type) {
    ix <- which(rownames(res.augur) == j)
    jx <- which(augur$AUC$cell_type == j)
    res.augur[ix,i] <- augur$AUC$auc[jx]
  }
}

res.augur <- as.data.frame(res.augur)
df.gg <- reshape2::melt(t(res.augur))

p <- ggplot(data=df.gg,aes(x=Var2,y=value))
p <- p + geom_boxplot(fill="lavender")
p <- p+geom_hline(yintercept=0.50,color="red",
                  linetype="dashed")
p <- p + theme_linedraw()
p <- p + ylab("Condition difference (Augur)")
p <- p + xlab("Cell type")
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

