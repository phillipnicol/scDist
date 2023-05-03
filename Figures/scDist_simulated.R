


set.seed(1)
library(scDist)

tau <- c(0.1,0.5,1)
dist_true <- seq(0,100,by=10)
reps <- 25

res <- array(0, c(length(dist_true),length(tau),reps))

for(i in 1:length(dist_true)) {
  for(j in 1:length(tau)) {
    cat(i, " ", j, "\n")
    for(k in 1:reps) {
      scd <- scDist::sim_scDist(nn=100,dist_true=dist_true[i],
                                N1=5,N2=5,d=20,G=1000,
                                J=50,tau=tau[j])
      res[i,j,k] <- scd$scDist_obj$results$Dist.
    }
  }
}


library(reshape2)
library(tidyverse)

df <- melt(res)

df <- df %>% group_by(Var1,Var2) %>% summarise(med=median(value),
                                               q3=quantile(value,0.75),
                                               q1=quantile(value,0.25))

p <- ggplot(data=df,aes(x=dist_true[Var1],
                        y=med,
                        group=c("Low", "Medium", "High")[Var2],
                        color=c("Low", "Medium", "High")[Var2]))
p <- p + geom_abline(slope=1,intercept=0,linetype="dashed")
p <- p + geom_point() + geom_line() + theme_bw() + labs(color="Sample-level variability")
p <- p + xlab("Ground truth distance") + ylab("Estimated distance (Posterior median)")
p
