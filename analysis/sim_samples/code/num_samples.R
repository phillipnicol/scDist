library(lme4)
library(lmerTest)
library(scDist)
library(Augur)

simCellType <- function(D,tau,G=1000,N1=5,N2=5,J=50,label="A",my.pi=0.9) {
  beta_true <- rep(0,G)

  print(my.pi)
  z <- sample(1:2,size=G,replace=TRUE, prob=c(my.pi,1-my.pi))
  beta_true[z == 1] <- rnorm(n=sum(z==1), mean=0, sd=0.1)
  beta_true[z == 2] <- rnorm(n=sum(z==2), mean=0, sd=1)

  beta_true <- D/sqrt(sum(beta_true^2))*beta_true

  print(J)
  y <- matrix(0, nrow=(N1+N2)*J,ncol=G)
  for(i in 1:G) {
    cntr <- 1
    for(j in 1:N1) {
      omega <- rnorm(1,mean=0,sd=tau)
      y[cntr:(cntr+J-1),i] <- omega+rnorm(J)
      cntr <- cntr+J
    }
    for(j in 1:N2) {
      omega <- rnorm(1,mean=0, sd=tau)
      y[cntr:(cntr+J-1),i] <- beta_true[i]+omega+rnorm(J)
      cntr <- cntr+J
    }
  }

  response <- rep(0,(N1+N2)*J)
  response[1:(N1*J)] <- 1

  samples <- c()
  for(i in 1:(N1+N2)) {
    samples <- c(samples, rep(i,J))
  }

  meta.data <- data.frame(response=response,patient=as.factor(samples),clusters=label)
  out <- list()
  out$Y <- t(y); out$meta.data <- meta.data
  return(out)
}

count_deg <- function(Y, meta.data, cts) {
  perturb <- rep(0, length(cts))
  for(k in 1:length(cts)) {
    ixs <- which(meta.data$clusters == cts[k])
    Y.sub <- Y[,ixs]
    meta.sub <- meta.data[ixs,]
    response <- meta.sub$response
    patient <- meta.sub$patient
    p.vals <- apply(Y.sub, 1, function(y) {
      fit <- lmer(y~response + (1|patient))
      summary(fit)$coefficients[2,5]
    })
    perturb[k] <- sum(p.vals < 0.05)
    print(perturb[k])
  }

  return(perturb)
}

set.seed(1)
nct <- 10; cts <- letters[1:10]
reps <- 10
samples <- c(2,3,5,10,20,50,100)
res <- array(0, dim=c(nct, length(samples), reps, 3))
D.true <- rexp(n=10,rate=0.05)
pi.true <- rbeta(n=10,shape1=1, shape2=1)
J <- 50
tau <- 1
G <- 1000

for(i in 1:length(samples)) {
  for(j in 1:reps) {
    Y <- matrix(0,nrow=G,ncol=0)
    meta.data <- data.frame(response=NULL,
                            patient=NULL,
                            clusters=NULL)
    for(k in 1:nct) {
      out <- simCellType(D=D.true[k],J=rpois(n=1,lambda=J),
                         N1=samples[i],N2=samples[i],
                         label=letters[k],tau=tau,
                         my.pi=pi.true[k])
      Y <- cbind(Y,out$Y)
      meta.data <- rbind(meta.data,out$meta.data)
    }

    out <- scDist(Y,meta.data,fixed.effects = "response",
                  random.effects="patient",
                  clusters="clusters",d=2*samples[i]+10)

    res[,i,j,1] <- out$results$Dist.
    res[,i,j,2] <- count_deg(Y, meta.data, cts)

    N1 <- samples[i]; N2 <- samples[i]
    meta.data.a <- data.frame(cell_type=meta.data$clusters,label=as.character(meta.data$response))
    rownames(Y) <- 1:G; colnames(Y) <- 1:ncol(Y)
    auc <- calculate_auc(input=Y,meta=meta.data.a)
    res[,i,j,3] <- as.vector(auc$AUC[order(auc$AUC$cell_type),]$auc)
  }
}

saveRDS(res, "../data/sim_results.RDS")

res <- readRDS("../data/sim_results.RDS")
samples <- c(2,3,5,10,20,50,100)

library(tidyverse)

#plot of cell types paths

df <- reshape2::melt(res)
df <- df |> group_by(Var1, Var2, Var4) |> summarize(mean = mean(value))

df$Var1 <- LETTERS[1:10][df$Var1]
df$Var2 <- samples[df$Var2]
df$Var4 <- c("scDist", "nDEG (lme4)", "Augur")[df$Var4]

p <- df |> ggplot(aes(x=Var2,y=mean,color=Var1)) +
  geom_point() +
  geom_line(linetype="dashed") +
  scale_x_log10() +
  facet_wrap(~Var4, nrow=3, scales="free_y") +
  theme_bw() +
  labs(color="Cell type")+
  xlab("Number of patients") +
  ylab("Perturbation")

ggsave(p,
       filename="../plots/many_samples_comparison.png")

## Root Mean squared error
df <- res[,,,1]
df <- reshape2::melt(df)
df <- df |> mutate(ground_truth=D.true[Var1]) |>
  mutate(se=(value-ground_truth)^2) |>
  group_by(Var2) |>
  summarize(rmse=sqrt(mean(se)))

df$Var2 <- samples[df$Var2]
p <- ggplot(df,aes(x=Var2, y=rmse)) +
  geom_point() + geom_line() +
  scale_x_log10() +
  xlab("Number of patients") +
  ylab("RMSE")+
  theme_bw()

ggsave(p,
       filename="../plots/many_samples_rmse.png")


## Pearson correlation

sp.cor <- apply(res, 2:4, function(x) {
  cor(x, D.true, method="pearson")
})
sp.cor[is.na(sp.cor)] <- 0
df <- reshape2::melt(sp.cor)

df <- df |> group_by(Var1, Var3) |> summarize(mean=mean(value))

df$Var1 <- samples[df$Var1]
df$Var3 <- c("scDist", "nDEG (lme4)", "Augur")[df$Var3]

p <- df |> ggplot(aes(x=Var1,y=mean,color=Var3)) +
  geom_point() +
  geom_line(linetype="dashed") +
  scale_x_log10() +
  theme_bw() +
  labs(color="Cell type")+
  xlab("Number of patients") +
  ylab("Pearson correlation")


ggsave(p,
       filename="../plots/many_samples_cor.png")
