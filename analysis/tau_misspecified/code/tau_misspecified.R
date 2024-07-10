library(scDist)

simCellType <- function(D,G=1000,N1=5,N2=5,J=50,label="A",my.pi=0.9,rate=1) {
  beta_true <- rep(0,G)

  beta_true <- rnorm(n=G, mean=0,sd=1)
  beta_true <- D/sqrt(sum(beta_true^2))*beta_true

  print(J)
  y <- matrix(0, nrow=(N1+N2)*J,ncol=G)
  for(i in 1:G) {
    sigma_g <- rgamma(n=1,shape=rate, rate=rate)
    tau_g <- rgamma(n=1,shape=0.5*rate, rate=rate)
    cntr <- 1
    for(j in 1:N1) {
      omega <- rnorm(1,mean=0,sd=tau_g)
      y[cntr:(cntr+J-1),i] <- omega+rnorm(J, sd=sigma_g)
      cntr <- cntr+J
    }
    for(j in 1:N2) {
      omega <- rnorm(1,mean=0, sd=tau_g)
      y[cntr:(cntr+J-1),i] <- beta_true[i]+omega+rnorm(J,sd=sigma_g)
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

simCellType2 <- function(D,G=1000,N1=5,N2=5,J=50,label="A",my.pi=0.9) {
  beta_true <- rep(0,G)

  beta_true <- rnorm(n=G, mean=0,sd=1)
  beta_true <- D/sqrt(sum(beta_true^2))*beta_true

  print(J)
  y <- matrix(0, nrow=(N1+N2)*J,ncol=G)
  for(i in 1:G) {
    cntr <- 1
    for(j in 1:N1) {
      omega <- rnorm(1,mean=0,sd=tau_g)
      y[cntr:(cntr+J-1),i] <- omega+rnorm(J, sd=0.5)
      cntr <- cntr+J
    }
    for(j in 1:N2) {
      omega <- rnorm(1,mean=0, sd=tau_g)
      y[cntr:(cntr+J-1),i] <- beta_true[i]+omega+rnorm(J,sd=1)
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


simData <- function(nct=10, J=50, N1, N2, G=1000, nn=100,rate=rate) {
  Y <- matrix(0,nrow=G,ncol=0)
  meta.data <- data.frame(response=NULL,
                          patient=NULL,
                          clusters=NULL)
  D.true <- rep(0,nct)
  z <- matrix(0,nrow=G,ncol=nct)
  for(i in 1:nct) {
    D.true[i] <- rexp(n=1,rate=0.05)
    out <- simCellType(D=D.true[i],J=rpois(n=1,lambda=J),N1=N1,N2=N2,label=letters[i],rate=rate)
    Y <- cbind(Y,out$Y)
    meta.data <- rbind(meta.data,out$meta.data)
  }

  out$Y <- Y
  out$meta.data <- meta.data
  out$D.true <- D.true
  return(out)
}

simData2 <- function(nct=10, J=50, N1, N2, G=1000, nn=100) {
  Y <- matrix(0,nrow=G,ncol=0)
  meta.data <- data.frame(response=NULL,
                          patient=NULL,
                          clusters=NULL)
  D.true <- rep(0,nct)
  z <- matrix(0,nrow=G,ncol=nct)
  for(i in 1:nct) {
    D.true[i] <- rexp(n=1,rate=0.05)
    out <- simCellType(D=D.true[i],J=rpois(n=1,lambda=J),N1=N1,N2=N2,label=letters[i])
    Y <- cbind(Y,out$Y)
    meta.data <- rbind(meta.data,out$meta.data)
  }

  out$Y <- Y
  out$meta.data <- meta.data
  out$D.true <- D.true
  return(out)
}


## Simulate Data
set.seed(1)
reps <- 20
rate.try <- sqrt(1:4)
res <- matrix(0, nrow = reps, ncol=length(rate.try) + 1)

for(i in 1:reps) {
  for(j in 1:length(rate.try)) {
    #InCorrect model
    sim <- simData(nct=1,J=100,N1=5,N2=5,rate=rate.try[j])

    #Correct
    out <- scDist(normalized_counts = sim$Y,
                  meta.data=sim$meta.data,
                  fixed.effects="response",
                  random.effects="patient",
                  clusters="clusters")
    res[i,j] <- abs(out$results$Dist. - sim$D.true)
  }
}

for(i in 1:reps) {
  sim <- simData2(nct=1,J=100,N1=5,N2=5)
  #Correct
  out <- scDist(normalized_counts = sim$Y,
                meta.data=sim$meta.data,
                fixed.effects="response",
                random.effects="patient",
                clusters="clusters")
  res[i,ncol(res)] <- abs(out$results$Dist. - sim$D.true)
}

saveRDS(res, file="../data/tau_misspecified.RDS")


res <- readRDS(file="../data/tau_misspecified.RDS")

rate.try <- sqrt(1:4)

var.tested <- c(1/(rate.try^2), 0)

library(tidyverse)

df <- reshape2::melt(res)
df$Var2 <- factor(round(var.tested[df$Var2],2))

p <- ggplot(data=df,aes(x=Var2, y=value)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, alpha = 0.5)  +
  theme_bw() +
  xlab("Variance of sigma_g") +
  ylab("Absolute error") +
  theme(legend.position = "none")



ggsave(p, filename="../plots/tau_misspecified.png")




