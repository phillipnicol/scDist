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



simData <- function(nct=10, J=50, N1, N2, G=1000, nn=100,tau=0.5) {
  Y <- matrix(0,nrow=G,ncol=0)
  meta.data <- data.frame(response=NULL,
                          patient=NULL,
                          clusters=NULL)
  D.true <- rep(0,nct)
  for(i in 1:nct) {
    D.true[i] <- rexp(n=1,rate=0.05)
    out <- simCellType(D=D.true[i],J=rpois(n=1,lambda=J),N1=N1,N2=N2,label=letters[i],tau=tau)
    Y <- cbind(Y,out$Y)
    meta.data <- rbind(meta.data,out$meta.data)
  }

  out$Y <- Y
  out$meta.data <- meta.data
  out$D.true <- D.true
  return(out)
}
