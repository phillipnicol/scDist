
#' @export
pcDiffPop <- function(normalized_counts,
                      meta.data,
                      fixed.effects,
                      random.effects=c(),
                      clusters,
                      d=50,
                      truncate=FALSE) {
  # Save relevant info to variables
  design <- makeDesign(fixed.effects,random.effects)
  print(design)
  design.null <- makeDesign(fixed.effects[-1],random.effects)
  print(design.null)

  RE <- TRUE
  if(length(random.effects)==0) {
    RE <- FALSE
  }

  #get relevant metadata
  meta.cols <- vapply(c(fixed.effects,random.effects),function(x) {
    which(colnames(meta.data)==x)
  }, integer(1))
  data <- meta.data[,meta.cols, drop=FALSE]
  print(head(data))
  data$y <- rep(0,nrow(data))

  clusters <- meta.data[[clusters]]
  all_clusters <- sort(unique(clusters))
  distances <- c()
  res <- matrix(0,nrow=0,ncol=6)
  out <- list()
  out$vals <- list()
  p <- c()
  for(i in all_clusters) {
    print(i)
    ix <- which(clusters==i)
    normalized_counts.sub <- normalized_counts[ix,]
    #pca.sub <- prcomp(normalized_counts.sub)
    pca.sub <- irlba::prcomp_irlba(x=normalized_counts.sub,n=d)
    data.sub <- data[ix,]
    vals <- pcDiff(pca.sub,data.sub,design,design.null,d,RE,truncate)
    vals$loadings <- pca.sub$rotation
    out$vals[[i]] <- vals
    res <- rbind(res,c(vals$D.hat,
                       vals$D.se,
                       vals$W,
                       vals$LRT.p,
                       vals$F.p,
                       sum(vals$BioVar)/sum(vals$lambda^2)))
  }
  res <- as.data.frame(res)
  rownames(res) <- all_clusters
  colnames(res) <- c("Dist.",
                     "S.e.",
                     "Stat.",
                     "p.LRT",
                     "p.F",
                     "BioVar")
  out$results <- res
  out$design <- design
  return(out)
}

pcDiff <- function(pca,
                   data,
                   design,
                   design.null,
                   d,
                   RE,
                   truncate) {
  beta <- rep(0,d)
  beta_sd <- rep(0,d)
  dfs <- rep(0,d)
  likelihood_full <- 0
  likelihood_null <- 0
  BioVar <- rep(0,d)
  CI_lb <- rep(0,d)
  CI_ub <- rep(0,d)
  for(i in 1:d) {
    data$y <- pca$x[,i]
    if(RE) {
      fit <- lmer(formula=design,data=data,REML=TRUE)
      sumfit <- summary(fit)
      dfs[i] <- sumfit$coefficients[2,3]
      dfs[i] <- sumfit$coefficients[2,3]

      X <- model.matrix(fit)

      a <- fit@beta[1]
      b <- fit@beta[2]

      numerator <- a + X[,2]*b
      numerator <- var(numerator)
      denominator <- numerator+fit@theta^2+fit@sigma^2
      R2 <- numerator/denominator
      BioVar[i] <- R2*pca$sdev[i]^2

      fit.full <- lmer(formula=design,data=data,REML=FALSE)
      fit.null <- lmer(formula=design.null,data=data,REML=FALSE)
      likelihood_full <- likelihood_full+logLik(fit.full)
      likelihood_null <- likelihood_null+logLik(fit.null)

    } else {
      fit <- lm(formula=design,data=data)
      dfs[i] <- nrow(data)-length(fit$coefficients)
      sumfit <- summary(fit)
    }
    beta[i] <- sumfit$coefficients[2,1]
    beta_sd[i] <- sumfit$coefficients[2,2]
  }

  LRT = -2*(likelihood_null-likelihood_full)

  #CI
  CI_total_lb <- sum(CI_lb^2)
  CI_total_ub <- sum(CI_ub^2)

  #Estimate D
  D.hat <- sum(beta^2-beta_sd^2)
  if(truncate) {
    D.hat <- max(D.hat,0)
  }

  #Estimate standard error
  #D.se <- sqrt(sum(3*beta_sd^4))
  B2 <- ifelse(beta^2 - beta_sd^2 < 0, 0, beta^2-beta_sd^2)
  D.se <- sqrt(sum(4*B2*beta_sd^2))

  CI_total_lb <- qnorm(0.025,mean=D.hat,sd=D.se)
  CI_total_ub <- qnorm(0.975,mean=D.hat,sd=D.se)

  #Wald stat
  W <- sum((beta/beta_sd)^2)

  #Monte carlo p-value
  mcreps <- 10^4
  mysum <- rep(0,mcreps)
  for(i in 1:d) {
    myf <- rf(mcreps,df1=1,df2=dfs[i])
    mysum <- mysum+myf
  }
  p <- (sum(mysum > W)+1)/(mcreps+1)

  out <- list()
  out$D.hat <- D.hat
  out$D.se <- D.se
  out$W <- W
  out$LRT.p <- pchisq(LRT,df=d,lower.tail=FALSE)
  out$beta <- beta
  out$beta_sd <- beta_sd
  out$lambda <- pca$sdev[1:d]
  out$scores <- pca$x[,1:d]
  out$data <- data
  out$dfs <- dfs
  out$CI <- c(CI_total_lb, CI_total_ub)
  out$F.p <- p
  out$BioVar <- BioVar
  out
}


makeDesign <- function(fixed.effects, random.effects) {
  design <- "y~"
  for(i in fixed.effects) {
    if(i != fixed.effects[1]) {
      design <- paste0(design,"+",i)
    } else {
      design <- paste0(design, i)
    }
  }
  for(i in random.effects) {
    if(length(fixed.effects) > 0) {
      design <- paste0(design,"+(1|",i,")")
    } else {
      design <- paste0(design,"(1|",i,")")
    }
  }
  if(sum(c(length(fixed.effects),length(random.effects)))==0) {
    design <- "y~1"
  }
  as.formula(design)
}
