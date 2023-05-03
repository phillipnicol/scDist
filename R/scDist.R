
#' @export
#'
#' @title scDist: Identify perturbed cell types in
#' single-cell RNA-seq data
#'
#' @description Estimate the distance between
#' condition means in gene expression space.
#'
#' @param normalized_counts A matrix containing
#' normalized data with genes on rows and cells
#' on columns
#' @param meta.data A data frame containing meta data for each cell.
#' @param fixed.effects The columns in meta.data corresponding to the fixed effects. In
#' a typical case, this would be the condition of interest.
#' @param random.effects The columns in meta.data corresponding to the random effects.
#' In a typical use case this would be the column containing patient ID.
#' @param clusters The column containing the cell-type annotation.
#' @param d The number of PCs to use.
#' @param truncate Whether or not to round negative distances to 0.
#' @param min.counts.per.cell The minimum number of cells per cluster to perform the estimation.
#'
#' @return A list with components
#' \itemize{
#' \item \code{results} - A data frame containing the cell
#' type, estimated distance, and other statistics such as p-value.
#' \item \code{vals} For each cell type a list of more detailed
#' information (such as raw data) and coefficients for each PC are
#' included.
#' }
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
scDist <- function(normalized_counts,
                      meta.data,
                      fixed.effects,
                      random.effects=c(),
                      clusters,
                      d=20,
                      truncate=FALSE,
                      min.counts.per.cell=20) {
  #Normalized counts currently in cells x genes
  normalized_counts <- t(normalized_counts)

  # Save relevant info to variables
  design <- makeDesign(fixed.effects,random.effects)
  design.null <- makeDesign(fixed.effects[-1],random.effects)

  RE <- TRUE
  if(length(random.effects)==0) {
    RE <- FALSE
  }

  #get relevant metadata
  meta.cols <- vapply(c(fixed.effects,random.effects),function(x) {
    which(colnames(meta.data)==x)
  }, integer(1))
  data <- meta.data[,meta.cols, drop=FALSE]
  data$y <- rep(0,nrow(data))

  clusters <- meta.data[[clusters]]
  all_clusters <- sort(unique(clusters))
  distances <- c()
  res <- matrix(0,nrow=0,ncol=4)
  out <- list()
  out$vals <- list()
  p <- c()
  bar <- txtProgressBar(min=0,max=length(all_clusters),initial = 0)
  cntr <- 1
  for(i in all_clusters) {
    setTxtProgressBar(bar,cntr)
    cntr <- cntr+1
    ix <- which(clusters==i)
    normalized_counts.sub <- normalized_counts[ix,]
    #pca.sub <- prcomp(normalized_counts.sub)
    if(length(ix) < min.counts.per.cell) {
      next
    }
    pca.sub <- irlba::prcomp_irlba(x=normalized_counts.sub,n=d)
    data.sub <- data[ix,]
    vals <- pcDiff(pca.sub,data.sub,design,design.null,d,RE,truncate)
    vals$loadings <- pca.sub$rotation
    out$vals[[i]] <- vals
    res <- rbind(res,c(vals$D.post.med,
                       vals$D.post.lb,
                       vals$D.post.ub,
                       vals$p.sum))
  }
  res <- as.data.frame(res)
  rownames(res) <- all_clusters
  colnames(res) <- c("Dist.",
                     "95% CI (low)",
                     "95% CI (upper)",
                     "p.val")
  out$results <- res
  out$design <- design

  close(bar)
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
  for(i in 1:d) {
    data$y <- pca$x[,i]
    if(RE) {
      fit <- lmerTest::lmer(formula=design,data=data,REML=TRUE,
                  control = lme4::lmerControl(check.conv.singular="ignore"))
      sumfit <- summary(fit)
      dfs[i] <- sumfit$coefficients[2,3]
      dfs[i] <- sumfit$coefficients[2,3]

    } else {
      fit <- lm(formula=design,data=data)
      dfs[i] <- nrow(data)-length(fit$coefficients)
      sumfit <- summary(fit)
    }
    beta[i] <- sumfit$coefficients[2,1]
    beta_sd[i] <- sumfit$coefficients[2,2]
  }

  beta.ash <- ashr::ash(betahat=beta,sebetahat = beta_sd)
  beta.post <- ashr::get_post_sample(beta.ash,10^4)
  D.post <- apply(beta.post,1,function(b) {
    sqrt(sum(b^2))
  })
  D.post.med <- median(D.post)
  D.post.lb <- quantile(D.post, 0.025)
  D.post.ub <- quantile(D.post,0.975)

  #Estimate D
  D.hat <- sum(beta^2-beta_sd^2)
  if(truncate) {
    D.hat <- max(D.hat,0)
  }

  #Estimate standard error
  beta2.var <- 4*(beta^2-beta_sd^2)*beta_sd^2+2*beta_sd^2
  D.se <- sum(beta2.var)
  D.se <- ifelse(D.se < 0, 0, sqrt(D.se))

  #Wald stat
  W <- sum((beta/beta_sd)^2)
  W.max <- max((beta/beta_sd)^2)

  #Monte carlo p-value
  mcreps <- 10^5
  mymax <- rep(0,mcreps)
  mysum <- rep(0,mcreps)
  for(i in 1:d) {
    myf <- rf(mcreps,df1=1,df2=dfs[i])
    mysum <- mysum+myf
    mymax <- ifelse(myf > mymax,myf,mymax)
  }
  p.sum <- (sum(mysum > W)+1)/(mcreps+1)
  p.max <- (sum(mymax > W.max)+1)/(mcreps+1)

  out <- list()
  out$D.hat <- D.hat
  out$D.se <- D.se
  out$W <- W
  out$beta <- beta
  out$beta_sd <- beta_sd
  out$lambda <- pca$sdev[1:d]
  out$scores <- pca$x[,1:d]
  out$data <- data
  out$dfs <- dfs
  out$p.sum <- p.sum
  out$p.max <- p.max
  out$D.post.med <- D.post.med
  out$D.post.lb <- D.post.lb
  out$D.post.ub <- D.post.ub
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
