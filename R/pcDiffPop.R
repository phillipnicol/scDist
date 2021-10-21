
#' @export
pcDiffPop <- function(pca,
                      meta.data,
                      fixed.effects,
                      random.effects,
                      clusters,
                      shrinkage=FALSE) {
  # Save relevant info to variables
  npcs <- ncol(pca)
  design <- makeDesign(fixed.effects,random.effects)

  #get relevant metadata
  meta.cols <- vapply(c(fixed.effects,random.effects),function(x) {
    which(colnames(meta.data)==x)
  }, integer(1))
  data <- meta.data[,meta.cols]
  data$y <- rep(0,nrow(data))

  clusters <- meta.data[[clusters]]
  all_clusters <- sort(unique(clusters))
  distances <- c()
  res <- matrix(0,nrow=0,ncol=3)
  out <- list()
  out$vals <- list()
  p <- c()
  for(i in all_clusters) {
    ix <- which(clusters==i)
    pca.sub <- pca[ix,]
    data.sub <- data[ix,]
    vals <- pcDiff(pca.sub,data.sub,design,shrinkage)
    out$vals[[i]] <- vals
    p <- c(p, vals$combinep)
    distances <- c(distances, sqrt(sum(vals$beta^2)))
  }
  res <- data.frame(row.names = all_clusters,
                    distance = distances,
                    p.val = p)
  out$results <- res
  out$design <- design
  return(out)
}

pcDiff <- function(pca, data, design, shrinkage) {
  d <- ncol(pca)
  beta <- rep(0,d); p <- rep(1,d)
  for(i in 1:d) {
    data$y <- pca[,i]
    try({
      fit <- lmer(formula=design,data=data)
      sumfit <- summary(fit)
      beta[i] <- sumfit$coefficients[2,1]
      p[i] <- sumfit$coefficients[2,5]
      if(shrinkage) {
        stderr <- sumfit$coefficients[2,2]
        beta[i] <- beta[i]/(stderr^2 + 1)
      }
    })
  }
  #Combine p values using empirical browns method
  combinep <- EmpiricalBrownsMethod::empiricalBrownsMethod(data_matrix=t(pca),
                                               p_values=p)
  out <- list()
  out$beta <- beta; out$p <- p; out$combinep <- combinep
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
    design <- paste0(design,"+(1|",i,")")
  }
  as.formula(design)
}
