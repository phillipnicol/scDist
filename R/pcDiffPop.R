

pcDiffPop <- function(sc.object,
                      fixed.effects,
                      random.effects,
                      clusters) {
  pca <- sc.object@reductions$pca@cell.embeddings
  npcs <- ncol(pca)
  design <- makeDesign(fixed.effects,random.effects)
  weights <- sc.object@reductions$pca@stdev

  #get data
  meta.cols <- vapply(c(fixed.effects,random.effects),function(x) {
    which(colnames(sc.object@meta.data)==x)
  }, integer(1))
  data <- sc.object@meta.data[,meta.cols]
  data$y <- rep(0,nrow(data))

  clusters <- as.vector(unlist(sc.object[[clusters]]))
  all_clusters <- sort(unique(clusters))
  print(all_clusters)
  distances <- c()
  res <- matrix(0,nrow=0,ncol=3)
  out <- list()
  out$vals <- list()
  p <- c()
  for(i in all_clusters) {
    ix <- which(clusters==i)
    pca.sub <- pca[ix,]
    data.sub <- data[ix,]
    vals <- pcDiff(pca.sub,data.sub,design)
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

pcDiff <- function(pca, data, design) {
  d <- ncol(pca)
  beta <- rep(0,d); p <- rep(1,d)
  for(i in 1:d) {
    data$y <- pca[,i]
    try({
      fit <- lmer(formula=design,data=data)
      sumfit <- summary(fit)
      beta[i] <- sumfit$coefficients[2,1]
      p[i] <- sumfit$coefficients[2,5]
    },silent=TRUE)
  }
  combinep <- EmpiricalBrownsMethod::empiricalBrownsMethod(data_matrix=t(pca),
                                               p_values=p)
  out <- list()
  out$beta <- beta; out$p <- p; out$combinep <- combinep
  out
}

combineP <- function(vals, npcs) {
  #Fisher's method
  p.val <- lapply(vals, function(x) {
    min(1,min(npcs*x$p))
    #stat <- -2*sum(log(x$p))
    #pchisq(stat,df=2*npcs,lower.tail = FALSE)
  })
  p.val
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
