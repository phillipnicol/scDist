

pcDiffPop <- function(sc.object,
                      fixed.effects,
                      random.effects,
                      clusters) {
  pca <- sc.object@reductions$pca@cell.embeddings
  npcs <- ncol(pca)
  design <- makeDesign(fixed.effects,random.effects)

  #get data
  meta.cols <- vapply(c(fixed.effects,random.effects),function(x) {
    which(colnames(sc.object@meta.data)==x)
  }, integer(1))
  print(meta.cols)
  data <- sc.object@meta.data[,meta.cols]
  data$y <- rep(0,nrow(data))

  all_clusters <- unique(clusters)
  res <- matrix(0,nrow=0,ncol=3)
  out <- list()
  out$vals <- list()
  for(i in all_clusters) {
    print(i)
    ix <- which(clusters==i)
    pca.sub <- pca[ix,]
    data.sub <- data[ix,]
    vals <- pcDiff(pca.sub,data.sub,design)
    out$vals[[i]] <- vals
    discov <- sum(vals$p < 0.05)
    pval <- pbinom(discov,size=npcs,prob=0.05,lower.tail=FALSE)
    #Holm's correction
    #pval <- min(p.adjust(vals$p, method="holm"))
    #Fisher's test
    #stat <- -2*sum(log(vals$p))
    #pval <- pchisq(stat,df=2*npcs,lower.tail = FALSE)
    row <- c(i,sqrt(sum(vals$beta^2)),pval)
    res <- rbind(res,row)
  }
  res <- data.frame(row.names=res[,1],
                    distance=res[,2],
                    p.val=p.adjust(res[,3],method="fdr"))
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
  out <- list()
  out$beta <- beta; out$p <- p
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
