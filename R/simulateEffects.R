

simulateEffects <- function(sc.object,
                            bio_effect=0,
                            sample_effect=0,
                            nsamples.0=5,
                            nsamples.1=5,
                            ncells=50,
                            augur=TRUE) {
  sc.object$response <- "0"
  sc.object$patient <- "P1"
  sc.object$cell_type <- "C1"

  sc.object <- sc.object[,1:ncells]
  sc.1 <- sc.object
  m <- nsamples.0; n <- nsamples.1
  for(i in 2:m) {
    ixs <- sample(1:ncells, size=ncells, replace=TRUE)
    sc.1 <- sc.1[,ixs]
    sc.1@meta.data$response <- "0"
    sc.1@meta.data$patient <- paste0("P",i)
    sc.1@meta.data$cell_type <- "C1"
    sc.object <- merge(sc.object,sc.1)
  }
  for(i in (m+1):(m+n)) {
    sc.1@meta.data$response <- "1"
    sc.1@meta.data$patient <- paste0("P",i)
    sc.object <- merge(sc.object,sc.1)
  }
  #Process data
  sc.object <- Seurat::NormalizeData(sc.object)
  sc.object <- Seurat::FindVariableFeatures(sc.object)
  sc.object <- Seurat::ScaleData(sc.object)
  sc.object <- Seurat::RunPCA(sc.object)

  results <- matrix(0,nrow=length(bio_effect),
                    ncol=length(sample_effect))
  rownames(results) <- bio_effect
  colnames(results) <- sample_effect

  betaMatrix <- results; pMatrix <- results
  augurMatrix <- results

  for(i in 1:nrow(results)) {
    for(j in 1:ncol(results)) {
      sim <- testSimulation(sc.object,
                            bio_effect[i],
                            sample_effect[j],
                            augur)
      betaMatrix[i,j] <- sim$beta
      pMatrix[i,j] <- sim$p
      if(augur) {
        augurMatrix[i,j] <- sim$augur
      }
    }
  }

  out <- list()
  out$betaMatrix <- betaMatrix
  out$pMatrix <- pMatrix
  out$augurMatrix <- augurMatrix
  return(out)
}

testSimulation <- function(sc.object, bio, pat, augur_flag) {
  pca <- sc.object@reductions$pca@cell.embeddings
  lambda <- apply(pca,2,var)
  d <- ncol(pca)

  #biological effect
  ixs <- which(sc.object@meta.data$response=="1")
  for(i in 1:d) {
    pca[ixs,i] <- pca[ixs,i]+bio*sqrt(lambda[i])
  }

  for(patient in unique(sc.object@meta.data$patient)) {
    ixs <- which(sc.object@meta.data$patient == patient)
    for(j in 1:d) {
      omega <- rnorm(1,mean=0,sd=pat*lambda[j])
      pca[ixs,j] <- pca[ixs,j]+omega
    }
  }

  sc.object@reductions$pca@cell.embeddings <- pca

  #expr <-sc.object@reductions$pca@feature.loadings %*% t(sc.object@reductions$pca@cell.embeddings)
  expr <- t(sc.object@reductions$pca@cell.embeddings)
  augur.md <- data.frame(cell_type="A",label=sc.object$response)
  if(augur_flag) {
    augur <- calculate_auc(expr,meta=augur.md,n_threads=4)
  }

  pcdp <- pcDiffPop(sc.object,fixed.effects="response",
                   random.effects="patient",
                   clusters="cell_type",
                   shrinkage=TRUE)
  out <- list()
  out$beta <- pcdp$results$distance[1]
  out$p <- pcdp$results$p.val[1]
  if(augur_flag) {
    out$augur <- augur$AUC$auc[1]
  }
  return(out)
}
