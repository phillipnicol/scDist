
makeCellTree <- function(sco) {
  cell_types <- unique(sco@meta.data$cellType)
  group <- rep(1, length(cell_types))
  TreeMat <- matrix(1, nrow = length(cell_types), ncol = 1)
  counter <- 1

  while(length(unique(group)) != length(cell_types)) {
    #Get group with most clusters
    ixs <- which(group == Mode(group))
    print(group)
    cat("Current group: ", cell_types[ixs], "\n")
    split <- splitGroup(sco[,sco$cellType %in% cell_types[ixs]], cell_types[ixs])

    scd.meta$cluster.new <- scd.meta$cluster
    for(j in 1:nrow(sco@meta.data)) {
      scd.meta$cluster.new[j] <- which(scd.meta$cluster[j] == cell_types[ixs])
    }

    if(length(unique(cell_types[ixs])) > 0) {
      out <- scDist(fixed.effects="response",
                    random.effects="patient",
                    cluster="cluster.new")
      print(out$results)
    }

    if(length(unique(split)) > 1) {
      counter <- counter + 1
      subgroup <- group[ixs]
      subgroup[split == 2] <- counter
      counter <- counter + 1
      subgroup[split == 1] <- counter
      group[ixs] <- subgroup

      TreeMat <- cbind(TreeMat, group)
    }
    else {
      for(j in ixs) {
        counter <- counter + 1
        group[j] <- counter
        TreeMat <- cbind(TreeMat, group)
      }
    }
  }

  Tree <- matrix(nrow = 0, ncol = 2)

  counter <- 1

  for(i in 2:ncol(TreeMat)) {
    newvals <- unique(TreeMat[TreeMat[,i] > counter,i])
    for(j in newvals) {
      counter <- counter + 1
      ixs <- which(TreeMat[,i] == counter)
      newvec <- c(counter, TreeMat[ixs[1], i-1])
      Tree <- rbind(Tree, newvec)
    }
  }

  ### Process Tree

  #Change names of leaves
  cntr <- 1
  for(ct in TreeMat[,ncol(TreeMat)]) {
    ix <- which(Tree[,1] == ct)
    Tree[ix,1] <- cell_types[cntr]
    cntr <- cntr+1
  }

  # Also give the columns a name
  colnames(Tree) <- c("Child", "Parent")
  rownames(Tree) <- c(1:nrow(Tree))

  return(Tree)
}

splitGroup <- function(sco_sub, ixs) {
  downsample <- min(100, nrow(sco_sub@meta.data))
  down_ixs <- c()
  for(i in ixs) {
    rxs <- which(sco_sub@meta.data$cellType == i)
    downsample <- min(100, length(rxs))
    rxs <- rxs[sample(length(rxs), downsample, replace = FALSE)]
    down_ixs <- c(down_ixs, rxs)
  }
  sco_sub <- sco_sub[,down_ixs]
  if(length(ixs) == 2) {
    return(c(1,2))
  }

  sco_sub <- Seurat::NormalizeData(sco_sub)

  #Find variable features
  sco_sub <- Seurat::FindVariableFeatures(sco_sub, verbose = FALSE)

  #Scale data
  sco_sub <- Seurat::ScaleData(sco_sub, verbose = FALSE)

  #Run PCA on the variable features. Get 50 dimensional embeddings
  sco_sub <- Seurat::RunPCA(sco_sub, verbose = FALSE)
  embeddings <- Seurat::Embeddings(object = sco_sub, reduction = 'pca')[,1:50]

  pseudobulk <- matrix(0, nrow = 0, ncol = 50)
  for(i in ixs) {
    rxs <- which(sco_sub@meta.data$cellType == i)
    pseudobulk <- rbind(pseudobulk, colMeans(embeddings[rxs,]))
  }
  #Cluster via k means
  km <- kmeans(pseudobulk, centers = 2, iter.max = 10)$cluster

  return(km)
}

Mode <- function(x) {
  xu <- unique(x)
  return(xu[which.max(tabulate(match(x, xu)))])
}

DFS <- function(Tree, node, cell_types) {
  if(node == 0) {
    return(cell_types)
  }

  leaves <- c()

  row <- Tree[node,]

  newparent <- row[1]

  allchild <- which(Tree[,2] == row[1])

  for(child in allchild) {
    leaves <- c(leaves, DFS(Tree, child))
  }

  if(length(allchild) == 0) {
    return(row[1])
  }

  return(leaves)
}

getPseudoBulkCounts <- function(sco, subtree, response) {
  counts <- sco@assays$RNA@counts
  pseudobulk <- matrix(0, nrow = nrow(counts), ncol = length(unique(sco@meta.data$patient)))
  rownames(pseudobulk) <- rownames(counts)
  response <- rep(0, ncol(pseudobulk))
  pseudoMetaData <- matrix(0, nrow = ncol(pseudobulk), ncol = ncol(sco@meta.data))
  pseudoMetaData <- as.data.frame(pseudoMetaData); colnames(pseudoMetaData) <- colnames(sco@meta.data)
  for(i in 1:ncol(pseudobulk)) {
    ixs <- which(sco@meta.data$patient == unique(sco@meta.data$patient)[i])
    response[i] <- sco@meta.data$response[ixs[1]]
    pseudoMetaData[i,] <- sco@meta.data[ixs[1],]
  }
  sco_sub <- sco[,sco$cellType %in% subtree]
  counts <- sco_sub@assays$RNA@counts
  for(i in 1:ncol(pseudobulk)) {
    ixs <- which(sco_sub@meta.data$patient == unique(sco@meta.data$patient)[i])
    if(length(ixs) == 0) {
      pseudobulk[,i] <- 1
    }
    else if(length(ixs) == 1) {
      pseudobulk[,i] <- counts[,ixs] + 1
    }
    else {
      pseudobulk[,i] <- Matrix::rowSums(counts[,ixs]) + 1
    }
  }

  colData <- data.frame(gene = c(1:length(response)), response = as.factor(response))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk, colData = colData, design = ~response)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  return(dds)
}

runDESeq <- function(counts, colData, response) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~response)
  dds <- DESeq2::DESeq(dds)
  results <- DESeq2::results(dds)
  return(results)
}

permutation_test <- function(x, nperm, dds, stat, ncores) {
  vec <- c(1:ncol(dds))
  print(dds)
  null_dist <- unlist(parallel::mclapply(c(1:nperm), function(i) {
    vec <- sample(vec, size = length(vec), replace = FALSE)
    dds$response <- dds$response[vec]
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    y <- res$stat; names(y) <- rownames(res); y <- y[names(y) %in% names(x)]; y <- y[names(x)]
    return(sum(x*y))
  }, mc.cores = ncores))
  print(null_dist)
  pval <- (sum(abs(null_dist) > abs(stat)) + 1)/(nperm + 1)
  print("PVAL:"); print(pval)
  return(pval)
}
