#' @import ggraph
NULL

#' @import Seurat
NULL


#' @export
#'
#' @title Tree-based analysis of scDist distances
#'
#' @description Finds a cluster tree and then uses
#' the scDist procedure to estimate the distance at
#' each node of the tree
#'
#' @param sco A `seurat` object containing the
#' following columns in meta data: `cellType`,
#' `condition`, and patient.
scDist_tree <- function(sco) {
  tree <- makeCellTree(sco)

  Tree <- tree$Tree

  sco$identity <- "A"
  nodes <- rep(0, length(unique(as.vector(Tree))))
  names(nodes) <- unique(as.vector(Tree))
  for(i in 1:nrow(Tree)) {
    subtree <- as.vector(DFS(Tree, i, cell_types))
    #print(subtree)

    ixs <- which(sco$cellType %in% subtree)
    out <- scDist(sco@assays$SCT@scale.data[,ixs],
                  meta.data=sco@meta.data[ixs,],
                  fixed.effects="condition",
                  random.effects="patient",
                  cluster="identity")
    #print(out$results)

    nodes[Tree[i,1]] <- out$results$Dist.
  }

  TreeIG <- cbind(Tree[,2], Tree[,1])
  TreeIG <- as.matrix(TreeIG)

  G <- igraph::graph_from_edgelist(TreeIG)

  xy <- layout_as_tree(G)
  V(G)$x <- xy[, 1]
  V(G)$y <- xy[, 2]

  graph_data <- data.frame(name = V(G)$name, x = V(G)$x, y = V(G)$y)
  graph_data$color <- nodes[V(G)$name]

  p <- ggraph(G, "manual", x= V(G)$x, y=V(G)$y) + geom_edge_link()
  #p <- p + geom_node_circle(aes(x0=x,y0=y, r=0.1), colour = "black", fill="white")
  p <- p + geom_node_circle(aes(x0=x,y0=y, fill=color, r=0.1), colour = "black",
                            data=graph_data)
  p <- p + scale_fill_gradient(low = "white", high = "turquoise")
  p <- p + geom_node_label(aes(label = name, angle = 90),
                           repel = FALSE, nudge_y = 0.25,
                           col = "midnightblue")
  p <- p + labs(fill="Distance")+theme_graph()+ theme(legend.position = "bottom")

  return(p)
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



makeCellTree <- function(sco) {
  cell_types <- unique(sco@meta.data$cellType)
  group <- rep(1, length(cell_types))
  TreeMat <- matrix(1, nrow = length(cell_types), ncol = 1)
  counter <- 1

  my.dist <- matrix(0,nrow=0,ncol=2)

  while(length(unique(group)) != length(cell_types)) {
    #Get group with most clusters
    ixs <- which(group == Mode(group))
    print(group)
    cat("Current group: ", cell_types[ixs], "\n")
    split <- splitGroup(sco[,sco$cellType %in% cell_types[ixs]], cell_types[ixs])

    cluster.new <- vapply(sco$cellType, FUN.VALUE = numeric(1), function(ct) {
      ix <- which(ct == cell_types)
      group[ix]
    })

    sco$cluster.new <- cluster.new

    ##out <- scDist(sco@assays$SCT@scale.data,
    #              meta.data=sco@meta.data,
    #              fixed.effects="Status",
    #              random.effects="Donor",
    #              cluster="cluster.new")

    #for(j in 1:nrow(out$results)) {
    #  if(rownames(out$results)[j] %in% my.dist[,1]) {
    #    next
    #  } else {
    #    new.row <- c(rownames(out$results)[j], out$results$Dist.[j])
    #    my.dist <- rbind(my.dist, new.row)
    #  }
    #}

    #print(out$results)


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
      print("WARNING ")
      for(j in ixs) {
        group[j] <- counter
        counter <- counter + 1
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

  out <- list()
  out$Tree <- Tree

  return(out)
}

