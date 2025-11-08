#'
#' @title Condition vs Cell-type distances
#'
#' @description Compare the condition perturbation to the estimated distance between
#' different cell types
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
#' @param weights An optional vector of length equal to the number of genes specifying the weight
#' to place on each gene in the distance estimate.
#'
#' @return A list with components
#' \itemize{
#' \item \code{p} - A boxplot comparing the condition perturbation to the
#' inter-celltype distances.
#' \item \code{D} - A distance matrix for the cell types.
#' }
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
condition_v_celltype_distance <- function(normalized_counts,
                                          meta.data,
                                          fixed.effects,
                                          random.effects=c(),
                                          clusters,
                                          d=20,
                                          truncate=FALSE,
                                          min.count.per.cell=20,
                                          weights=NULL) {
  require(tidyverse)

  out.full <- scDist(normalized_counts=normalized_counts,
                     meta.data = meta.data,
                     fixed.effects = fixed.effects,
                     random.effects = random.effects,
                     clusters = clusters,
                     d=d,
                     truncate=truncate,
                     min.count.per.cell=min.count.per.cell,
                     weights = weights)

  meta.data$cell.type <- meta.data[[clusters]]
  meta.data$cluster <- "A"

  cell.types <- unique(meta.data$cell.type)
  nc <- length(cell.types)
  D <- matrix(0, nrow=nc, ncol=nc)
  for(i in 1:nc) {
    for(j in 1:nc) {
      if(i < j) {
        next
      } else if(i == j) {
        next
      }
      ixs <- which(meta.data$cell.type %in% c(cell.types[i], cell.types[j]))
      expr.sub <- as.matrix(normalized_counts[,ixs])
      meta.sub <- meta.data[ixs,]
      cat(i, " ", j, "\n")
      try({
        out <- scDist(expr.sub,meta.sub,fixed.effects="cell.type",
                      random.effects=random.effects,
                      clusters="cluster")
        D[i,j] <- out$results$Dist.
      })
    }
  }

  rownames(D) <- cell.types
  colnames(D) <- rownames(D)
  out_names <- rownames(out.full[["results"]])
  # D <- D[rownames(out.full$results),rownames(out.full$results)]
  D <- D[out_names,out_names]
  D <- D + t(D)
  D2 <- cbind(D, "condition"=out.full$results[out_names,"Dist."])
  # colnames(D)[14] <- "condition"  ## this may cause error for fixed num

  df <- reshape2::melt(D2)

  df.sub <- df |> filter(value > 10^{-10}) |>
    filter(Var1 != "condition") |>
    filter(Var2 != "condition")

  df.sub2 <- df |> filter(Var2 == "condition")

  p <- ggplot(data=df.sub,aes(x=Var1, y = value)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    geom_point(data=df.sub2, aes(x=Var1, y=value),#pch="*",
               color="red", size=2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Cell type") + ylab("Distance")
  out <- list()
  out$D <- D
  out$plot <- p
  return(out)
}
