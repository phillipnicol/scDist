

condition_v_celltype_distance <- function(normalized_counts,
                                          meta.data,
                                          fixed.effects,
                                          random.effects=c(),
                                          clusters,
                                          d=20,
                                          truncate=FALSE,
                                          min.counts.per.cell=20,
                                          weights=NULL) {
  out.full <- scDist(normalized_counts=normalized_counts,
                     meta.data = meta.data,
                     fixed.effects = fixed.effects,
                     random.effects = random.effects,
                     clusters = clusters,
                     d=d,
                     truncate=truncate,
                     min.counts.per.cell=min.counts.per.cell,
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
  D <- D[rownames(out.full$results),rownames(out.full$results)]

  D <- D + t(D)

  D <- cbind(D, out.full$results$Dist.)
  colnames(D)[14] <- "condition"

  df <- reshape2::melt(D)

  library(tidyverse)

  df.sub <- df |> filter(value > 10^{-10}) |>
    filter(Var1 != "condition") |>
    filter(Var2 != "condition")

  df.sub2 <- df |> filter(Var2 == "condition")

  p <- ggplot(data=df.sub,aes(x=Var1, y = value)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    geom_point(data=df.sub2, aes(x=Var1, y=value), pch="*",
               color="red", size=5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Cell type") + ylab("Distance")

  out <- list()
  out$D <- D
  out$plot <- plot
  return(out)
}
