
#' @export
nullPermutation <- function(sc.object,
                            response,
                            samples,
                            clusters,
                            iters=100,
                            level=seq(0.01,0.1,by=0.01)) {
  orig_clusters <- clusters
  clusters <- as.vector(unlist(sc.object[[clusters]]))
  all_clusters <- unique(clusters)
  vals <- matrix(0, nrow=length(all_clusters),ncol=length(level))
  rc <- which(colnames(sc.object@meta.data) == response)
  sc <- which(colnames(sc.object@meta.data) == samples)
  sc.sub <- sc.object[,sc.object@meta.data[,rc] ==
                        sc.object@meta.data[1,rc]]
  m <- length(unique(sc.sub@meta.data[,sc]))
  false_discov <- matrix(0, nrow=iters,ncol=length(level))
  for(i in 1:iters) {
    sample_perm <- sample(unique(sc.object@meta.data[,sc]),
                          length(unique(sc.object@meta.data[,sc])),
                          replace=FALSE)
    sc.object@meta.data[,rc] <- ifelse(sc.object@meta.data[,sc] %in% sample_perm[1:m],
                                       "1", "0")
    out <- pcDiffPop(sc.object, fixed.effects = response,
                     random.effects = samples,
                     clusters=orig_clusters)
    out$results$p.val <- p.adjust(out$results$p.val, method="fdr")
    print(out$results)
    for(j in 1:length(level)) {
      vals[,j] <- vals[,j] + ifelse(out$results$p.val < level[j], 1, 0)
      false_discov[i,j] <- sum(out$results$p.val < level[j])
    }
  }
  vals <- vals/iters
  rownames(vals) <- rownames(out$results); colnames(vals) <- level
  rownames(false_discov) <- c(1:iters); colnames(false_discov) <- level
  out <- list()
  out$vals <- vals
  out$false_discov <- false_discov/length(all_clusters)
  out
}
