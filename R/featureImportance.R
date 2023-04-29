
#' @export
featureImportance <- function(expr,
                              meta.data,
                              scd,
                              cluster,
                              component) {

  ixs <- which(meta.data$cell_type == cluster)
  expr <- expr[,ixs]
  expr_diff <- rep(0, nrow(expr))
  data.sub <- scd$vals[[cluster]]$data
  expr_diff <- vapply(1:nrow(expr),FUN.VALUE=numeric(1),function(i) {
    print(i)
    data.sub$y <- expr[i,]
    fit <- lmerTest::lmer(formula=scd$design,
                data=data.sub)
    fit@beta[2]
  })

  loadings <- scd$vals[[cluster]]$loadings[,component]
  FI <- abs(loadings)*expr_diff
  names(FI) <- rownames(expr)

  df <- data.frame(x=loadings,y=expr_diff)
  up5 <- order(expr_diff)[1:5]
  down5 <- order(expr_diff,decreasing=TRUE)[1:5]
  df$label <- ifelse(1:nrow(expr) %in% c(up5,down5), names(FI), "")
  df$color <- ifelse(1:nrow(expr) %in% c(up5,down5), "red", "grey")
  df$alpha <- ifelse(1:nrow(expr) %in% c(up5,down5), 1, 0.75)

  p <- ggplot(df,aes(x=x,y=y,color=color,alpha=alpha,
                     label=label))
  p <- p+scale_color_manual(values=c("grey","red"))
  p <- p + geom_point()
  p <- p + ggrepel::geom_text_repel(max.overlaps=Inf)
  p <- p + geom_hline(yintercept=0,color="blue",
                      linetype="dashed")
  p <- p + geom_vline(xintercept=0,color="blue",
                      linetype="dashed")
  p <- p + theme_linedraw()
  p <- p + xlab("Loading")+ylab("Normalized expr. difference")
  p <- p + guides(color="none",alpha="none")
  p

  out <- list()
  out$FI <- FI
  out$plot <- p
  return(out)
}
