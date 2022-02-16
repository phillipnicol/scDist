#' @export
pcDiffPlot <- function(pcdp.object) {
  results <- pcdp.object$results
  results <- results[order(results$Dist,decreasing=FALSE),]

  df <- data.frame(cell_type=rownames(results),
                   dist=results$Dist,
                   se_up=results$Dist+results$S.e.,
                   se_down=ifelse(results$Dist-results$S.e. < 0,
                                  0,
                                  results$Dist-results$S.e.))
  df$cell_type <- factor(df$cell_type,levels=df$cell_type)

  p <- ggplot(data=df,aes(x=cell_type,y=dist))
  p <- p + geom_errorbar(aes(ymin=se_down,ymax=se_up))
  p <- p + geom_point(color="red")
  p <- p + xlab("Cell type")+ylab("Dist.")
  p <- p + theme_linedraw()
  p

}





#' @export
plotBetas <- function(pcdp, cluster) {
  ix <- which(names(pcdp$vals) == cluster)
  df <- data.frame(beta=pcdp$vals[[ix]]$beta,
                   p=pcdp$vals[[ix]]$F.p)
  df$dim <- c(1:length(df$p))
  df$signif <- apply(df,1,function(x) {
    if(x["p"] < 0.001/nrow(df)) {
      "***"
    }
    else if(x["p"] < 0.01/nrow(df)) {
      "**"
    }
    else if(x["p"] < 0.05/nrow(df)) {
      "*"
    }
    else {" "}
  })
  p <- ggplot2::ggplot(df,ggplot2::aes(x=dim,y=beta))
  p <- p + geom_segment(ggplot2::aes(x=dim,xend=dim,y=0,yend=beta))
  p <- p + ggplot2::geom_point(size=2,color="purple")
  df$beta <- ifelse(df$beta > 0, df$beta+0.35, df$beta-0.35)
  p <- p + ggplot2::geom_text(data=df, ggplot2::aes(label=signif,x=dim,y=beta))
  p <- p + ggplot2::labs(x="PC Dimension", y=expression(beta[k]))
  p <- p + ggplot2::ggtitle(label=cluster)
  p <- p + ggplot2::theme_linedraw()
  p
}


