#' @import ggplot2
NULL

#' @importFrom ggrepel geom_text_repel
NULL

#' @export
#'
#' @title Plot the results of scDist
#'
#' @description Plot the distance estimates and corresponding standard
#' errors
#'
#' @param scd.object A list obtained from applying the main function \code{\link{scDist}}.
#'
#' @return If return.plot is true, a ggplot object is returned. Otherwise,
#' nothing is returned.
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>

DistPlot <- function(
    scd.object,
    return.plot=FALSE
) {
  results <- scd.object$results
  results <- results[order(results$Dist.,decreasing=FALSE),]

  df <- data.frame(cell_type=rownames(results),
                   dist=results$Dist.,
                   se_up=results$`95% CI (upper)`,
                   se_down=results$`95% CI (low)`)
  df$cell_type <- factor(df$cell_type,levels=df$cell_type)

  p <- ggplot(data=df,aes(x=cell_type,y=dist))
  p <- p + geom_errorbar(aes(ymin=se_down,ymax=se_up))
  p <- p + geom_point(color="red")
  p <- p + xlab("Cell type")+ylab("Dist.")
  p <- p + theme_linedraw()
  p <- p + coord_flip()

  return(p)
}

#' FDRDistPlot
#'
#' ADD DESCRIPTION HERE
#'
#' @param scd.object A list obtained from applying the main function \code{\link{scDist}}.
#'
#' @export

FDRDistPlot <- function(
    scd.object
) {
  p.value <- p.adjust(scd.object$results$p.val,method="fdr")
  dist <- scd.object$results$Dist.

  df <- data.frame(x=dist,y=-log10(p.value))
  df$color <- ifelse(df$y > 1, "orange", "grey")
  df$label <- rownames(scd.object$results)

  p <- ggplot(df,aes(x=x,y=y,color=color,label=label))
  p <- p + scale_color_manual(values=c("grey","orange"))
  p <- p + geom_point()
  p <- p + ggrepel::geom_text_repel(max.overlaps=Inf)
  p <- p + geom_hline(yintercept=1,color="blue",
                      linetype="dashed")
  p <- p + theme_linedraw()
  p <- p + xlab("Estimated distance")
  p <- p + ylab("-log10 FDR")
  p <- p + guides(color="none",alpha="none")
  print(p)

  out <- list()
  out$plot <- p

  return(out)
}

#' plotBetas
#'
#' ADD DESCRIPTION HERE
#'
#' @param scd.object A list obtained from applying the main function \code{\link{scDist}}.
#' @param cluster The cluster to make the plot for
#'
#' @export

plotBetas <- function(
    scd.object,
    cluster
) {
  ix <- which(names(scd.object$vals) == cluster)
  df <- data.frame(beta=scd.object$vals[[ix]]$beta,
                   p=scd.object$vals[[ix]]$p.sum)
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
  #df$beta <- ifelse(df$beta > 0, df$beta+0.35, df$beta-0.35)
  #p <- p + ggplot2::geom_text(data=df, ggplot2::aes(label=signif,x=dim,y=beta))
  p <- p + ggplot2::labs(x="PC Dimension", y=expression(U*beta[k]))
  p <- p + ggplot2::ggtitle(label=cluster)
  p <- p + ggplot2::theme_linedraw()
  p
  return(p)
}

# Create a sample dataset
set.seed(123)
data <- data.frame(
  value = rnorm(20, mean = 10, sd = 2)
)

#' @export
#'
#' @title Plot the direction of the perturbation
#'
#' @description Plot the distance estimates and corresponding standard
#' errors
#'
#' @param scd.object A list obtained from applying the main function \code{\link{scDist}}.
#' @param cluster The cluster to make the plot for
#'
#' @return A `ggplot2` object containing the plot
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>

distGenes <- function(
    scd.object,
    cluster
) {
  # Create the stripchart
  G <- nrow(scd.object$vals[[cluster]]$loadings)
  up5 <- order(scd.object$vals[[cluster]]$beta.hat)[1:5]
  down5 <- order(scd.object$vals[[cluster]]$beta.hat,decreasing=TRUE)[1:5]
  color <- rep("normal", G)
  color[up5] <- "up"; color[down5] <- "down"

  label <- ifelse(1:G %in% c(up5,down5), scd.object$gene.names, "")

  data <- data.frame(value=scd.object$vals[[cluster]]$beta.hat,
                     color=color,
                     label=label)
  p <- ggplot(data, aes(x = "", y = value, color=color,
                   label=label)) +
    geom_jitter(width = 0.01, alpha = 0.5) +  # Add jittered points
    scale_color_manual(values=c("red", "grey90", "blue")) +
    ggrepel::geom_text_repel(max.overlaps=Inf) +
    labs(x = NULL, y = "Condition difference") +
    coord_flip()+
    guides(color="none") +
    theme_bw()
  return(p)
}
