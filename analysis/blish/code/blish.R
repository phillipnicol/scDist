
Sco <- readRDS("../../data/blish_covid.seu.rds")

library(scDist)

out <- scDist(normalized_counts = Sco@assays$SCT@scale.data,
              meta.data=Sco@meta.data,
              fixed.effects = "Status",
              random.effects="Donor",
              clusters="cell.type.coarse")


saveRDS(out$results, "../data/scDist_results.RDS")


out2 <- condition_v_celltype_distance(normalized_counts = Sco@assays$SCT@scale.data,
                                      meta.data=Sco@meta.data,
                                      fixed.effects = "Status",
                                      random.effects = "Donor",
                                      clusters = "cell.type.coarse")
saveRDS(out2$D, "../data/celltype_v_condition.RDS")


p.dist.plot <- scDist::DistPlot(out, return.plot = TRUE)
p.fdr.plot <- scDist::FDRDistPlot(out)$plot
p.new.plot <- scDist::distGenes(out,"CD14 Monocyte") + ggtitle("CD14 Monocytes")

p <- ggarrange(ggarrange(p.dist.plot, p.fdr.plot, nrow=1,ncol=2),
          p.new.plot, nrow=2, ncol=1, heights=c(1.5,1))

ggsave(p, filename="../plots/blish_example.png",
       width=8.64, height=7.1)

