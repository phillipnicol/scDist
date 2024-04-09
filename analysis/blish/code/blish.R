
Sco <- readRDS("../../data/blish_covid.seu.rds")

library(scDist)

out <- scDist(normalized_counts = Sco@assays$SCT@scale.data,
              meta.data=Sco@meta.data,
              fixed.effects = "Status",
              random.effects="Donor",
              clusters="cell.type.coarse")


saveRDS(out$results, "../data/scDist_results.RDS")
