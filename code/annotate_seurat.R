library(tidyverse)
library(Seurat)

sc <- readRDS("../data/seurat.processed.rds")
#sc2 <- readRDS("../data/seurat.processed2.rds")

# rename samples
sc$sample <- str_replace(sc$sample, "day", "")

# filter out clusters w <= 5k cells
small_clusts <- as_tibble(sc$seurat_clusters, rownames="cell") %>% 
  group_by(value) %>%
  count %>%
  filter(n<=5000) %>%
  .$value
keep_cells <- as_tibble(sc$seurat_clusters, rownames="cell") %>%
  filter(!value %in% small_clusts) %>%
  .$cell

sc <- sc[,keep_cells]
sc$seurat_clusters <- sc$seurat_clusters %>%
  as.integer %>%
  factor(levels=seq(1,7))

# annotate cell types
clust2type <- function(c) {
  if (c == 1) {
    type = "CF"
  } else if (c == 2) {
    type = "MES"
  } else if (c == 3) {
    type = "IPSC"
  } else if (c == 4) {
    type = "PROG"
  } else if (c == 5) {
    type = "CMES"
  } else if (c == 6) {
    type = "CM"
  } else if (c == 7) {
    type = "UNK"
  }
}

cell.types <- c("IPSC", "MES", "CMES", "PROG", "CM", "CF", "UNK")
sc$type <- factor(sapply(sc$seurat_clusters, clust2type), levels=cell.types)
saveRDS(sc, "../data/seurat.annotated.rds")

# subset to the two lineages present, rescale and reduce dimensionality
cm_cells <- as_tibble(sc$type, rownames="cell") %>%
  filter(!value %in% c("UNK", "CF")) %>%
  .$cell
sc_cm <- sc[,cm_cells]
sc_cm <- sc_cm %>%
  ScaleData %>%
  RunPCA(npcs=50)
cell.types <- setdiff(levels(sc_cm$type), c("UNK", "CF"))
sc_cm$type <- factor(sc_cm$type, cell.types)
saveRDS(sc_cm, "../data/seurat.cm.rds")

cf_cells <- as_tibble(sc$type, rownames="cell") %>%
  filter(!value %in% c("UNK", "CM")) %>%
  .$cell
sc_cf <- sc[,cf_cells]
sc_cf <- sc_cf %>%
  ScaleData %>%
  RunPCA(npcs=50)
cell.types <- setdiff(levels(sc_cf$type), c("UNK", "CM"))
sc_cf$type <- factor(sc_cf$type, cell.types)
saveRDS(sc_cf, "../data/seurat.cf.rds")

# include both lineages, just excluding unknown (likely non-cardiac cells)
known_cells <- as_tibble(sc$type, rownames="cell") %>%
  filter(!value %in% c("UNK")) %>%
  .$cell
sc_cardiac <- sc[,known_cells]
sc_cardiac <- sc_cardiac %>%
  ScaleData %>%
  RunPCA(npcs=50)
cell.types <- setdiff(levels(sc_cardiac$type), "UNK")
sc_cardiac$type <- factor(sc_cardiac$type, cell.types)
saveRDS(sc_cardiac, "../data/seurat.cardiac.rds")
