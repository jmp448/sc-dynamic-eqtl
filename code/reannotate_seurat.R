library(Seurat)
library(tidyverse)

sc <- readRDS("../data/seurat.normalized.rds")

# add cluster assignments
leiden <- read_tsv("../data/marcc_transfer3/leiden.tsv")
leiden_md <- leiden$leiden
names(leiden_md) <- leiden$X1
sc <- AddMetaData(sc, leiden_md, col.name="leiden")

clust2type <- function(c) {
  if (c == 0) {
    type = "MES"
  } else if (c == 1) {
    type = "CF"
  } else if (c == 2) {
    type = "IPSC"
  } else if (c == 3) {
    type = "PROG"
  } else if (c == 4) {
    type = "CM"
  } else if (c == 5) {
    type = "CMES"
  } else {
    type = "UNK"
  }
}

cell.types <- c("IPSC", "MES", "CMES", "PROG", "CM", "CF", "UNK")
sc$type <- factor(sapply(sc$leiden, clust2type), levels=cell.types)

umap <- read_tsv("../data/marcc_transfer3/umap.tsv") %>%
  rename(UMAP_1=`0`, UMAP_2=`1`) %>% 
  column_to_rownames("X1") %>%
  as.matrix()
sc[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = "SCT")

saveRDS(sc, "../data/seurat.clustered.rds")

# subset to cardiac cells 
sc <- sc[,sc$leiden %in% seq(0,5)]

pseudotime <- read_tsv("../data/marcc_transfer3/pseudotime.tsv")
pseudo_md <- pseudotime$dpt_pseudotime
names(pseudo_md) <- pseudotime$X1
sc <- AddMetaData(sc, pseudo_md, col.name="pseudotime")

fa <- read_tsv("../data/marcc_transfer3/graph_embedding.tsv") %>%
  rename(FA_1=`0`, FA_2=`1`) %>% 
  column_to_rownames("X1") %>%
  as.matrix()
sc[["fa"]] <- CreateDimReducObject(embeddings = fa, key = "FA_", assay = "SCT")
saveRDS(sc, "../data/seurat.annotated.rds")

# subset to the two lineages present
cm_cells <- as_tibble(sc$type, rownames="cell") %>%
  filter(!value %in% c("CF")) %>%
  .$cell
sc_cm <- sc[,cm_cells]
cell.types <- setdiff(levels(sc_cm$type), c("CF"))
sc_cm$type <- factor(sc_cm$type, cell.types)
saveRDS(sc_cm, "../data/seurat.cm.rds")

cf_cells <- as_tibble(sc$type, rownames="cell") %>%
  filter(!value %in% c("CM")) %>%
  .$cell
sc_cf <- sc[,cf_cells]
cell.types <- setdiff(levels(sc_cf$type), c("CM"))
sc_cf$type <- factor(sc_cf$type, cell.types)
saveRDS(sc_cf, "../data/seurat.cf.rds")
