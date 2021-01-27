library(Seurat)
library(tidyverse)
library(Matrix)

# pointer to the original seurat file, no filtering
full_file_loc <- "/project2/gilad/reem/singlecellCM/Sobjs/xfinalmerged/Sobj_allrounds_finalmerged_withNA_demuxlabelsastry2.rds"

# set thresholds
prb_dbl_cutoff <- 0.3  # max demuxlet-assigned doublet posterior prob
pct_mito_cutoff <- 25  # max percent mito
n_gene_cutoff <- 300  # min genes per cell
n_cell_cutoff <- 10  # min cells per gene
n_sd_reads_cutoff <- 4  # max standard deviations above the median (counts)
n_sd_genes_cutoff <- 4  # max standard deviations above the median (genes)

sc <- readRDS(full_file_loc)

quality_genes <- sc@assays$RNA@counts %>%
  (function(x)(x>0)) %>%
  rowSums %>%
  as_tibble(rownames="gene") %>%
  filter(value>=n_cell_cutoff) %>%
  .$gene

metadata <- as_tibble(sc@meta.data, rownames="cell")

quality_cells <- metadata %>%
  filter(PRB.DBL <= prb_dbl_cutoff) %>%
  filter(individual != "doublet_ambiguous") %>%
  filter(percent.mito <= pct_mito_cutoff) %>%
  filter(nFeature_RNA >= n_gene_cutoff) 

counts_cut <- median(quality_cells$nCount_RNA) + n_sd_reads_cutoff*sd(quality_cells$nCount_RNA)
feats_cut <- median(quality_cells$nFeature_RNA) + n_sd_reads_cutoff*sd(quality_cells$nFeature_RNA)

quality_cells <- quality_cells %>%
  filter(nCount_RNA <= counts_cut) %>%
  filter(nFeature_RNA <= feats_cut) %>%
  .$cell

sc <- sc[quality_genes, quality_cells]

# clean up metadata
sc$diffday <- str_replace(sc$diffday_try2, "Day ", "day") %>%
  as_factor() %>%
  fct_relevel(paste0("day", c(0, 1, 3, 5, 7, 11, 15)))
sc$individual <- str_replace(sc$individual_try2, "NA", "") %>%
  as_factor()
sc$sample <- sc$sample_try2 %>% 
  str_replace("NA", "") %>%
  str_replace(".Day ", "_") %>%
  as_factor()

sample_levels <- c()
for (i in levels(sc$individual)) {
  for (d in c(0, 1, 3, 5, 7, 11, 15)) {
    sample_levels <- c(sample_levels, paste(i, d, sep="_"))
  }
}

sc$sample <- sc$sample %>%
  fct_relevel(sample_levels)

saveRDS(sc, "../data/seurat.filtered.rds")
