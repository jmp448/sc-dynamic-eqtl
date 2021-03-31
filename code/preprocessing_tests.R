library(tidyverse)

cis.dist <- "50k"
all_tests <- read_tsv(paste0("../data/all_tests.", cis.dist, ".tsv"), col_names=c("gene", "snp")) %>%
  filter(!duplicated(.))

get_filters <- function(dataset, agg, ...) {
  # at least 10 samples with cpm >= 0.1
  cpm_tpm <- if_else(dataset %in% c("bulk", "bulk7"), "tpm", "cpm")
  cpm <- read_tsv(paste0("../data/", dataset, "/", agg, "/", cpm_tpm, ".tsv")) %>%
    column_to_rownames("gene")
  filter1 <- rownames(cpm)[rowSums(cpm>=0.1)>=10]
  
  # at least 10 samples w (corrected) counts >= 6
  counts <- read_tsv(paste0("../data/", dataset, "/", agg, "/counts.tsv")) %>%
    column_to_rownames("gene")
  filter2 <- rownames(counts)[rowSums(counts>=6)>=10]
  
  filtered_tests <- all_tests %>%
    filter(gene %in% intersect(filter1, filter2)) %>%
    write_tsv(paste0("../data/", dataset, "/", agg, "/filtered_tests.50k.tsv"))
}

# aggregation_schemes <- tibble("dataset"=c("bulk", "bulk7", "pseudobulk", "pseudobulk-cm", "pseudobulk-cm", "pseudobulk-cf", "pseudobulk-cf"),
#                               "agg"=c("day", "day", "day", "day", "bin15", "day", "bin15"))
aggregation_schemes <- tibble("dataset"=rep(c("pseudobulk-cm", "pseudobulk-cf"), each=4),
                              "agg"=paste0("bin", rep(c(7, 16, 25, 30), times=2)))

pmap(aggregation_schemes, get_filters)