library(tidyverse)

cis.dist <- "50k"
all_tests <- read_tsv(paste0("../data/all_tests.", cis.dist, ".tsv"), col_names=c("gene", "snp")) %>%
  filter(!duplicated(.))

get_groups <- function(dataset, agg) {
  if (dataset %in% c("bulk", "bulk7")) {
    paste0("day", seq(0, 15))
  } else if ((dataset == "pseudobulk") & (agg == "day")) {
    paste0("day", c(0, 1, 3, 5, 7, 11, 15)) 
  } else if ((dataset == "pseudobulk") & (agg == "type")) {
    c("IPSC", "MES", "CMES", "PROG", "CF", "CM", "UNK")
  }
}

get_filters <- function(dataset, agg, ...) {
  cpm_tpm <- if_else(dataset %in% c("bulk", "bulk7"), "tpm", "cpm")
  for (g in get_groups(dataset, agg)) {
    # at least 5 samples with cpm >= 0.1
    cpm <- vroom(paste0("../data/static/", dataset, "/", agg, "/", g, "/", cpm_tpm, ".tsv")) %>%
      column_to_rownames("gene")
    filter1g <- rownames(cpm)[rowSums(cpm>=0.1)>=5]
    
    # at least 5 samples w (corrected) counts >= 5
    counts <- read_tsv(paste0("../data/", dataset, "/", agg, "/counts.tsv")) %>%
      column_to_rownames("gene")
    filter2g <- rownames(counts)[rowSums(counts>=5)>=5]
    
    # keep track of the intersection across all groups
    if (g == get_groups(dataset, agg)[1]) {
      filter1 <- filter1g
      filter2 <- filter2g
    } else {
      filter1 <- intersect(filter1, filter1g)
      filter2 <- intersect(filter2, filter2g)
    }
    
    filtered_tests_g <- all_tests %>%
      filter(gene %in% intersect(filter1g, filter2g)) %>%
      write_tsv(paste0("../data/static/", dataset, "/", agg, "/", g, "/filtered_tests.50k.tsv"))
  }
  
  filtered_tests <- all_tests %>%
    filter(gene %in% intersect(filter1, filter2)) %>%
    write_tsv(paste0("../data/static/", dataset, "/", agg, "/overlap_filtered_tests.50k.tsv"))
}

aggregation_schemes <- tibble("dataset"=c("bulk", "pseudobulk", "pseudobulk"),
                              "agg"=c("day", "day", "type"))

pmap(aggregation_schemes, get_filters)