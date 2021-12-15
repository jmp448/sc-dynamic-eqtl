suppressPackageStartupMessages(
  c(library("tidyverse"),
    library("vroom"),
    library("patchwork"))
)

cm_results_permute <- vroom("../results/eqtl_dynamic/linear_dQTL/permuted_cells/pseudobulk-cm/bin16/50k-5clpcs-0pcs-notypes.mtc.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)
cf_results_permute <- vroom("../results/eqtl_dynamic/linear_dQTL/permuted_cells/pseudobulk-cf/bin16/50k-5clpcs-0pcs-notypes.mtc.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)

cm_pcomp <- ggplot(sample_n(cm_results_permute, 100000), aes(x=-log10(p.unadj), y=-log10(emp.pvals))) +
  geom_point() +
  theme_classic(base_size=14) +
  xlab("-log10(p.nominal)") + ylab("-log10(p.empirical)") +
  ggtitle("Cardiomyocyte lineage")
cm_hist1 <- ggplot(cm_results_permute, aes(x=emp.pvals)) +
  geom_histogram(bins=100) +
  theme_classic(base_size=14) +
  xlab("p.empirical")
cm_hist2 <- ggplot(cm_results_permute, aes(x=p.unadj)) +
  geom_histogram(bins=100) +
  theme_classic(base_size=14) +
  xlab("p.nominal")

cf_pcomp <- ggplot(sample_n(cf_results_permute, 100000), aes(x=-log10(p.unadj), y=-log10(emp.pvals))) +
  geom_point() +
  theme_classic(base_size=14) +
  xlab("-log10(p.nominal)") + ylab("-log10(p.empirical)") +
  ggtitle("Cardiac fibroblast lineage")
cf_hist1 <- ggplot(cf_results_permute, aes(x=emp.pvals)) +
  geom_histogram(bins=100) +
  theme_classic(base_size=14) +
  xlab("p.empirical")
cf_hist2 <- ggplot(cf_results_permute, aes(x=p.unadj)) +
  geom_histogram(bins=100) +
  theme_classic(base_size=14) +
  xlab("p.nominal")
