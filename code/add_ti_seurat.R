library(tidyverse)
library(Seurat)
library(ggpubr)
library(scales)

pseudotime <- read_tsv("../data/pseudotimes.hacky.tsv")
sc_cm <- readRDS("../data/seurat.cm.rds")
cells <- tibble("cell"=colnames(sc_cm)) %>%
  left_join(pseudotime, by="cell")
sc_cm$pseudotime <- cells$dpt_pseudotime
saveRDS(sc_cm, "../data/seurat.cm.rds")

# sc_cf <- readRDS("../data/sc_cf.rds")
cells <- tibble("cell"=colnames(sc_cf)) %>%
  left_join(pseudotime, by="cell")
sc_cf$pseudotime <- cells$dpt_pseudotime
saveRDS(sc_cf, "../data/seurat.cf.rds")

sc <- readRDS("../data/seurat.annotated.rds")
proportions <- read_tsv("../results/cibersort/pseudobulk_indday.full.true.tsv")
cells <- tibble("cell"=colnames(sc)) %>%
  left_join(pseudotime, by="cell") %>%
  mutate(type=sc$type) %>%
  mutate(diffday=factor(str_replace(sc$diffday, "day", ""),
                    levels=c("0", "1", "3", "5", "7", "11", "15"))) %>%
  mutate(day=as.numeric(str_replace(sc$diffday, "day", ""))) %>%
  mutate(ind=sc$individual) %>%
  mutate(sample=paste(ind, day, sep="_")) %>%
  left_join(proportions, by="sample")
# now we set both CM and CF pseudotime
ct.pseudotime <- function(dpt_pseudotime, type, lineage, ...) {
  if (lineage == "CM") {
    nulltypes <- c("CF", "UNK")
  } else {
    nulltypes <- c("CM", "UNK")
  }
  if (type %in% nulltypes) {
    t <- NA
  } else {
    t <- dpt_pseudotime
  }
  t
}
cells <- cells %>%
  mutate(cm_time=pmap_dbl(.l=., .f=ct.pseudotime, lineage="CM")) %>%
  mutate(cf_time=pmap_dbl(.l=., .f=ct.pseudotime, lineage="CF"))
my_cols_type <- c(hue_pal()(7)[c(5,4,3,2,1,6)], "#949494")
my_cols_day <- c(hue_pal()(7))
ggplot(filter(cells, type!="UNK"), aes(x=type, y=dpt_pseudotime, fill=type)) +
  geom_violin() +
  scale_fill_manual(values=my_cols_type[-c(7)]) +
  theme_classic() +
  xlab("Cell Type") +
  ylab("Pseudotime")
ggplot(filter(cells, type!="UNK"), aes(x=day, y=dpt_pseudotime, fill=diffday, group=diffday)) +
  geom_violin() +
  scale_fill_manual(values=my_cols_day) +
  theme_classic() +
  xlab("Differentiation Day") +
  ylab("Pseudotime") 
cells2 <- cells %>%
  gather(cm_time, cf_time, key="lineage", value="pseudotime") %>%
  filter(ind %in% inds[1:6])
ggplot(cells2, aes(x=day, y=pseudotime, fill=diffday, group=diffday)) +
  geom_violin() +
  scale_fill_manual(values=my_cols_day) +
  theme_classic() +
  facet_grid(cols=vars(lineage))
ggplot(cells2, aes(x=CM, y=pseudotime, group=factor(CM))) +
  geom_violin() +
  scale_fill_manual(values=my_cols_day) +
  theme_classic() +
  facet_grid(cols=vars(lineage))

ggplot(filter(cells, type!="UNK"), aes(x=day, y=dpt_pseudotime)) +
  geom_point() +
  theme_classic() +
  xlab("Differentiation Day") +
  ylab("Pseudotime") +
  stat_cor(aes(x=day, y=dpt_pseudotime, label = ..r.label..))


inds <- unique(cell$ind)
proportions.jr <- proportions %>%
  filter(ind %in% inds[1:6]) %>%
  group_by(sample)
ggplot(proportions.jr, aes(x=day, y=dpt_pseudotime, fill=day)) +
  geom_point() +
  # facet_grid(rows=vars(ind)) +
  stat_cor(aes(x=day, y=dpt_pseudotime, label = ..r.label..))
ggplot(proportions.jr, aes(x=CM, y=dpt_pseudotime, fill=day)) +
  geom_point() +
  # facet_grid(rows=vars(ind)) +
  stat_cor(aes(x=CM, y=dpt_pseudotime, label = ..r.label..))
ggplot(days, aes(x=value, y=dpt_pseudotime, fill=value)) +
  geom_violin()

sc$pseudotime <- cells$dpt_pseudotime
saveRDS(sc, "../data/seurat.annotated.rds")

# sc_cardiac <- readRDS("data/seurat.cardiac.rds")
cells <- tibble("cell"=colnames(sc_cardiac)) %>%
  left_join(pseudotime, by="cell")
sc_cardiac$pseudotime <- cells$dpt_pseudotime
saveRDS(sc_cardiac, "../data/seurat.cardiac.rds")

