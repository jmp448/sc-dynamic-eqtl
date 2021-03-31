library(tidyverse)
library(vroom)
library(patchwork)
library(snpStats)
source("viz_helpers.R")

gwas_metadata <- read_tsv("../data/gwas/gwas_metadata.txt")
gwas_cardiometabolic <- filter(gwas_metadata, Category == "Cardiometabolic")
gwas_anthropometric <- filter(gwas_metadata, Category == "Anthropometric")

hypertension <- "../data/gwas/imputed_gwas_hg38_1.1/imputed_UKB_20002_1065_self_reported_hypertension.txt.gz"
hypertension_gwas <- vroom(hypertension, col_select=c(chromosome, position, pvalue))
hypertension_hits <- hypertension_gwas %>%
  filter(pvalue<=0.1/nrow(.))

cholesterol <- "../data/gwas/imputed_gwas_hg38_1.1/imputed_UKB_20002_1473_self_reported_high_cholesterol.txt.gz"
cholesterol_gwas <- vroom(cholesterol, col_select=c(chromosome, position, pvalue)) 
cholesterol_hits <- cholesterol_gwas %>%
  filter(pvalue<=0.1/nrow(.))

bmi <- "../data/gwas/imputed_gwas_hg38_1.1/imputed_UKB_21001_Body_mass_index_BMI.txt.gz"
bmi_gwas <- vroom(bmi, col_select=c(chromosome, position, pvalue)) 
bmi_hits <- bmi_gwas %>%
  filter(pvalue<=0.1/nrow(.))

bfp <- "../data/gwas/imputed_gwas_hg38_1.1/imputed_UKB_23099_Body_fat_percentage.txt.gz"
bfp_gwas <- vroom(bmi, col_select=c(chromosome, position, pvalue)) %>%
  filter(pvalue<=0.1/nrow(.))

ha <- "../data/gwas/imputed_gwas_hg38_1.1/imputed_UKB_6150_1_Vascular_or_heart_problems_diagnosed_by_doctor_Heart_attack.txt.gz"
ha_gwas <- vroom(ha, col_select=c(chromosome, position, pvalue)) %>%
  filter(pvalue<=0.1/nrow(.))

t2d <- "../data/gwas/imputed_gwas_hg38_1.1/imputed_UKB_20002_1223_self_reported_type_2_diabetes.txt.gz"
t2d_gwas <- vroom(t2d, col_select=c(chromosome, position, pvalue)) 
t2d_hits <- t2d_gwas %>%
  filter(pvalue<=0.1/nrow(.))

qtl_day <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk/day/50k-5clpcs-0pcs-notypes.mtc.tsv")
qtl_cm <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin16/50k-5clpcs-0pcs-notypes.mtc.tsv")
sig_qtl_cm <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin16/50k-5clpcs-0pcs-notypes.tophits.tsv") %>% 
  filter(qval.unadj<=0.05)

qtl_cf <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-5clpcs-0pcs-notypes.mtc.tsv")
sig_qtl_cf <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-5clpcs-0pcs-notypes.tophits.tsv") %>% 
  filter(qval.unadj<=0.05)

ieqtl_cm <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-0pcs-regtypes.mtc.tsv")
sig_ieqtl_cm <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-0pcs-regtypes.tophits.tsv") %>% 
  filter(qval.unadj<=0.05)

ieqtl_cf <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-regtypes.mtc.tsv")
sig_ieqtl_cf <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-regtypes.tophits.tsv") %>% 
  filter(qval.unadj<=0.05)

ieqtl_cmes <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CMES/50k-5clpcs-0pcs-notypes.mtc.tsv")
sig_ieqtl_cmes <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CMES/50k-5clpcs-0pcs-notypes.tophits.tsv") %>% 
  filter(qval.unadj<=0.05)

nonlinear_bulk <- vroom("../results/eqtl_dynamic/nonlinear_dQTL/bulk/day/t_test/50k-5clpcs-0pcs-notypes.mtc.tsv")
sig_nl_qtl_bulk <- vroom("../results/eqtl_dynamic/nonlinear_dQTL/bulk/day/t_test/50k-5clpcs-0pcs-notypes.tophits.tsv") %>% 
  filter(qval.unadj<=0.05)

nl_qtl_cm <- vroom("../results/eqtl_dynamic/nonlinear_dQTL/pseudobulk-cm/bin16/t_test/50k-5clpcs-0pcs-notypes.mtc.tsv")
sig_nl_qtl_cm <- vroom("../results/eqtl_dynamic/nonlinear_dQTL/pseudobulk-cm/bin16/t_test/50k-5clpcs-0pcs-notypes.tophits.tsv") %>% 
  filter(qval.unadj<=0.05)

nl_qtl_cf <- vroom("../results/eqtl_dynamic/nonlinear_dQTL/pseudobulk-cf/bin16/t_test/50k-5clpcs-0pcs-notypes.mtc.tsv")
sig_nl_qtl_cf <- vroom("../results/eqtl_dynamic/nonlinear_dQTL/pseudobulk-cf/bin16/t_test/50k-5clpcs-0pcs-notypes.tophits.tsv") %>% 
  filter(qval.unadj<=0.05)


hypertension_cm_cands <- filter(sig_qtl_cm, snp %in% paste(hypertension_gwas$chromosome, hypertension_gwas$position, sep="_"))$gene
hypertension_cf_cands <- filter(sig_qtl_cf, snp %in% paste(hypertension_gwas$chromosome, hypertension_gwas$position, sep="_"))$gene
hypertension_ie_cm_cands <- filter(sig_ieqtl_cm, snp %in% paste(hypertension_gwas$chromosome, hypertension_gwas$position, sep="_"))$gene
hypertension_ie_cf_cands <- filter(sig_ieqtl_cf, snp %in% paste(hypertension_gwas$chromosome, hypertension_gwas$position, sep="_"))$gene
hypertension_ie_cmes_cands <- filter(sig_ieqtl_cmes, snp %in% paste(hypertension_gwas$chromosome, hypertension_gwas$position, sep="_"))$gene
hypertension_nl_cm_cands <- filter(sig_nl_qtl_cm, snp %in% paste(hypertension_gwas$chromosome, hypertension_gwas$position, sep="_"))$gene
hypertension_nl_cf_cands <- filter(sig_nl_qtl_cf, snp %in% paste(hypertension_gwas$chromosome, hypertension_gwas$position, sep="_"))$gene

bmi_nl_cm_cands <- filter(sig_nl_qtl_cm, snp %in% paste(bmi_gwas$chromosome, bmi_gwas$position, sep="_"))$gene
bmi_nl_cf_cands <- filter(sig_nl_qtl_cf, snp %in% paste(bmi_gwas$chromosome, bmi_gwas$position, sep="_"))$gene
bmi_ie_cm_cands <- filter(sig_ieqtl_cm, snp %in% paste(bmi_hits$chromosome, bmi_hits$position, sep="_"))$gene
bmi_ie_cf_cands <- filter(sig_ieqtl_cf, snp %in% paste(bmi_hits$chromosome, bmi_hits$position, sep="_"))$gene

bfp_nl_cm_cands <- filter(sig_nl_qtl_cm, snp %in% paste(bfp_gwas$chromosome, bfp_gwas$position, sep="_"))$gene
bfp_nl_cf_cands <- filter(sig_nl_qtl_cf, snp %in% paste(bfp_gwas$chromosome, bfp_gwas$position, sep="_"))$gene

chol_nl_cm_cands <- filter(sig_nl_qtl_cm, snp %in% paste(cholesterol_gwas$chromosome, cholesterol_gwas$position, sep="_"))$gene
chol_nl_cf_cands <- filter(sig_nl_qtl_cf, snp %in% paste(cholesterol_gwas$chromosome, cholesterol_gwas$position, sep="_"))$gene
chol_ie_cm_cands <- filter(sig_ieqtl_cm, snp %in% paste(cholesterol_hits$chromosome, cholesterol_hits$position, sep="_"))$gene
chol_ie_cf_cands <- filter(sig_ieqtl_cf, snp %in% paste(cholesterol_hits$chromosome, cholesterol_hits$position, sep="_"))$gene

ha_nl_cm_cands <- filter(sig_nl_qtl_cm, snp %in% paste(ha_gwas$chromosome, ha_gwas$position, sep="_"))$gene
ha_nl_cf_cands <- filter(sig_nl_qtl_cf, snp %in% paste(ha_gwas$chromosome, ha_gwas$position, sep="_"))$gene
ha_ie_cm_cands <- filter(sig_ieqtl_cm, snp %in% paste(ha_gwas$chromosome, ha_gwas$position, sep="_"))$gene
ha_ie_cf_cands <- filter(sig_ieqtl_cf, snp %in% paste(ha_gwas$chromosome, ha_gwas$position, sep="_"))$gene

t2d_nl_cm_cands <- filter(sig_nl_qtl_cm, snp %in% paste(t2d_gwas$chromosome, t2d_gwas$position, sep="_"))$gene
t2d_nl_cf_cands <- filter(sig_nl_qtl_cf, snp %in% paste(t2d_gwas$chromosome, t2d_gwas$position, sep="_"))$gene
t2d_ie_cm_cands <- filter(sig_ieqtl_cm, snp %in% paste(t2d_hits$chromosome, t2d_hits$position, sep="_"))$gene
t2d_ie_cf_cands <- filter(sig_ieqtl_cf, snp %in% paste(t2d_hits$chromosome, t2d_hits$position, sep="_"))$gene


# visualize dynamic eQTL
genodict <- vroom("../data/YRI.hg38.filtered.tsv", col_select=c("#CHROM", "POS", "REF", "ALT", "ID")) %>%
  mutate(snp=paste(`#CHROM`, POS, sep="_"), .keep="unused")
egene <- "ARHGAP42"
evar <- filter(sig_qtl_cf, gene==egene)$snp
evar.ref <- filter(genodict, snp==evar)$REF
evar.alt <- filter(genodict, snp==evar)$ALT
evar.rsid <- filter(genodict, snp==evar)$ID
p1 <- viz_residuals(egene, evar, "pseudobulk-cf", "bin16", nclpcs=1, ref=evar.ref, alt=evar.alt, rsid=evar.rsid)
p2 <- viz_residuals(egene, evar, "pseudobulk-cm", "bin16", nclpcs=1)
print(p1 + p2 + plot_layout(ncol=1))

# look for epigenetics
snp.cand <- filter(sig_qtl_cf, snp==evar) %>%
  select(snp) %>%
  mutate(chr=str_extract(snp, "[^_]+")) %>%
  mutate(stop=as.numeric(str_extract(snp, "[^_]+$"))) %>%
  select(-snp) %>%
  mutate(start=stop-1, .after=1) %>%
  write_tsv("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.bed", col_names = F)
  
declare -a ipsc_epis=( "E018" "E019" "E020" "E021" "E022" )
for epi in ${ipsc_epis[@]}; do
  bedtools intersect -a ../data/epigenomes/iPSC/${epi}_15state.bed \
  -b ../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.bed \
  >> ../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.ipsc.bed
done

ipsc.annotation <- read_tsv("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.ipsc.bed", col_names=F)$X4

declare -a heart_epis=( "E083" "E104" "E095" "E105" "E065" )
for epi in ${heart_epis[@]}; do
bedtools intersect -a ../data/epigenomes/heart/${epi}_15state.bed \
-b ../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.bed \
>> ../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.heart.bed
done

heart.annotation <- read_tsv("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.heart.bed", col_names=F)$X4

declare -a muscle_epis=( "E076" "E078" "E103" "E111" )
for epi in ${muscle_epis[@]}; do
bedtools intersect -a ../data/epigenomes/muscle/${epi}_15state.bed \
-b ../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.bed \
>> ../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.muscle.bed
done

muscle.annotation <- read_tsv("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/gwas_candidates.muscle.bed", col_names=F)$X4


egene <- "APOE"
evar <- filter(sig_nl_qtl_cf, gene==egene)$snp
p1 <- viz_residuals_nonlinear(egene, evar, "pseudobulk-cf", "bin16", nclpcs=5, ref="C", alt="G")
p2 <- viz_residuals_nonlinear(egene, evar, "pseudobulk-cm", "bin16", nclpcs=5)
print(p1 + p2 + plot_layout(ncol=1))

# Manhattan plot
manhattan_qtl <- qtl_cf %>%
  filter(gene == "ARHGAP42") %>%
  mutate(pos=as.numeric(str_extract(snp, "[^_]+$"))) %>%
  mutate(log10p=-log10(p.unadj))
my.snp <- filter(manhattan_qtl, log10p==max(log10p))$pos

manhattan_gwas <- hypertension_gwas %>%
  filter((chromosome == str_extract(manhattan_qtl$snp[1], "[^_]+")) & (position >= min(manhattan_qtl$pos)) & (position <= max(manhattan_qtl$pos))) %>%
  mutate(log10p=-log10(pvalue))
p3 <- ggplot(manhattan_gwas, aes(x=position, y=log10p)) +
  geom_point() + 
  theme_classic()
manhattan_combined <- manhattan_gwas %>%
  select(position, log10p) %>%
  rename(pos=position) %>%
  mutate(assay="Hypertension GWAS") %>%
  bind_rows(mutate(select(manhattan_qtl, c(pos, log10p)), assay="Dynamic eQTL"))

ggplot(manhattan_combined, aes(x=pos, y=log10p, color=assay)) +
  geom_point() +
  facet_grid(rows=vars(assay), scales="free_y") +
  geom_line(aes(x=my.snp)) +
  theme_classic() +
  theme(legend.position="none") +
  ylab("-Log10(P-value)") +
  xlab("Chromosome 11 Position")

# compute LD in the Yoruba population
genotypes <- vroom("../data/genotypes.tsv") %>%
  filter(snp %in% manhattan_qtl$snp) %>% 
  column_to_rownames("snp") %>%
  t
ld <- cor(genotypes)^2
lead.hit <- arrange(manhattan_qtl, desc(log10p))$snp[1]
test <- manhattan_qtl %>% 
  mutate(ld.lead=ld[snp, lead.hit])
ggplot(test, aes(x=pos, y=log10p, color=ld.lead)) +
  geom_point() +
  theme_classic() +
  scale_color_distiller(palette="RdYlBu")

# compute LD on GWAS window
snps <- manhattan_gwas %>%
  select(chromosome, position) %>%
  mutate(chromosome=str_replace(chromosome, "chr", "")) %>%
  rename(`#CHROM`=chromosome) %>%
  rename(POS=position) %>%
  write_tsv(paste0("../data/gwas/ref1000g/snp_subset_", egene, ".tsv"))
# ! vcftools --vcf ../data/gwas/ref1000g/chr11_1000genomes.GRCh38.vcf --positions ../data/gwas/ref1000g/snp_subset_ARHGAP42.tsv --out ARHGAP42_window
# ! vcftools --vcf ../data/gwas/ref1000g/chr11_1000genomes.GRCh38.vcf --positions ../data/gwas/ref1000g/test.tsv --out test

test <- vroom("../data/gwas/ref1000g/chr11_1000genomes.GRCh38.vcf", col_select=c(POS), skip=19)
test2 <- test %>% mutate(`#CHROM`=11, .before=1) %>%
  write_tsv("../data/gwas/ref1000g/test.tsv")
test3 <- vroom("../data/gwas/ref1000g/snp_subset_ARHGAP42.tsv")
test4 <- qtl_cm %>%
  mutate(`#CHROM`=str_extract(snp, "[^chr_]+")) %>%
  mutate(POS=str_extract(snp, "[^_]+$")) %>%
  select(`#CHROM`, POS) %>%
  write_tsv("../data/gwas/ref1000g/test2.tsv")
# run coloc
dqtl.data <- qtl_cm %>% filter(gene == "NPPA")
gwas.data <- hypertension_gwas %>%
    filter((chromosome == "chr1") & (position >= min(manhattan_qtl$pos)) & (position <= max(manhattan_qtl$pos)))

# save gwas for locuszoom 
manhattan_print <- vroom(hypertension) %>%
  filter((chromosome == "chr1") & (position >= min(manhattan_qtl$pos)) & (position <= max(manhattan_qtl$pos))) %>%
  select(c(chromosome, position, effect_allele, non_effect_allele, pvalue, effect_size, standard_error, frequency)) %>%
  mutate(log10p=-log10(pvalue), .keep="unused") %>%
  arrange(position) %>%
  write_tsv("../data/locuszoom/hypertension_nppa.tsv.gz")
  
# ieQTL visualization
# Manhattan plot
egene <- "TBC1D10B"
manhattan_qtl <- ieqtl_cm %>%
  filter(gene == egene) %>%
  mutate(pos=as.numeric(str_extract(snp, "[^_]+$"))) %>%
  mutate(log10p=-log10(p.unadj))
my.snp <- as.numeric(str_extract(filter(sig_ieqtl_cm, gene==egene)$snp, "[^_]+$"))

manhattan_gwas <- bmi_gwas %>%
  filter((chromosome == str_extract(manhattan_qtl$snp[1], "[^_]+")) & (position >= min(manhattan_qtl$pos)) & (position <= max(manhattan_qtl$pos))) %>%
  mutate(log10p=-log10(pvalue))
# ggplot(manhattan_gwas, aes(x=position, y=log10p)) +
#   geom_point() + 
#   theme_classic()
manhattan_combined <- manhattan_gwas %>%
  select(position, log10p) %>%
  rename(pos=position) %>%
  mutate(assay="Hypertension GWAS") %>%
  bind_rows(mutate(select(manhattan_qtl, c(pos, log10p)), assay="Dynamic eQTL"))
ggplot(manhattan_combined, aes(x=pos, y=log10p, color=assay)) +
  geom_point() +
  facet_grid(rows=vars(assay), scales="free_y") +
  geom_line(aes(x=my.snp)) +
  theme_classic() +
  theme(legend.position="none") +
  ylab("-Log10(P-value)") +
  xlab("Chromosome 11 Position")

