library(splatter)
library(scater)
library(limma)
library(vroom)
library(VariantAnnotation)
library(qvalue)
library(Seurat)
library(tidyverse)
library(Matrix)
library(Matrix.utils)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")

counts2cpm <- function(c, lognorm=F) {
  g <- c$gene
  c <- select(c, !gene)
  
  cpm <- DGEList(counts=c) %>% calcNormFactors(method="TMMwsp") %>% cpm(log=lognorm) %>%
    as_tibble %>% mutate(gene=g, .before=1) %>% arrange(gene)
  cpm
}

# load some real data
exp <- vroom("../data/pseudobulk-cm/bin15/logcpm.tsv")
all_tests <- read_tsv("../data/pseudobulk-cm/bin15/filtered_tests.50k.tsv")

# take the first 850 of 8500 genes
tests.sim <- all_tests %>% filter(gene %in% .$gene[1:850]) 
genes.sim <- all_tests %>% .$gene %>% unique %>% .[1:850]

all_tests <- all_tests %>% filter(gene %in% genes.sim)
genotypes <- vroom("../data/genotypes/standardized_genotypes.tsv") %>%
  filter(snp %in% all_tests$snp)
snps <- unique(all_tests$snp)

### SIMULATION 1 - PSEUDOBULK DATA (NO MAIN EFFECTS)
# get the actual variance in their log normalized counts
sim.vars <- exp %>% filter(gene %in% genes.sim) %>% dplyr::select(!gene) %>%
  apply(1, var)
sim.means <- exp %>% filter(gene %in% genes.sim) %>% dplyr::select(!gene) %>%
  apply(1, mean)

# for each gene, simulate logcpm
sim <- tibble("mean"=sim.means, "var"=sim.vars) %>% mutate(sd=sqrt(var))
sim.logcpm <- mapply(rnorm, mean=sim$mean, sd=sim$sd, ncol(exp)-1) %>%
  `rownames<-`(colnames(exp)[-c(1)]) %>%
  `colnames<-`(genes.sim)

# run PCA for the time variable
pcs <- sim.logcpm %>% t %>% as_tibble(rownames="gene") %>% regular.pca
days <- pcs$u %>% dplyr::select(c("sample", "PC1")) %>% dplyr::rename(t=PC1)


### SIMULATION 2 - SINGLE CELL DATA (NO MAIN EFFECTS)
params <- newSplatPopParams(eqtl.n=0, similarity.scale=50)
geno19 <- mockVCF(n.samples=19)
gff <- mockGFF(n.genes=850)

sim.means <- splatPopSimulateMeans(vcf=geno19, gff=gff, 
                                   params=params)
sim.sc <- splatPopSimulateSC(params=params,
                             key=sim.means$key,
                             sim.means=sim.means$means,
                             batchCells=10000)

# replace sample with one of our individuals
inds <- colnames(exp)[-c(1)] %>% str_extract("[^_]+") %>% unique
stopifnot(length(inds) == length(unique(colData(sim.sc)$Sample)))
colData(sim.sc)$Sample <- tibble("old"=colData(sim.sc)$Sample) %>% left_join(tibble("old"=unique(colData(sim.sc)$Sample), "new"=inds), by="old") %>% .$new
sim.sc <- logNormCounts(sim.sc)
sim.sc <- runPCA(sim.sc, ncomponents = 10)
plotPCA(sim.sc, colour_by = "Sample")

# aggregate pseudobulk into 15 PC1 bins
counts <- assay(sim.sc, "counts")
colData(sim.sc)$CMBin <- as_tibble(reducedDim(sim.sc, "PCA")[,1]) %>%
  `colnames<-`(c("PC1")) %>% 
  rowid_to_column(var="orig.order") %>%
  arrange(PC1) %>%
  rowid_to_column("time_order") %>% 
  mutate(bin=floor(15*time_order/nrow(.)), .keep="unused") %>%
  mutate(bin=sapply(bin, function(x){min(x, 14)})) %>%
  arrange(orig.order) %>% .$bin
colData(sim.sc)$BinInd <- paste(colData(sim.sc)$Sample, colData(sim.sc)$CMBin, sep="_")
plotPCA(sim.sc, colour_by = "BinInd")

sim.logcpm <- t(counts) %>% aggregate.Matrix(colData(sim.sc)$BinInd, fun="sum") %>%
  `colnames<-`(genes.sim) 
days <- tibble("sample"=colData(sim.sc)$BinInd, "PC1"=reducedDim(sim.sc, "PCA")[,1]) %>%
  group_by(sample) %>% summarise(t=median(PC1))


### SIMULATION 3 - SINGLE CELL DATA (1% EQTL EFFECTS, PSEUDOTIME)
sim.snps <- vroom("../data/genotypes/snp_locs.bed", col_names=c("chr", "start", "stop", "ID")) %>% 
  mutate(snp=paste(chr, stop, sep="_")) %>% 
  filter(snp %in% snps) %>%
  select(ID) %>%
  write_tsv("../data/genotypes/simulation_snps.txt", col_names=F)
# get new vcf file with ../data/genotypes/simulation_filter.sh
genes <- vroom("../data/gencode.hg38.filtered.gff") %>%
  dplyr::mutate(hgnc=gsub(".*gene_name=([^;]+)[;].*", "\\1", attribute)) %>% 
  filter(hgnc %in% genes.sim) %>%
  dplyr::select(!hgnc) %>%
  write_tsv("../data/gencode.hg38.simulation.gff")

params <- newSplatPopParams(eqtl.n=0.05, 
                            similarity.scale=5,
                            de.facLoc=0.5,
                            de.facScale=0.5,
                            de.prob=0.5)
geno19 <- readVcf("../data/genotypes/simulated.recode.vcf")
gff <- genes

sim.means <- splatPopSimulateMeans(vcf=geno19, gff=gff, 
                                   params=params)
sim.sc <- splatPopSimulateSC(params=params,
                             key=sim.means$key,
                             sim.means=sim.means$means,
                             batchCells=1000,
                             method="paths")

# replace sample with one of our individuals
inds <- colnames(exp)[-c(1)] %>% str_extract("[^_]+") %>% unique
stopifnot(length(inds) == length(unique(colData(sim.sc)$Sample)))
colData(sim.sc)$Sample <- tibble("old"=colData(sim.sc)$Sample) %>% left_join(tibble("old"=unique(colData(sim.sc)$Sample), "new"=inds), by="old") %>% .$new
sim.sc <- logNormCounts(sim.sc)
sim.sc <- runPCA(sim.sc, ncomponents = 10)
plotPCA(sim.sc, colour_by = "Sample")

# aggregate pseudobulk into 15 PC1 bins
counts <- assay(sim.sc, "counts")
colData(sim.sc)$CMBin <- as_tibble(reducedDim(sim.sc, "PCA")[,1]) %>%
  `colnames<-`(c("PC1")) %>% 
  rowid_to_column(var="orig.order") %>%
  arrange(PC1) %>%
  rowid_to_column("time_order") %>% 
  mutate(bin=floor(15*time_order/nrow(.)), .keep="unused") %>%
  mutate(bin=sapply(bin, function(x){min(x, 14)})) %>%
  arrange(orig.order) %>% .$bin
colData(sim.sc)$BinInd <- paste(colData(sim.sc)$Sample, colData(sim.sc)$CMBin, sep="_")
plotPCA(sim.sc, colour_by = "CMBin")

sim.logcpm <- t(counts) %>% aggregate.Matrix(colData(sim.sc)$BinInd, fun="sum") %>%
  `colnames<-`(genes.sim) 
days <- tibble("sample"=colData(sim.sc)$BinInd, "PC1"=reducedDim(sim.sc, "PCA")[,1]) %>%
  group_by(sample) %>% summarise(t=median(PC1))


### SIMULATION ANALYSIS - needs `sim.logcpm` and `days`
genotypes <- genotypes %>% column_to_rownames("snp") %>% t %>% 
  as_tibble(rownames="ind") %>% arrange(ind)
inds <- tibble("sample"=rownames(sim.logcpm)) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
  arrange(ind)
ind.indices <- as.numeric(factor(inds$ind))
geno.mat <- genotypes %>% column_to_rownames("ind") %>% as.matrix
genotypes <- geno.mat[ind.indices,] %>% `rownames<-`(inds$sample) %>%
  as_tibble(rownames="sample")

snps <- colnames(genotypes)[-c(1)]

# call eQTLs
eqtl.call <- function(i) {
  snp.genes = all_tests$gene[all_tests$snp==snps[i]]
  exp.snp = t(sim.logcpm[,all_of(snp.genes)])
  design = left_join(genotypes[,c(1,i+1)], days, by="sample") %>% 
    dplyr::rename(g=starts_with("chr"))
  design.snp = model.matrix(formula("~g*t"), design)
  model=lmFit(exp.snp, design.snp)
  beta = model$coefficients[,"g:t"]
  sd = model$stdev.unscaled[,"g:t"]*model$sigma
  tibble("gene"=snp.genes, "snp"=snps[i], "beta"=beta, "sd"=sd)
}

qtl_results <- bind_rows(lapply(1:length(snps), eqtl.call))

# bonferroni correction
qtl_results <- qtl_results %>%
  mutate(t.unadj = beta / sd) %>%
  mutate(p.unadj = 2*pt(-abs(t.unadj), df=276))
gene_ct <- table(qtl_results$gene)
qtl_results$bonf.p <- sapply(gene_ct[qtl_results$gene] * qtl_results$p.unadj, function(x){min(1,x)})
qtl_results <- qtl_results %>% write_tsv("../results/simulation/sim3_results.tsv")

tophits <- qtl_results %>%
  arrange(bonf.p) %>%
  filter(!duplicated(gene)) %>%
  mutate(qval=qvalue(bonf.p)$q) %>%
  write_tsv("../results/simulation/sim3_tophits.tsv")

ggplot(qtl_results, aes(x=t.unadj)) + geom_density()

