library(tidyverse)
library(vroom)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")

viz_expression <- function(egene, evar, dataset, agg) {
  
  logexp <- if_else(dataset=="bulk", "logtpm", "logcpm")
  
  expr <- vroom(paste0("../data/", dataset, "/", agg, "/", logexp, ".tsv")) %>% 
    filter(gene==!!egene) %>% 
    column_to_rownames("gene") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("sample", "logcpm")) %>%
    mutate(day=str_extract(sample, "[^_]+$")) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  
  geno <- vroom("../data/genotypes.tsv") %>% 
    filter(snp==evar) %>%
    column_to_rownames("snp") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("ind", "genotype")) %>%
    mutate(ind=str_replace(ind, "NA", ""))
  
  dynvis <- left_join(expr, geno, by="ind") %>%
    mutate(logtpm=as.numeric(logcpm)) %>%
    mutate(genotype=factor(round(genotype), levels=c("0","1","2"))) %>%
    mutate(ind=factor(ind)) %>%
    mutate(day=factor(day, levels=seq(0, 15)))
  
  p <- ggplot(dynvis, aes(x=day, y=logcpm, fill=genotype)) + 
    geom_boxplot() +
    ylab(egene) +
    theme_classic()
  
  p
  
}

selective.append <- function(s) {
  if ("" %in% s) {
    s
  } else {
    append(s, "")
  }
}

dosage2allele <- function(d, ref.allele, alt.allele) {
  if (d == "0") {
    geno = paste0(ref.allele, ref.allele)
  } else if (d == "1") {
    geno = paste0(ref.allele, alt.allele)
  } else {
    geno = paste0(alt.allele, alt.allele)
  }
}


viz_residuals <- function(egene, evar, dataset, agg, 
                          nspcs=0, nclpcs=5, ctreg=F,
                          ref=NA, alt=NA, rsid=NA) {
  
  logexp <- if_else(dataset=="bulk", "logtpm", "logcpm")
  
  expr <- vroom(paste0("../data/", dataset, "/", agg, "/", logexp, ".tsv")) %>% 
    filter(gene==!!egene) %>% 
    column_to_rownames("gene") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("sample", logexp)) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  
  inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
    arrange(ind)
  
  geno <- vroom("../data/genotypes.tsv") %>% 
    filter(snp==evar) %>%
    column_to_rownames("snp") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("ind", "genotype")) %>%
    mutate(ind=str_replace(ind, "NA", "")) %>%
    right_join(inds, by="ind")
  
  if (agg == "day") {
    days <- expr %>% dplyr::select(sample) %>% 
      mutate(day=as.numeric(str_sub(sample, 7))) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  } else {
    medians <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/bin_medians.tsv")) %>%
      rename(sample=binind)
    days <- expr %>% dplyr::select(sample) %>% 
      left_join(medians, by="sample") %>%
      rename(day=t) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  }
  
  cvrts <- read_tsv(paste0("../data/", dataset, "/", agg, "/clpcs.tsv")) %>%
    dplyr::select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% dplyr::select(!ind) %>%
    arrange(sample) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (nspcs > 0) {
    samp.pcs <- read_tsv(paste0("../data/", dataset, "/", agg, "/pcs.tsv")) %>%
      dplyr::select(1:(nspcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
  if (ctreg) {
    if (dataset == "bulk") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK)
    } else if (dataset == "bulk7") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK) %>% filter(sample %in% expr$sample)
    } else if (dataset %in% c("pseudobulk-cm", "pseudobulk-cf")) {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% dplyr::select(!IPSC)
    } else {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% dplyr::select(!UNK)
    }
    cell.types <- setdiff(colnames(ct.props), "sample")
    cvrts <- left_join(cvrts, ct.props, by="sample")
  }
  
  clpc.str <- cvrts %>%
    dplyr::select(starts_with("cl_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  samppc.str <- cvrts %>%
    dplyr::select(starts_with("samp_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  if (ctreg) {
    ctprop.str <- cell.types %>%
      selective.append %>%
      paste(collapse="+")
  } else {
    ctprop.str <- ""
  }
  
  formula.str <- paste0(logexp, "~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day-g:day")
  
  design = left_join(geno, days, by="sample") %>% 
    rename(g=genotype) %>%
    left_join(cvrts, by="sample") %>%
    left_join(dplyr::select(expr, c(sample, !!logexp)), by="sample")
  
  residuals <- tibble("res"=residuals(lm(formula.str, design)), 
                      "sample"=design$sample, "genotype"=factor(round(design$g), levels=c("0", "1", "2")),
                      "day"=factor(days$day.ind, levels=seq(0, 15)))
  
  if (!is.na(ref)) {
    residuals <- residuals %>% rowwise() %>% 
      mutate(genotype=dosage2allele(genotype, ref.allele=ref, alt.allele=alt))
    # residuals$genotype[residuals$genotype=="AA"] <- "AA (n=14)"
    # residuals$genotype[residuals$genotype=="GA"] <- "GA (n=4)"
    # residuals$genotype[residuals$genotype=="GG"] <- "GG (n=1)"
  }
  
  p <- ggplot(residuals, aes(x=day, y=res, fill=genotype)) + 
    geom_boxplot() +
    ylab(egene) +
    theme_classic(base_size=20)
  if (str_extract(dataset, "[^-]+") == "pseudobulk") {
    p <- p + xlab("Pseudotime Bin")
  } else {
    p <- p + xlab("Collection Day")
  }
  if (!is.na(rsid)) {
    p <- p + ggtitle(rsid)
  }
  p
  
}

viz_residuals_nonlinear <- function(egene, evar, dataset, agg, 
                          nspcs=0, nclpcs=5, ctreg=F,
                          ref=NA, alt=NA, rsid=NA) {
  
  logexp <- if_else(dataset=="bulk", "logtpm", "logcpm")
  
  expr <- vroom(paste0("../data/", dataset, "/", agg, "/", logexp, ".tsv")) %>% 
    filter(gene==!!egene) %>% 
    column_to_rownames("gene") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("sample", logexp)) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  
  inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
    arrange(ind)
  
  geno <- vroom("../data/genotypes.tsv") %>% 
    filter(snp==evar) %>%
    column_to_rownames("snp") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("ind", "genotype")) %>%
    mutate(ind=str_replace(ind, "NA", "")) %>%
    right_join(inds, by="ind")
  
  if (agg == "day") {
    days <- expr %>% dplyr::select(sample) %>% 
      mutate(day=as.numeric(str_sub(sample, 7))) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  } else {
    medians <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/bin_medians.tsv")) %>%
      rename(sample=binind)
    days <- expr %>% select(sample) %>% 
      left_join(medians, by="sample") %>%
      rename(day=t) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  }
  days <- days %>% mutate(day.sq=day^2)
  
  cvrts <- read_tsv(paste0("../data/", dataset, "/", agg, "/clpcs.tsv")) %>%
    select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    arrange(sample) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (nspcs > 0) {
    samp.pcs <- read_tsv(paste0("../data/", dataset, "/", agg, "/pcs.tsv")) %>%
      select(1:(nspcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
  if (ctreg) {
    if (dataset == "bulk") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK)
    } else if (dataset == "bulk7") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK) %>% filter(sample %in% expr$sample)
    } else if (dataset %in% c("pseudobulk-cm", "pseudobulk-cf")) {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!IPSC)
    } else {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!UNK)
    }
    cell.types <- setdiff(colnames(ct.props), "sample")
    cvrts <- left_join(cvrts, ct.props, by="sample")
  }
  
  clpc.str <- cvrts %>%
    select(starts_with("cl_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  samppc.str <- cvrts %>%
    select(starts_with("samp_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  if (ctreg) {
    ctprop.str <- cell.types %>%
      selective.append %>%
      paste(collapse="+")
  } else {
    ctprop.str <- ""
  }
  
  formula.str <- paste0(logexp, "~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day+(", clpc.str, "g)*day.sq-g:day-g:day.sq")
  
  design = left_join(geno, days, by="sample") %>% 
    rename(g=genotype) %>%
    left_join(cvrts, by="sample") %>%
    left_join(select(expr, c(sample, !!logexp)), by="sample")
  
  residuals <- tibble("res"=residuals(lm(formula.str, design)), 
                      "sample"=design$sample, "genotype"=factor(round(design$g), levels=c("0", "1", "2")),
                      "day"=factor(days$day.ind, levels=seq(0, 15)))
  
  if (!is.na(ref)) {
    residuals <- residuals %>% rowwise() %>% 
      mutate(genotype=dosage2allele(genotype, ref.allele=ref, alt.allele=alt))
    # residuals$genotype[residuals$genotype=="AA"] <- "AA (n=4)"
    # residuals$genotype[residuals$genotype=="TA"] <- "TA (n=8)"
    # residuals$genotype[residuals$genotype=="TT"] <- "TT (n=7)"
  }
  
  p <- ggplot(residuals, aes(x=day, y=res, fill=genotype)) + 
    geom_boxplot() +
    ylab(egene) +
    theme_classic(base_size=20)
  
  if (str_extract(dataset, "[^-]+") == "pseudobulk") {
    p <- p + xlab("Pseudotime Bin")
  } else {
    p <- p + xlab("Collection Day")
  }
  if (!is.na(rsid)) {
    p <- p + ggtitle(rsid)
  }
  p
  
}

viz_residuals_ieqtl_notypes <- function(egene, evar, celltype, 
                                nspcs=0, nclpcs=5, ctreg=F,
                                ref=NA, alt=NA, rsid=NA) {
  
  expr <- vroom(paste0("../data/bulk/day/logtpm.tsv")) %>% 
    filter(gene==!!egene) %>% 
    column_to_rownames("gene") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("sample", "logtpm")) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  
  inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
    arrange(ind)
  
  geno <- vroom("../data/genotypes.tsv") %>% 
    filter(snp==evar) %>%
    column_to_rownames("snp") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("ind", "genotype")) %>%
    mutate(ind=str_replace(ind, "NA", "")) %>%
    right_join(inds, by="ind")
  
  cvrts <- read_tsv(paste0("../data/bulk/day/clpcs.tsv")) %>%
    select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (nspcs > 0) {
    samp.pcs <- read_tsv(paste0("../data/bulk/day/pcs.tsv")) %>%
      select(1:(nspcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")}) %>%
      mutate(ind=str_extract(sample, "[^_]+"))
    cvrts <- right_join(cvrts, samp.pcs, by="sample") 
  }
  ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
    select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
    arrange(sample) %>% select(!UNK)
  cell.types <- setdiff(colnames(ct.props), "sample")
  cvrts <- left_join(cvrts, ct.props, by="sample")
  
  stopifnot(all.equal(expr$sample, geno$sample))
  stopifnot(all.equal(expr$sample, cvrts$sample))
  
  # for model specification, need cell line PC and sample PC terms in a convenient format
  selective.append <- function(s) {
    if ("" %in% s) {
      s
    } else {
      append(s, "")
    }
  }
  clpc.str <- cvrts %>%
    select(starts_with("cl_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  samppc.str <- cvrts %>%
    select(starts_with("samp_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  if (ctreg) {
    ctprop.str <- setdiff(c("IPSC", "MES", "CMES", "PROG", "CM", "CF"), celltype) %>%
      selective.append %>%
      paste(collapse="+")
  } else {
    ctprop.str <- ""
  }
  
  formula.str <- paste0("logtpm~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*", celltype, "-g:", celltype)
  
  design = geno %>% 
    rename(g=genotype) %>%
    left_join(cvrts, by="sample") %>%
    left_join(select(expr, c(sample, logtpm)), by="sample")
  
  residuals <- tibble("res"=residuals(lm(formula.str, design)), 
                      "sample"=design$sample, "genotype"=factor(round(design$g), levels=c("0", "1", "2")),
                      "ctprop"=design[[celltype]])
  
  if (!is.na(ref)) {
    residuals <- residuals %>% rowwise() %>% 
      mutate(genotype=dosage2allele(genotype, ref.allele=ref, alt.allele=alt))
  }
  
  p <- ggplot(residuals, aes(x=as.numeric(ctprop), y=res, color=genotype)) + 
    geom_point() +
    ylab(egene) +
    theme_classic() +
    xlab(celltype) +
    geom_smooth(method="lm")
  if (!is.na(rsid)) {
    p <- p + ggtitle(rsid)
  }
  p
  
}

viz_residuals_ieqtl_regtypes <- function(egene, evar, celltype, 
                          nspcs=0, nclpcs=5, ctreg=T,
                          ref=NA, alt=NA, rsid=NA) {
  
  expr <- vroom(paste0("../data/bulk/day/logtpm.tsv")) %>% 
    filter(gene==!!egene) %>% 
    column_to_rownames("gene") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("sample", logexp)) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  
  inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
    arrange(ind)
  
  geno <- vroom("../data/genotypes.tsv") %>% 
    filter(snp==evar) %>%
    column_to_rownames("snp") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("ind", "genotype")) %>%
    mutate(ind=str_replace(ind, "NA", "")) %>%
    right_join(inds, by="ind")
  
  cvrts <- read_tsv(paste0("../data/", dataset, "/day/clpcs.tsv")) %>%
    select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (nspcs > 0) {
    samp.pcs <- read_tsv(paste0("../data/bulk/day/pcs.tsv")) %>%
      select(1:(nspcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")}) %>%
      mutate(ind=str_extract(sample, "[^_]+"))
    cvrts <- right_join(cvrts, samp.pcs, by="sample") 
  }
  ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
    select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
    arrange(sample) %>% select(!UNK)
  cell.types <- setdiff(colnames(ct.props), "sample")
  cvrts <- left_join(cvrts, ct.props, by="sample")
  
  stopifnot(all.equal(expr$sample, geno$sample))
  stopifnot(all.equal(expr$sample, cvrts$sample))
  
  # for model specification, need cell line PC and sample PC terms in a convenient format
  selective.append <- function(s) {
    if ("" %in% s) {
      s
    } else {
      append(s, "")
    }
  }
  clpc.str <- cvrts %>%
    select(starts_with("cl_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  samppc.str <- cvrts %>%
    select(starts_with("samp_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  if (ctreg) {
    ctprop.str <- setdiff(c("IPSC", "MES", "CMES", "PROG", "CM", "CF"), celltype) %>%
      selective.append %>%
      paste(collapse="+")
  } else {
    ctprop.str <- ""
  }
  
  formula.str <- paste0("logcpm~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*", celltype, "-g:", celltype)
  
  design = geno %>% 
    rename(g=genotype) %>%
    left_join(cvrts, by="sample") %>%
    left_join(select(expr, c(sample, !!logexp)), by="sample")
  
  residuals <- tibble("res"=residuals(lm(formula.str, design)), 
                      "sample"=design$sample, "genotype"=factor(round(design$g), levels=c("0", "1", "2")),
                      "ctprop"=design[[celltype]])
  
  if (!is.na(ref)) {
    residuals <- residuals %>% rowwise() %>% 
      mutate(genotype=dosage2allele(genotype, ref.allele=ref, alt.allele=alt))
  }
  
  p <- ggplot(residuals, aes(x=as.numeric(ctprop), y=res, color=genotype)) + 
    geom_point() +
    ylab(egene) +
    theme_classic() +
    xlab(celltype) +
    geom_smooth(method="lm")
  if (!is.na(rsid)) {
    p <- p + ggtitle(rsid)
  }
  p
  
}

viz_partial_regression <- function(egene, evar, dataset, agg, 
                                  nspcs=0, nclpcs=5, ctreg=F,
                                  ref=NA, alt=NA, rsid=NA,
                                  color.geno=F) {
  
  logexp <- if_else(dataset=="bulk", "logtpm", "logcpm")
  
  expr <- vroom(paste0("../data/", dataset, "/", agg, "/", logexp, ".tsv")) %>% 
    filter(gene==!!egene) %>% 
    column_to_rownames("gene") %>%
    apply(1, center.scale) %>%
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("sample", logexp)) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  
  inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
    arrange(ind)
  
  geno <- vroom("../data/genotypes.tsv") %>% 
    filter(snp==evar) %>%
    column_to_rownames("snp") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("ind", "genotype")) %>%
    mutate(ind=str_replace(ind, "NA", "")) %>%
    right_join(inds, by="ind")
  
  if (agg == "day") {
    days <- expr %>% select(sample) %>% 
      mutate(day=as.numeric(str_sub(sample, 7))) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  } else {
    medians <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/bin_medians.tsv")) %>%
      rename(sample=binind)
    days <- expr %>% select(sample) %>% 
      left_join(medians, by="sample") %>%
      rename(day=t) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  }
  
  cvrts <- read_tsv(paste0("../data/", dataset, "/", agg, "/clpcs.tsv")) %>%
    select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    arrange(sample) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (nspcs > 0) {
    samp.pcs <- read_tsv(paste0("../data/", dataset, "/", agg, "/pcs.tsv")) %>%
      select(1:(nspcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
  if (ctreg) {
    if (dataset == "bulk") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK)
    } else if (dataset == "bulk7") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK) %>% filter(sample %in% expr$sample)
    } else if (dataset %in% c("pseudobulk-cm", "pseudobulk-cf")) {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!IPSC)
    } else {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!UNK)
    }
    cell.types <- setdiff(colnames(ct.props), "sample")
    cvrts <- left_join(cvrts, ct.props, by="sample")
  }
  
  clpc.str <- cvrts %>%
    select(starts_with("cl_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  samppc.str <- cvrts %>%
    select(starts_with("samp_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  if (ctreg) {
    ctprop.str <- cell.types %>%
      selective.append %>%
      paste(collapse="+")
  } else {
    ctprop.str <- ""
  }
  
  intx.y.str <- paste0(logexp, "~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day-g:day")
  intx.x.str <- paste0(" intx ~ 1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day-g:day")
  
  design = left_join(geno, days, by="sample") %>% 
    rename(g=genotype) %>%
    left_join(cvrts, by="sample") %>%
    left_join(select(expr, c(sample, !!logexp)), by="sample")
  resid.y <- residuals(lm(intx.y.str, design))
  
  design.expanded <- model.matrix(as.formula(paste0("~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day")), design) %>%
    as_tibble() %>% select(-`(Intercept)`) %>%
    rename(intx="g:day")
  resid.x <- residuals(lm(intx.x.str, design.expanded))
  
  residuals <- tibble("res"=resid.x, 
                      "partial.res"=resid.y,
                      "sample"=design$sample, "genotype"=factor(round(design$g), levels=c("0", "1", "2")),
                      "day"=factor(days$day.ind, levels=seq(0, 15)))
  
  if (!is.na(ref)) {
    residuals <- residuals %>% rowwise() %>% 
      mutate(genotype=dosage2allele(genotype, ref.allele=ref, alt.allele=alt))
  }
  
  if (color.geno) {
    p <- ggplot(residuals, aes(x=res, y=partial.res, color=genotype)) + 
      geom_point() +
      ylab(expression(resid(y%~%x[-i]))) +
      xlab(expression(resid(x[i]%~%x[-i]))) +
      theme_classic(base_size=14) +
      scale_color_manual(values=c("#0077BB", "#E08EF3", "#BF0202"))
  } else {
    p <- ggplot(residuals, aes(x=res, y=partial.res)) + 
      geom_point() +
      ylab(expression(Y[""%.%group("[", g%*%t, "]")])) +
      xlab(expression(X[g%*%t%.%group("[", g%*%t, "]")])) +
      theme_classic(base_size=14)
  }
  
  if (!is.na(rsid)) {
    p <- p + ggtitle(paste0(egene, "--", rsid))
  }
  p
}

# viz_partial_residuals <- function(egene, evar, dataset, agg, 
#                                   nspcs=0, nclpcs=5, ctreg=F,
#                                   ref=NA, alt=NA, rsid=NA,
#                                   color.geno=F) {
#   
#   logexp <- if_else(dataset=="bulk", "logtpm", "logcpm")
#   
#   expr <- vroom(paste0("../data/", dataset, "/", agg, "/", logexp, ".tsv")) %>% 
#     filter(gene==!!egene) %>% 
#     column_to_rownames("gene") %>% t %>% 
#     as_tibble(rownames="ind") %>% 
#     `colnames<-`(c("sample", logexp)) %>%
#     mutate(ind=str_extract(sample, "[^_]+"))
#   
#   inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
#     arrange(ind)
#   
#   geno <- vroom("../data/genotypes.tsv") %>% 
#     filter(snp==evar) %>%
#     column_to_rownames("snp") %>% t %>% 
#     as_tibble(rownames="ind") %>% 
#     `colnames<-`(c("ind", "genotype")) %>%
#     mutate(ind=str_replace(ind, "NA", "")) %>%
#     right_join(inds, by="ind")
#   
#   if (agg == "day") {
#     days <- expr %>% select(sample) %>% 
#       mutate(day=as.numeric(str_sub(sample, 7))) %>%
#       mutate(day.ind=as.numeric(str_sub(sample, 7)))
#   } else {
#     medians <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/bin_medians.tsv")) %>%
#       rename(sample=binind)
#     days <- expr %>% select(sample) %>% 
#       left_join(medians, by="sample") %>%
#       rename(day=t) %>%
#       mutate(day.ind=as.numeric(str_sub(sample, 7)))
#   }
#   
#   cvrts <- read_tsv(paste0("../data/", dataset, "/", agg, "/clpcs.tsv")) %>%
#     select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
#     right_join(inds, by="ind") %>% select(!ind) %>%
#     arrange(sample) %>%
#     rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
#   if (nspcs > 0) {
#     samp.pcs <- read_tsv(paste0("../data/", dataset, "/", agg, "/pcs.tsv")) %>%
#       select(1:(nspcs+1)) %>%
#       rename_with(function(s){str_replace(s,"PC","samp_PC_")})
#     cvrts <- left_join(cvrts, samp.pcs, by="sample")
#   }
#   if (ctreg) {
#     if (dataset == "bulk") {
#       ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
#         select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
#         arrange(sample) %>% select(!UNK)
#     } else if (dataset == "bulk7") {
#       ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
#         select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
#         arrange(sample) %>% select(!UNK) %>% filter(sample %in% expr$sample)
#     } else if (dataset %in% c("pseudobulk-cm", "pseudobulk-cf")) {
#       ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
#         arrange(sample) %>% select(!IPSC)
#     } else {
#       ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
#         arrange(sample) %>% select(!UNK)
#     }
#     cell.types <- setdiff(colnames(ct.props), "sample")
#     cvrts <- left_join(cvrts, ct.props, by="sample")
#   }
#   
#   clpc.str <- cvrts %>%
#     select(starts_with("cl_PC")) %>%
#     colnames() %>%
#     selective.append %>%
#     paste(collapse="+")
#   samppc.str <- cvrts %>%
#     select(starts_with("samp_PC")) %>%
#     colnames() %>%
#     selective.append %>%
#     paste(collapse="+")
#   if (ctreg) {
#     ctprop.str <- cell.types %>%
#       selective.append %>%
#       paste(collapse="+")
#   } else {
#     ctprop.str <- ""
#   }
#   
#   intx.y.str <- paste0(logexp, "~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day-g:day")
#   intx.x.str <- paste0(" intx ~ 1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day-g:day")
#   
#   design = left_join(geno, days, by="sample") %>% 
#     rename(g=genotype) %>%
#     left_join(cvrts, by="sample") %>%
#     left_join(select(expr, c(sample, !!logexp)), by="sample")
#   resid.y <- residuals(lm(intx.y.str, design))
#   
#   design.expanded <- model.matrix(as.formula(paste0("~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day")), design) %>%
#     as_tibble() %>% select(-`(Intercept)`) %>%
#     rename(intx="g:day")
#   
#   residuals <- tibble("res"=design.expanded$intx, 
#                       "partial.res"=resid.y,
#                       "sample"=design$sample, "genotype"=factor(round(design$g), levels=c("0", "1", "2")),
#                       "day"=factor(days$day.ind, levels=seq(0, 15)))
#   
#   if (!is.na(ref)) {
#     residuals <- residuals %>% rowwise() %>% 
#       mutate(genotype=dosage2allele(genotype, ref.allele=ref, alt.allele=alt))
#   }
#   
#   if (color.geno) {
#     p <- ggplot(residuals, aes(x=res, y=partial.res, color=genotype)) + 
#       geom_point() +
#       ylab(expression(resid(y%~%x[-i]))) +
#       xlab(expression(resid(x[i]%~%x[-i]))) +
#       theme_classic(base_size=14) +
#       scale_color_manual(values=c("#0077BB", "#E08EF3", "#BF0202"))
#   } else {
#     p <- ggplot(residuals, aes(x=res, y=partial.res)) + 
#       geom_point() +
#       ylab(expression(resid(y%~%x[-g%*%t]))) +
#       xlab(expression(x[g%*%t])) +
#       theme_classic(base_size=14)
#   }
#   
#   if (!is.na(rsid)) {
#     p <- p + ggtitle(paste0(egene, "--", rsid))
#   }
#   p
# }

viz_partial_residuals <- function(egene, evar, dataset, agg, 
                                  nspcs=0, nclpcs=5, ctreg=F,
                                  ref=NA, alt=NA, rsid=NA,
                                  color.geno=F) {
  
  logexp <- if_else(dataset=="bulk", "logtpm", "logcpm")
  
  expr <- vroom(paste0("../data/", dataset, "/", agg, "/", logexp, ".tsv")) %>% 
    filter(gene==!!egene) %>% 
    column_to_rownames("gene") %>% 
    apply(1, center.scale) %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("sample", logexp)) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  
  inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
    arrange(ind)
  
  geno <- vroom("../data/genotypes.tsv") %>% 
    filter(snp==evar) %>%
    column_to_rownames("snp") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("ind", "genotype")) %>%
    mutate(ind=str_replace(ind, "NA", "")) %>%
    right_join(inds, by="ind")
  
  if (agg == "day") {
    days <- expr %>% select(sample) %>% 
      mutate(day=as.numeric(str_sub(sample, 7))) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  } else {
    medians <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/bin_medians.tsv")) %>%
      rename(sample=binind)
    days <- expr %>% select(sample) %>% 
      left_join(medians, by="sample") %>%
      rename(day=t) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  }
  
  cvrts <- read_tsv(paste0("../data/", dataset, "/", agg, "/clpcs.tsv")) %>%
    select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    arrange(sample) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (nspcs > 0) {
    samp.pcs <- read_tsv(paste0("../data/", dataset, "/", agg, "/pcs.tsv")) %>%
      select(1:(nspcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
  if (ctreg) {
    if (dataset == "bulk") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK)
    } else if (dataset == "bulk7") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK) %>% filter(sample %in% expr$sample)
    } else if (dataset %in% c("pseudobulk-cm", "pseudobulk-cf")) {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!IPSC)
    } else {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!UNK)
    }
    cell.types <- setdiff(colnames(ct.props), "sample")
    cvrts <- left_join(cvrts, ct.props, by="sample")
  }
  
  clpc.str <- cvrts %>%
    select(starts_with("cl_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  samppc.str <- cvrts %>%
    select(starts_with("samp_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  if (ctreg) {
    ctprop.str <- cell.types %>%
      selective.append %>%
      paste(collapse="+")
  } else {
    ctprop.str <- ""
  }
  
  model.str <- paste0(logexp, "~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day")
  design.str <- paste0("~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day")
  
  reg.df = left_join(geno, days, by="sample") %>% 
    rename(g=genotype) %>%
    left_join(cvrts, by="sample") %>%
    left_join(select(expr, c(sample, !!logexp)), by="sample")
  
  m <- lm(model.str, reg.df)
  
  intx = model.matrix(as.formula(design.str), reg.df)[,"g:day"]
  res = residuals(m)
  pres = res + intx * coefficients(m)['g:day']
  
  partial_res <- tibble("partial.res"=pres, "intx"=intx)
  
  p <- ggplot(partial_res, aes(x=intx, y=partial.res)) + 
        geom_point() +
        ylab(expression(hat(beta)[g%*%t]*X[g%*%t]+res)) +
        xlab(expression(X[g%*%t])) +
        theme_classic(base_size=14)
  
  if (!is.na(rsid)) {
    p <- p + ggtitle(paste0(egene, "--", rsid))
  }
  p
}

viz_partial_regression_nonlinear <- function(egene, evar, dataset, agg, 
                                            nspcs=0, nclpcs=5, ctreg=F,
                                            ref=NA, alt=NA, rsid=NA,
                                            color.geno=F) {
  
  logexp <- if_else(dataset=="bulk", "logtpm", "logcpm")
  
  expr <- vroom(paste0("../data/", dataset, "/", agg, "/", logexp, ".tsv")) %>% 
    filter(gene==!!egene) %>% 
    column_to_rownames("gene") %>% 
    apply(1, center.scale) %>%
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("sample", logexp)) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  
  inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
    arrange(ind)
  
  geno <- vroom("../data/genotypes.tsv") %>% 
    filter(snp==evar) %>%
    column_to_rownames("snp") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("ind", "genotype")) %>%
    mutate(ind=str_replace(ind, "NA", "")) %>%
    right_join(inds, by="ind")
  
  if (agg == "day") {
    days <- expr %>% select(sample) %>% 
      mutate(day=as.numeric(str_sub(sample, 7))) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  } else {
    medians <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/bin_medians.tsv")) %>%
      rename(sample=binind)
    days <- expr %>% select(sample) %>% 
      left_join(medians, by="sample") %>%
      rename(day=t) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  }
  days <- days %>% mutate(day.sq=day^2)
  
  cvrts <- read_tsv(paste0("../data/", dataset, "/", agg, "/clpcs.tsv")) %>%
    select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    arrange(sample) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (nspcs > 0) {
    samp.pcs <- read_tsv(paste0("../data/", dataset, "/", agg, "/pcs.tsv")) %>%
      select(1:(nspcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
  if (ctreg) {
    if (dataset == "bulk") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK)
    } else if (dataset == "bulk7") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK) %>% filter(sample %in% expr$sample)
    } else if (dataset %in% c("pseudobulk-cm", "pseudobulk-cf")) {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!IPSC)
    } else {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!UNK)
    }
    cell.types <- setdiff(colnames(ct.props), "sample")
    cvrts <- left_join(cvrts, ct.props, by="sample")
  }
  
  clpc.str <- cvrts %>%
    select(starts_with("cl_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  samppc.str <- cvrts %>%
    select(starts_with("samp_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  if (ctreg) {
    ctprop.str <- cell.types %>%
      selective.append %>%
      paste(collapse="+")
  } else {
    ctprop.str <- ""
  }
  
  intx.y.str <- paste0(logexp, "~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day+(", clpc.str, "g)*day.sq-g:day.sq")
  intx.x.str <- paste0("intx ~ 1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day+(", clpc.str, "g)*day.sq-g:day.sq")
  
  design = left_join(geno, days, by="sample") %>% 
    rename(g=genotype) %>%
    left_join(cvrts, by="sample") %>%
    left_join(select(expr, c(sample, !!logexp)), by="sample")
  resid.y <- residuals(lm(intx.y.str, design))
  
  design.expanded <- model.matrix(as.formula(paste0("~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day+(", clpc.str, "g)*day.sq")), design) %>%
    as_tibble() %>% select(-`(Intercept)`) %>%
    rename(intx="g:day.sq")
  resid.x <- residuals(lm(intx.x.str, design.expanded))
  
  residuals <- tibble("res"=resid.x, 
                      "partial.res"=resid.y,
                      "sample"=design$sample, "genotype"=factor(round(design$g), levels=c("0", "1", "2")),
                      "day"=factor(days$day.ind, levels=seq(0, 15)))
  
  if (!is.na(ref)) {
    residuals <- residuals %>% rowwise() %>% 
      mutate(genotype=dosage2allele(genotype, ref.allele=ref, alt.allele=alt))
  }
  
  if (color.geno) {
    p <- ggplot(residuals, aes(x=res, y=partial.res, color=genotype)) + 
      geom_point() +
      ylab(expression(resid(y%~%x[-g%*%t^2]))) +
      xlab(expression(resid(x[g%*%t^2]%~%x[-g%*%t^2]))) +
      theme_classic(base_size=14) +
      scale_color_manual(values=c("#0077BB", "#E08EF3", "#BF0202"))
  } else {
    p <- ggplot(residuals, aes(x=res, y=partial.res)) + 
      geom_point() +
      ylab(expression(Y[""%.%group("[", g%*%t^2, "]")])) +
      xlab(expression(X[g%*%t^2%.%group("[", g%*%t^2, "]")])) +
      theme_classic(base_size=14)
  }
  
  if (!is.na(rsid)) {
    p <- p + ggtitle(paste0(egene, "--", rsid))
  }
  p
}

# viz_partial_residuals_nonlinear <- function(egene, evar, dataset, agg, 
#                                              nspcs=0, nclpcs=5, ctreg=F,
#                                              ref=NA, alt=NA, rsid=NA,
#                                              color.geno=F) {
#   
#   logexp <- if_else(dataset=="bulk", "logtpm", "logcpm")
#   
#   expr <- vroom(paste0("../data/", dataset, "/", agg, "/", logexp, ".tsv")) %>% 
#     filter(gene==!!egene) %>% 
#     column_to_rownames("gene") %>% t %>% 
#     as_tibble(rownames="ind") %>% 
#     `colnames<-`(c("sample", logexp)) %>%
#     mutate(ind=str_extract(sample, "[^_]+"))
#   
#   inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
#     arrange(ind)
#   
#   geno <- vroom("../data/genotypes.tsv") %>% 
#     filter(snp==evar) %>%
#     column_to_rownames("snp") %>% t %>% 
#     as_tibble(rownames="ind") %>% 
#     `colnames<-`(c("ind", "genotype")) %>%
#     mutate(ind=str_replace(ind, "NA", "")) %>%
#     right_join(inds, by="ind")
#   
#   if (agg == "day") {
#     days <- expr %>% select(sample) %>% 
#       mutate(day=as.numeric(str_sub(sample, 7))) %>%
#       mutate(day.ind=as.numeric(str_sub(sample, 7)))
#   } else {
#     medians <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/bin_medians.tsv")) %>%
#       rename(sample=binind)
#     days <- expr %>% select(sample) %>% 
#       left_join(medians, by="sample") %>%
#       rename(day=t) %>%
#       mutate(day.ind=as.numeric(str_sub(sample, 7)))
#   }
#   days <- days %>% mutate(day.sq=day^2)
#   
#   cvrts <- read_tsv(paste0("../data/", dataset, "/", agg, "/clpcs.tsv")) %>%
#     select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
#     right_join(inds, by="ind") %>% select(!ind) %>%
#     arrange(sample) %>%
#     rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
#   if (nspcs > 0) {
#     samp.pcs <- read_tsv(paste0("../data/", dataset, "/", agg, "/pcs.tsv")) %>%
#       select(1:(nspcs+1)) %>%
#       rename_with(function(s){str_replace(s,"PC","samp_PC_")})
#     cvrts <- left_join(cvrts, samp.pcs, by="sample")
#   }
#   if (ctreg) {
#     if (dataset == "bulk") {
#       ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
#         select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
#         arrange(sample) %>% select(!UNK)
#     } else if (dataset == "bulk7") {
#       ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
#         select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
#         arrange(sample) %>% select(!UNK) %>% filter(sample %in% expr$sample)
#     } else if (dataset %in% c("pseudobulk-cm", "pseudobulk-cf")) {
#       ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
#         arrange(sample) %>% select(!IPSC)
#     } else {
#       ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
#         arrange(sample) %>% select(!UNK)
#     }
#     cell.types <- setdiff(colnames(ct.props), "sample")
#     cvrts <- left_join(cvrts, ct.props, by="sample")
#   }
#   
#   clpc.str <- cvrts %>%
#     select(starts_with("cl_PC")) %>%
#     colnames() %>%
#     selective.append %>%
#     paste(collapse="+")
#   samppc.str <- cvrts %>%
#     select(starts_with("samp_PC")) %>%
#     colnames() %>%
#     selective.append %>%
#     paste(collapse="+")
#   if (ctreg) {
#     ctprop.str <- cell.types %>%
#       selective.append %>%
#       paste(collapse="+")
#   } else {
#     ctprop.str <- ""
#   }
#   
#   intx.y.str <- paste0(logexp, "~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day+(", clpc.str, "g)*day.sq-g:day.sq")
#   intx.x.str <- paste0("intx ~ 1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day+(", clpc.str, "g)*day.sq-g:day.sq")
#   
#   design = left_join(geno, days, by="sample") %>% 
#     rename(g=genotype) %>%
#     left_join(cvrts, by="sample") %>%
#     left_join(select(expr, c(sample, !!logexp)), by="sample")
#   resid.y <- residuals(lm(intx.y.str, design))
#   
#   design.expanded <- model.matrix(as.formula(paste0("~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day+(", clpc.str, "g)*day.sq")), design) %>%
#     as_tibble() %>% select(-`(Intercept)`) %>%
#     rename(intx="g:day.sq")
#   
#   residuals <- tibble("xi"=design.expanded$intx, 
#                       "partial.res"=resid.y,
#                       "sample"=design$sample, "genotype"=factor(round(design$g), levels=c("0", "1", "2")),
#                       "day"=factor(days$day.ind, levels=seq(0, 15)))
#   
#   if (!is.na(ref)) {
#     residuals <- residuals %>% rowwise() %>% 
#       mutate(genotype=dosage2allele(genotype, ref.allele=ref, alt.allele=alt))
#   }
#   
#   if (color.geno) {
#     p <- ggplot(residuals, aes(x=xi, y=partial.res, color=genotype)) + 
#       geom_point() +
#       ylab(expression(resid(y%~%x[-g%*%t^2]))) +
#       xlab(expression(x[g%*%t^2])) +
#       theme_classic(base_size=14) +
#       scale_color_manual(values=c("#0077BB", "#E08EF3", "#BF0202"))
#   } else {
#     p <- ggplot(residuals, aes(x=xi, y=partial.res)) + 
#       geom_point() +
#       ylab(expression(resid(y%~%x[-g%*%t^2]))) +
#       xlab(expression(x[g%*%t^2])) +
#       theme_classic(base_size=14)
#   }
#   
#   if (!is.na(rsid)) {
#     p <- p + ggtitle(paste0(egene, "--", rsid))
#   }
#   p
# }

viz_partial_residuals_nonlinear <- function(egene, evar, dataset, agg, 
                                            nspcs=0, nclpcs=5, ctreg=F,
                                            ref=NA, alt=NA, rsid=NA,
                                            color.geno=F) {
  
  logexp <- if_else(dataset=="bulk", "logtpm", "logcpm")
  
  expr <- vroom(paste0("../data/", dataset, "/", agg, "/", logexp, ".tsv")) %>% 
    filter(gene==!!egene) %>% 
    column_to_rownames("gene") %>% 
    apply(1, center.scale) %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("sample", logexp)) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  
  inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
    arrange(ind)
  
  geno <- vroom("../data/genotypes.tsv") %>% 
    filter(snp==evar) %>%
    column_to_rownames("snp") %>% t %>% 
    as_tibble(rownames="ind") %>% 
    `colnames<-`(c("ind", "genotype")) %>%
    mutate(ind=str_replace(ind, "NA", "")) %>%
    right_join(inds, by="ind")
  
  if (agg == "day") {
    days <- expr %>% select(sample) %>% 
      mutate(day=as.numeric(str_sub(sample, 7))) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  } else {
    medians <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/bin_medians.tsv")) %>%
      rename(sample=binind)
    days <- expr %>% select(sample) %>% 
      left_join(medians, by="sample") %>%
      rename(day=t) %>%
      mutate(day.ind=as.numeric(str_sub(sample, 7)))
  }
  days <- days %>% mutate(day.sq=day^2)
  
  cvrts <- read_tsv(paste0("../data/", dataset, "/", agg, "/clpcs.tsv")) %>%
    select(1:(nclpcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    arrange(sample) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (nspcs > 0) {
    samp.pcs <- read_tsv(paste0("../data/", dataset, "/", agg, "/pcs.tsv")) %>%
      select(1:(nspcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
  if (ctreg) {
    if (dataset == "bulk") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK)
    } else if (dataset == "bulk7") {
      ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
        select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
        arrange(sample) %>% select(!UNK) %>% filter(sample %in% expr$sample)
    } else if (dataset %in% c("pseudobulk-cm", "pseudobulk-cf")) {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!IPSC)
    } else {
      ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
        arrange(sample) %>% select(!UNK)
    }
    cell.types <- setdiff(colnames(ct.props), "sample")
    cvrts <- left_join(cvrts, ct.props, by="sample")
  }
  
  clpc.str <- cvrts %>%
    select(starts_with("cl_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  samppc.str <- cvrts %>%
    select(starts_with("samp_PC")) %>%
    colnames() %>%
    selective.append %>%
    paste(collapse="+")
  if (ctreg) {
    ctprop.str <- cell.types %>%
      selective.append %>%
      paste(collapse="+")
  } else {
    ctprop.str <- ""
  }
  
  model.str <- paste0(logexp, "~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day+(", clpc.str, "g)*day.sq")
  design.str <- paste0("~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day+(", clpc.str, "g)*day.sq")
  
  reg.df = left_join(geno, days, by="sample") %>% 
    rename(g=genotype) %>%
    left_join(cvrts, by="sample") %>%
    left_join(select(expr, c(sample, !!logexp)), by="sample")
  
  m <- lm(model.str, reg.df)
  
  intx = model.matrix(as.formula(design.str), reg.df)[,"g:day.sq"]
  res = residuals(m)
  pres = res + intx * coefficients(m)['g:day.sq']
  
  partial_res <- tibble("partial.res"=pres, "intx"=intx)
  
  p <- ggplot(partial_res, aes(x=intx, y=partial.res)) + 
    geom_point() +
    ylab(expression(hat(beta)[g%*%t^2]*X[g%*%t^2]+res)) +
    xlab(expression(X[g%*%t^2])) +
    theme_classic(base_size=14)
  
  if (!is.na(rsid)) {
    p <- p + ggtitle(paste0(egene, "--", rsid))
  }
  p
}