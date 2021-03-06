library(tidyverse)
library(vroom)

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
  
  formula.str <- paste0(logexp, "~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*day-g:day")
  
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
  }
  
  p <- ggplot(residuals, aes(x=day, y=res, fill=genotype)) + 
    geom_boxplot() +
    ylab(egene) +
    theme_classic()
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
  }
  
  p <- ggplot(residuals, aes(x=day, y=res, fill=genotype)) + 
    geom_boxplot() +
    ylab(egene) +
    theme_classic()
  
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