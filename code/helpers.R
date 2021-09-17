### HELPER FUNCTIONS ###

make.normal <- function(day.str) {
  str_replace(day.str, "Day ", "day")
}

a2ma.row <- function(gt) {
  if (sum(gt) > length(gt)) {
    gt <- 2 - gt
  }
  gt
}

a2ma.mat <- function(mat) {
  # input/ output is snp x individual matrix
  mat = apply(mat, 1, a2ma.row)
  t(mat)
}

ac2maf <- function(genotypes) {
  genotypes = a2ma.mat(genotypes)
  maf = rowSums(genotypes)/(2*ncol(genotypes))
}

threshold <- function(x, thresh, op="ge") {
  if (op=="ge") {
    pass.bool <- (x >= thresh)
  } else if (op=="gt") {
    pass.bool <- (x > thresh)
  } else if (op=="lt") {
    pass.bool <- (x < thresh)
  }
  pass.bool
}

lognorm <- function(c, scale.factor=1e6, pseudocounts=1) {
  # log2(CPM+1) for a (individual, cell type) pair
  log2(scale.factor*c/sum(c)+pseudocounts)
}

center.scale <- function(g) {
  # center and scale within a gene's expression vector
  if (sd(g) == 0) {
    g - g
  } else {
    (g-mean(g))/sd(g)
  }
}

leiden2type <- function(c) {
  if (c == 0) {
    "EPDC"
  } else if (c == 1) {
    "iPSC"
  } else if (c == 2) {
    "meso"
  } else if (c == 3) {
    "prog"
  } else if (c == 4) {
    "CM"
  } else if (c == 5) {
    "cardiomes"
  } else if (c == 6) {
    "EMT"
  }
}

ct2 <- function(x) {
  if (x %in% c("EMT", "meso", "cardiomes", "prog")) {
    "mes"
  } else {
    x
  }
}