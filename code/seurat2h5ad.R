library(Seurat)
library(loomR)
library(sceasy)
library(reticulate)
use_condaenv('singlecell')
loompy <- reticulate::import('loompy')

# all data 
sc <- readRDS("../data/seurat.normalized.rds")
sceasy::convertFormat(sc, from="seurat", to="anndata",
                      assay="SCT", main_layer='data', transfer_layers=c('scale.data'),
                      outFile='../data/seurat.normalized.h5ad')

# all data (including unknown cells)
# sc <- readRDS("../data/seurat.annotated.rds")
# sceasy::convertFormat(sc, from="seurat", to="anndata",
#                        outFile='../data/seurat.annotated.h5ad')

# # cardiac data (excluding unknown cells)
# sc_cardiac <- readRDS("../data/seurat.cardiac.rds")
# sceasy::convertFormat(sc_cardiac, from="seurat", to="anndata",
#                        outFile='../data/seurat.cardiac.h5ad')
