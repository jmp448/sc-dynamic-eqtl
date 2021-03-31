library(Seurat)
library(tidyverse)
library(loomR)
library(sceasy)
library(reticulate)
use_condaenv('singlecell')
loompy <- reticulate::import('loompy')

set.seed(2021)

sc <- readRDS("../data/seurat.annotated.rds")

sc@assays$SCT@data <- sc@assays$SCT@data[VariableFeatures(object = sc)[1:1000],]
sc@assays$SCT@scale.data <- sc@assays$SCT@scale.data[VariableFeatures(object = sc)[1:1000],]
sc <- RunPCA(sc, npcs=100, features=VariableFeatures(object = sc)[1:1000])
sc@assays$SCT@meta.features <- sc@assays$SCT@meta.features[VariableFeatures(object = sc)[1:1000],]
sceasy::convertFormat(sc, from="seurat", to="anndata",
                      assay="SCT", main_layer='data', transfer_layers=c('scale.data'),
                      outFile='../data/seurat.annotated.sct.h5ad')


sc <- readRDS("../data/seurat.conservative.rds")
sc@assays$SCT@data <- sc@assays$SCT@data[VariableFeatures(object = sc)[1:1000],]
sc@assays$SCT@scale.data <- sc@assays$SCT@scale.data[VariableFeatures(object = sc)[1:1000],]
sc@assays$SCT@meta.features <- sc@assays$SCT@meta.features[VariableFeatures(object = sc)[1:1000],]
sceasy::convertFormat(sc, from="seurat", to="anndata",
                      assay="SCT", main_layer='data', transfer_layers=c('scale.data'),
                      outFile='../data/seurat.conservative.h5ad')

sc <- readRDS("../data/seurat.liberal.rds")
sc@assays$SCT@data <- sc@assays$SCT@data[VariableFeatures(object = sc)[1:1000],]
sc@assays$SCT@scale.data <- sc@assays$SCT@scale.data[VariableFeatures(object = sc)[1:1000],]
sc@assays$SCT@meta.features <- sc@assays$SCT@meta.features[VariableFeatures(object = sc)[1:1000],]
sceasy::convertFormat(sc, from="seurat", to="anndata",
                      assay="SCT", main_layer='data', transfer_layers=c('scale.data'),
                      outFile='../data/seurat.liberal.h5ad')

sc <- readRDS("../data/seurat.normalized.rds")
sc@assays$SCT@data <- sc@assays$SCT@data[VariableFeatures(object = sc)[1:1000],]
sc@assays$SCT@scale.data <- sc@assays$SCT@scale.data[VariableFeatures(object = sc)[1:1000],]
sc@assays$SCT@meta.features <- sc@assays$SCT@meta.features[VariableFeatures(object = sc)[1:1000],]
sceasy::convertFormat(sc, from="seurat", to="anndata",
                      assay="SCT", main_layer='data', transfer_layers=c('scale.data'),
                      outFile='../data/seurat.sctransformed.h5ad')

sc <- readRDS("../data/seurat.mitonormalized.rds")
sc@assays$SCT@data <- sc@assays$SCT@data[VariableFeatures(object = sc)[1:1000],]
sc@assays$SCT@scale.data <- sc@assays$SCT@scale.data[VariableFeatures(object = sc)[1:1000],]
sc@assays$SCT@meta.features <- sc@assays$SCT@meta.features[VariableFeatures(object = sc)[1:1000],]
sceasy::convertFormat(sc, from="seurat", to="anndata",
                      assay="SCT", main_layer='data', transfer_layers=c('scale.data'),
                      outFile='../data/seurat.sctransformed.mito.h5ad')
