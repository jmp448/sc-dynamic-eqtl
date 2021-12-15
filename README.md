# Single-Cell Sequencing Reveals Lineage-Specific Dynamic Genetic Regulation of Gene Expression During Human Cardiomyocyte Differentiation
Reem Elorbany,  Joshua M Popp, Katherine Rhodes, Benjamin J Strober, Kenneth Barr, Guanghao Qi, Yoav Gilad, and Alexis Battle

A [workflowr](https://github.com/jdblischak/workflowr) project.

Preprint available online at [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.03.446970v1)

## Preprocessing of Single Cell Data
- Quality control filtering of single cell data is done at `code/create_seurat.R`
  - The (cells x genes) matrix of UMI counts produced by this file is available online at GEO Series GSE175634 as a .mtx file (`GSE175634_cell_counts.mtx.gz`, cell and gene indices at `GSE175634_cell_indices.tsv.gz` and `GSE175634_gene_indices.tsv.gz` respectively)
- Normalization of the single cell counts data is done at `code/normalize_seurat.R`
  - The (cells x genes) matrix of *corrected* UMI counts produced by this file is available as a .mtx file (`GSE175634_cell_counts_sctransform.mtx.gz`, same cell and gene indices as the uncorrected counts)
  - The (cells x highly variable genes) matrix of pearson residuals produced by this file is available as a .tsv file (`GSE175634_pearson_residuals_sctransform.mtx.gz`)
  
## Trajectory Inference
- Trajectory inference is performed at `analysis/ti.ipynb`

## Dynamic eQTL Calling
- Pseudobulk aggregation, normalization, and computation of PC covariates for both pseudobulk and bulk data are done at `code/preprocessing.R`)
- Dynamic eQTL calling is done at `code/eqtl_dynamic.R` for linear eQTLs, `code/eqtl_dynamic_nonlinear.R` for nonlinear dynamic eQTLs, and multiple testing correction is done at `code/mtc_dynamic.R` for linear eQTLs, `code/mtc_dynamic_ttest.R` for nonlinear dynamic eQTLs

## Cell Type Interaction eQTL Calling
- Cell type interaction eQTL calling is performed at `code/eqtl_interaction.R`, with multiple testing at `code/mtc_interaction.R`

## Further Analyses
- Overlap of dynamic eQTLs with GWAS hits is done at `code/gwas.R`
