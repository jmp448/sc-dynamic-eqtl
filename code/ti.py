import scanpy as sc
import numpy as np
import sys

# adata = sc.read_h5ad("../data/seurat.annotated.sct.h5ad")
# sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)
# adata.write_h5ad("../data/scanpy.neighbors.h5ad")
adata = sc.read_h5ad("../data/scanpy.neighbors.h5ad")
sc.tl.leiden(adata, resolution=0.45, key_added="leiden_045")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_05")
sc.tl.leiden(adata, resolution=0.55, key_added="leiden_055")
sc.tl.leiden(adata, resolution=0.6, key_added="leiden_06")
sc.tl.leiden(adata, resolution=0.65, key_added="leiden_065")
adata.write_h5ad("../data/scanpy.clustered.h5ad")

# sc.tl.paga(adata, groups='leiden')
# adata.write_h5ad("../data/scanpy.paga.h5ad")
# sc.tl.draw_graph(adata, init_pos='paga')
# sc.tl.diffmap(adata)
# adata.write_h5ad("../data/scanpy.pseudotime.h5ad")

# adata = sc.read_h5ad("../data/scanpy.paga.h5ad")
# adata = adata[adata.obs['leiden']!="14",:]
# sc.pp.scale(adata)
# sc.tl.pca(adata)
# sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)
# adata.write_h5ad("../data/scanpy.subbed.neighbors.h5ad")
# adata = sc.read_h5ad("../data/scanpy.subbed.neighbors.h5ad")
# sc.tl.leiden(adata, resolution=0.4)
# sc.tl.paga(adata, groups='leiden')
# adata.write_h5ad("../data/scanpy.subbed.sct.h5ad")

# adata = sc.read_h5ad("../data/scanpy.paga.h5ad")
# sc.tl.draw_graph(adata, init_pos='paga')
# adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden']  == '2')[0]
# sc.tl.dpt(adata)
# sc.tl.umap(adata, init_pos='paga')
# adata.write_h5ad("../data/scanpy.graphed.h5ad")
