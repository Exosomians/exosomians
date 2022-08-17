import scanpy as sc
sc.settings.set_figure_params(dpi=300)

def plot_umap(data, labels=None, show=False, save=None):
    adata = sc.AnnData(data)
    adata.obs['label'] = labels

    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='label', show=show, frameon=False, save=save)
