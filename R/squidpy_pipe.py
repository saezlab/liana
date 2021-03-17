import pandas as pd
import scanpy as sc
import squidpy as sp
import re


def get_squidpy_res(op_resource, adata_seurat, ident, **kwargs):
    """ Helper function to map over each resource and call the CellPhoneDB implementation of SquidPy
     Parameters
        ----------
        op_resource
            OmniPath database resource to use in LR inference
        adata_seurat
            Anndata/ScanPy clustering object
        ident
            Cluster identity column
        Returns
            A LigRec Results Object"""
    try:
        res = sp.gr.ligrec(
            adata_seurat,
            ident,
            fdr_method=None, copy=True,
            interactions=op_resource,
            corr_method = "fdr_bh",
            threshold=0.1, seed=1004,
            n_perms=10000, n_jobs=1 
            ) # should replace with kwargs and elipses
        return res
    except ValueError as e:
        print(e)
        print(op_resource)
        return None


def get_ligrec(intercell_resources, adata_seurat, ident):
    return list(map(lambda x: (x, get_squidpy_res(x, adata_seurat, ident)), intercell_resources))


def reformat(x):
    """Reformat SquidPy LigRec results to non-multiindex DFs"""
    # flatten multiindex columns (i.e. cell types)
    x_reform = pd.DataFrame(x.to_records())
    # fix column formatting
    x_reform = x_reform. \
        rename(columns=lambda x: re.sub(r', ', '_', x)). \
        rename(columns=lambda x: re.sub(r"'", '', x))
    x_reform.columns = list(map(lambda x: x[1:-1] if x.startswith(r"(") else x, x_reform.columns))
    return x_reform


def convert_anndata(exprs,
                    meta,
                    feature_meta,
                    embedding):
    """Convert R Seurat to Python Anndata
     Parameters
        ----------
        exprs
            Expression Matrix
        meta
            Clustering metadata
        feature_meta
            Feature metadata (e.g. vst parameters from Seurat)
        embedding
            Dimensionality reduction data (e.g. UMAP, tSNE)

        Returns
            A converted AnnData Object
    """
    adata_seurat = sc.AnnData(X=exprs.T, obs=meta, var=feature_meta)
    adata_seurat.obsm['umap'] = embedding
    adata_seurat.raw = adata_seurat.copy()
    return adata_seurat


def call_squidpy(intercell_resources,
                 exprs,
                 meta,
                 feature_meta,
                 embedding,
                 ident):
    """Call Squidpy
         Parameters
        ----------
        intercell_resources
            List with OmniPath resources

        exprs
            Expression Matrix
        meta
            Clustering metadata
        feature_meta
            Feature metadata (e.g. vst parameters from Seurat)
        embedding
            Dimensionality reduction data (e.g. UMAP, tSNE)
        ident
            Cluster identity column

        Returns
            Two lists: One with LR interaction pvalue results for each resource, and one with means.
    """
    
    adata_seurat = convert_anndata(exprs, meta, feature_meta, embedding)
    # call squidpy
    squidpy_res = get_ligrec(intercell_resources, adata_seurat, ident)
    
    # reformat results
    squidpy_pvalues = list(map(lambda x: reformat(x[1]["pvalues"]), squidpy_res))
    squidpy_means = list(map(lambda x: reformat(x[1]["means"]), squidpy_res))

    return {"pvalues": squidpy_pvalues, "means": squidpy_means}
