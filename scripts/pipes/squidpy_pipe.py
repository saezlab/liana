import pandas as pd
import scanpy as sc
import squidpy as sp
import re


def get_omni_resources(op_resource, adata_seurat):
    """ Helper function to map over each resource and call the CellPhoneDB implementation of SquidPy
     Parameters
        ----------
        op_resource
            OmniPath database to use in LR inference
        adata_seurat
            Anndata/ScanPy clustering object
        Returns
            A LigRec Results Object"""
    try:
        res_filt = None if op_resource == "OmniPath" else op_resource

        y = sp.gr.ligrec(adata_seurat, "seurat_annotations",
                         fdr_method=None, copy=True,
                         interactions_params={"resources": res_filt},
                         transmitter_params={"categories": "ligand", "resources": res_filt},
                         receiver_params={"categories": "receptor", "resources": res_filt},
                         threshold=0.1, seed=1004, n_perms=10000, n_jobs=1)
        return y
    except ValueError:
        print(op_resource)
        return None


def get_ligrec(intercell_resources, adata_seurat):
    return dict(map(lambda x: (x, get_omni_resources(x, adata_seurat)), intercell_resources))


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
                 embedding):
    """Call Squidpy
         Parameters
        ----------
        intercell_resources,

        exprs
            Expression Matrix
        meta
            Clustering metadata
        feature_meta
            Feature metadata (e.g. vst parameters from Seurat)
        embedding
            Dimensionality reduction data (e.g. UMAP, tSNE)

        Returns
            Two lists: One with LR interaction pvalue results for each resource, and one with means.
    """
    intercell_resources = list(intercell_resources)
    adata_seurat = convert_anndata(exprs, meta, feature_meta, embedding)
    # call squidpy
    squidpy_res = get_ligrec(intercell_resources, adata_seurat)

    squidpy_pvalues = list(map(lambda x: reformat(squidpy_res[x].pvalues), intercell_resources))
    squidpy_means = list(map(lambda x: reformat(squidpy_res[x].means), intercell_resources))

    return {"pvalues": squidpy_pvalues, "means": squidpy_means}
