import numpy as np
import pandas as pd
import scanpy as sc

def scRNA_recipe_czz(adata,norm=None,
                     doublet_detection=None,
                     log=None,
                     scale=None,
                     pca=None,
                     harmony=None,
                     neighbors=None,
                     leiden=None
                     help=False
                     ):
    if help:
        
        return
    sc.settings.verbosity = 3
    if norm is None:
        sc.pp.normalize_total(adata)
    else:
        sc.pp.normalize_total(adata, **norm)
    return