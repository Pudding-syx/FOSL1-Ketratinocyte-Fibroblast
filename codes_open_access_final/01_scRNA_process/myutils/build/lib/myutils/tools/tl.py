import numpy as np
import pandas as pd
import scanpy as sc

def calculate_cell_fraction(adata,group_by,sample_key=None,condition_key=None,sample_meta_key=None,melt=False,melt_key=None):
    if sample_key is not None:
        fractions = np.zeros((len(np.unique(adata.obs[sample_key])), len(np.unique(adata.obs[group_by]))))
        cells = np.zeros((len(np.unique(adata.obs[sample_key])), len(np.unique(adata.obs[group_by]))))
        for i, p in enumerate(np.unique(adata.obs[sample_key])):
            tmp = adata.obs.iloc[np.where(adata.obs[sample_key] == p)[0],:]
            n_cells = tmp.shape[0]
            for j, c in enumerate(np.unique(adata.obs[group_by])):
                fractions[i,j] = tmp.iloc[np.where(tmp[group_by] == c)[0],:].shape[0]/n_cells
                cells[i,j] = tmp.iloc[np.where(tmp[group_by] == c)[0],:].shape[0]
        patient_meta = pd.DataFrame(fractions)
        patient_meta.columns = np.unique(adata.obs[group_by])
        patient_meta[sample_key] = np.unique(adata.obs[sample_key])
        if condition_key is not None:
            patient_meta[condition_key] = [adata.obs.iloc[np.where(adata.obs[sample_key] == p)[0],:][condition_key].iloc[0] for p in patient_meta[sample_key]]
        if sample_meta_key is not None:
            for name in sample_meta_key:
                patient_meta[name] = [adata.obs.iloc[np.where(adata.obs[sample_key] == p)[0],:][name].iloc[0] for p in patient_meta[sample_key]]
        patient_meta['n_cells'] = [adata.obs.iloc[np.where(adata.obs[sample_key] == p)[0],:].shape[0] for p in patient_meta[sample_key]]
    else:
        fractions = np.zeros((1, len(np.unique(adata.obs[group_by]))))
        cells = np.zeros((1, len(np.unique(adata.obs[group_by]))))
        tmp = adata.obs
        n_cells = tmp.shape[0]
        for j, c in enumerate(np.unique(adata.obs[group_by])):
            fractions[0, j] = tmp.iloc[np.where(tmp[group_by] == c)[0], :].shape[0] / n_cells
            cells[0, j] = tmp.iloc[np.where(tmp[group_by] == c)[0], :].shape[0]
            patient_meta = pd.DataFrame(fractions)
            patient_meta.columns = np.unique(adata.obs[group_by])
            patient_meta['n_cells'] = [adata.shape[0]]
    if melt:
        if melt_key is not None:
            patient_meta = pd.melt(tmp, id_vars=melt_key)
        else:
            patient_meta = pd.melt(tmp, id_vars=list(sample_key)+list(condition_key)+["n_cells"])
    return patient_meta

def scanpy_ranking_genes(adata,group_by,save=None,calculate_average_expression=False,return_result=True,**kwargs):
    try:
        tmp = adata.uns['log1p']["base"]
        sc.tl.rank_genes_groups(adata=adata, groupby=group_by, **kwargs)
    except:
        adata.uns['log1p']["base"] = None
        sc.tl.rank_genes_groups(adata=adata, groupby=group_by, **kwargs)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    cell_markers = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'scores', 'pvals_adj', "logfoldchanges"]})
    if calculate_average_expression:
        for group in groups:
            tmp = adata[adata.obs[group_by] == group].copy()
            tmp = tmp[:,list(cell_markers[group + "_names"])]
            tmp = tmp.X.mean(0).T
            cell_markers.insert(np.where(cell_markers.columns == group+"_logfoldchanges")[0][0]+1, group + "_expr", tmp)
    if save is not None:
        if save.endswith(".tsv"):
            cell_markers.to_csv(save, sep="\t")
        elif save.endswith(".xlxs"):
            cell_markers.to_excel(save)
        elif save.endswith(".json"):
            cell_markers.to_json(save)
        elif save.endswith(".html"):
            cell_markers.to_html(save)
        elif save.endswith(".pkl"):
            cell_markers.to_pickle(save)
        else:
            cell_markers.to_csv(save)
    if return_result:
        return cell_markers
    else:
        return

def drawing_marker_genes(database,**kwargs):
    return pd.read_csv("../data/"+database+".tsv",index_col=None,header=0)