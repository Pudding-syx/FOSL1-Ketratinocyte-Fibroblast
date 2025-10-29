import numpy as np
import pandas as pd
import scanpy as sc

import scipy.stats as stats

import os

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
    current_path = os.path.dirname(__file__)
    return pd.read_csv(current_path+"/../_data/_"+database+".py",index_col=None,header=0,sep="\t")

def cell_interaction(adata, group_by, database, method="fast", all_or_none=False,interaction_pairs=None, **kwargs):
    methods = ['cellphone', "fast"]
    if method not in methods:
        print("ERROR:considering method in :", end="")
        print(methods)
        return
    current_path = os.path.dirname(__file__)
    var_names = list(adata.var_names)
    if database == "CellChatDB_human":
        db = pd.read_csv(current_path+"/../_data/_"+database+".py",index_col=None,header=0,sep="\t")
        ligand_list = db['ligand.symbol'].to_numpy()
        receptor_list = db['receptor.symbol'].to_numpy()
        index = np.zeros(db.shape[0])
        for i in range(db.shape[0]):
            a = 1
            tmp = ligand_list[i].split(", ")
            for gene in tmp:
                if gene not in var_names:
                    a = 0
            tmp = receptor_list[i].split(", ")
            for gene in tmp:
                if gene not in var_names:
                    a = 0
            index[i] = a
        ligand_list = ligand_list[index > 0]
        receptor_list = receptor_list[index > 0]
        genes = []
        for i in ligand_list:
            genes += i.split(", ")
        for i in receptor_list:
            genes += i.split(", ")
        genes = np.unique(genes)
        index = [i in adata.var_names for i in genes]
        genes = genes[index]
        adata = adata[:, genes].copy()

    if all_or_none:
        adata.X[adata.X>0.] = 1.
        adata.X[adata.X<0.] = 0.
    if method == "fast":
        genes = np.array(list(adata.var_names))

        clusters_names = [str(i) for i in np.unique(adata.obs[group_by])]
        if interaction_pairs is None:
            inter_names = []
            for l in clusters_names:
                inter_names += [l + "|" + r for r in clusters_names]
        else:
            inter_names = interaction_pairs

        ligand_strength = np.zeros((len(clusters_names), len(ligand_list)), float)
        receptor_strength = np.zeros((len(clusters_names), len(receptor_list)), float)
        cluster_expr = np.zeros((len(clusters_names), adata.shape[1]), float)

        tmp = np.array([str(i) for i in adata.obs[group_by].to_numpy().flatten()])
        for i, c in enumerate(clusters_names):
            cluster_expr[i, :] = np.array(adata[tmp == c].X.mean(0))[0, :]
        for i in range(ligand_strength.shape[1]):
            tmp = ligand_list[i].split(", ")
            index = [i in tmp for i in genes]
            ligand_strength[:, i] = cluster_expr[:, index].mean(1)
        for i in range(receptor_strength.shape[1]):
            tmp = receptor_list[i].split(", ")
            index = [i in tmp for i in genes]
            receptor_strength[:, i] = cluster_expr[:, index].min(1)

        inter_strength = np.zeros((len(inter_names), len(ligand_list)), float)
        for i in range(inter_strength.shape[0]):
            l = inter_names[i].split("|")[0]
            l = [c == l for c in clusters_names]
            l = ligand_strength[l, :]
            r = inter_names[i].split("|")[1]
            r = [c == r for c in clusters_names]
            r = receptor_strength[r, :]
            inter_strength[i, :] = (l * r)[0, :]
        index = inter_strength.mean(0) > 0
        ligand_list = ligand_list[index]
        receptor_list = receptor_list[index]
        inter_strength = inter_strength[:, index]
        inter_strength = np.sqrt(inter_strength)
        fold_changes = inter_strength / inter_strength.mean(0)

        result = pd.DataFrame()
        indexes = np.argsort(fold_changes, axis=1)[:, ::-1]
        for i, name in enumerate(inter_names):
            index = indexes[i, :]
            result[name + "_ligand"] = ligand_list[index]
            result[name + "_receptor"] = receptor_list[index]
            result[name + "_strength"] = inter_strength[i, index]
            result[name + "_foldchange"] = fold_changes[i, index]
        return result

    if method == "cellphone":
        ligand_strength = np.zeros((adata.shape[0], len(ligand_list)), float)
        receptor_strength = np.zeros((adata.shape[0], len(receptor_list)), float)
        for i in range(ligand_strength.shape[1]):
            print(i, end="\r")
            tmp = ligand_list[i].split(", ")
            ligand_strength[:, i] = np.array(adata[:, tmp].X.mean(1))[:, 0]
        for i in range(receptor_strength.shape[1]):
            print(i, end="\r")
            tmp = receptor_list[i].split(", ")
            receptor_strength[:, i] = np.array(adata[:, tmp].X.min(1).todense())[:, 0]
        return

def annotate_cell(adata,key_column,copy=False, **kwargs):
    if copy:
        adata = adata.copy()
    for new_column in kwargs.keys():
        adata.obs[new_column] = "Others"
        for annotation in kwargs[new_column].keys():
            index = [i in kwargs[new_column][annotation] for i in adata.obs[key_column]]
            adata.obs[new_column][index] = annotation
    if copy:
        return adata
    else:
        return


def hypergeometric_test(setA, setB, M=None):
    """
    Accepts to lists
    M is the population size (previously N)
    n is the number of successes in the population
    N is the sample size (previously n)
    x is still the number of drawn “successes”
    """

    if M is None:
        M = 24000
    n = len(setA)
    N = len(setB)
    x = len(setA.intersection(setB))

    return stats.hypergeom.cdf(x, M, n, N), stats.hypergeom.sf(x - 1, M, n, N)