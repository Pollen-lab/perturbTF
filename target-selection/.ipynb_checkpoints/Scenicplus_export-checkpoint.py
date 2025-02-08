#!/usr/bin/env python
# coding: utf-8



import warnings
import os
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

import dill
work_dir = sys.argv[1]
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))


from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
apply_std_filtering_to_eRegulons(scplus_obj)


from scenicplus.eregulon_enrichment import score_eRegulons
region_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/region_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
gene_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/gene_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 12)
score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 12)
from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks

generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'GEX_supervised_name',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Gene_based')
generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'GEX_supervised_name',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Region_based')

TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'GEX_supervised_name',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Gene_based',
            out_key = 'filtered_gene_based')
TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'GEX_supervised_name',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Region_based',
            out_key = 'filtered_region_based')



from matplotlib import pyplot as plt
import numpy as np
n_targets = [int(x.split('(')[1].replace('r)', '')) for x in scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
rho = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
adj_pval = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

thresholds = {
        'rho': [-0.50, 0.50],
        'n_targets': 0
}
import seaborn as sns
fig, ax = plt.subplots(figsize = (10, 5))
sc = ax.scatter(rho, n_targets, c = -np.log10(adj_pval), s = 5)
ax.set_xlabel('Correlation coefficient')
ax.set_ylabel('nr. target regions')
#ax.hlines(y = thresholds['n_targets'], xmin = min(rho), xmax = max(rho), color = 'black', ls = 'dashed', lw = 1)
ax.vlines(x = thresholds['rho'], ymin = 0, ymax = max(n_targets), color = 'black', ls = 'dashed', lw = 1)
ax.text(x = thresholds['rho'][0], y = max(n_targets), s = str(thresholds['rho'][0]))
ax.text(x = thresholds['rho'][1], y = max(n_targets), s = str(thresholds['rho'][1]))
sns.despine(ax = ax)
fig.colorbar(sc, label = '-log10(adjusted_pvalue)', ax = ax)
plt.show()

selected_cistromes = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].loc[
        np.logical_or(
                scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] > thresholds['rho'][1],
                scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] < thresholds['rho'][0]
        )]['Cistrome'].to_list()
selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
selected_eRegulons_gene_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()
        if x.split('_(')[0] in selected_eRegulons]
selected_eRegulons_region_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()
        if x.split('_(')[0] in selected_eRegulons]
#save the results in the scenicplus object
scplus_obj.uns['selected_eRegulon'] = {'Gene_based': selected_eRegulons_gene_sig, 'Region_based': selected_eRegulons_region_sig}
print(f'selected: {len(selected_eRegulons_gene_sig)} eRegulons')


#hvg/hvr
from pycisTopic.diff_features import find_highly_variable_features

hvg = find_highly_variable_features(scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Gene']))].T, n_top_features=1000, plot = False)
file = open(os.path.join(work_dir, 'scenicplus/hvg-1000.txt'),'w')
for item in hvg:
    file.write(item+"\n")
file.close()

hvg = find_highly_variable_features(scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Gene']))].T, n_top_features=1500, plot = False)
file = open(os.path.join(work_dir, 'scenicplus/hvg-1500.txt'),'w')
for item in hvg:
    file.write(item+"\n")
file.close()

hvr = find_highly_variable_features(scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Region']))], n_top_features=1000, plot = False)
file = open(os.path.join(work_dir, 'scenicplus/hvr-1000.txt'),'w')
for item in hvr:
    file.write(item+"\n")
file.close()

hvr = find_highly_variable_features(scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Region']))], n_top_features=1500, plot = False)
file = open(os.path.join(work_dir, 'scenicplus/hvr-1500.txt'),'w')
for item in hvr:
    file.write(item+"\n")
file.close()




#gene-based
scplus_obj.uns['TF2G_adj'].to_csv(os.path.join(work_dir, 'scenicplus/TF2G_adj_all.csv'))
import pandas as pd
values = [x.split('_')[-1][1:-2] for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()]
directions = [x.split('_')[-2] for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()]
keys = [x.split('_')[0] for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()]
TFs = [x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()]
map = pd.DataFrame({'TF': TFs, 'GRNs': keys, 'regions':values,'directions':directions})
map.to_csv(os.path.join(work_dir, 'scenicplus/GRN_signatures_gene.csv'))
map.head()

if not os.path.exists(os.path.join(work_dir, 'scenicplus/GRNs_gene')):
    os.makedirs(os.path.join(work_dir, 'scenicplus/GRNs_gene'))
for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based']:
    name= x.replace('(','')
    name= name.replace(')','')
    file = open(os.path.join(work_dir, 'scenicplus/GRNs_gene',name+'.txt'),'w')
    for item in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'][x]:
        file.write(item+"\n")
    file.close()


#region-based
import pandas as pd
values = [x.split('_')[-1][1:-2] for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()]
directions = [x.split('_')[-2] for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()]
keys = [x.split('_')[0] for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()]
TFs = [x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()]
map = pd.DataFrame({'TF': TFs, 'GRNs': keys, 'regions':values,'directions':directions})
map.to_csv(os.path.join(work_dir, 'scenicplus/GRN_signatures_region.csv'))
map.head()

if not os.path.exists(os.path.join(work_dir, 'scenicplus/GRNs_region')):
    os.makedirs(os.path.join(work_dir, 'scenicplus/GRNs_region'))
for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based']:
    name= x.replace('(','')
    name= name.replace(')','')
    file = open(os.path.join(work_dir, 'scenicplus/GRNs_region',name+'.txt'),'w')
    for item in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'][x]:
        file.write(item+"\n")
    file.close()


#TF to gene
keys = [x.split('_')[0] for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys() if x.split('_')[1] == '+']
df = scplus_obj.uns['TF2G_adj']
df = df[df['TF'].isin(keys)]
df = df[df['target'].isin(keys)]
df.to_csv(os.path.join(work_dir, 'scenicplus/TF2G_adj.csv'))
df.head()

scplus_obj.uns['TF2G_adj'].to_csv(os.path.join(work_dir, 'scenicplus/TF2G_adj_all.csv'))
scplus_obj.uns['region_to_gene'].to_csv(os.path.join(work_dir, 'scenicplus/R2G_adj_all.csv'))

#metadata
scplus_obj.uns['eRegulon_metadata_filtered'].to_csv(os.path.join(work_dir, 'scenicplus/eRegulon_metadata_filtered.csv'))


#RSS
from scenicplus.RSS import *
regulon_specificity_scores(
        scplus_obj,
        variable = 'GEX_supervised_name',
        auc_key = 'eRegulon_AUC_filtered',
        signature_keys = ['Region_based'],
        selected_regulons = [x for x in scplus_obj.uns['eRegulon_metadata_filtered']['Region_signature_name'].unique()],
        out_key_suffix = '_filtered')

scplus_obj.uns['RSS']['GEX_supervised_name_filtered'].to_csv(os.path.join(work_dir, 'scenicplus/RSS_region.csv'))


from scenicplus.RSS import *
regulon_specificity_scores(
        scplus_obj,
        variable = 'GEX_supervised_name',
        auc_key = 'eRegulon_AUC_filtered',
        signature_keys = ['Gene_based'],
        selected_regulons = [x for x in scplus_obj.uns['eRegulon_metadata_filtered']['Gene_signature_name'].unique()],
        out_key_suffix = '_filtered')

scplus_obj.uns['RSS']['GEX_supervised_name_filtered'].to_csv(os.path.join(work_dir, 'scenicplus/RSS_gene.csv'))


#repressors
import pandas as pd
import sklearn

#TF expression
dgem = pd.DataFrame(scplus_obj.X_EXP, index=scplus_obj.cell_names,
                        columns=scplus_obj.gene_names).copy()
dgem = dgem.T

#region accessibility    
auc_key= 'eRegulon_AUC_filtered'
signature_keys= ['Region_based']
#scaled
region_mat = pd.concat([pd.DataFrame(sklearn.preprocessing.StandardScaler().fit_transform(
            scplus_obj.uns[auc_key][x].T), index=scplus_obj.uns[auc_key][x].T.index.to_list(), columns=scplus_obj.uns[auc_key][x].T.columns) for x in signature_keys])
selected_regulons = [x for x in region_mat.index if x.split('_')[0] in dgem.index]
#subset to only repressors
selected_regulons = [x for x in selected_regulons if x.split('_')[-2] == '-']
region_mat = region_mat.loc[region_mat.index.isin(selected_regulons)]

#gene expression 
auc_key= 'eRegulon_AUC_filtered'
signature_keys= ['Gene_based']
#scaled
exp_mat = pd.concat([pd.DataFrame(sklearn.preprocessing.StandardScaler().fit_transform(
            scplus_obj.uns[auc_key][x].T), index=scplus_obj.uns[auc_key][x].T.index.to_list(), columns=scplus_obj.uns[auc_key][x].T.columns) for x in signature_keys])
selected_regulons = [x for x in exp_mat.index if x.split('_')[0] in dgem.index]
#subset to only repressors
selected_regulons = [x for x in selected_regulons if x.split('_')[-2] == '-']
exp_mat = exp_mat.loc[exp_mat.index.isin(selected_regulons)]

selected_genes = [x.split('_')[0] for x in region_mat.index]
dgem = dgem.loc[dgem.index.isin(selected_genes)]

#create anndata
import anndata 
adata = anndata.AnnData(layers= {'exp': dgem.sort_index(ascending=True).T,'access': region_mat.sort_index(ascending=True).T,'target': exp_mat.sort_index(ascending=True).T},
                        obs= exp_mat.columns.tolist(),
                        var=[x.rsplit('_', 1)[0] for x in exp_mat.sort_index(ascending=True).index])


adata.obs['Cellname'] = adata.obs[0]
adata.obs.index = adata.obs['Cellname']
del adata.obs[0]

adata.var['Gene'] = adata.var[0]
adata.var.index = [x.split('_')[0] for x in adata.var['Gene']]
del adata.var[0]

#save raw data
adata.X = adata.layers['exp']
adata.raw = adata
dgem = dgem.T / dgem.T.sum(0) * 10**6 #skip normalization
dgem = np.log1p(dgem).T
adata.X = dgem.T

adata.write(os.path.join(work_dir, 'scenicplus/TF_rep.h5ad'))