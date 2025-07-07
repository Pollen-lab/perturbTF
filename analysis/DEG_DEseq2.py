#!/usr/bin/env python
# coding: utf-8

# # Pseudo-bulk functional analysis




import scanpy as sc
import decoupler as dc

# Only needed for processing
import numpy as np
import pandas as pd

# Needed for some plotting
import matplotlib.pyplot as plt

# Plotting options, change to your liking
sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))


# ## Loading the data

'''
import scanpy as sc
import pandas as pd
import os
adata = sc.read('/wynton/group/pollen/jding/brainchromatin/HM2D/HM2D-filtered.h5ad', compression='gzip')
adata.obs['Gene_target'] = adata.obs['gene_NKS'] 
targets = [x for x in adata.obs['Gene_target'].unique() if x in adata.var_names]
targets.append('non-targeting')
adata = adata[adata.obs['Gene_target'].isin(targets)]
adata.obs['sex'] = [x.split('-')[-1] for x in adata.obs['individual']]
adata.obs['stage'] = [x.split('#')[0] for x in adata.obs['individual']]
adata.obs['perturbation'] = ['NT' if x in ['WT','non-targeting'] else 'Perturbed' for x in adata.obs['gene_NKS']]



## define species, cluster, output directory to use for this analysis

species = 'human'
groups_col = 'supervised_name'

dir='/wynton/group/pollen/jding/brainchromatin/HM2D/DEseq2'
if not os.path.exists(os.path.join(dir, species + '-class')):
    os.makedirs(os.path.join(dir, species + '-class'))    
dir = os.path.join(dir, species + '-class')


old_to_new = {
    'IN_dLGE-CGE':'InhibitoryNeurons',
    'RG_DIV':'RadialGlia',
    'RG-Astro':'RadialGlia',
    'IPC_IN':'InhibitoryNeurons',
    'EN_ImN':'ExcitatoryNeurons',
    'EN':'ExcitatoryNeurons',
    'IPC_EN':'ExcitatoryNeurons',
    'OPC-Oligo':'Glia',
    'Technical':'Technical',
    'MG':'Microglia',
    'Vasc':'Vascular'
}
adata.obs['class'] = (
    adata.obs['supervised_name']
    .map(old_to_new)
    .astype('category')
)
'''

import os
import pandas as pd
dir='/wynton/group/pollen/jding/brainchromatin/HM2D/DEseq2/2D_guides_bulk'
df = pd.read_csv(os.path.join(dir, 'min_log2FC.csv'),index_col=0)

#load adata
adata = sc.read('/wynton/group/pollen/jding/brainchromatin/HM2D/D7-filtered_guides.h5ad', compression='gzip')
adata.obs['sex'] = [x.split('-')[-1] for x in adata.obs['individual']]
adata.obs['stage'] = [x.split('#')[0] for x in adata.obs['individual']]

#define 'sgRNA_effecitve' as sgRNAs that are not non-targeting
adata.obs['sgRNA_effecitve'] = [','.join(list(set(i for i in x.split(',') if not i.startswith('Non-Targeting')))) for x in adata.obs['sgRNA'].astype(str)]
adata.obs['num_sgRNA_effective'] = [len(list(set(x.split(',')))) for x in adata.obs['sgRNA_effecitve'].astype(str)]
adata.obs['sgRNA_effecitve'] = ['non-targeting' if x == '' else x for x in adata.obs['sgRNA_effecitve']]
adata = adata[adata.obs['num_sgRNA_effective'] == 1]

#removev guides that had <25% KD (log2FC > -0.4)
adata = adata[~adata.obs['sgRNA_effecitve'].isin(df[df['log2FC'] > -0.4]['guide'].to_list())]

adata.obs['Gene_target'] = adata.obs['gene_NKS'] 
targets = [x for x in adata.obs['Gene_target'].unique() if x in adata.var_names]
targets.append('non-targeting')
adata = adata[adata.obs['Gene_target'].isin(targets)]
adata.obs['perturbation'] = ['NT' if x=='non-targeting' else 'Perturbed' for x in adata.obs['Gene_target'] ]


species = 'human'
groups_col = 'supervised_name'

dir='/wynton/group/pollen/jding/brainchromatin/HM2D/DEseq2/2D_guides_cleaned_min5'



# ## Processing
# 

adata = adata[adata.obs['species'] == species]


# Store raw counts in layers
#adata.X = np.round(adata.X)
adata.raw = adata.copy()
adata.layers['counts'] = adata.X.copy()


# ## Generation of pseudo-bulk profiles

adata.obs['individual-Gene_target'] = adata.obs['individual'].astype(str) + adata.obs['Gene_target'].astype(str)


# In[ ]:


adata.obs['supervised_name'].value_counts()



# In[14]:


# Get pseudo-bulk profile
#remove conditions that have less than 5 cells
pdata = dc.get_pseudobulk(
    adata,
    sample_col='individual-Gene_target',
    groups_col= groups_col,
    layer='counts',
    mode='sum',
    min_cells=5,
    min_counts=0
)

'''
for x in pdata.obs.columns.tolist():
    pdata.obs[x]=pdata.obs[x].astype(str)
for x in pdata.var.columns.tolist():
    pdata.var[x]=pdata.var[x].astype(str)

pdata.write(os.path.join(dir,'pdata.h5ad'))
'''


# ### Exploration of pseudobulk profiles


# Store raw counts in layers
pdata.layers['counts'] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

# Return raw counts to X
dc.swap_layer(pdata, 'counts', X_layer_key=None, inplace=True)


# In[9]:


sc.pl.pca(pdata, color=['individual', groups_col], ncols=1, size=300)
sc.pl.pca_variance_ratio(pdata)




# Import DESeq2
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats


targets=[x for x in adata.obs['Gene_target'].unique() if x in adata.var_names]
inference = DefaultInference(n_cpus=1)

for celltype in pdata.obs[groups_col].unique():
    print(celltype)
    sub = pdata[pdata.obs[groups_col] == celltype]
    for gene in targets:
        print(gene)
        if gene in sub.obs['Gene_target'].unique():
            if (pd.crosstab(sub.obs[groups_col],sub.obs['Gene_target'])[gene] > 1).bool(): #remove conditions that have less than 2 biological replicates
                sub2 = sub[sub.obs['Gene_target'].isin([gene,'non-targeting'])]
                #print(sub2.obs.head(5))
                try:
                    dds = DeseqDataSet(adata= sub2, design_factors='perturbation',ref_level=['perturbation', 'NT'],
                                       refit_cooks=True,inference=inference)
                except KeyError:
                    dds = DeseqDataSet(adata= sub2, design_factors='perturbation',ref_level=['perturbation', 'NT'],
                                   refit_cooks=True,inference=inference)
                dds.deseq2()
                stat_res = DeseqStats(dds,contrast=["perturbation", 'Perturbed', 'NT'],inference=inference)
                stat_res.summary()
                stat_res.results_df.to_csv(os.path.join(dir,gene + '.' + celltype +'.csv'))
        else:
            continue
