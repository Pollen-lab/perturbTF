#!/usr/bin/env python
# coding: utf-8

#script to process and generate figures using annotated scanpy object


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import anndata
import numpy as np
import os
import scanpy as sc
#import bbknn
warnings.simplefilter("ignore")

sc.set_figure_params(figsize=(4, 4))


#subset data for figure1
import scanpy as sc 
adata = sc.read('/wynton/group/pollen/jding/brainchromatin/Li/query_HM2Dall_guides.h5ad')
print(adata)

#define 'sgRNA_effecitve' as sgRNAs that are not non-targeting
feature_ref = df= pd.read_csv('/wynton/home/pollenlab/jding/BrainChromatin/Perturb/cellranger/feature_ref.csv',index_col=0)
feature_ref['guide_name'] = feature_ref['target_gene_name'] + '_' +  feature_ref['guideID']
dic = dict(map(lambda i,j : (i,j) , feature_ref.index, feature_ref['guide_name']))
adata.obs['sgRNA'] = [','.join(list(set([dic[i] if i in dic.keys() else i for i in x.split(',')]))) for x in adata.obs['sgRNA_NKS'].astype(str)]
adata.obs['num_sgRNA'] = [len(list(set(x.split(',')))) for x in adata.obs['sgRNA'].astype(str)]

adata.obs['sgRNA_effecitve'] = [','.join(list(set(i for i in x.split(',') if not i.lower().startswith('non-targeting')))) for x in adata.obs['sgRNA'].astype(str)]
adata.obs['sgRNA_effecitve'] = ['non-targeting' if x == '' else x for x in adata.obs['sgRNA_effecitve']]
adata.obs['num_sgRNA_effective'] = [len(list(set(x.split(',')))) for x in adata.obs['sgRNA_effecitve'].astype(str)]


adata.obs['sex'] = [x.split('-')[-1] for x in adata.obs['individual']]
adata.obs['stage'] = [x.split('#')[0] for x in adata.obs['individual']]
adata.obs['perturbation'] = ['NT' if x in ['WT','non-targeting'] else 'Perturbed' for x in adata.obs['gene_NKS']]




#adata = adata[adata.obs['individual'].isin(['GW18#140-M','GW16#172-F','GW18#151-F','GW16#130-F'])]
#adata = adata[adata.obs['timepoint'] == 'D0']
adata = adata[adata.obs['batch'].isin(['HMD7','HM2ndD7'])]
#adata = adata[adata.obs['species'] == 'human']


#preprocessing
sc.pp.normalize_total(adata, target_sum=1e4,exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
#sc.pp.highly_variable_genes(adata, n_top_genes=2500, batch_key="batch_name", subset= False)

from matplotlib import rcParams
rcParams["figure.figsize"] = (2, 2)
sc.pl.umap(adata, color=['AQP4','GLI3','FOXG1','DCX','SATB2','PPP1R17','OLIG2','TBR1','LHX6','PAX6','SCGN','CALB2',
                         ],use_raw=False,
           frameon=False, legend_fontsize=13, legend_fontoutline=3,ncols=4 ,save='12genes-HM2D')


sc.settings.file_format_figs='pdf'
# Inital setting for plot size
sc.pl.violin(adata,groupby='individual',keys=['XIST'], rotation=90, multi_panel=True,use_raw=False,save='HM2D_indivdual')


from matplotlib import rcParams
rcParams["figure.figsize"] = (2, 2)
sc.pl.umap(adata, color=['sex','stage','perturbation'],frameon=False, legend_fontsize=13, legend_fontoutline=3,
           palette='Set3',ncols=1,save='distribution-HM2D')

'''
from matplotlib import rcParams
rcParams["figure.figsize"] = (2, 2)
sc.pl.umap(adata, color=['AQP4','GLI3','DCX','SATB2','PPP1R17','OLIG2','TBR1','LHX6','PAX6','SCGN','ST18','CALB2',
                         'TSHZ1','FOXG1','NR2F2','SP8'],use_raw=False,
           frameon=False, legend_fontsize=13, legend_fontoutline=3,ncols=4 ,save='12genes-2D')


sc.pl.umap(adata, color=['AQP4','GLI3','DCX','SATB2','PPP1R17','OLIG2','TBR1','LHX6','PAX6','SCGN','ST18','CALB2'],use_raw=False,
           frameon=False, legend_fontsize=13, legend_fontoutline=3,ncols=4 ,save='12genes-2D')

adata.obs['class'] = ['Oligodendrocytes' if x == 'Glia' else x for x in adata.obs['class']]
sc.pl.umap(adata, color=['timepoint','perturbation','class','individual','sex'],frameon=False, legend_fontsize=13, legend_fontoutline=3,
           palette='Set3',ncols=1,save='distribution-2D')

'''
sc.settings.file_format_figs='svg'
sc.pl.umap(adata, color=['class'], 
           frameon=False, legend_fontsize=13, legend_fontoutline=3,save='HM2D_class')
sc.pl.umap(adata, color=['supervised_name'], 
           frameon=False, legend_fontsize=13, legend_fontoutline=3,save='HM2D_leiden')

sc.settings.file_format_figs='pdf'
# Inital setting for plot size
from matplotlib import rcParams
rcParams["figure.figsize"] = (2, 2)
sc.pl.umap(adata[adata.obs['species'] == 'macaque'], color=['supervised_name'],frameon=False, legend_fontsize=13, legend_fontoutline=3,legend_loc='none',
           ncols=1,save='celltype-HM2DQ')
sc.pl.umap(adata[adata.obs['species'] == 'human'], color=['supervised_name'],frameon=False, legend_fontsize=13, legend_fontoutline=3,legend_loc='none',
           ncols=1,save='celltype-HM2DH')

adata.obs['species'] = pd.Categorical(adata.obs['species'], 
                      categories=["macaque","human"],
                      ordered=True)
adata.obs['size'] = [120000 / 130000 if x == 'human' else 3* 120000 / 13000 for x in adata.obs['species']]
sc.pl.umap(adata, color=['species'],frameon=False, legend_fontsize=13, legend_fontoutline=3,size = adata.obs['size'],
           palette='Set3',ncols=1,save='species-HM2D')


adata.obs['supervised_name'] = adata.obs['supervised_name'].astype('category')
#sc.tl.dendrogram(adata, 'supervised_name')
dp = sc.pl.dotplot(adata, ['CDC20','HMGB2','GLIS3','VIM','EOMES','NEUROD6','EPHA3','MEF2C','ST18','PAX6','SCGN','BCAS1','EBF1','AIF1'], 
                   'supervised_name', categories_order = ['RG_DIV','RG-Astro','IPC_EN','EN_ImN','EN','IPC_IN','IN_dLGE-CGE','OPC-Oligo','Vasc','MG','Technical'],
                                 standard_scale='var',dot_max=0.6, return_fig=True)
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5,cmap='Reds' ).savefig("./figures/dotplot-HM2D.pdf")

sc.pl.violin(adata,groupby='individual',keys=['XIST'], rotation=90, multi_panel=True,use_raw=False,save='HM2D_indivdual')


adata.obs['sex'] = [x.split('-')[-1] for x in adata.obs['individual']]
adata.obs['stage'] = [x.split('#')[0] for x in adata.obs['individual']]
adata.obs['perturbation'] = ['NT' if x in ['WT','non-targeting'] else 'Perturbed' for x in adata.obs['gene_NKS']]


sc.pl.violin(adata,groupby='individual',keys=['XIST'], rotation=90, multi_panel=True,use_raw=False,save='HM2D_indivdual')

sc.settings.file_format_figs='svg'
sc.pl.umap(adata, color=['class'], 
           frameon=False, legend_fontsize=13, legend_fontoutline=3,save='HM2D_class')
sc.pl.umap(adata, color=['supervised_name'], 
           frameon=False, legend_fontsize=13, legend_fontoutline=3,save='HM2D_leiden')
sc.pl.umap(adata[adata.obs['species'] == 'human'], color=['supervised_name'], 
           frameon=False, legend_fontsize=13, legend_fontoutline=3,save='HM2DH_leiden')
