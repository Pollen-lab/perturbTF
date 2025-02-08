#!/usr/bin/env python
# coding: utf-8

# # Tutorial: Perturbation simulation
# 
# In this tutorial we illustrate how the predictions from SCENIC+ can be utilized to simulate the effect of transcription factor perturbations.
# 
# Here, the predictions of SCENIC+ serve as a feature selection method. We will use the expression of transcription factors (TFs) as predictors for their target gene expression. 
# For this a random forest regression model will be fitted for each gene with the expression of TFs which are predicted to regulate them by SCENIC+ as predictor for their target gene expression.
# After fitting the models we can alter the expression of a TF of choice and we can simulate a new gene expression matrix. This simulation is repeated for several iterations to simulate indirect effects.
# The simulated cells in this new matrix can be projected in an embedding of choice to visualize the effect of the perturbation.
# 
# For this tutorial we will continue with the SCENIC+ analysis which he have done in the "Mix of melanoma cell lines" tutorial.

# <div class="alert alert-warning">
# 
# **Warning:**
# 
# In order to continue you need to have the python package `velocyto` installed. 
# </div>

# In[ ]:


#%pip install velocyto


# ## Setting up work environment

# In[1]:


#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
work_dir = '/wynton/group/pollen/jding/brainchromatin/Li/scenicplus/seacells'
tmp_dir = '/wynton/scratch/jding/ray'

if not os.path.exists(work_dir):
    os.makedirs(work_dir)


# Load the scenicplus object from previous tutorial.

# In[2]:


import dill
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
        signature_key = 'Gene_based',
        nr_cells=6)
generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'GEX_supervised_name',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Region_based',
        nr_cells=6)

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


thresholds = {
        'rho': [-0.65, 0.60],
        'n_targets': 0
}

import numpy as np
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



# Let's calculate a target gene based PCA embedding to plot the perturbations on

# In[3]:


from scenicplus.dimensionality_reduction import run_eRegulons_pca
run_eRegulons_pca(
        scplus_obj,
        auc_key = 'eRegulon_AUC_filtered',
        reduction_name = 'eRegulons_PCA_gene_based',
        selected_regulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
        n_pcs = 40)


# ## Plotting perturbation effect on an embedding
# 
# Let's simulate the effect of SOX10 knockdown. In this example of the melanoma cell lines it is known, from previous studies it is known that this perturbation can cause a phenotype swith of melanocytic states towards a more mesenchymal like state.
# 
# For the sake of computational time we will only simulate the expression of the top 200 highly variable genes.
# 
# You might opt to select more or all features, depending on your analysis

# In[4]:


from pycisTopic.diff_features import find_highly_variable_features
hvg = find_highly_variable_features(scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Gene']))].T, n_top_features = 200, plot = False)



from scenicplus.simulation import plot_perturbation_effect_in_embedding
import seaborn as sns
plot_perturbation_effect_in_embedding(
        scplus_obj = scplus_obj,
        reduction_name = 'eRegulons_PCA_gene_based',
        n_cpu = 5,
        perturbation = {'NR2E1': 0}, #specifies that we want to set the expression of SOX10 to 0 in all cells.
        variable = 'GEX_supervised_name',
        genes_to_use = hvg,
        figsize = (4, 4),
        save = '/wynton/home/pollenlab/jding/BrainChromatin/Li/figures/NR2E1.pdf')

plot_perturbation_effect_in_embedding(
        scplus_obj = scplus_obj,
        reduction_name = 'eRegulons_PCA_gene_based',
        n_cpu = 5,
        perturbation = {'ARX': 0}, #specifies that we want to set the expression of SOX10 to 0 in all cells.
        variable = 'GEX_supervised_name',
        genes_to_use = hvg,
        figsize = (4, 4),
        save = '/wynton/home/pollenlab/jding/BrainChromatin/Li/figures/ARX.pdf')

plot_perturbation_effect_in_embedding(
        scplus_obj = scplus_obj,
        reduction_name = 'eRegulons_PCA_gene_based',
        n_cpu = 5,
        perturbation = {'SOX2': 0}, #specifies that we want to set the expression of SOX10 to 0 in all cells.
        variable = 'GEX_supervised_name',
        genes_to_use = hvg,
        figsize = (4, 4),
        save = '/wynton/home/pollenlab/jding/BrainChromatin/Li/figures/SOX2.pdf')

plot_perturbation_effect_in_embedding(
        scplus_obj = scplus_obj,
        reduction_name = 'eRegulons_PCA_gene_based',
        n_cpu = 5,
        perturbation = {'ZNF219': 0}, #specifies that we want to set the expression of SOX10 to 0 in all cells.
        variable = 'GEX_supervised_name',
        genes_to_use = hvg,
        figsize = (4, 4),
        save = '/wynton/home/pollenlab/jding/BrainChromatin/Li/figures/ZNF219.pdf')

plot_perturbation_effect_in_embedding(
        scplus_obj = scplus_obj,
        reduction_name = 'eRegulons_PCA_gene_based',
        n_cpu = 5,
        perturbation = {'NEUROD2': 0}, #specifies that we want to set the expression of SOX10 to 0 in all cells.
        variable = 'GEX_supervised_name',
        genes_to_use = hvg,
        figsize = (4, 4),
        save = '/wynton/home/pollenlab/jding/BrainChromatin/Li/figures/NEUROD2.pdf')

# We predict that MM074, MM057 and MM087 (intermidate melanocytic cell state) will move to the right along principle component 0 which in this case corresponds to MEL-MES variation.

# ## Prioritizing TFs based on their preturbation effect.
# 
# We can also make use of the perturbation simulation to prioritize TFs for a certain effect. For example here we want to find TFs which will move cells along principle component 0 (i.e. change their state to a more melanocytic or mesenchymal one).
# 
# To do this we need some custom code.

# In[5]:


from scenicplus.simulation import *
from scenicplus.simulation import _make_rankings


# In[6]:

#from pyscenic.diff_features import *
from scenicplus.diff_features import *
get_differential_features(scplus_obj, 'GEX_supervised_name', use_hvg = True, contrast_type = ['DARs', 'DEGs'])
scplus_obj.uns['DEGs'].keys()


# First we train the gene expression random forrest models (this was previously done in the `plot_perturbation_effect_in_embedding` function).
# 
# We will train models for all differentially expressed genes (+ all identified TFs)

# In[7]:


flatten_list = lambda t: [item for sublist in t for item in sublist]
DEGs = list(set(flatten_list([list(scplus_obj.uns['DEGs']['GEX_supervised_name'][k].index) for k in scplus_obj.uns['DEGs']['GEX_supervised_name'].keys()])))
genes_to_use = list(set([*DEGs, * [x.split('_')[0] for x in scplus_obj.uns['selected_eRegulon']['Gene_based']]]))


# In[ ]:


regressors = train_gene_expression_models(
        scplus_obj,
        eRegulon_metadata_key = 'eRegulon_metadata',
        genes = genes_to_use)


# Next we will simulate a knockdown of all identified TFs and recalculate eRegulon enrichment values on the perturbed matrices using AUCell.

# In[ ]:


len(scplus_obj.uns['selected_eRegulon']['Gene_based'])


# In[ ]:


#parameters
TFs_of_interest = list(set([x.split('_')[0] for x in scplus_obj.uns['selected_eRegulon']['Gene_based']]))
TFs_of_interest = ['NR4A2', 'POU2F1', 'SOX5', 'VEZF1', 'EMX1', 'BHLHE22', 'ZBTB20', 'SOX2', 'TCF7L1', 'TFAP2C', 'MEF2C',  'POU3F2', 'ARX', 'TCF3', 'ZBTB18', 'SOX6', 'SOX9', 'NR2E1', 'CTCF', 'NEUROD6', 'JUN', 'ZNF148', 'NFIB', 
 'ZNF219', 'NEUROD2', 'ZNF281', 'NFIA', 'E2F1', 'ZNF441', 'FOS', 'TCF12', 'POU3F1', 'TBR1', 'KLF11', 'KLF10', 'SATB2', 'KAT7', 'ATF7', 'RFX2', 'KLF3', 'MEIS2', 'PHF21A', 'DEAF1', 'ASCL1']
eRegulon = [x.split('_')[0] for x in scplus_obj.uns['eRegulon_metadata_filtered'].Region_signature_name.unique()]
TFs_of_interest = [x for x in TFs_of_interest if x in eRegulon]

#all detected TFs
#TFs_of_interest = [x.split('_')[0] for x in scplus_obj.uns['eRegulon_metadata_filtered'].Region_signature_name.unique()]

TFs_of_interest = list(set(TFs_of_interest) | set([x.split('_')[0] for x in scplus_obj.uns['selected_eRegulon']['Gene_based']]))


n_iter = 5

from tqdm.notebook import tqdm
import logging
logging.basicConfig(level=logging.CRITICAL)
import warnings
from scenicplus.eregulon_enrichment import score_eRegulons
from scenicplus.dimensionality_reduction import run_eRegulons_pca
for TF in tqdm(TFs_of_interest, total = len(TFs_of_interest)):
        perturbed_matrices = simulate_perturbation(
                scplus_obj,
                perturbation = {TF: 0},
                regressors = regressors,
                keep_intermediate = True,
                n_iter = n_iter)
        perturbed_rankings = {k: _make_rankings(perturbed_matrices[k]) for k in perturbed_matrices.keys()}
        for k in perturbed_rankings.keys():
                with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        score_eRegulons(
                                scplus_obj,
                                ranking = perturbed_rankings[k],
                                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                                key_added = f'{TF}_KD_sim_eRegulon_AUC_iter_{k}',
                                enrichment_type = 'gene',
                                n_cpu = 8)



# Calculate shift along principle components

# In[ ]:


from scenicplus.simulation import _project_perturbation_in_embedding
shifts_PC0 = {}
shifts_PC1 = {}
import sys
sys.stderr = open(os.devnull, "w")  # silence stderr    
for TF in TFs_of_interest:
        delta_embedding = _project_perturbation_in_embedding(
                scplus_obj,
                original_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_0']['Gene_based'],
                perturbed_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_4']['Gene_based'],
                reduction_name = 'eRegulons_PCA_gene_based')
        mean_shift = pd.DataFrame(delta_embedding).groupby(scplus_obj.metadata_cell['GEX_supervised_name'].to_numpy()).mean()
        shifts_PC0[TF] = mean_shift[0]
        shifts_PC1[TF] = mean_shift[1]
#sys.stderr = sys.__stderr__  # unsilence stderr


#export shift tables
shift_df = pd.DataFrame(shifts_PC0).T
shift_df.to_csv(os.path.join(work_dir, 'scenicplus/shifts_PCO.csv'))
shift_df = pd.DataFrame(shifts_PC1).T
shift_df.to_csv(os.path.join(work_dir, 'scenicplus/shifts_PC1.csv'))


'''
# Plot factors with a strong effect in a heatmap


# In[17]:


shift_df = pd.DataFrame(shifts_PC1).T
factors_to_plot = [
        *shift_df.max(1).sort_values(ascending = False).index,
        *reversed(shift_df.min(1).sort_values(ascending = True).index)]
index_order = [ 'Astro', 'Oligo','OPC','RG1', 'RG2', 'IPC', 'EN_NB', 'EN',  'IN_CGE','IN_MGE', 'Vascular','MG']
import seaborn as sns
fig, ax = plt.subplots(figsize = (10, 5))
sns.heatmap(
        shift_df.loc[factors_to_plot, line_order].T,
        yticklabels=True,vmin = -0.3, vmax = 0.3, ax = ax, cmap = 'bwr')
#for ytick in ax.get_yticklabels():
#        ytick.set_color(color_dict_line[ytick.get_text()])
fig.tight_layout()
plt.savefig('/wynton/home/pollenlab/jding/BrainChromatin/scenicplus/scenicplus/PC0shift.pdf')



# In[51]:


shift_df = pd.DataFrame(shifts_PC1).T
factors_to_plot = [
        *shift_df.max(1).sort_values(ascending = False).head(40).index,
        *reversed(shift_df.min(1).sort_values(ascending = True).head(40).index)]
index_order = [ 'Astro', 'Oligo','OPC','RG1', 'RG2', 'IPC', 'EN_NB', 'EN',  'IN_CGE','IN_MGE', 'Vascular','MG']
import seaborn as sns
fig, ax = plt.subplots(figsize = (10, 5))
sns.heatmap(
        shift_df.loc[factors_to_plot, line_order].T,
        yticklabels=True,vmin = -0.3, vmax = 0.3, ax = ax, cmap = 'bwr')
#for ytick in ax.get_yticklabels():
#        ytick.set_color(color_dict_line[ytick.get_text()])
fig.tight_layout()
plt.savefig('/wynton/home/pollenlab/jding/BrainChromatin/scenicplus/scenicplus/PC0shift_40.pdf')


#export shift tables
shift_df = pd.DataFrame(shifts_PC0).T
shift_df.to_csv(os.path.join(work_dir, 'scenicplus/shifts_PCO.csv'))
shift_df = pd.DataFrame(shifts_PC1).T
shift_df.to_csv(os.path.join(work_dir, 'scenicplus/shifts_PC1.csv'))

'''


