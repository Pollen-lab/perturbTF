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
work_dir = '/wynton/group/pollen/jding/brainchromatin/multiome/seacells/atac_meta_allfrag/'
tmp_dir = '/wynton/scratch/jding/scenicplus'



# Load the scenicplus object from previous tutorial.

# In[2]:


import dill
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))


#EN = scplus_obj.cell_data[scplus_obj.cell_data.GEX_celltype.isin(['RG','IPC','EN1','EN2','EN3','EN4','EN5','EN6'])].index.tolist()
#scplus_obj = scplus_obj.subset(EN, copy=True)

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






# We predict that MM074, MM057 and MM087 (intermidate melanocytic cell state) will move to the right along principle component 0 which in this case corresponds to MEL-MES variation.

# ## Prioritizing TFs based on their preturbation effect.
# 
# We can also make use of the perturbation simulation to prioritize TFs for a certain effect. For example here we want to find TFs which will move cells along principle component 0 (i.e. change their state to a more melanocytic or mesenchymal one).
# 
# To do this we need some custom code.

# In[5]:


from scenicplus.simulation import *
from scenicplus.simulation import _make_rankings
#from scenicplus.diff_features import get_differential_features

# In[6]:

#variable = ['GEX_celltype']
#get_differential_features(scplus_obj, variable, use_hvg = True, contrast_type = ['DEGs', 'DARs'])
print(scplus_obj.uns['DEGs'].keys())


# First we train the gene expression random forrest models (this was previously done in the `plot_perturbation_effect_in_embedding` function).
# 
# We will train models for all differentially expressed genes (+ all identified TFs)

# In[7]:


flatten_list = lambda t: [item for sublist in t for item in sublist]
DEGs = list(set(flatten_list([list(scplus_obj.uns['DEGs']['GEX_celltype'][k].index) for k in scplus_obj.uns['DEGs']['GEX_celltype'].keys()])))
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
                                n_cpu = 12)



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
        mean_shift = pd.DataFrame(delta_embedding).groupby(scplus_obj.metadata_cell['GEX_celltype'].to_numpy()).mean()
        shifts_PC0[TF] = mean_shift[0]
        shifts_PC1[TF] = mean_shift[1]
#sys.stderr = sys.__stderr__  # unsilence stderr


#export shift tables
shift_df = pd.DataFrame(shifts_PC0).T
shift_df.to_csv(os.path.join(work_dir, 'scenicplus/shifts_PC0.csv'))
shift_df = pd.DataFrame(shifts_PC1).T
shift_df.to_csv(os.path.join(work_dir, 'scenicplus/shifts_PC1.csv'))


# Plot factors with a strong effect in a heatmap


# In[17]:

index_order = [ 'RG',  'IPC', 'EN1', 'EN2', 'EN3', 'EN4', 'EN6', 'IN1', 'IN2','IN3']
shift_df = pd.DataFrame(shifts_PC1).T
factors_to_plot = [
        *shift_df.max(1).sort_values(ascending = False).index,
        *reversed(shift_df.min(1).sort_values(ascending = True).index)]
line_order = index_order
import seaborn as sns
fig, ax = plt.subplots(figsize = (10, 5))
sns.heatmap(
        shift_df.loc[factors_to_plot, line_order].T,
        yticklabels=True,vmin = -0.3, vmax = 0.3, ax = ax, cmap = 'bwr')
#for ytick in ax.get_yticklabels():
#        ytick.set_color(color_dict_line[ytick.get_text()])
fig.tight_layout()
plt.savefig('/wynton/home/pollenlab/jding/BrainChromatin/scenicplus/figures/PC0shift.pdf')



# In[51]:


shift_df = pd.DataFrame(shifts_PC1).T
factors_to_plot = [
        *shift_df.max(1).sort_values(ascending = False).head(40).index,
        *reversed(shift_df.min(1).sort_values(ascending = True).head(40).index)]
line_order = index_order
import seaborn as sns
fig, ax = plt.subplots(figsize = (10, 5))
sns.heatmap(
        shift_df.loc[factors_to_plot, line_order].T,
        yticklabels=True,vmin = -0.3, vmax = 0.3, ax = ax, cmap = 'bwr')
#for ytick in ax.get_yticklabels():
#        ytick.set_color(color_dict_line[ytick.get_text()])
fig.tight_layout()
plt.savefig('/wynton/home/pollenlab/jding/BrainChromatin/scenicplus/figures/PC0shift_40.pdf')





