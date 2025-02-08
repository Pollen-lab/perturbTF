#!/usr/bin/env python
# coding: utf-8

# 
# # Kreigstein Lab Li Seacelled Multiome 
# 
# 
# The data consists of *PBMC from a Healthy Donor - Granulocytes Removed Through Cell Sorting (3k)* which is freely available from 10x Genomics (click [here](https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0), some personal information needs to be provided before you can gain access to the data). This is a multi-ome dataset.
# 
# <div class="alert alert-info">
# 
# **Note:**
# 
# In this notebook we will only show the minimal steps needed for running the SCENIC+ analysis. For more information on analysing scRNA-seq data and scATAC-seq data we refer the reader to other tutorials (e.g. [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/index.html) and [pycisTopic](https://pycistopic.readthedocs.io/en/latest/) in python or [Seurat](https://satijalab.org/seurat/) and [cisTopic](https://github.com/aertslab/cisTopic) or [Signac](https://satijalab.org/signac/) in R).
# 
# </div>
# 
# 

# ## Set-up environment and download data 
# We will first create a directory to store the data and results

# In[1]:


#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
_stderr = sys.stderr                                                         
null = open(os.devnull,'wb')


# In[2]:


import os
work_dir = '/wynton/group/pollen/jding/brainchromatin/Li/scenicplus/seacells'
tmp_dir = '/wynton/scratch/jding/'


# ## Loading preprocessed @BrainChromatin/Li/SEACell_computation_ATAC_allfrag.py
# First we preprocess the scRNA-seq side of the mutliome datasets. Most importantly we will use this side of the data to annotate celltypes. 
# 
# For this we will make use of [Scanpy](https://scanpy.readthedocs.io/en/stable/). 
# 
# 
# <div class="alert alert-info">
# 
# **Note:**
# 
# You may also use [Seurat](https://satijalab.org/seurat/) (or any other tool in fact) to preprocess your data, however this will require some extra steps to import the data in python.
# </div>
# 
# <div class="alert alert-info">
# 
# **Note:**
# 
# Further on in the actual SCENIC+ analysis the raw count matrix will be used.
# </div>
# 

# In[6]:


import scanpy as sc
#set some figure parameters for nice display inside jupyternotebooks.
#get_ipython().run_line_magic('matplotlib', 'inline')
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(5, 5), facecolor='white')

#make a directory for to store the processed scRNA-seq data.
if not os.path.exists(os.path.join(work_dir, 'scRNA')):
    os.makedirs(os.path.join(work_dir, 'scRNA'))


# We now have preprocessed the scRNA-seq side of the multiome data.
# 
# In particular we have:
# 
# 1. fitlered the data to only contain high quality cells.
# 2. annotated cells to cell types.
# 
# We also did some preliminary visualization of the data for which we needed to normalize the gene expression counts. Note that SCENIC+ uses the raw gene expression counts (i.e. without normalization and scaling). We have kept this raw data in `adata.raw`.
# 
# Now that we have clusters of annotated cells we can continue with preprocessing the scATAC-seq data. There we will use the annotated clusters of cells to generate pseudobulk ATAC-seq profiles per cell type which will be used for peak calling.

# ## scATAC-seq preprocessing using pycisTopic
# 
# Now we will preprocess the scATAC-seq side of the multiome data.
# 
# Most importantly we will:
# 
# 1. generate pseudobulk ATAC-seq profiles per cell type and call peaks
# 2. merge these peaks into a consensus peak-set
# 3. do quality control on the scATAC-seq barcodes
# 4. run topic modeling to find sets of co-accessible regions and to impute chromatin accessibility resolving the issue of drop outs
# 
# For this we will use the python package [pycisTopic](https://pycistopic.readthedocs.io/en/latest/). 
# 
# <div class="alert alert-info">
# 
# **Note:**
# 
# pycisTopic can also be used for independent analysis of scATAC-seq data and has many more features which will not be demonstrated here. For more information see the read the docs page of [pycisTopic](https://pycistopic.readthedocs.io/en/latest/)
# </div>

# In[3]:


import os
work_dir = '/wynton/group/pollen/jding/brainchromatin/Li/seacells'
import pycisTopic
#set some figure parameters for nice display inside jupyternotebooks.
#get_ipython().run_line_magic('matplotlib', 'inline')

#make a directory for to store the processed scRNA-seq data.
if not os.path.exists(os.path.join(work_dir, 'scATAC')):
    os.makedirs(os.path.join(work_dir, 'scATAC'))
tmp_dir = '/scratch/leuven/330/vsc33053/'

#read seacelled RNA/ATAC
adata = sc.read(os.path.join(work_dir,'rna_meta_ad.h5ad'))
ad = sc.read(os.path.join(work_dir,'atac_meta_ad.h5ad'))
ad

# ### Creating a cisTopic object and topic modeling
# 
# Now that we have good quality barcodes we will generate a binary count matrix of ATAC-seq fragments over consensus peaks. This matrix, along with metadata, will be stored in a cisTopic object and be used for topic modeling.
# 
# We will start by reading cell metadata from the scRNA-seq side of the analysis. For SCENIC+ we will only keep cells which passed quality metrics in both assays.
# 
# <div class="alert alert-info">
# 
# **Note:**
# 
# For independent scATAC-seq analysis you probably want to keep all cells (not only the cells passing the scRNA-seq filters). 
# </div>

# Load scATAC-seq data

# In[25]:


import pandas as pd
fragment_matrix = pd.DataFrame(ad.X.todense(), index=ad.obs_names,
                        columns=ad.var_names).copy()


# Create cisTopic object

# In[29]:


from pycisTopic.cistopic_class import create_cistopic_object
cistopic_obj = create_cistopic_object(fragment_matrix=fragment_matrix.T)
print(cistopic_obj)


adata.obs.index = adata.obs_names.astype(str) + '___cisTopic'
cistopic_obj.add_cell_data(adata.obs[['supervised_name']])


# In[44]:

cistopic_obj.cell_data['supervised_name'].value_counts()
index = cistopic_obj.cell_data.supervised_name.dropna().index.tolist()
cistopic_obj = cistopic_obj.subset(index, copy=True)
print(cistopic_obj)

# Save object.

# In[46]:


import pickle
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))


# Run topic modeling. The purpose of this is twofold:
# 
# 1. To find sets of co-accessible regions (topics), this will be used downstream as candidate enhancers (together with Differentially Accessible Regions (DARs)).
# 2. To impute dropouts. 
# 
# <div class="alert alert-info">
# 
# **Note:**
# 
# scATAC-seq data is *very* sparse. This is because of the nature of the technique. It profiles chromatin accessibility in single cells and each cell only has two copies (two alleles) of each genomic region, also the genome is realtively large (so the number of regions which can be measured is large). For this reason drop outs are a *real* problem, the chance of measuring a specific genomic region (out of *many* potential regions) which only has two copies per cell is relatively small. Compare this to scRNA-seq analysis, here the number of genes which can be measure is relatively small (compared to the number of genomic regions) and each cell has potentially hundres of copies of each transcript. To account for drop outs in the scATAC-seq assay imputation techniques are often used, in this case we make use of topic modeling, making use of the fact that the data often contains cells which are quite similar to each other but might have slightly different features measured. 
# </div>
# 
# 
# Before running the topic modeling we are not sure what the best number of topics will be to explain the data. Analog to PCA where you also often don't know before hand what the best number of principle components is. For this reason we will generate models with increasing numbers of topics and after the fact choose the model with the optimal amount of topics. For demonstration purposes we will only try a few amount of models, you might want to explore a larger number of topics.
# 
# <div class="alert alert-warning">
# 
# **Warning:**
# 
# Topic modeling can be computationaly intense!
# </div>

# In[47]:

tmp_dir = '/wynton/scratch/jding/'
import pickle
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
from pycisTopic.cistopic_class import *


os.environ['MALLET_MEMORY'] = '800G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path="/wynton/home/pollenlab/jding/Mallet-202108/bin/mallet"
# Run models
models=run_cgs_models_mallet(
    mallet_path, cistopic_obj,
    n_topics=[2,5,10,15,30,45,60],
    n_cpu=16,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=os.path.join(tmp_dir + 'ray_spill'),
    save_path=None
)
# Save results

# In[ ]:


if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))

pickle.dump(models,
            open(os.path.join(work_dir, 'scATAC/models/models_500_iter_mallet.pkl'), 'wb'))



# Analyze models.
# 
# We will make use of four quality metrics to select the model with the optimal amount of topics:
# 1. [Arun *et al.* 2010](http://link.springer.com/10.1007/978-3-642-13657-3_43)
# 2. [Cao & Juan *et al.* 2009](https://linkinghub.elsevier.com/retrieve/pii/S092523120800372X)
# 3. [Mimno *et al.* 2011](http://dirichlet.net/pdf/mimno11optimizing.pdf)
# 4. Log likelihood
# 
# For more information on these metrics see publications (linked above) and the [read the docs](https://pycistopic.readthedocs.io/en/latest/) page of pycisTopic.

# In[15]:

import pickle
models = pickle.load(open(os.path.join(work_dir, 'scATAC/models/models_500_iter_mallet.pkl'), 'rb'))
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
from pycisTopic.lda_models import *
model = evaluate_models(models,
                       select_model=30, 
                       return_model=True, 
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False)


# The metrics seem to stabelise with a model using 16 topics, so let's choose that model. 

# In[23]:


cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))



# For further analysis see the [pycisTopic read the docs page](https://pycistopic.readthedocs.io/en/latest/)

# ### Inferring candidate enhancer regions
# 
# Next we will infer candidate enhancer regions by:
# 
# 1. binarization of region-topic probabilites.
# 2. calculation differentially accessibile regions (DARs) per cell type.
# 
# These regions will be used as input for the next step, [pycistarget](https://pycistarget.readthedocs.io/en/latest/), in which we will look which motifs are enriched in these regions.

# First we will binarize the topics using the [otsu](http://ieeexplore.ieee.org/document/4310076/) method and by taking the top 3k regions per topic.

# In[19]:
import pickle
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))

from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)


# Next we will calculate DARs per cell type

# In[20]:


from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='supervised_name', var_features=variable_regions, split_pattern = '-')


# Save results

# In[25]:


if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))
import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'wb'))


# We now completed all the mininal scATAC-seq preprocessing steps. 
# 
# In particular we:
# 
# 1. generated a set of consensus peaks
# 2. performed quality control steps, only keeping cell barcods which passed QC metrics in both the scRNA-seq and scATAC-seq assay
# 3. performed topic modeling
# 4. inferred candidate enhancer regions by binarizing the region-topic probabilities and DARs per cell type
# 
# In the next step we will perform motif enrichment analysis on these candidate enhancer regions using the python package, [pycistarget](phttps://pycistarget.readthedocs.io/en/latest/). For this a precomputed motif-score database is needed. A sample specific database can be generated by scoring the consensus peaks with motifs or a general pre-scored database can be used as well.

# ## Motif enrichment analysis using pycistarget
# 
# After having identified candidate enhancer regions we will use [pycistarget](https://pycistarget.readthedocs.io/en/latest/) to find which motifs are enriched in these regions. 

# ### Cistarget databases
# 
# In order to run pycistarget one needs a precomputed database containing motif scores for genomic regions.
# 
# You can choose to compute this database yourself by scoring the consensus peaks generated in the scATAC-seq analysis using a set of motifs. The advantage of creating a sample specific database is that you can potentially pick up more target regions, given that only regions included/overlappig with regions in the cistarget database will be used for the SCENIC+ analysis. For more information checkout the [create_cisTarget_databases repo on github](https://github.com/aertslab/create_cisTarget_databases). 
# 
# We also provide several precomputed databases containing regions covering many experimentally defined candidate cis-regulatory elements. These databases are available on: [https://resources.aertslab.org/cistarget/](https://resources.aertslab.org/cistarget/).
# 
# For this analysis we will use a precomputed database using [screen regions](https://screen.encodeproject.org/).
# 
# Next to a precomputed motif database we also need a motif-to-tf annotation database. This is also available on [https://resources.aertslab.org/cistarget/](https://resources.aertslab.org/cistarget/).

# Load candidate enhancer regions identified in previous step.

# In[2]:


import pickle
import os
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'rb'))


# Convert to dictionary of pyranges objects.

# In[3]:


import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))


# In[22]:


for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')


# Define rankings, score and motif annotation database.
# 
# The ranking database is used for running the cistarget analysis and the scores database is used for running the DEM analysis. For more information see [the pycistarget read the docs page](https://pycistarget.readthedocs.io/en/latest/)
# 

# In[5]:


fpath = "/wynton/group/pollen/jding/cisTarget"


# In[6]:


rankings_db = os.path.join(fpath, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(fpath, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(fpath, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')


# Next we will run pycistarget using the `run_pycistarget` wrapper function.
# 
# This function will run cistarget based and DEM based motif enrichment analysis with or without promoter regions.
# 

# In[31]:


if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))


# In[42]:


import pandas as pd
annot = pd.read_csv('/wynton/home/pollenlab/jding/scenicplus/hg38_annot.csv',index_col=0)
annot


import pandas as pd
transcripts ='/wynton/group/pollen/genomes/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/regions/transcripts.bed'
transcripts = pd.read_csv(transcripts, sep='\t',  header = None)
transcripts.columns=['Chromosome', 'Start','End','GeneID','4','Strand']
transcripts.index = transcripts['GeneID']
transcripts = transcripts.drop(columns=['GeneID','4'])


#assign gene names from geneInfo
geneInfo ='/wynton/group/pollen/genomes/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star/geneInfo.tab'
geneInfo = pd.read_csv(geneInfo, sep='\t',  skiprows = 1, header=None, index_col=0)
geneInfo.columns=['Gene', 'Transcript_type']
geneInfo = geneInfo[geneInfo.Transcript_type == 'protein_coding']
transcripts = transcripts[transcripts.index.isin(geneInfo.index)]
transcripts['Gene'] = 'NA'
transcripts['Transcript_type'] = 'protein_coding'
transcripts['Strand'] = [1 if x == '+' else -1 for x in transcripts['Strand']]


#only keep genes from annotated chromosomes
filter = transcripts['Chromosome'].str.startswith('chr')
transcripts = transcripts[filter]

for x in transcripts.index:
    transcripts.loc[x,'Gene'] = geneInfo['Gene'][x]
    
transcripts


# In[ ]:


from scenicplus.wrappers.run_pycistarget import run_pycistarget
tmp_dir = '/scratch/leuven/330/vsc33053/'
run_pycistarget(
    region_sets = region_sets,
    species = 'custom',
    custom_annot =  transcripts,  
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 16,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10nr_clust',
    ignore_reinit_error=True,
    include_dashboard=False,
    )

'''
from scenicplus.wrappers.run_pycistarget import run_pycistarget
tmp_dir = '/scratch/leuven/330/vsc33053/'
run_pycistarget(
    region_sets = region_sets,
    species = 'homo_sapiens',
    #custom_annot =  transcripts,  
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 16,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10nr_clust',
    #ignore_reinit_error=True,
    #include_dashboard=False,
    )


run_pycistarget(
    region_sets = region_sets,
    species = 'custom',
    custom_annot =  transcripts,  
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 8,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10nr_clust',
    ignore_reinit_error=True,
    include_dashboard=False,
    annotation=['Direct_annot', 'Orthology_annot','Motif_similarity_annot','Motif_similarity_and_Orthology_annot'],
    ctx_auc_threshold = 0.004,
    ctx_nes_threshold = 3.0,
    ctx_rank_threshold = 0.04,
    dem_log2fc_thr = 0.4,
    dem_motif_hit_thr = 3.0,
    dem_max_bg_regions = 400,
    motif_similarity_fdr = 0.00001,
    )

'''
# Let's explore some of the results. Below we show the motifs found for topic 8 (specific to B-cells) using DEM.

# In[33]:


import dill
menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))


# In[34]:


menr['DEM_topics_otsu_All'].DEM_results('Topic8')


# We now have completed all the steps necessary for starting the SCENIC+ analysis 😅.
# 
# In particalular, we have
# 
# 1. preprocessed the scRNA-seq side of the data, selecting high quality cells and annotation these cells.
# 2. preprocessed the scATAC-seq side of the data, selecting high quality cells, performing topic modeling and identifying candidate enhacer regions.
# 3. looked for enriched motifs in candidate enhancer regions.
# 
# In the next section we will combine all these analysis and run SCENIC+
# 

# ## inferring enhancer-driven Gene Regulatory Networks (eGRNs) using SCENIC+
# 
# We now have completed all the steps for running the SCENIC+ analysis. 
# 
# We will start by creating a scenicplus object containing all the analysis we have done up to this point.
# 
# For this we will need to load:
# 
# 1. the AnnData object containing the scRNA-seq side of the analysis.
# 2. the cisTopic object containing the scATAC-seq side of the analysis.
# 3. the motif enrichment dictionary containing the motif enrichment results.

# In[44]:


import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
import numpy as np
_stderr = sys.stderr                                                         
null = open(os.devnull,'wb')
tmp_dir = 'wynton/scratch/jding/leuven/330/vsc33053/'
adata= sc.read(os.path.join(work_dir, 'scRNA/adata2.h5ad'))
cistopic_obj = dill.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))


# In[45]:


adata.obs.index = adata.obs_names.astype(str) + '___cisTopic'
adata.obs.head()


# In[47]:


cistopic_obj.cell_data


# Create the SCENIC+ object. It will store both the gene expression and chromatin accessibility along with motif enrichment results and cell/region/gene metadata.
# 
# Cell metadata comming from the cistopic_obj will be prefixed with the string `ACC_` and metadata comming from the adata object will be prefixed with the string `GEX_`.

# In[49]:

adata.X = adata.layers['raw']

from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr,
    multi_ome_mode = True,
    key_to_group_by = 'supervised_name'
)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
scplus_obj




# <div class="alert alert-info">
# 
# **Note:**
# 
# the scenicplus package contains many function, if you need help with any of them just run the `help()` function. For example see below:
# </div>

# In[7]:


from scenicplus.scenicplus_class import create_SCENICPLUS_object
help(create_SCENICPLUS_object)


# Before running SCENIC+ it is important to check with which biomart host the gene names used in your analysis match best. Biomart will be used to find transcription starting sites of each gene. The names of genes (symbols) change quite often, so it is important to select the biomart host with the largest overlap, otherwise a lot of genes can potentially be lost. 
# 


# Once the object is created we can run SCENIC+. The whole analysis workflow can be run using a single wrapper function or each step can be run individually.
# 
# <div class="alert alert-info">
# 
# **Note:**
# 
# Here we will use the wrapper function. For more information on the different steps, see the following tutorial: TO BE ADDED. 
# </div>

# In[41]:


biomart_host = "http://sep2019.archive.ensembl.org/"



# In[4]:


import dill
import scanpy as sc


# Now we are ready to run the analysis. 
# 
# <div class="alert alert-warning">
# 
# **Warning:**
# 
# Running SCENIC+ can be computationaly expensive. We don't recommend to run it on your local machine.
# </div>

# In[52]:


import pandas as pd
import pyranges as pr
f_fasta_index = '/wynton/group/pollen/genomes/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa.fai' #provde the path to the fasta index file, often this file is named simething like: genome.fa.fai
chromsizes = pd.read_csv(
    f_fasta_index, sep = '\t', names = ['Chromosome', 'Length', 'Offset', 'Linebases', 'Linewidth'])[['Chromosome', 'Length']] #we only care about the chromosome name and the length of the chromosomes
chromsizes['Start'] = 0 #for PyRanges the dataframe should be formatted as Chromosome, Start, End. We can add an extra column defining the start of each chromosome to 0.
chromsizes = chromsizes[['Chromosome', 'Start', 'Length']] #set correct order of columns.
chromsizes.columns = ['Chromosome', 'Start', 'End'] #rename Length to End
filter = chromsizes['Chromosome'].str.startswith('chr')
chromsizes = chromsizes[filter]
pr_chromsizes = pr.PyRanges(chromsizes) #convert pandas DataFrame to pyranges object.
pr_chromsizes.head()


# In[70]:


import pandas as pd
transcripts ='/wynton/group/pollen/genomes/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/regions/transcripts.bed'
transcripts = pd.read_csv(transcripts, sep='\t',  header = None)
transcripts.columns=['Chromosome', 'Start','End','GeneID','4','Strand']
transcripts.index = transcripts['GeneID']
transcripts = transcripts.drop(columns=['GeneID','4'])


#assign gene names from geneInfo
geneInfo ='/wynton/group/pollen/genomes/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star/geneInfo.tab'
geneInfo = pd.read_csv(geneInfo, sep='\t',  skiprows = 1, header=None, index_col=0)
geneInfo.columns=['Gene', 'Transcript_type']
geneInfo = geneInfo[geneInfo.Transcript_type == 'protein_coding']
transcripts = transcripts[transcripts.index.isin(geneInfo.index)]
transcripts['Gene'] = 'NA'
transcripts['Transcript_type'] = 'protein_coding'

#only keep genes from annotated chromosomes
filter = transcripts['Chromosome'].str.startswith('chr')
transcripts = transcripts[filter]

for x in transcripts.index:
    transcripts.loc[x,'Gene'] = geneInfo['Gene'][x]
    
transcripts.head()


# In[71]:


transcripts['Transcription_Start_Site'] = transcripts['Start'] 

pr_annot = pr.PyRanges(transcripts.dropna(axis=0))

if not any(['chr' in c for c in scplus_obj.region_names]):
    pr_annot.Chromosome = pr_annot.Chromosome.str.replace('chr', '')
    
pr_annot.head() 


# In[53]:


#for species different from human, mouse or fruit fly
from scenicplus.enhancer_to_gene import get_search_space
get_search_space(
    scplus_obj,
    pr_annot = pr_annot, #see above
    pr_chromsizes = pr_chromsizes, #see above
    upstream = [1000, 150000], downstream = [1000, 150000],
    #biomart_host = 'http://apr2020.archive.ensembl.org/'
)


# In[50]:


from scenicplus.wrappers.run_scenicplus import run_scenicplus
tmp_dir = 'wynton/scratch/jding/'
biomart_host = "http://sep2019.archive.ensembl.org/"
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['supervised_name'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = '/wynton/home/pollenlab/jding/scenicplus/utoronto_human_tfs_v_1.01.txt',
        save_path = os.path.join(work_dir, 'scenicplus'),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = True,
        export_to_UCSC_file = True,
        path_bedToBigBed = '/wynton/home/pollenlab/jding/Scenicplus/',
        n_cpu = 16,
        _temp_dir = os.path.join(tmp_dir, 'ray'),
        #address='127.0.0.1'
    )
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)

'''


# In[33]:


scplus_obj.uns.keys()


# In[ ]:


dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)


# ## Note on the output of SCENIC+
# 
# After running the SCENIC+ analysis the `scplus_obj` will be populated with the results of several analysis. Below we will describe how you can explore these.
# 
# For structuring this data we took inspiration from the [AnnData](https://anndata.readthedocs.io/en/latest/) format.

# In[ ]:


import os
os._exit(00)


# In[14]:


import warnings
import os
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')


# In[15]:


import dill
work_dir = '/wynton/group/pollen/jding/brainchromatin/multiome/seacells/atac_meta_allfrag/'
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))
#scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj2.pkl'), 'rb'))


# In[3]:


scplus_obj


# ### Gene expression and chromatin accessibility data
# 
# Both the raw gene expression counts and chromatin accessibility data are stored in the `scplus_obj` and can be accessed by running, `scplus_obj.to_df('EXP')` or `scplus_obj.to_df('ACC')`.

# In[4]:


scplus_obj.to_df('EXP').head()


# In[5]:


scplus_obj.to_df('ACC').head()


# ### Cell, region and gene metadata
# 
# Cell metatdata is stored in the `.metadata_cell` slot. Data comming from the gene expression side is prefixed with the string `GEX_` ad data comming from the chromatin accessibility side is prefixed with the string `ACC_`.
# 
# Region metadata is stored in the `.metadata_region` slot.
# 
# Gene metadata is stored in the `.metadata_gene` slot.

# In[6]:


scplus_obj.metadata_cell.head()


# In[7]:


scplus_obj.metadata_regions.head()


# In[8]:


scplus_obj.metadata_genes.head()


# ### Motif enrichment data
# 
# Motif enrichment data is stored in the `.menr` slot.

# In[9]:


scplus_obj.menr.keys()


# ### Dimensionality reduction data
# 
# Dimensionality reductions of the cells are stored in the `.dr_cell` slot. Reductions comming from the gene expression side of the data are prefixed with the string `GEX_` and those comming from the chromatin accessibility side with `ACC_`.
# 
# IF region based dimensionality reductions were calculated then those will be stored in the `.dr_region` slote (not the case in this example).

# In[10]:


scplus_obj.dr_cell.keys()


# ### Unstructured data
# 
# Additional unstructured data will be stored in the `.uns` slot. After running the standard scenicplus wrapper function this slot will contain the following entries:
# 
# 1. `Cistromes`: this contains TFs together with target regions based on the motif enrichment analysis (i.e. **prior to running SCENIC+**)
# 2. `search_space`: this is a dataframe containing the search space for each gene.
# 3. `region_to_gene`: this is a dataframe containing region to gene links **prior to running SCENIC+** (i.e unfiltered/raw region to gene importance scores and correlation coefficients).
# 4. `TF2G_adj`: this is a datafram containing TF to gene links **prior to running SCENIC+** (i.e unfiltered/raw TF to gene importance scores and correlation coefficients).
# 
# These four slots contain all the data necessary for running the SCENIC+ analysis. The following entries are produced by the SCENIC+ analysis:
# 
# 1. `eRegulons`: this is the raw output from the SCENIC+ analysis. We will go into a bit more detail for these below.
# 2. `eRegulon_metadata`: this is a dataframe containing the same information as `eRegulons` bit in a format which is a bit easier to parse for a human.
# 3. `eRegulon_signatures`: this is a dictionary with target regions and genes for each eRegulon
# 4. `eRegulon_AUC`: this slot contains dataframes with eRegulon enrichment scores calculated using AUCell (see below).
# 5. `pseudobulk`: contains pseudobulked gene expression and chromatin accessibility data, this is used to calculated TF to eRegulon correlation values.
# 6. `TF_cistrome_correlation`: contains correlation values between TF expression and eRegulon enrichment scores (seperate entries for target gene and target region based scores).
# 7. `eRegulon_AUC_thresholds`: contains thresholds on the AUC values (eRegulon enrichment scores), this is necessary to be able to visualize the results in [SCope](https://scope.aertslab.org/)
# 8. `RSS`: contains eRegulon Specificity Scores (RSS), a measure on how cell type specific an eRegulon is.
# 9. `DEGs`: contains Differentially Expressed Genes.
# 10. `DARs`: contains Differentially Accessibile Regions. 

# In[51]:


scplus_obj.uns.keys()


# #### The eRegulons entry
# 
# The main output of SCENIC+ are eRegulons.
# 
# This is initially stored in a list of `eRegulon` classes as depicted below.

# In[52]:


scplus_obj.uns['eRegulons'][0:5]


# each eRegulon has the following information (attributes):
# 
# 1. `cistrome_name`: name of the cistrome (from `scenicplus.uns['Cistromes']`) from which this eRegulon was created.
# 2. `context`: specifies the binarization method(s) used for binarizing region to gene relationships and wether positive/negative region-to-gene and TF-to-gene relationships were used.
# 3. `gsea_adj_pval`/`gsea_enrichment_score`/`gsea_pval`/`in_leading_edge`: are internal parameters used for generating the eRegulons. The values are lost when generating the final eRegulons because results from several analysis (different binarization methods) are combined.
# 4. `is_extended`: specifies wether extended (i.e. non-direct) motif-to-TF annotations were used.
# 5. `n_target_genes`: number of target genes.
# 6. `n_target_regions`: number of target regions.
# 7. `regions2genes`: region to gene links **after running SCENIC+**.
# 8. `target_genes`: target genes of the eRegulon
# 9. `target_regions`: target regions of the eRegulon
# 10. `transcription_factor`: TF name

# In[40]:


for attr in dir(scplus_obj.uns['eRegulons'][0]):
    if not attr.startswith('_'):
        print(f"{attr}: {getattr(scplus_obj.uns['eRegulons'][0], attr) if not type(getattr(scplus_obj.uns['eRegulons'][0], attr)) == list else getattr(scplus_obj.uns['eRegulons'][0], attr)[0:5]}")


# The information of all eRegulons is combined in the `eRegulon_metadata` dataframe.

# In[24]:


scplus_obj.uns['eRegulon_metadata'].head()


# For the eRegulon names we use the following convetion:
# 
# `<TF NAME>_<TF-TO-GENE RELATIONSHIP (+/-)>_<REGION-TO-GENE RELATIONSHIP (+/-)>_(NUMBER OF TARGET REGIONS(r)/GENES(g))`
# 
# For example the name: `ARID3A_+_+_(364r)` and `ARID3A_+_+_(278g)` indicates that the we found an eRegulon for the TF `ARID3A` which has `364` target regions and `278` target genes, that expression of the TF correlates positively with the expression of all the target genes (first `+` sign) and that the accessibility of all target regions correlates positively with the expression of all target regions (seconf `+` sign).
# 
# Given this convention we can have at maximum 4 (actually 8, see note below) different eRegulons for each TF, however not all four combinations are always found.
# 
# The table below describes a possible biological interpretation of each combination:
# 
# 
# combinations of signs | TF-to-gene correlation | region-to-gene correlation | interpretation
# ----------------------|------------------------|----------------------------|---------------
# "+ +"                 | positive               | positive                   | This eRegulon is an activator which opens the chromatin of the target regions and induces gene expression of the target genes.
# "+ -"                 | positive               | negative                   | This eRegulon is an activator which closes the chromatin of the target regions and induces gene expression of the target genes (this is biologically quite implausible).
# "- +"                 | negative               | positive                   | This eRegulon is a repressor which closes the chromatin of the target regions and represses the expression of the target genes.
# "- -"                 | negative               | negative                   | This eRegulon is a repressor which opens the chromatin of the target regions and represses the expression of the target genes (this is biologically quite implausible).
# 
# *table of possible eRegulon names*
# 
# <div class="alert alert-warning">
# 
# **Warning:**
# 
# We suggest to mainly focus on activator eRegulons. Repressor eRegulons contain more false predicitons because it is more difficult to find negative correlations (it relies on the abscence of data instead of the prescence) and they are often found because of the prescence of a different TF from the same familly (i.e. one which has the same/similar DNA binding domain and thus DNA motif) for which the expression is anti-correlated. You can never trust predictions from SCENIC+ blindly but certainly the repressor eRegulons you should handle with caution.
# </div>
# 
# <div class="alert alert-info">
# 
# **Note:**
# 
# Each eRegulon can be derived from motifs which are annotated *directly* to the TF (for example the annotation comes from a ChIP-seq experiment in a cell line of the same species) or the annotation can be *extended* (the annotation is infered based on orthology or motif similarity).
# The direct annotations are of higher quality, for that reason the we add the suffix `_extended` to eRegulons derived from extended annotations. Because both annotations can exists for the same TF we can have a maximum of 8 eRegulons per TF (4 from all the combinations described above each of which can be either extended or direct). We will handle all these types below and simplify the output. 
# </div>

# ## Downstream analysis
# 
# We have finally finished the SCENIC+ analysis 🎉. Now we can start analysing the results, we will show some example analysis below but really the sky is the limit so feel free to be creative with your newly acquired eRegulons.

# In[11]:


import warnings
import os
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')


# In[12]:


import warnings
import os
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')
import dill
import pyranges
import os
#from pyranges import pyranges
work_dir = '/wynton/group/pollen/jding/brainchromatin/multiome/seacells/atac_meta_allfrag/'
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))


# In[89]:


import dill 
#from pyranges import .pyranges
work_dir = '/wynton/group/pollen/jding/brainchromatin/multiome/seacells/atac_meta_allfrag_old/'
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))


# In[93]:


import dill 
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))


# In[9]:


import dill
dill.__version__


# In[10]:


import scenicplus
scenicplus.__version__


# ### Simplifying and filtering SCENIC+ output
# 
# Given the multitude of eRegulons that can be generated for each TF (see above) we will first simplify the result by:
# 1. Only keeping eRegulons with an extended annotation if there is no direct annotation available (given that the confidence of direct motif annotations is in genral higher).
# 2. Discarding eRegulons for which the region-to-gene correlation is negative (these are often noisy).
# 3. Renaming the eRegulons so that eRegulons with the suffix `TF_+_+` become `TF_+` and those with `TF_-_+` become `TF_-`.

# In[13]:


from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
apply_std_filtering_to_eRegulons(scplus_obj)


# This will create two new entries in the scenicplus object: `scplus_obj.uns['eRegulon_metadata_filtered']` and `scplus_obj.uns['eRegulon_signatures_filtered']` containing the simplified results.
# We will use these for downstream analysis.

# In[12]:


scplus_obj.uns['eRegulon_metadata_filtered'].head()


# In[14]:


scplus_obj.uns['eRegulon_metadata_filtered'].head()


# In[16]:


scplus_obj.uns['eRegulon_metadata_filtered']['TF'].unique()


# ### eRegulon enrichment scores
# We can score the enrichment of eRegulons using the AUCell function. This function takes as input a gene or region based ranking (ranking of genes/regions based on the expression/accessibility per cell) and a list of eRegulons.
# 
# These values were already calculated in the wrapper function but let's recalculate them using the filtered output.

# In[17]:


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


# ### eRegulon dimensionality reduction
# Based on the enrichment scores calculated above we can generate dimensionality reductions (e.g. tSNE and UMAP).
# 
# To calculate these dimensionality reductions we use both the regions **and** gene based enrichment scores.

# In[26]:


from scenicplus.dimensionality_reduction import run_eRegulons_tsne, run_eRegulons_umap
run_eRegulons_umap(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_UMAP', #overwrite previously calculated UMAP
)
run_eRegulons_tsne(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_tSNE', #overwrite previously calculated tSNE
)


# In[15]:


from scenicplus.dimensionality_reduction import harmony
harmony(
    scplus_obj = scplus_obj,
    variable = 'GEX_sample',
    auc_key = 'eRegulon_AUC_filtered',
    out_key = 'eRegulon_AUC_harmony'
)


# Let's visualize the UMAP and tSNE stored respectively in `eRegulons_UMAP` and `eRegulons_tSNE`, these are calculated based on the combined region and gene AUC values described above. 
# 
# Let's also add some nice colours by specifying a color_dictionary.

# In[16]:


from scenicplus.dimensionality_reduction import plot_metadata_given_ax
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

#specify color_dictionary
'''
color_dict = {
    'B_cells_1': "#065143",
    'B_cells_2': "#70B77E",
    'CD4_T_cells': "#E0A890",
    'CD8_T_cells': "#F56476",
    'NK_cells': "#CE1483",
    'Dendritic_cells': "#053C5E" ,
    'FCGR3A+_Monocytes': "#38A3A5",
    'CD14+_Monocytes': "#80ED99"
}
'''

fig, axs = plt.subplots(ncols=2, figsize = (16, 8))
plot_metadata_given_ax(
    scplus_obj=scplus_obj,
    ax = axs[0],
    reduction_name = 'eRegulons_UMAP',
    variable = 'GEX_celltype', #note the GEX_ prefix, this metadata originated from the gene expression metadata (on which we did the cell type annotation before)
    dot_size = 50,
    #color_dictionary={'GEX_celltype': color_dict}
)

fig.tight_layout()
sns.despine(ax = axs[0]) #remove top and right edge of axis border
sns.despine(ax = axs[1]) #remove top and right edge of axis border
plt.show()


# ### plot the activity / expression of an eRegulon on the dimensionality reduction
# 
# Nex we visualize the gene expression and target gene and region activity of some eRegulons on the tSNE.

# In[87]:


from scenicplus.dimensionality_reduction import plot_eRegulon
plot_eRegulon(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_UMAP',
    selected_regulons = ['EOMES_+', 'INSM1_+','PAX6_+','SOX4_+','NR2E1_+'],
    scale = True,
    dot_size = 50,
    auc_key = 'eRegulon_AUC_filtered')


# In[5]:


from scenicplus.dimensionality_reduction import plot_eRegulon
plot_eRegulon(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_UMAP',
    selected_regulons = ['NFIA_+','TBR1_+', 'CTCF_+','MEIS2_+'],
    scale = True,
    dot_size = 50,
    auc_key = 'eRegulon_AUC_filtered')


# In[85]:


from scenicplus.dimensionality_reduction import plot_eRegulon
plot_eRegulon(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_UMAP',
    selected_regulons = ['NEUROD6_extended_+','NEUROD2_+','SOX6_+','SOX2_+','SOX1_+','EMX2_+'],
    scale = True,
    dot_size = 50,
    auc_key = 'eRegulon_AUC_filtered')


# In[23]:


from scenicplus.dimensionality_reduction import plot_eRegulon
plot_eRegulon(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_UMAP',
    selected_regulons = ['TCF3_-', 'RFX2_-', 'ZBTB18_-', 'ELF1_-', 'ZBTB20_-', 'KLF3_-', 'NFIB_-', 'EMX1_-'],
    scale = True,
    dot_size = 50,
    auc_key = 'eRegulon_AUC_filtered')


# #### Riverplot- TF expression/accessibility/target gene expression

# In[4]:


import warnings
import os
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import pandas as pd
import numpy as np
import sklearn
_stderr = sys.stderr
null = open(os.devnull,'wb')


# In[3]:


import dill
work_dir = '/wynton/group/pollen/jding/brainchromatin/multiome'
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))


# In[6]:


#TF expression
dgem = pd.DataFrame(scplus_obj.X_EXP, index=scplus_obj.cell_names,
                        columns=scplus_obj.gene_names).copy()
dgem = dgem.T / dgem.T.sum(0) * 10**6 #normalize
dgem = np.log1p(dgem).T
dgem = dgem.T


# In[7]:


#region accessibility    
auc_key= 'eRegulon_AUC_filtered'
signature_keys= ['Region_based']
#scaled
region_mat = pd.concat([pd.DataFrame(sklearn.preprocessing.StandardScaler().fit_transform(
            scplus_obj.uns[auc_key][x].T), index=scplus_obj.uns[auc_key][x].T.index.to_list(), columns=scplus_obj.uns[auc_key][x].T.columns) for x in signature_keys])
selected_regulons = [x for x in region_mat.index if x.split('_')[0] in dgem.index]
#subset to only actiavtors
selected_regulons = [x for x in selected_regulons if x.split('_')[-2] == '+']
region_mat = region_mat.loc[region_mat.index.isin(selected_regulons)]


# In[8]:


#gene expression 
auc_key= 'eRegulon_AUC_filtered'
signature_keys= ['Gene_based']
#scaled
exp_mat = pd.concat([pd.DataFrame(sklearn.preprocessing.StandardScaler().fit_transform(
            scplus_obj.uns[auc_key][x].T), index=scplus_obj.uns[auc_key][x].T.index.to_list(), columns=scplus_obj.uns[auc_key][x].T.columns) for x in signature_keys])
selected_regulons = [x for x in exp_mat.index if x.split('_')[0] in dgem.index]
selected_regulons = [x for x in selected_regulons if x.split('_')[-2] == '+']
exp_mat = exp_mat.loc[exp_mat.index.isin(selected_regulons)]


# In[10]:


selected_genes = [x.split('_')[0] for x in region_mat.index]
dgem = dgem.loc[dgem.index.isin(selected_genes)]


# In[44]:


import anndata 
adata = anndata.AnnData(layers= {'exp': dgem.sort_index(ascending=True).T,'access': region_mat.sort_index(ascending=True).T,'target': exp_mat.sort_index(ascending=True).T},
                        obs= exp_mat.columns.tolist(),
                        var=[x.rsplit('_', 1)[0] for x in exp_mat.sort_index(ascending=True).index])


# In[46]:


adata.obs.index = adata.obs[0]


# In[47]:


adata.var.index = [x.split('_')[0] for x in adata.var[0]]


# In[48]:


import scanpy as sc
EN = anndata.read_h5ad('/wynton/home/pollenlab/jding/BrainChromatin/mira/data/RNA_EN_mira.h5ad')
EN = EN[EN.obs.index.isin(adata.obs.index)]
adata = adata[adata.obs.index.isin(EN.obs.index)]
adata.obs['mira_pseudotime'] = EN.obs['mira_pseudotime'] 


# In[49]:


adata


# In[267]:


scplus_obj.uns['eRegulon_metadata'][scplus_obj.uns['eRegulon_metadata']['TF'] == 'SOX4'].sort_values(by = 'TF2G_rho', ascending=False).head(50)


# We can also plot only the activity of an eRegulon

# In[104]:


from scenicplus.dimensionality_reduction import plot_AUC_given_ax

fig, ax = plt.subplots(figsize = (8,8))
plot_AUC_given_ax(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_tSNE',
    feature = 'PAX5_+_(119g)',
    ax = ax,
    auc_key = 'eRegulon_AUC_filtered',
    signature_key = 'Gene_based')
sns.despine(ax = ax)
plt.show()


# ### dotplot-heatmap
# 
# For eRegulons it is often usefull to visualize both information on the TF/target genes expression and region accessibility at the same time.
# 
# A dotplot-heatmap is a useful way to visualize this. Here the color of the heatmap can be used to visualize one aspect of the eRegulon (for example TF expression) and the size of the dot can be used to visualize another aspect (for example the enrichment (AUC value) of eRegulon target regions).
# 
# Before we plot the the dotplot-heatmap let's first select some high quality eRegulons to limit the amount of space we need for the plot. One metric which can be used for selecting eRegulons is the correlation between TF expression and target region enrichment scores (AUC values). Let's (re)calculate this value based on the simplified eRegulons
# 
# We first generate pseudobulk gene expression and region accessibility data, per celltype, to limit the amount of noise for the correlation calculation.

# In[18]:


scplus_obj.metadata_cell['GEX_celltype'] = scplus_obj.metadata_cell['GEX_celltype'].cat.remove_unused_categories()
scplus_obj.metadata_cell['GEX_celltype'].value_counts()


# In[22]:


from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks

generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'GEX_celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Gene_based',
        nr_cells=6)
generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'GEX_celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Region_based',
        nr_cells=6)

TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'GEX_celltype',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Gene_based',
            out_key = 'filtered_gene_based')

TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'GEX_celltype',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Region_based',
            out_key = 'filtered_region_based')


# In[20]:


scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].head()


# In[21]:


scplus_obj.uns['TF_cistrome_correlation']['filtered_gene_based'].sort_values(by=[ 'Rho'],ascending = False).head(120).sort_values(by=['TF'],ascending = True)


# In[22]:


df = scplus_obj.uns['TF_cistrome_correlation']['filtered_gene_based']
df.index= df['TF']
df.loc['EOMES']


# In[23]:


df = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']
df.index= df['TF']
df.loc['EOMES']


# In[4]:


df = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']
df.index= df['TF']
df.loc['ZBTB18']


# In[24]:


df = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']
df.index= df['TF']
df.loc['INSM1']


# In[25]:


df = scplus_obj.uns['TF_cistrome_correlation']['filtered_gene_based']
df.index= df['TF']
df.loc['INSM1']


# In[49]:


df = scplus_obj.uns['TF_cistrome_correlation']['filtered_gene_based']
df.index= df['TF']
df.loc['FEZF2']


# Let's visualize these correlations in a scatter plot and select eRegulons for which the correlaiton coefficient is above 0.70 or below -0.75

# In[26]:


import numpy as np
n_targets = [int(x.split('(')[1].replace('r)', '')) for x in scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
rho = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
adj_pval = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

thresholds = {
        'rho': [-0.75, 0.70],
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


# In[24]:


import numpy as np
import matplotlib.pyplot as plt

n_targets = [int(x.split('(')[1].replace('r)', '')) for x in scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
rho = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
adj_pval = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

thresholds = {
        'rho': [-0.75, 0.70],
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


# In[8]:


selected_eRegulons_gene_sig


# In[14]:


[x for x in selected_eRegulons_gene_sig if x.split('_')[-2] == '+']


# In[22]:


[x.split('_')[0] for x in selected_eRegulons_gene_sig if x.split('_')[-2] == '+']


# In[27]:


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


# In[49]:


print(scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based']['SOX2_+_(305g)'])


# In[52]:


print(scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based']['NR2E1_+_(257g)'])


# In[27]:


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


# Let's save these changes we have made to the scenicplus_obj

# In[28]:


dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)


# Let's plot the heatmap-dotplot

# In[68]:


index_order = [ 'RG',  'IPC', 'EN1', 'EN2', 'EN3', 'EN4', 'EN6', 'IN1', 'IN2','IN3']


# In[30]:


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
#del scplus_obj.uns['selected_eRegulon'] 
scplus_obj.uns['selected_eRegulon'] = {'Gene_based': selected_eRegulons_gene_sig, 'Region_based': selected_eRegulons_region_sig}
print(f'selected: {len(selected_eRegulons_gene_sig)} eRegulons')


from scenicplus.plotting.dotplot import heatmap_dotplot
heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'], #specify what to plot as dot sizes, target region enrichment in this case
        color_matrix = scplus_obj.to_df('EXP'), #specify  what to plot as colors, TF expression in this case
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'GEX_celltype',
        subset_eRegulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
        index_order = index_order,
        figsize = (10, 15),
        orientation = 'vertical',
        save = '/wynton/home/pollenlab/jding/BrainChromatin/seacells/heatmap.pdf'
)


# In[39]:


from scenicplus.plotting.dotplot import heatmap_dotplot
heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'], #specify what to plot as dot sizes, target region enrichment in this case
        color_matrix = scplus_obj.to_df('EXP'), #specify  what to plot as colors, TF expression in this case
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'GEX_celltype',
        subset_eRegulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
        index_order = index_order,
        figsize = (10, 15),
        orientation = 'vertical',
        save = '/wynton/home/pollenlab/jding/BrainChromatin/seacells/heatmap.pdf')


# In[71]:


from scenicplus.plotting.dotplot import heatmap_dotplot
heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'], #specify what to plot as dot sizes, target region enrichment in this case
        color_matrix = scplus_obj.to_df('EXP'), #specify  what to plot as colors, TF expression in this case
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'GEX_celltype',
        subset_eRegulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
        index_order = index_order,
        figsize = (10, 15),
        orientation = 'vertical',
        #save = '/wynton/home/pollenlab/jding/BrainChromatin/seacells/heatmap.pdf'
)


# ### Integrated Multiome Plot
# We can also generate plots showing the chromatin profiles per group, region-to-gene relationships and TF and gene expression.

# In[ ]:


import pickle
work_dir = '/wynton/group/pollen/jding/brainchromatin/multiome'
infile = open('/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/10x_multiome_brain/output/scenicplus/scplus_obj.pkl', 'rb')
scplus_obj = pickle.load(infile)
infile.close()


# In[74]:


# Generate interaction and annotation pyranges
import matplotlib.pyplot as plt
import os
from scenicplus.utils import get_interaction_pr
import pyranges as pr
bigwig_dir = os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bw_files/')
bw_dict = {x.replace('.bw', ''): os.path.join(bigwig_dir, x) for x in os.listdir(bigwig_dir) if '.bw' in x}
pr_consensus_bed = pr.read_bed(os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'))
pr_interact = get_interaction_pr(scplus_obj, 'hsapiens', 'hg38', inplace = False, subset_for_eRegulons_regions = True, eRegulons_key = 'eRegulons')


# In[88]:


del bw_dict['Astro_Oligo']
del bw_dict['Endo_Peri']
bw_dict['Astro/Oligo']= bw_dict.pop('AstroOligo')
bw_dict['Endo/Peri']= bw_dict.pop('EndoPeri')
bw_dict


# In[81]:


gtf_file = "/wynton/group/pollen/jding/ref/hg38/gencode.v41.annotation.gtf.gz"
pr_gtf = pr.read_gtf(gtf_file)


# In[318]:


# Plot
from importlib import reload
from scenicplus.plotting.coverageplot import *
fig = coverage_plot(
        SCENICPLUS_obj = scplus_obj,
        bw_dict = bw_dict,
        region = 'chr3:27990337-28004260',
        figsize = (10,20),
        pr_gtf = pr_gtf,
        color_dict = None,
        plot_order = None,
        pr_interact = pr_interact,
        genes_violin_plot = ['EOMES'],
        meta_data_key = 'GEX_celltype',
        pr_consensus_bed = pr_consensus_bed,
        #fontsize_dict={'bigwig_label': 12, 'gene_label': 0, 'violinplots_xlabel': 10, 'title': 12, 'bigwig_tick_label': 0, 'violinplots_ylabel': 3},
        #height_ratios_dict = {'bigwig_violin': 1, 'genes': 0.5, 'arcs': 10, 'custom_ax': 5}
)
plt.tight_layout()
    


# In[322]:


# Plot
from importlib import reload
from scenicplus.plotting.coverageplot import *
fig = coverage_plot(
        SCENICPLUS_obj = scplus_obj,
        bw_dict = bw_dict,
        region = 'chr3:27950337-28004260',
        figsize = (10,20),
        pr_gtf = pr_gtf,
        color_dict = None,
        plot_order = None,
        pr_interact = pr_interact,
        #genes_violin_plot = ['EOMES'],
        meta_data_key = 'GEX_celltype',
        pr_consensus_bed = pr_consensus_bed,
        #fontsize_dict={'bigwig_label': 12, 'gene_label': 0, 'violinplots_xlabel': 10, 'title': 12, 'bigwig_tick_label': 0, 'violinplots_ylabel': 3},
        #height_ratios_dict = {'bigwig_violin': 1, 'genes': 0.5, 'arcs': 10, 'custom_ax': 5}
)
plt.tight_layout()
    


# In[103]:


# Plot
from importlib import reload
from scenicplus.plotting.coverageplot import *
fig = coverage_plot(
        SCENICPLUS_obj = scplus_obj,
        bw_dict = bw_dict,
        region = 'chr5:126755164-126781991',#ZNF43 regulon'chr5:126758449-126758949'
        figsize = (10,20),
        pr_gtf = pr_gtf,
        color_dict = None,
        plot_order = None,
        pr_interact = pr_interact,
        genes_violin_plot = ['LMNB1','ZNF43'],
        meta_data_key = 'GEX_celltype',
        pr_consensus_bed = pr_consensus_bed,
        #fontsize_dict={'bigwig_label': 12, 'gene_label': 0, 'violinplots_xlabel': 10, 'title': 12, 'bigwig_tick_label': 0, 'violinplots_ylabel': 3},
        #height_ratios_dict = {'bigwig_violin': 1, 'genes': 0.5, 'arcs': 10, 'custom_ax': 5}
)
plt.tight_layout()
    


# In[112]:


# Plot
from importlib import reload
from scenicplus.plotting.coverageplot import *
fig = coverage_plot(
        SCENICPLUS_obj = scplus_obj,
        bw_dict = bw_dict,
        region = 'chr2:172076751-172098430',#ZNF536 regulon chr2:172093567-172094067
        figsize = (10,20),
        pr_gtf = pr_gtf,
        color_dict = None,
        plot_order = None,
        pr_interact = pr_interact,
        genes_violin_plot = ['DLX1','ZNF536'], 
        meta_data_key = 'GEX_celltype',
        pr_consensus_bed = pr_consensus_bed,
        #fontsize_dict={'bigwig_label': 12, 'gene_label': 0, 'violinplots_xlabel': 10, 'title': 12, 'bigwig_tick_label': 0, 'violinplots_ylabel': 3},
        #height_ratios_dict = {'bigwig_violin': 1, 'genes': 0.5, 'arcs': 10, 'custom_ax': 5}
)
plt.tight_layout()
    


# In[116]:


# Plot
from importlib import reload
from scenicplus.plotting.coverageplot import *
fig = coverage_plot(
        SCENICPLUS_obj = scplus_obj,
        bw_dict = bw_dict,
        region = 'chr11:115436826-115660864', #chr11:115645634-115646134
        figsize = (10,20),
        pr_gtf = pr_gtf,
        color_dict = None,
        plot_order = None,
        pr_interact = pr_interact,
        genes_violin_plot = ['CADM1','ZNF536'], 
        meta_data_key = 'GEX_celltype',
        pr_consensus_bed = pr_consensus_bed,
        #fontsize_dict={'bigwig_label': 12, 'gene_label': 0, 'violinplots_xlabel': 10, 'title': 12, 'bigwig_tick_label': 0, 'violinplots_ylabel': 3},
        #height_ratios_dict = {'bigwig_violin': 1, 'genes': 0.5, 'arcs': 10, 'custom_ax': 5}
)
plt.tight_layout()
    


# ### overlap of predicted target regions
# An interesting aspect of gene regulation is transcription factor cooperativity (i.e. multiple TFs cobinding the same enhancer together driving gene expression). 
# 
# By looking at the overlap of predicted target regions of TFs we can infer potential cooperativity events.
# 
# Let's look at the overlap of target regions of the top 5 TFs per cell type based on the Regulon Specificity Score (RSS). 
# 
# First we calculate the RSS for the target regions of the selected eRegulons.

# In[28]:


from scenicplus.RSS import *
regulon_specificity_scores(
        scplus_obj,
        variable = 'GEX_celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_keys = ['Region_based'],
        selected_regulons = [x for x in scplus_obj.uns['selected_eRegulon']['Region_based'] if '-' not in x],
        out_key_suffix = '_filtered')


# Let's visualize the RSS values using a scatter plot

# In[32]:


plot_rss(scplus_obj, 'GEX_celltype_filtered', num_columns=4, top_n=10, figsize = (20, 10))


# Next we select the top 10 eRegulons per cell type

# In[37]:


df = scplus_obj.uns['TF2G_adj']
df[df['target'] =='NR2E1'].sort_values(by = 'importance_x_rho',ascending = False).head(20)


# In[40]:


df = scplus_obj.uns['TF2G_adj']
df[df['target'] =='SOX2'].sort_values(by = 'importance_x_rho',ascending = False).head(50)


# In[28]:


highlight_features = ['EN2', 'NEUROD2', 'BHLHE22', 'ZBTB18', 'MEF2C', 'ETS1', 'EMX1', 'NFIB', 'JUN', 'POU3F1', 'SOX2', 'POU3F2', 'POU2F1', 'TCF12', 'ZNF441', 'ARX', 'ZNF219', 'ZBTB20', 'ASCL1', 'TCF7L1', 'TCF3', 'ZNF148', 'NR2E1', 'VEZF1', 'SOX9', 'NEUROD6', 'TFAP2C', 'PHF21A', 'KLF11', 'RFX2', 'ATF7', 'KAT7', 'ZNF444', 'KLF3', 'KLF10']
len([x for x in scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'].columns.to_list() if x.split('_')[0] in highlight_features])


# In[34]:


highlight_features = ['EN2', 'NEUROD2', 'BHLHE22', 'ZBTB18', 'MEF2C', 'ETS1', 'EMX1', 'NFIB', 'JUN', 'POU3F1', 'SOX2', 'POU3F2', 'POU2F1', 'TCF12', 'ZNF441', 'ARX', 'ZNF219', 'ZBTB20', 'ASCL1', 'TCF7L1', 'TCF3', 'ZNF148', 'NR2E1', 'VEZF1', 'SOX9', 'NEUROD6', 'TFAP2C', 'PHF21A', 'KLF11', 'RFX2', 'ATF7', 'KAT7', 'ZNF444', 'KLF3', 'KLF10']
selected_markers = [x for x in scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'].columns.to_list() if x.split('_')[0] in highlight_features]
selected_markers = [x for x in selected_markers if x.split('_')[1] == '+' or  x.split('_')[2] == '+']
len(selected_markers)


# In[7]:


flat_list = lambda t: [item for sublist in t for item in sublist]
selected_markers = list(set(flat_list(
    [scplus_obj.uns['RSS']['GEX_celltype_filtered'].loc[celltype].sort_values(ascending = False).head(10).index.to_list() 
    for celltype in scplus_obj.uns['RSS']['GEX_celltype_filtered'].index])))


# In[38]:


from scenicplus.plotting.correlation_plot import *

region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
        use_plotly = False,
        selected_regulons = selected_markers,
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (10, 10), return_data = True, vmax = 0.5, cmap = 'plasma')


# In[39]:


red = ['TFAP2C', 'ARX', 'SOX2', 'ASCL1', 'JUN', 'ZNF441', 'TCF3', 'ZBTB20', 'NR2E1', 'KLF11', 'ZNF148', 'ZNF219', 'POU2F1', 'NFIB', 'PHF21A', 'TCF12', 'KAT7', 'RFX2', 'KLF11', 'ZBTB20']
selected_markers = [x for x in scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'].columns.to_list() if x.split('_')[0] in red]
selected_markers = [x for x in selected_markers if x.split('_')[1] == '+' or  x.split('_')[2] == '+']

from scenicplus.plotting.correlation_plot import *

region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
        use_plotly = False,
        selected_regulons = selected_markers,
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (10, 10), return_data = True, vmax = 0.5, cmap = 'plasma')


# In[40]:


blue = ['KLF10', 'ZNF444', 'VEZF1', 'KLF3', 'ATF7', 'POU3F1', 'NEUROD2', 'POU3F2', 'BHLHE22', 'NEUROD6','ETS1', 'EN2', 'KLF10', 'ZBTB18', 'SOX9', 'EMX1', 'NEUROD2', 'MEF2C', 'TCF7L1', 'NEUROD6']
selected_markers = [x for x in scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'].columns.to_list() if x.split('_')[0] in blue]
selected_markers = [x for x in selected_markers if x.split('_')[1] == '+' or  x.split('_')[2] == '+']

from scenicplus.plotting.correlation_plot import *

region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
        use_plotly = False,
        selected_regulons = selected_markers,
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (10, 10), return_data = True, vmax = 0.5, cmap = 'plasma')


# In[34]:


from scenicplus.plotting.correlation_plot import *

region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
        use_plotly = False,
        selected_regulons = selected_markers,
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (10, 10), return_data = True, vmax = 0.5, cmap = 'plasma')


# In[110]:


import dill
work_dir = '/wynton/group/pollen/jding/brainchromatin/multiome/seacells/atac_meta_allfrag/'
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))


# In[111]:


scplus_obj.to_df('EXP').head()


# In[112]:


#TF expression
dgem = pd.DataFrame(scplus_obj.X_EXP, index=scplus_obj.cell_names,
                        columns=scplus_obj.gene_names).copy()
#dgem = dgem.T / dgem.T.sum(0) * 10**6 #skip normalization
#dgem = np.log1p(dgem).T
dgem = dgem.T


# In[113]:


#region accessibility    
auc_key= 'eRegulon_AUC_filtered'
signature_keys= ['Region_based']
#scaled
region_mat = pd.concat([pd.DataFrame(sklearn.preprocessing.StandardScaler().fit_transform(
            scplus_obj.uns[auc_key][x].T), index=scplus_obj.uns[auc_key][x].T.index.to_list(), columns=scplus_obj.uns[auc_key][x].T.columns) for x in signature_keys])
selected_regulons = [x for x in region_mat.index if x.split('_')[0] in dgem.index]
#subset to only actiavtors
selected_regulons = [x for x in selected_regulons if x.split('_')[-2] == '+']
region_mat = region_mat.loc[region_mat.index.isin(selected_regulons)]


# In[114]:


#gene expression 
auc_key= 'eRegulon_AUC_filtered'
signature_keys= ['Gene_based']
#scaled
exp_mat = pd.concat([pd.DataFrame(sklearn.preprocessing.StandardScaler().fit_transform(
            scplus_obj.uns[auc_key][x].T), index=scplus_obj.uns[auc_key][x].T.index.to_list(), columns=scplus_obj.uns[auc_key][x].T.columns) for x in signature_keys])
selected_regulons = [x for x in exp_mat.index if x.split('_')[0] in dgem.index]
selected_regulons = [x for x in selected_regulons if x.split('_')[-2] == '+']
exp_mat = exp_mat.loc[exp_mat.index.isin(selected_regulons)]


# In[115]:


selected_genes = [x.split('_')[0] for x in region_mat.index]
dgem = dgem.loc[dgem.index.isin(selected_genes)]


# In[126]:


import anndata 
adata = anndata.AnnData(layers= {'exp': dgem.sort_index(ascending=True).T,'access': region_mat.sort_index(ascending=True).T,'target': exp_mat.sort_index(ascending=True).T},
                        obs= exp_mat.columns.tolist(),
                        var=[x.rsplit('_', 1)[0] for x in exp_mat.sort_index(ascending=True).index])


# In[127]:


adata.obs['Cellname'] = adata.obs[0]
adata.obs.index = adata.obs['Cellname']
del adata.obs[0]


# In[128]:


adata.var['Gene'] = adata.var[0]
adata.var.index = [x.split('_')[0] for x in adata.var['Gene']]
del adata.var[0]


# In[130]:


adata.write(os.path.join(work_dir, 'scenicplus/TF.h5ad'))


# ### Plotting a network
# 
# eRegulons can also be visualized in a network. Simple plots can be made using python. For more complicated plots (i.e. containing many nodes and edges) we suggest exporting your network to cytoscape.
# 
# Let's create a very simple network for B cells. We will use the top 1000 highly variable regions and genes in this plot. If you want to use more feautures please export your nework to cytoscape.

# In[61]:


from pycisTopic.diff_features import find_highly_variable_features
hvr = find_highly_variable_features(scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Region']))], n_top_features=2300, plot = False)
hvg = find_highly_variable_features(scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Gene']))].T, n_top_features=2300, plot = False)


# In[46]:


from pycisTopic.diff_features import find_highly_variable_features
hvr = find_highly_variable_features(scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Region']))], n_top_features=1000, plot = False)
hvg = find_highly_variable_features(scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Gene']))].T, n_top_features=1000, plot = False)


# First we format the eRegulons into a table which can be used to create a network using the package [networkx](https://networkx.org/)

# In[ ]:


genes = ['PAX6', 'SOX9','SOX2', 'FOXG1', 'ASCL1']


# In[47]:


genes = ['NEUROD2', 'POU3F1','ZBTB18', 'BLHLE22', 'POU3F2']


# In[62]:


genes = ['SOX2', 'SOX9','TFAP2C', 'NR2E1', 'ARX']


# In[63]:


from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx, export_to_cytoscape
nx_tables = create_nx_tables(
    scplus_obj = scplus_obj,
    eRegulon_metadata_key ='eRegulon_metadata_filtered',
    subset_eRegulons = genes,
    subset_regions = hvr,
    subset_genes = hvg,
    add_differential_gene_expression = True,
    add_differential_region_accessibility = True,
    differential_variable = ['GEX_celltype'])


# In[72]:


nx_tables['Edge']['R2G']['Gene_signature_name'].unique()


# Next we layout the graph.

# In[64]:


###IPC

G, pos, edge_tables, node_tables = create_nx_graph(nx_tables, 
                   use_edge_tables = ['TF2R','R2G'],
                   color_edge_by = {'TF2R': {'variable' : 'TF', 'category_color' : {genes[0]: 'Orange', 
                                                                                    genes[1]: 'Purple', 
                                                                                    genes[2]: 'Red',
                                                                                    genes[3]: 'Green',
                                                                                    genes[4]: 'Blue'}},
                                    'R2G': {'variable' : 'R2G_rho', 'continuous_color' : 'viridis', 'v_min': -1, 'v_max': 1}},
                   transparency_edge_by =  {'R2G': {'variable' : 'R2G_importance', 'min_alpha': 0.1, 'v_min': 0}},
                   width_edge_by = {'R2G': {'variable' : 'R2G_importance', 'max_size' :  1.5, 'min_size' : 1}},
                   color_node_by = {'TF': {'variable': 'TF', 'category_color' : {genes[0]: 'Orange', 
                                                                                    genes[1]: 'Purple', 
                                                                                    genes[2]: 'Red',
                                                                                    genes[3]: 'Green',
                                                                                    genes[4]: 'Blue'}},
                                    'Gene': {'variable': 'GEX_celltype_Log2FC_IPC', 'continuous_color' : 'bwr'},
                                    'Region': {'variable': 'GEX_celltype_Log2FC_IPC', 'continuous_color' : 'viridis'}},
                   transparency_node_by =  {'Region': {'variable' : 'GEX_celltype_Log2FC_IPC', 'min_alpha': 0.1},
                                    'Gene': {'variable' : 'GEX_celltype_Log2FC_IPC', 'min_alpha': 0.1}},
                   size_node_by = {'TF': {'variable': 'fixed_size', 'fixed_size': 30},
                                    'Gene': {'variable': 'fixed_size', 'fixed_size': 15},
                                    'Region': {'variable': 'fixed_size', 'fixed_size': 10}},
                   shape_node_by = {'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}},
                   label_size_by = {'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 30.0},
                                    'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 13.0},
                                    'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}},
                   layout='kamada_kawai_layout',
                   scale_position_by=250)


# In[49]:


###RG

G, pos, edge_tables, node_tables = create_nx_graph(nx_tables, 
                   use_edge_tables = ['TF2R','R2G'],
                   color_edge_by = {'TF2R': {'variable' : 'TF', 'category_color' : {genes[0]: 'Orange', 
                                                                                    genes[1]: 'Purple', 
                                                                                    genes[2]: 'Red',
                                                                                    genes[3]: 'Green',
                                                                                    genes[4]: 'Blue'}},
                                    'R2G': {'variable' : 'R2G_rho', 'continuous_color' : 'viridis', 'v_min': -1, 'v_max': 1}},
                   transparency_edge_by =  {'R2G': {'variable' : 'R2G_importance', 'min_alpha': 0.1, 'v_min': 0}},
                   width_edge_by = {'R2G': {'variable' : 'R2G_importance', 'max_size' :  1.5, 'min_size' : 1}},
                   color_node_by = {'TF': {'variable': 'TF', 'category_color' : {genes[0]: 'Orange', 
                                                                                    genes[1]: 'Purple', 
                                                                                    genes[2]: 'Red',
                                                                                    genes[3]: 'Green',
                                                                                    genes[4]: 'Blue'}},
                                    'Gene': {'variable': 'GEX_celltype_Log2FC_RG', 'continuous_color' : 'bwr'},
                                    'Region': {'variable': 'GEX_celltype_Log2FC_RG', 'continuous_color' : 'viridis'}},
                   transparency_node_by =  {'Region': {'variable' : 'GEX_celltype_Log2FC_RG', 'min_alpha': 0.1},
                                    'Gene': {'variable' : 'GEX_celltype_Log2FC_RG', 'min_alpha': 0.1}},
                   size_node_by = {'TF': {'variable': 'fixed_size', 'fixed_size': 30},
                                    'Gene': {'variable': 'fixed_size', 'fixed_size': 15},
                                    'Region': {'variable': 'fixed_size', 'fixed_size': 10}},
                   shape_node_by = {'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}},
                   label_size_by = {'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 30.0},
                                    'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 13.0},
                                    'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}},
                   layout='kamada_kawai_layout',
                   scale_position_by=250)


# Finally we can visualize the network.
# 
# In this network diamond shapes represent regions and they are color coded by their log2fc value in B cells target genes and TFs are visualized using circles and are labeled.

# In[60]:


from scenicplus.dimensionality_reduction import plot_metadata_given_ax
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

plt.figure(figsize=(15,15))
plot_networkx(G, pos)


# In[50]:


from scenicplus.dimensionality_reduction import plot_metadata_given_ax
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

plt.figure(figsize=(15,15))
plot_networkx(G, pos)


# In[44]:


from scenicplus.dimensionality_reduction import plot_metadata_given_ax
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

plt.figure(figsize=(15,15))
plot_networkx(G, pos)



'''