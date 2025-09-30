#!/usr/bin/env python
# coding: utf-8

# # Preprocessing 

# In[3]:


import os
os.chdir("/Users/gclyu07/Desktop/XMAS_analysis/")


# In[4]:


import pycisTopic
pycisTopic.__version__


# In[5]:


import os
out_dir = "outputs"
os.makedirs(out_dir, exist_ok = True)


# In[4]:


fragments_dict = {
    "D11": "/Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D11/atac_fragments.tsv.gz",
    "D16": "/Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D16/atac_fragments.tsv.gz",
    "D28": "/Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D28/atac_fragments.tsv.gz",
    "D42": "/Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D42/atac_fragments.tsv.gz",
    "D56": "/Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D56/atac_fragments.tsv.gz"
}

fragments_dict = {
    "XMAS": "/Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_cellranger_ARC_aggr_02/atac_fragments.tsv.gz",
}
# In[6]:


import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as plt


# In[7]:


adata = sc.read_h5ad('/Users/gclyu07/Desktop/XMAS_analysis/outputs/s04_XMAS_velocyto_spliced_unspliced_corrected.h5ad')


# In[8]:


cell_data = adata.obs


# In[9]:


cell_data['barcode'] = [s.split('-')[0] + '-1' for s in cell_data.index]
cell_data['sampleid'] = 'XMAS'


# In[9]:


cell_data


# In[10]:


cell_data.index = cell_data['barcode']


# In[12]:


cell_data


# # Getting pseudobulk profiles from cell annotations

# In[17]:


chromsizes = pd.read_table(
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
chromsizes.head()


# In[20]:


from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)


# In[22]:


bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data,
    variable = "bimod_oc_des",
    sample_id_col = "orig.ident",
    chromsizes = chromsizes,
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragments_dict,
    normalize_bigwig = True,
    temp_dir = "/tmp",
    split_pattern = "-"
)


# In[23]:


with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")


# In[24]:


with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")


# # Getting pseudobulk profiles from cell annotations

# In[37]:


bw_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bw_paths.update({v: p})


# In[38]:


bed_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})


# In[39]:


bw_paths


# In[40]:


bed_paths


# In[27]:


from pycisTopic.pseudobulk_peak_calling import peak_calling


# In[28]:


macs_path = "/Users/gclyu07/anaconda3/envs/SC_v4/bin/macs2"


# In[41]:


os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok = True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(out_dir, "consensus_peak_calling/MACS")),
    genome_size = 'hs',
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = '/tmp/peakcalling'
)


# In[44]:


from pycisTopic.iterative_peak_calling import get_consensus_peaks
# Other param
peak_half_width=250
path_to_blacklist="/Users/gclyu07/pycisTopic/blacklist/hg38-blacklist.v2.bed"
# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist)


# In[45]:


consensus_peaks.to_bed(
    path = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep =True,
    compression = 'infer',
    chain = False)


# # Quality control

# In[46]:


get_ipython().system('pycistopic tss gene_annotation_list | grep Human')


# In[47]:


get_ipython().system('mkdir -p outputs/qc')
get_ipython().system('pycistopic tss get_tss      --output outputs/qc/tss.bed      --name "hsapiens_gene_ensembl"      --to-chrom-source ucsc      --ucsc hg38')


# In[48]:


get_ipython().system('head outputs/qc/tss.bed | column -t')


# In[54]:


fragments_dict.items()


# In[55]:


regions_bed_filename = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
tss_bed_filename = os.path.join(out_dir, "qc", "tss.bed")

pycistopic_qc_commands_filename = "pycistopic_qc_commands.txt"

# Create text file with all pycistopic qc command lines.
with open(pycistopic_qc_commands_filename, "w") as fh:
    for sample, fragment_filename in fragments_dict.items():
        print(
            "pycistopic qc",
            f"--fragments {fragment_filename}",
            f"--regions {regions_bed_filename}",
            f"--tss {tss_bed_filename}",
            f"--output {os.path.join(out_dir, 'qc')}/{sample}",
            sep=" ",
            file=fh,
        )

cat pycistopic_qc_commands.txt | parallel -j 4 {}~/anaconda3/envs/scenicplus/bin/pycistopic qc --fragments /Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D11/atac_fragments.tsv.gz --regions outputs/consensus_peak_calling/consensus_regions.bed --tss outputs/qc/tss.bed --output outputs/qc/D11
~/anaconda3/envs/scenicplus/bin/pycistopic qc --fragments /Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D16/atac_fragments.tsv.gz --regions outputs/consensus_peak_calling/consensus_regions.bed --tss outputs/qc/tss.bed --output outputs/qc/D16
~/anaconda3/envs/scenicplus/bin/pycistopic qc --fragments /Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D28/atac_fragments.tsv.gz --regions outputs/consensus_peak_calling/consensus_regions.bed --tss outputs/qc/tss.bed --output outputs/qc/D28
~/anaconda3/envs/scenicplus/bin/pycistopic qc --fragments /Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D42/atac_fragments.tsv.gz --regions outputs/consensus_peak_calling/consensus_regions.bed --tss outputs/qc/tss.bed --output outputs/qc/D42
~/anaconda3/envs/scenicplus/bin/pycistopic qc --fragments /Users/gclyu07/Desktop/XMAS_analysis/dataset/XMAS_D56/atac_fragments.tsv.gz --regions outputs/consensus_peak_calling/consensus_regions.bed --tss outputs/qc/tss.bed --output outputs/qc/D56
# In[56]:


from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
import matplotlib.pyplot as plt


# In[57]:


for sample_id in fragments_dict:
    fig = plot_sample_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = "outputs/qc"
    )


# In[58]:


from pycisTopic.qc import get_barcodes_passing_qc_for_sample
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id = sample_id,
            pycistopic_qc_output_dir = "outputs/qc",
            unique_fragments_threshold = None, # use automatic thresholding
            tss_enrichment_threshold = None, # use automatic thresholding
            frip_threshold = 0,
            use_automatic_thresholds = True,
    )


# In[59]:


for sample_id in fragments_dict:
    fig = plot_barcode_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = "outputs/qc",
        bc_passing_filters = sample_id_to_barcodes_passing_filters[sample_id],
        detailed_title = False,
        **sample_id_to_thresholds[sample_id]
    )


# # Creating a cisTopic object 

# In[190]:


from pycisTopic.cistopic_class import *


# In[61]:


path_to_regions = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
path_to_blacklist = "/Users/gclyu07/pycisTopic/blacklist/hg38-blacklist.v2.bed"
pycistopic_qc_output_dir = "outputs/qc"

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl

cistopic_obj_list = []
for sample_id in fragments_dict:
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = 1,
        project = sample_id,
        split_pattern = '-'
    )
    cistopic_obj_list.append(cistopic_obj)


# In[71]:


cistopic_obj = cistopic_obj_list[0]
print(cistopic_obj)


# In[133]:


cistopic_obj_list[1].cell_data


# In[143]:


cistopic_obj_list[0].cell_data['ori_barcode'] =  [s.split('-')[0] + '-1' for s in cistopic_obj_list[0].cell_data.index]
cistopic_obj_list[1].cell_data['ori_barcode'] =  [s.split('-')[0] + '-2' for s in cistopic_obj_list[1].cell_data.index]
cistopic_obj_list[2].cell_data['ori_barcode'] =  [s.split('-')[0] + '-3' for s in cistopic_obj_list[2].cell_data.index]
cistopic_obj_list[3].cell_data['ori_barcode'] =  [s.split('-')[0] + '-4' for s in cistopic_obj_list[3].cell_data.index]
cistopic_obj_list[4].cell_data['ori_barcode'] =  [s.split('-')[0] + '-5' for s in cistopic_obj_list[4].cell_data.index]


# In[148]:


cell_data


# In[168]:


cistopic_obj_list[1].cell_data


# In[253]:


cistopic_obj = merge(cistopic_obj_list)


# In[210]:


import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_merged.pkl"), "wb")
)


# In[254]:


cistopic_obj.cell_data


# #  Adding metadata to cisTop objects
cistopic_obj.cell_data.index = cistopic_obj.cell_data['ori_barcode']
# In[255]:


# Extract barcodes from both DataFrames
cistopic_barcodes = cistopic_obj.cell_data['ori_barcode']
cell_data_barcodes = cell_data.index

# Find barcodes that are in both
common_barcodes = cistopic_barcodes[cistopic_barcodes.isin(cell_data_barcodes)]

# Find barcodes that are in cistopic_obj but not in cell_data
only_in_cistopic = cistopic_barcodes[~cistopic_barcodes.isin(cell_data_barcodes)]

# Find barcodes that are in cell_data but not in cistopic_obj
only_in_cell_data = cell_data_barcodes[~cell_data_barcodes.isin(cistopic_barcodes)]

# Display results
print(f"Common barcodes: {len(common_barcodes)}")
print(f"Barcodes only in cistopic_obj: {len(only_in_cistopic)}")
print(f"Barcodes only in cell_data: {len(only_in_cell_data)}")


# In[257]:


len(common_barcodes.index.to_list())


# In[259]:


cistopic_obj = cistopic_obj.subset(cells=common_barcodes.index.to_list())


# In[262]:


len(cistopic_obj.cell_names)


# In[264]:


import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_merged_matched.pkl"), "wb")
)


# In[267]:


cistopic_obj.cell_names = cistopic_obj.cell_data['ori_barcode']


# In[271]:


cistopic_obj.cell_data.index = cistopic_obj.cell_data['ori_barcode']


# In[280]:


cistopic_obj.cell_names


# In[281]:


cell_data


# In[282]:


import pandas as pd
cell_data.head()
cistopic_obj.add_cell_data(cell_data, split_pattern='-')
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)


# In[283]:


cistopic_obj.cell_data


# # Running scrublet (Skipped, using AMULET)

# In[284]:


import scrublet as scr
scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=0.1)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram();
scrub.call_doublets(threshold=0.22)
scrub.plot_histogram();
scrublet = pd.DataFrame([scrub.doublet_scores_obs_, scrub.predicted_doublets_], columns=cistopic_obj.cell_names, index=['Doublet_scores_fragments', 'Predicted_doublets_fragments']).T


# In[285]:


cistopic_obj.add_cell_data(scrublet, split_pattern = '-')
sum(cistopic_obj.cell_data.Predicted_doublets_fragments == True)


# In[286]:


pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_scrublet.pkl"), "wb")
)

# Remove doublets
singlets = cistopic_obj.cell_data[cistopic_obj.cell_data.Predicted_doublets_fragments == False].index.tolist()
# Subset cisTopic object
cistopic_obj_noDBL = cistopic_obj.subset(singlets, copy=True, split_pattern='-')
print(cistopic_obj_noDBL)pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)
# In[107]:


adata_cis = adata[adata.obs_names.isin(cistopic_obj.cell_names)]


# In[113]:


adata_cis


# In[111]:


adata_cis.write_h5ad("/Users/gclyu07/Desktop/XMAS_analysis/outputs/adata_pycisTopic_matched.h5ad")


# # Run models 

# In[288]:


get_ipython().system('wget https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz')
get_ipython().system('tar -xf Mallet-202108-bin.tar.gz')

#!/bin/bash -l

#SBATCH -A naiss2024-22-807
#SBATCH -t 12:00:00
#SBATCH -J 081524_XMAS_topic
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=180gb
#SBATCH --requeue
#SBATCH --output=081524_XMAS_topic.out
#SBATCH --error=081524_XMAS_topic.err
#SBATCH --mail-type=END
#SBATCH --mail-user=guochang.lyu@ki.se

cd /proj/gclyu07d/scenicplus

source activate scenicplus
python ./scripts/XMAS_topic.pyimport os
import pycisTopic
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as plt
import pickle

os.chdir("/proj/gclyu07d/scenicplus")
out_dir = "outputs"
os.makedirs(out_dir, exist_ok = True)

pycisTopic.__version__

with open('/proj/gclyu07d/scenicplus/cistopic_obj_scrublet.pkl', "rb") as file:
    cistopic_obj = pickle.load(file)

os.environ['MALLET_MEMORY'] = '200G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path="/proj/gclyu07d/scenicplus/Mallet-202108/bin/mallet"
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=10,
    n_iter=500,
    random_state=1121,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path="/proj/gclyu07d/scenicplus/tmp",
    save_path="/proj/gclyu07d/scenicplus/outputs",
    mallet_path=mallet_path,
)

pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)
# # Model selection 

# In[11]:


import pickle


# In[12]:


with open('/Users/gclyu07/Desktop/XMAS_analysis/outputs/cistopic_obj.pkl', "rb") as file:
    cistopic_obj = pickle.load(file)


# In[13]:


with open('/Users/gclyu07/Desktop/XMAS_analysis/outputs/models.pkl', "rb") as file:
    models = pickle.load(file)


# In[14]:


from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models,
    select_model = 40,
    return_model = True
)


# In[15]:


cistopic_obj.add_LDA_model(model)


# In[16]:


pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_modeled.pkl"), "wb")
)


# # Clustering and visualization 

# In[16]:


from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)


# In[17]:


cistopic_obj.cell_names = cistopic_obj.cell_data['ori_barcode'].to_list()


# In[18]:


cistopic_obj.cell_names


# In[19]:


find_clusters(
    cistopic_obj,
    target  = 'cell',
    k = 10,
    res = [0.6, 1.2, 3],
    prefix = 'pycisTopic_',
    scale = True,
    split_pattern = '-'
)


# In[20]:


cistopic_obj.projections["cell"] = {}


# In[21]:


run_umap(
    cistopic_obj,
    target  = 'cell', scale=True)


# In[22]:


run_tsne(
    cistopic_obj,
    target  = 'cell', scale=True)


# In[23]:


plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['orig.ident', 'bimod_oc','bimod_oc_des',
               'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
    target='cell', num_columns=3,
    text_size=10,
    dot_size=5)


# In[24]:


plot_topic(
    cistopic_obj,
    reduction_name = 'UMAP',
    target = 'cell',
    num_columns=5
)

cell_topic_heatmap(
    cistopic_obj,
    variables = ['orig.ident'],
    color_dictionary={"D11":"#F0AD4E", "D16":"#D9534F", 
                      "D28":"#428BCA", "D42":"#9933CC", "D56":"#66CCCC"},
    scale = False,
    legend_loc_x = 1.0,
    legend_loc_y = -1.2,
    legend_dist_y = -1,
    figsize = (10, 10)
)
# # Topic binarization & QC

# In[25]:


from pycisTopic.topic_binarization import binarize_topics


# In[26]:


region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, num_columns=5
)


# In[27]:


region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)


# In[28]:


binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)


# In[29]:


from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img


# In[30]:


topic_qc_metrics = compute_topic_metrics(cistopic_obj)


# In[31]:


fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)


# In[32]:


# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.show()


# In[33]:


topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='bimod_oc_des',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)


# In[34]:


topic_annot


# # Differentially Accessible Regions (DARs) 

# In[35]:


from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np


# In[36]:


imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)


# In[37]:


normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)


# In[38]:


variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)


# In[39]:


len(variable_regions)


# In[40]:


markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='bimod_oc_des',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir='/tmp/',
    split_pattern = '-'
)


# In[41]:


from pycisTopic.clust_vis import plot_imputed_features


# In[42]:


plot_imputed_features(
    cistopic_obj,
    reduction_name='UMAP',
    imputed_data=imputed_acc_obj,
    features=[markers_dict[x].index.tolist()[0] for x in ['ProgFP Development', 'NbM', 'DA Neuron projection', 
                                                          'DA Neurotransmitter release','DA Synaptic assembly','DA Synaptic modulation']],
    scale=False,
    num_columns=3
)


# In[43]:


print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")


# # Save region sets

# In[44]:


os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok = True)


# In[45]:


from pycisTopic.utils import region_names_to_coordinates


# In[46]:


for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )


# In[47]:


for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )


# In[48]:


for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )

pickle.dump(
    imputed_acc_obj,
    open(os.path.join(out_dir, "cistopic_obj_imputed.pkl"), "wb")
)
# # Gene activity

# In[49]:


import pyranges as pr
from pycisTopic.gene_activity import get_gene_activity


# In[50]:


chromsizes = pd.read_table(os.path.join(out_dir, "qc", "hg38.chrom_sizes_and_alias.tsv"))
chromsizes


# In[51]:


chromsizes.rename({"# ucsc": "Chromosome", "length": "End"}, axis = 1, inplace = True)
chromsizes["Start"] = 0
chromsizes = pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]])


# In[52]:


chromsizes


# In[53]:


pr_annotation = pd.read_table(
        os.path.join(out_dir, "qc", "tss.bed")
    ).rename(
        {"Name": "Gene", "# Chromosome": "Chromosome"}, axis = 1)
pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
pr_annotation = pr.PyRanges(pr_annotation)
pr_annotation


# In[54]:


gene_act, weigths = get_gene_activity(
    imputed_acc_obj,
    pr_annotation,
    chromsizes,
    use_gene_boundaries=True, # Whether to use the whole search space or stop when encountering another gene
    upstream=[1000, 100000], # Search space upstream. The minimum means that even if there is a gene right next to it
                             # these bp will be taken (1kbp here)
    downstream=[1000,100000], # Search space downstream
    distance_weight=True, # Whether to add a distance weight (an exponential function, the weight will decrease with distance)
    decay_rate=1, # Exponent for the distance exponential funciton (the higher the faster will be the decrease)
    extend_gene_body_upstream=10000, # Number of bp upstream immune to the distance weight (their value will be maximum for
                          #this weight)
    extend_gene_body_downstream=500, # Number of bp downstream immune to the distance weight
    gene_size_weight=False, # Whether to add a weights based on the length of the gene
    gene_size_scale_factor='median', # Dividend to calculate the gene size weigth. Default is the median value of all genes
                          #in the genome
    remove_promoters=False, # Whether to remove promoters when computing gene activity scores
    average_scores=True, # Whether to divide by the total number of region assigned to a gene when calculating the gene
                          #activity score
    scale_factor=1, # Value to multiply for the final gene activity matrix
    extend_tss=[10,10], # Space to consider a promoter
    gini_weight = True, # Whether to add a gini index weigth. The more unique the region is, the higher this weight will be
    return_weights= True, # Whether to return the final weights
    project='Gene_activity') # Project name for the gene activity object


# In[55]:


DAG_markers_dict= find_diff_features(
    cistopic_obj,
    gene_act,
    variable='bimod_oc_des',
    var_features=None,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir='/tmp',
    split_pattern = '-')


# In[56]:


print("Number of DAGs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(DAG_markers_dict[x])}")


# In[57]:


plot_imputed_features(
    cistopic_obj,
    reduction_name='UMAP',
    imputed_data=gene_act,
    features=['FOXA2', 'LMX1A', 'OTX2', 'EN1', 
              'SOX6', 'TH', 'DDC',
              'NR4A2', 'NEUROD1', 'CALB1', 'LMO3', 
              'BNC2', 'PITX3', 'PBX1'],
    scale=True,
    num_columns=3
)


# In[58]:


pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_final.pkl"), "wb")
)


# In[59]:


pickle.dump(
    imputed_acc_obj,
    open(os.path.join(out_dir, "cistopic_obj_imputed_final.pkl"), "wb")
)


# In[ ]:





# # Label transfer (Skipped) 

# # Exporting to loom

# In[99]:


from pycisTopic.loom import export_region_accessibility_to_loom, export_gene_activity_to_loom


# In[100]:


cluster_markers = {'bimod_oc_des': markers_dict}


# In[101]:


os.makedirs(os.path.join(out_dir, "loom"), exist_ok=True)


# In[102]:


export_region_accessibility_to_loom(
    accessibility_matrix = imputed_acc_obj,
    cistopic_obj = cistopic_obj,
    binarized_topic_region = region_bin_topics_otsu,
    binarized_cell_topic = binarized_cell_topic,
    selected_cells = cistopic_obj.projections['cell']['UMAP'].index.tolist(),
    out_fname = os.path.join(out_dir, "loom", "XMAS_pycisTopic_region_accessibility.loom"),
    cluster_annotation = ['bimod_oc_des'],
    cluster_markers = cluster_markers,
    tree_structure = ('XMAS', 'pycisTopic', 'noDBL_all'),
    title = 'Region accessibility all',
    nomenclature = "hg38",
    split_pattern = '-'
)


# In[103]:


export_gene_activity_to_loom(
    gene_activity_matrix = gene_act,
    cistopic_obj = cistopic_obj,
    out_fname = os.path.join(out_dir, "loom", "XMAS_pycisTopic_gene_activity.loom"),
    cluster_annotation = ['bimod_oc_des'],
    cluster_markers = cluster_markers,
    tree_structure = ('XMAS', 'pycisTopic', 'ATAC'),
    title = 'Gene activity',
    nomenclature = "hg38",
    split_pattern = '-'
)

