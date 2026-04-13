import pandas as pd 
import numpy as np 
import scipy.stats as stats
#from statsmodels.stats.multitest import multipletests 
from matplotlib import pyplot as plt
import os
import seaborn as sns
import math 
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import scale
import anndata as ad
import scanpy as sc
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx


# plt.rcParams(font='Arial')

# DATA_PATH = r'C:\Users\AngelLab\Documents\GitHub\spatial-proteomics-analyzer\fake_data_for_testing.xlsx'
DATA_PATH = '/Users/bryngerding/Documents/GitHub/spatial-proteomics-analyzer/fake_data_for_testing.xlsx'
data = pd.read_excel(DATA_PATH, sheet_name=None)

roi_labels = {
    'ROI_101':'DCIS', 
    'ROI_102':'DCIS', 
    'ROI_103':'DCIS', 
    'ROI_104':'IBC', 
    'ROI_105':'DCIS', 
    'ROI_106':'IBC', 
    'ROI_107':'DCIS', 
    'ROI_108':'IBC',
    'ROI_109':'IBC', 
    'ROI_110':'IBC',
    'ROI_128':'Normal'  
}

good_peptides = [758.4519, 797.4264, 976.4517]
good_peptides = [758.4519, 797.4264, 976.4517, 999.5946, 1060.5269, 
1068.5069, 1082.6317, 1084.5018, 1089.4847, 1098.5062, 
1098.5902, 1102.564, 1115.544, 1125.5283, 1128.528, 
1128.532, 1142.5073, 1154.5073, 1166.5324, 1175.5473] 


# ###### load_and_process ######
# print('Filtering ROIs...')
# to_remove=[]
# for region, label in roi_labels.items():
#     data = data[region]
#     intensities = data.iloc[:, 4:]      # Skip first 4 columns 
#     print(f'Processing {region}...')
#     intensities = intensities.fillna(0)        # Change all NaN values to 0 
#     zeros = (intensities==0).astype(int).sum().sum()      # Number of 0 intensities
#     total = intensities.shape[0]*data.shape[1]     # Total number of intensities
            
#     if zeros > (0.25*total):       # Remove ROI from list if >25% of peptides are 0 
#         to_remove.append(region)
#         print(f'Removed {region} with {label} label from list: >25% zero intensities.')
#     else: 
#         print(f'{region} accepted')
#         # print(f'Normalizing intensities in {region}')
#         # data.iloc[:, 4:] = np.log(data.iloc[:, 4:])  # Natural log intensities
#     for region in to_remove:
#         del roi_labels[region]
# print('ROI filtering complete.')


# ##### compute_peptide_sparsity #####
# print('Filtering peptides...')
# starting_length = len(good_peptides)
# zero_counts_per_peptide = pd.Series(dtype=int)      # Will track how many ROIs each peptide 'failed' in 

# for region in roi_labels: 
#     intensities = pd.read_excel(DATA_PATH, sheet_name=region).iloc[:, 4:]
#     zero_pct = (intensities==0).mean()      # Percentage of 0 intensities for each peptide 
#     zero_counts_per_peptide = zero_counts_per_peptide.add((zero_pct > 0.20).astype(int), fill_value=0)      # Label peptide for region with 0 if > 20% of values are 0

# # Keep peptides above 10% threshold
# threshold = len(roi_labels)*0.1
# good_peptides = zero_counts_per_peptide[zero_counts_per_peptide <= threshold].index.tolist()

# ending_length = len(good_peptides)

# print(f'Peptide filtering complete. {starting_length - ending_length} peptides filtered out.')


# ##### differential expression test #####
# results_list = []

# print('Identifying differentially expressed peptides...')
# for peptide in good_peptides:
#     group1 = [data[region][peptide].mean() for region in roi_labels if roi_labels[region] == 'DCIS']     # Mean intensities of peptide in DCIS (class 1) rois
#     group2 = [data[region][peptide].mean() for region in roi_labels if roi_labels[region] == 'IBC']      # Mean Intensities of peptide in IBC (class 2) rois
#     # Run KW
#     H_statistic, p_value = stats.kruskal(group1, group2)
#     # Add to data frame 
#     results_list.append({
#     'peptide': peptide, 
#     'pvalue_ols': p_value
#     })
# results_df = pd.DataFrame(results_list) 
# # Run BH FDR, add q val to df 
# reject, q_values, _, _ = multipletests(results_df['pvalue_ols'], alpha=0.05, method='fdr_bh')
# results_df['qvalue_bh'] = q_values
# # Save peptides whose q value is <= 0.05
# good_peptides = good_peptides = results_df[results_df['qvalue_bh'] <= 0.05]['peptide'].tolist()
# print('Differential analysis complete. Significant peptides: ', good_peptides)



##### plot spatial heatmaps ######
# print('Generating heatmaps for each peptide...')
# # Calculate spatial centroid of each roi
#     # Concatenate all ROI sheets, tagging each row with its ROI label
#     combined = pd.concat(
#         [df.assign(roi=roi) for roi, df in data.items() if roi in roi_labels],
#         ignore_index=True
#     )
#     # Groupby computes centroid AND all peptide means in one pass (vector operation)
#     agg_dict = {'x': 'mean', 'y': 'mean', **{p: 'median' for p in good_peptides}}
#     roi_stats = combined.groupby('roi').agg(agg_dict).reset_index()     # roi_stats columns: ['roi', 'x', 'y', peptide_1, peptide_2, ...]
#     roi_stats['class'] = roi_stats['roi'].map(roi_labels)

#     print('ROI calculations successful!')
#     return roi_stats

# roi_stats = get_roi_stats()



# def generate_boxplots_new(peptide):
#     plot_data = (
#         [{'Intensity': val, 'Class':'DCIS'} for val in [data[region][peptide].mean() for region in roi_labels if roi_labels[region] == 'DCIS' ]]
#         +
#         [{'Intensity': val, 'Class':'IBC'} for val in [data[region][peptide].mean() for region in roi_labels if roi_labels[region] == 'IBC' ]]
#     )
#     plot_df = pd.DataFrame(plot_data)
    
#     box = sns.catplot(
#         data=plot_df,
#         x='Class', y='Intensity',
#         hue='Class', 
#         palette={'DCIS':'dodgerblue', 'IBC':'orange'},
#         kind='box'
#     )
#     sns.swarmplot(
#         data=plot_df,
#         x='Class', y='Intensity', 
#         color='black'
#     )
#     plt.plot(figsize=(2.5,4))
#     plt.show()



# from sklearn.inspection import permutation_importance
# def get_random_forest_ranking(roi_stats):
#     good_peptides_str = [str(p) for p in good_peptides]
#     roi_stats_str = roi_stats.copy()
#     roi_stats_str.columns = [str(c) for c in roi_stats_str.columns]

#     X = roi_stats_str[good_peptides_str]
#     y = roi_stats_str['class']
#     X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    
#     print('Fitting forest model...')
#     model = RandomForestClassifier(n_estimators=1000, random_state=42)
#     model.fit(X_train, y_train)

#     # Feature importance based on mean decrease in impurity 
#     importances = model.feature_importances_
#     variable_importances = pd.Series(importances, index=good_peptides_str)

#     # Feature importance based on feature permutation
#     print('Calculating permutation importance per variable...')
#     # result = permutation_importance(model, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2)
#     # variable_importances = pd.Series(result.importances_mean, index=good_peptides_str)
#     variable_importances_sorted = variable_importances.sort_values()

#     # std_sorted = pd.Series(result.importances_std, index=good_peptides_str).loc[
#     #     variable_importances_sorted.index
#     # ]

#     print('Plotting...')
#     fig, ax = plt.subplots(figsize=(8, len(good_peptides) * 0.4 + 1))
#     variable_importances_sorted.plot.barh(
#         # xerr=std_sorted,
#         ax=ax,
#         color='mediumorchid',
#         ecolor='gray'
#     )

#     ax.set_title('Ranking of peptide importance in random forest classification model')
#     ax.set_xlabel('Feature importance')
#     ax.set_ylabel('Peptide m/z')
#     fig.tight_layout()
#     plt.show()



# print('Generating random forest models per peptide...')
# oob_list = []
# good_peptides = [str(p) for p in good_peptides]
# roi_stats.columns = [str(c) for c in roi_stats.columns]
# for peptide in good_peptides: 
#     X= roi_stats[[peptide]]         # data
#     y= roi_stats['class']         # labels
# # Randomly split the dataset so 20% is reserved for model validation
#     X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state=42)

#     forest_model = RandomForestClassifier(n_estimators=100, oob_score=True, random_state=42) 
#     forest_model.fit(X_train, y_train)
#     oob_list.append(round(forest_model.oob_score_, 3))

# # Plot each peptide oob on bar plot
# print('Constructing bar plot with each peptide oob score...')
#     # Sort them so they will be in ascending order in the plot
# sorted_pairs = sorted(zip(oob_list, good_peptides), key=lambda x: x[0])
# sorted_scores, sorted_peptides = zip(*sorted_pairs)
# plt.barh(sorted_peptides, sorted_scores, color='mediumorchid')

# plt.title('Peptides ranked by out-of-bag score in a random forest model')
# plt.ylabel('Peptide')
# plt.xlabel('OOB score')
# plt.xlim(0,1)
# plt.show()

########### pixel-level analysis ##############
print('Constructing anndata object...')


 # Concat ROI sheets into a single flat df
combined = pd.concat([df.assign(ROI=roi) for roi, df in data.items() if roi in roi_labels], ignore_index=True)
        
# Extract and scale peptide intensity matrix 
combined_intensities = combined[good_peptides].values.astype(float)
combined_intensities = scale(combined_intensities)      # z-score so peptides contribute equally to PCA
print('data keys:', list(data.keys()))
print('roi_labels keys:', list(roi_labels.keys()))
print('data is empty:', len(data)==0)

matches = [roi for roi in data.keys() if roi in roi_labels]
print("Matching keys:", matches)
print("Number of matches:", len(matches))

combined = pd.concat(
    [df.assign(ROI=roi) for roi, df in data.items() if roi in roi_labels],
    ignore_index=True
)

combined_intensities = combined[good_peptides].values.astype(float)
combined_intensities = scale(combined_intensities)      # z-score scale intentities so peptide contributions are equalized for PCA

pixel_metadata =  pd.DataFrame({
    'sample': combined['ROI'].values,
    'class': combined['ROI'].map(roi_labels).values, 
    'x': combined['x'].values.astype(float),
    'y': combined['y'].values.astype(float), 
    })
           
peptide_metadata = pd.DataFrame(
    {'peptide' : good_peptides}, 
    index=good_peptides
)          

# Anndata object (observations x variables) created -- this is the assumed format for scanpy function arguments
ann_obj = ad.AnnData(
    X = combined_intensities,   # shape: (n_pixels, n_peptides)
    obs = pixel_metadata,       # df w columns: sample, group, x, y (one row per pixel)
    var = peptide_metadata      # df with peptide names as index (one row per peptide)--just column titles
)

# NUM_PC = 5      # Number of principal components to search for in pca and passed to nearest neighbor graph 
# # NOTE should use some sort of heuristic here rather than a constant integer

# print('Running PCA...')
# sc.pp.pca(ann_obj, n_comps = NUM_PC)     # preprocess data for pca, n_comps give number of PCs to look for

# print('Running UMAP analysis...')
# sc.pp.neighbors(ann_obj, use_rep='X_pca', n_pcs=NUM_PC)       # generate nearest neighbor graph needed for umap
# sc.tl.umap(ann_obj)    

# # Do actual plotting
# print('Generating PCA and UMAP figures...')
# sc.pl.pca(ann_obj, color='group')
# sc.pl.umap(ann_obj, color='group')







### Run SCIMITAR ###
import scimitar.models 
import scimitar.plotting
import scimitar.morphing_mixture as mm 
import scimitar.differential_analysis
from collections import defaultdict



    ###### Metastable state graph #######
# metastable_graph : graph object whose nodes are GMM components/metastable states fit to full dataset. Each node carries the Gaussian parameters (mean vector and covariance) for that state
# bootstrap_replicates : raw record of all bootstrapping runs 
# edge_fractions : derived from replicates, for each pair of states, gives fraction of runs that found connection sig (i.e. confidence score of edge, also sig if assessed by Cohen's d)
# NOTE GMM is fit once to the full data to define states, then bootstrapping is done to get edge confidence 
print('Constructing metastable graph...')
metastable_graph, bootstrap_replicates, edge_fractions = scimitar.models.get_gmm_bootstrapped_metastable_graph(
    ann_obj.X, 
    n_boot = 20, 
    covariance_type = 'diag'
)
metastable_graph.edge_weights = edge_fractions     # store confidence scores on graph object--done so plotting function can later use thickness to show how supported the edge is

# Plot and capture outputs
print('Plotting metastable graph...')
plt.figure()
state_colors, embedding = scimitar.plotting.plot_metastable_graph(
    ann_obj.X, 
    metastable_graph, 
    edge_weights=edge_fractions
)

import numpy as np 
from matplotlib import pyplot as plt
from pyslingshot import Slingshot
from anndata import AnnData

load = True
num_cells = 1000
num_dims_reduced = 2
num_branches = 1

K=10    # cluster labels
filename = r'C:\Users\AngelLab\Downloads\fakedata-1branch.npy'
start_node = 4

if load:
    data = np.load(filename, allow_pickle=True).item()
    cluster_labels = data['cluster_labels']
    data = data['data']

#plt.scatter(data[:,0], data[:,1], c=cluster_labels)
#plt.show()

num_genes = 500
ad = AnnData(np.zeros((num_cells, num_genes)))
ad.obsm['X_umap'] = data
ad.obs['celltype'] = cluster_labels
ad

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,10))
custom_xlim = (-12,12)
custom_ylim = (-12,12)

slingshot = Slingshot(ad, celltype_key='celltype', obsm_key='X_umap', start_node=start_node, is_debugging='verbose')
slingshot.fit(num_epochs=1, debug_axes=axes)

fig,axes = plt.subplots(ncols=2, figsize=(12,4))
axes[0].set_title('Cluster')
axes[1].set_title('Pseudotime')
slingshot.plotter.curves(axes[0],slingshot.curves)
slingshot.plotter.clusters(axes[0], labels=np.arange(slingshot.num_clusters), s=4, alpha=0.5)
slingshot.plotter.clusters(axes[1], color_mode='pseudotime', s=5 )
plt.show()




    
# Do actual plotting
# print('Generating PCA and UMAP figures...')
# sc.pl.pca(ann_obj, color='class')
# sc.pl.umap(ann_obj, color='class')


##### generate MST ######
# print('Computing MST...')
# # Compute roi centroids in pca space
# n_comps = ann_obj.obsm['X_pca'].shape[1]
# pca_df = pd.DataFrame(
#     ann_obj.obsm['X_pca'], 
#     index = ann_obj.obs_names, 
#     columns=[f'PC{i+1}' for i in range(n_comps)]
# )
# pca_df['sample']=ann_obj.obs['sample'].values
# centroids_pca = pca_df.groupby('sample').mean()
        
# # Build mst on centroids
# dist_matrix = cdist(centroids_pca.values, centroids_pca.values, metric='euclidean')     # Gets Eucl. distance between every ROI centroid pair in PCA space
# mst_sparse = minimum_spanning_tree(dist_matrix)         # Actual MST algorithm -- finding edge subset that minimizes total distance, outputs sparse matrix identifying those edges
        
# graph = nx.from_scipy_sparse_array(mst_sparse)      # Converts sparse matrix into graph object 
# roi_names = list(centroids_pca.index)
# mst_graph = nx.relabel_nodes(graph, {i: name for i, name in enumerate(roi_names)})      # Swaps default node names for actual ROIs
        
# print(f'MST built with {mst_graph.number_of_nodes()} nodes and {mst_graph.number_of_edges()} edges')
        
# # Add MST data to ann_obj
# ann_obj.uns['mst'] = {
#     'graph': mst_graph, 
#     'centroids' : centroids_pca
# }


# ##### plot UMAP w MST overlaid #####
# mst_graph = ann_obj.uns['mst']['graph']
# umap_coords = ann_obj.obsm['X_umap']

# # Project each roi centroid into umap space
# umap_df = pd.DataFrame(umap_coords, columns=['UMAP1', 'UMAP2'])
# umap_df['sample'] = ann_obj.obs['sample'].values
# umap_df['class'] = ann_obj.obs['class'].values

# centroids_umap = umap_df.groupby('sample')[['UMAP1', 'UMAP2']].mean()
# centroids_umap['class'] = umap_df.groupby('sample')['class'].first()


# # Build edge lines from mst graph
# edge_lines = []
# for u, v in mst_graph.edges():
#     if u in centroids_umap.index and v in centroids_umap.index:
#         edge_lines.append({
#             'x': [centroids_umap.loc[u, 'UMAP1'], centroids_umap.loc[v, 'UMAP1']], 
#             'y': [centroids_umap.loc[u, 'UMAP2'], centroids_umap.loc[v, 'UMAP2']],
#         })

# # Plot
# color_map = {'DCIS': 'dodgerblue', 'IBC': 'orange', 'Normal': 'green'}
# fig, ax = plt.subplots(figsize=(8,6))

# #Pixel cloud
# for cls, grp in umap_df.groupby('class'):
#     ax.scatter(grp['UMAP1'], grp['UMAP2'], 
#                c=color_map[cls], s=1, alpha=0.3, label=cls, rasterized=True)
# # MST edges
# for edge in edge_lines:
#     ax.plot(edge['x'], edge['y'], 'k-', linewidth=1.5, zorder=3)
# # ROI centroids
# for roi, row in centroids_umap.iterrows():
#     ax.scatter(row['UMAP1'], row['UMAP2'], 
#                c=color_map[row['class']], s=80, edgecolors='black', 
#                linewidths=0.8, zorder=4)
#     ax.annotate(roi.replace('ROI_', ''), (row['UMAP1'], row['UMAP2']), 
#                 fontsize=6, ha='center', va='bottom', 
#                 xytext=(0,5), textcoords='offset points')
# # Formatting
# ax.legend(markerscale=5, framealpha=0.7)
# ax.set_xlabel('UMAP 1')
# ax.set_ylabel('UMAP 2')
# ax.set_title('Pixel UMAP with MST trajectory')
# plt.tight_layout()
# #plt.savefig(os.path.join.dirname(data_path), 'umap_mst.png')
# plt.show()


