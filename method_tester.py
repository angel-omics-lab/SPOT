import pandas as pd 
import numpy as np 
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests 
from matplotlib import pyplot as plt
import os
import seaborn as sns
import math 
from scipy.cluster.hierarchy import dendrogram, linkage 

# plt.rcParams(font='Arial')

DATA_PATH = r'C:\Users\AngelLab\Documents\GitHub\spatial-proteomics-analyzer\fake_data_for_testing.xlsx'
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

###### load_and_process ######
# print('Filtering ROIs...')
# for region, label in roi_labels.items():
#     data = pd.read_excel(DATA_PATH, sheet_name= region)
#     intensities = data.iloc[:, 4:]      # Skip first 4 columns 
#     print(f'Processing {region}...')
#     intensities = intensities.fillna(0)        # Change all NaN values to 0 
#     zeros = (intensities==0).astype(int).sum().sum()      # Number of 0 intensities
#     total = intensities.shape[0]*data.shape[1]     # Total number of intensities
            
#     if zeros > (0.25*total):       # Remove ROI from list if >25% of peptides are 0 
#         del roi_labels[region]
#         print(f'Removed {region} with {label} label from list: >25% zero intensities.')
#     else: 
#         print(f'{region} accepted')
#         print(f'Normalizing intensities in {region}')
#         data.iloc[:, 4:] = np.log(data.iloc[:, 4:])  # Natural log intensities
# print('ROI filtering complete. Filtered list: ', roi_labels.keys())


##### compute_peptide_sparsity #####
# print('Filtering peptides...')
# peptide_list = [758.4519, 797.4264, 976.4517]
# zero_counts_per_peptide = pd.Series(dtype=int)      # Will track how many ROIs each peptide 'failed' in 

# for region in roi_labels: 
#     intensities = pd.read_excel(DATA_PATH, sheet_name=region).iloc[:, 4:]
#     zero_pct = (intensities==0).mean()      # Percentage of 0 intensities for each peptide 
#     zero_counts_per_peptide = zero_counts_per_peptide.add((zero_pct > 0.20).astype(int), fill_value=0)      # Label peptide for region with 0 if > 20% of values are 0

# # Keep peptides above 10% threshold
# threshold = len(roi_labels)*0.1
# good_peptides = zero_counts_per_peptide[zero_counts_per_peptide <= threshold].index.tolist()


# print('Peptide filtering complete. Filtered list: ', good_peptides)


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
print('Generating heatmaps for each peptide...')
# Calculate spatial centroid of each roi
    # Concatenate all ROI sheets, tagging each row with its ROI label
combined = pd.concat(
    [df.assign(roi=roi) for roi, df in data.items() if roi in roi_labels],
    ignore_index=True
)
    # Groupby computes centroid AND all peptide means in one pass (vector operation)
agg_dict = {'x': 'mean', 'y': 'mean', **{p: 'mean' for p in good_peptides}}
roi_stats = combined.groupby('roi').agg(agg_dict).reset_index()     # roi_stats columns: ['roi', 'x', 'y', peptide_1, peptide_2, ...]
roi_stats['class'] = roi_stats['roi'].map(roi_labels)
print(roi_stats)
print(roi_stats.dtype())
        
# # Generate heatmaps for each peptide
# for peptide in good_peptides:
#     heatmap = plt.figure()
#     plt.scatter(roi_stats['x'], roi_stats['y'], 
#                 c=roi_stats[peptide], cmap='RdBu_r', s=100)
#     plt.colorbar()
#     plt.xticks([])
#     plt.yticks([])
#     plt.title(peptide, fontsize=16)
#     plt.show()
# print('Spatial heatmap generation successful. Saved to: ', os.path.dirname(DATA_PATH))
# print('ROI stats:', roi_stats)


# ### Generate boxplots 
# print('Generating box plots for significant peptides...')
# n_peptides = len(good_peptides)
# n_cols = 5
# n_rows = math.ceil(n_peptides/n_cols)

# fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols*4, n_rows*4))        # Grid of box plots
# axs = axs.flatten()     # Flatten so indexing is easier
# # Give overall plot name, axis, legend 
# for i, peptide in enumerate(good_peptides):
#     plot_data = (
#         [{'Intensity': val, 'Class':'DCIS'}
#          for region in roi_labels  
#          if roi_labels[region] == 'DCIS'
#          for val in [data[region][peptide].mean() for region in roi_labels if roi_labels[region] == 'DCIS' ]]
#         +
#         [{'Intensity': val, 'Class':'IBC'} for val in [data[region][peptide].mean() for region in roi_labels if roi_labels[region] == 'IBC' ]]
#     )
#     plot_df = pd.DataFrame(plot_data)

#     sns.boxplot(data=plot_df, x='Class', y='Intensity', 
#                 hue='Class', palette={'DCIS':'blue', 'IBC':'orange'},
#                 ax=axs[i], legend=False
#             )
#     axs[i].set_title(peptide)

# # Hide unused subplots
# for j in range(i+1, len(axs)): axs[j].set_visible(False)
            
# # Shared legend
# labels = [plt.Rectangle((0,0),1,1, color=c) for c in ['blue', 'orange']]
# fig.legend(labels, ['DCIS', 'IBC'], loc='upper right')

# # Shared title, yaxislabel
# fig.suptitle('Differentially Expressed Peptides (DCIS v IBC)')
# fig.supylabel('ln(Mean Intensity)')
# fig.tight_layout()

# plt.savefig(os.path.join(os.path.dirname(DATA_PATH), 'sig_peptide_boxplots.png'))
# plt.show()

# print(good_peptides)


# #### plot spatial dot plot for region types #####
# print('Generating spatial dot plot for regions...')
# # Calculate spatial centroid of each roi
#     # Concatenate all ROI sheets, tagging each row with its ROI label
# combined = pd.concat(
#     [df.assign(roi=roi) for roi, df in data.items() if roi in roi_labels],
#     ignore_index=True
# )
#     # Groupby computes centroid (vector operation)
# agg_dict = {'x': 'mean', 'y': 'mean'}
# roi_stats = combined.groupby('roi').agg(agg_dict).reset_index()     # roi_stats columns: ['roi', 'x', 'y', peptide_1, peptide_2, ...]
# roi_stats['class'] = roi_stats['roi'].map(roi_labels)
# roi_stats['roi'] = roi_stats['roi'].str.removeprefix('ROI_')    # remove prefix so label can fit inside dot
        
#     # Generate colored scatter plot
# sns.scatterplot(data=roi_stats, x='x', y='y', 
#                  hue='class', palette={'DCIS':'dodgerblue', 'IBC':'orange', 'Normal':'green'},
#                  s=250
# )
# for x, y, roi in zip(roi_stats['x'], roi_stats['y'], roi_stats['roi']):
#     plt.text(x, y, str(roi),
#              ha='center', va='center', 
#              fontsize=8, color='black')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.tight_layout()
# plt.xticks([])
# plt.yticks([])
# plt.xlabel(None)
# plt.ylabel(None)
# plt.title('Spatial ROI Plot')

# #plt.show()

# print('ROI stats:', roi_stats)


plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['mediumorchid', 'lightcoral', 'lightblue'])
print('Generating hierarchical clusters...')
# Reorder roi_stats based on class so they can be colored properly in the dendrogram
roi_stats = roi_stats.sort_values('class')
# Extract peptide intensities per roi 
feature_matrix = roi_stats[good_peptides].values
# Create dendrogram
print('Constructing dendrogram for visualization...')
linkage_data = linkage(feature_matrix, method='ward', metric='euclidean')
dend = dendrogram(linkage_data, labels=roi_stats['roi'].values, leaf_rotation=90)

# Formatting
color_map = {'DCIS': 'dodgerblue', 'IBC': 'orange', 'Normal': 'green'}
label_colors = [color_map[c] for c in roi_stats.set_index('roi').loc[roi_stats['roi']]['class']]
leaf_order = dend['ivl']  # ROIs in dendrogram order
label_colors = [color_map[roi_stats.set_index('roi').loc[roi, 'class']] for roi in leaf_order]
ax = plt.gca()
for tick, color in zip(ax.get_xticklabels(), label_colors):
    tick.set_color(color)
plt.yticks([])
# Create legend
class_labels = [plt.Rectangle((0,0),1,1, color=c) for c in ['dodgerblue', 'orange', 'green']]
plt.legend(class_labels, ['DCIS', 'IBC', 'Normal'], bbox_to_anchor=(1.05, -0.25), loc= 'lower left')

plt.title('ROI clusters based on peptide profile similarity')
plt.show()
