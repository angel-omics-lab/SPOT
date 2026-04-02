## Tester script for pca and umap generation 

import numpy as np 
import pandas as pd 
import anndata as ad 
import scanpy as sc
import json 
from sklearn.preprocessing import scale
import matplotlib.pyplot as plt 

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

good_peptides = [758.4519, 797.4264, 976.4517, 999.5946, 1060.5269, 
1068.5069, 1082.6317, 1084.5018, 1089.4847, 1098.5062, 
1098.5902, 1102.564, 1115.544, 1125.5283, 1128.528, 
1128.532, 1142.5073, 1154.5073, 1166.5324, 1175.5473]



### Extract region classes from json 
# try:
#     json_dict = json.load(DATA_PATH)
#     labels = json_dict['roi_labels']
#     print(type(json_dict))
#     print(labels)
# except FileNotFoundError:
#     print('Error: the .json file was not found. Check file path.')
# except json.JSONDecodeError:
#     print('Error: Failed to decode json from the file. Check file formatting.')




# Construct anndata object
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
    obs = pixel_metadata,       # df w columns: sample, class, x, y (one row per pixel)
    var = peptide_metadata      # df with peptide names as index (one row per peptide)--just column titles
)


# Generate pca and umap 
NUM_PC = 5      # Number of principal components to search for in pca and passed to nearest neighbor graph 
# # NOTE should use some sort of heuristic here rather than a constant integer

print('Running PCA...')
sc.pp.pca(ann_obj, n_comps = NUM_PC)     # preprocess data for pca, n_comps give number of PCs to look for

print('Running UMAP analysis...')
sc.pp.neighbors(ann_obj, use_rep='X_pca', n_pcs=NUM_PC)       # generate nearest neighbor graph needed for umap
sc.tl.umap(ann_obj)    

# # Do actual plotting
# print('Generating PCA and UMAP figures...')
# sc.pl.pca(ann_obj, color='class')
# sc.pl.umap(ann_obj, color='class')




# Slingshot pseudotime
from pyslingshot import Slingshot
print('Running pseudotime reconstruction. This may take a while...')
# Step 1: Clustering 
sc.tl.leiden(ann_obj, resolution=0.15)       # resolution up - more clusters V down = fewer clusters, ideally want # clusters to = # disease stages

# Step 2: Set root (identify cluster that Normal cells belong to)
normal_mask = ann_obj.obs['class']=='Normal'        
normal_clusters = ann_obj.obs.loc[normal_mask, 'leiden']
root_cluster = int(normal_clusters.mode()[0])    

# Step 3: Run 
ann_obj.uns['iroot'] = np.flatnonzero(ann_obj.obs['class']=='Normal')[0]    # Assign root
slingshot = Slingshot(
    ann_obj, 
    celltype_key='leiden',
    obsm_key='X_umap', 
    start_node=root_cluster
)
slingshot.fit()     # pseudotime values are stored in ann_obj.obs
print('Pseudotime fitting successful. Now proceeding to plotting...')

# Step 4: Plot 
ann_obj.obs['pseudotime'] = slingshot.unified_pseudotime    

# 4a: Color UMAP based on pseudotime
sc.pl.umap(ann_obj, color='pseudotime', show=True, color_map='RdBu')

# 4b: Overlay slingshot curves on UMAP
# Map classes to colors so scatterplot UMAP can be colored like original  
class_colors = {
    'Normal': 'green',
    'DCIS': 'dodgerblue',
    'IBC': 'orange'
}
point_colors = ann_obj.obs['class'].map(class_colors).values

# Plot UMAP as scatterplot
fig, ax = plt.subplots()
scatter = ax.scatter(
    ann_obj.obsm['X_umap'][:,0], 
    ann_obj.obsm['X_umap'][:,1],
    c=point_colors, 
    s=5
)
plt.colorbar(scatter, ax=ax)

# Overlay slingshot curves
for curve in slingshot.curves:
    points = curve.points
    ax.plot(points[:,0], points[:,1], color='black', linewidth=2)
ax.set_xlabel('UMAP1')
ax.set_ylabel('UMAP2')
plt.show()


# 4c: Add ROI centroids to umap with curves


print('Pseudotime reconstruction complete.')

