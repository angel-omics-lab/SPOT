# import pandas as pd 
# from matplotlib import pyplot as plt
# import seaborn as sns
# from sklearn.model_selection import train_test_split
# from sklearn.ensemble import RandomForestClassifier




# #plt.rcParams(font='Arial')
# DATA_PATH = r'F:\Bryn\DCIS_Souvik\dcis_data_full.xlsx'
# # DATA_PATH = r'C:\Users\AngelLab\Documents\GitHub\spatial-proteomics-analyzer\fake_data_for_testing.xlsx'
# # DATA_PATH = '/Users/bryngerding/Documents/GitHub/spatial-proteomics-analyzer/fake_data_for_testing.xlsx'
# data = pd.read_excel(DATA_PATH, sheet_name=None)

# # roi_labels = {
# #     'ROI_101':'DCIS', 
# #     'ROI_102':'DCIS', 
# #     'ROI_103':'DCIS', 
# #     'ROI_104':'IBC', 
# #     'ROI_105':'DCIS', 
# #     'ROI_106':'IBC', 
# #     'ROI_107':'DCIS', 
# #     'ROI_108':'IBC',
# #     'ROI_109':'IBC', 
# #     'ROI_110':'IBC',
# #     'ROI_128':'Normal'  
# # }

# roi_labels = {
#     'ROI_101':'DCIS', 
#     'ROI_102':'DCIS', 
#     'ROI_103':'DCIS', 
#     'ROI_104':'IBC', 
#     'ROI_105':'DCIS', 
#     'ROI_106':'IBC', 
#     'ROI_107':'DCIS', 
#     'ROI_108':'IBC',
#     'ROI_109':'IBC', 
#     'ROI_110':'IBC',
#     'ROI_111':'IBC', 
#     'ROI_112':'IBC',
#     'ROI_113':'IBC', 
#     'ROI_114':'IBC',
#     'ROI_115':'IBC', 
#     'ROI_116':'IBC',
#     'ROI_117':'IBC', 
#     'ROI_118':'IBC',
#     'ROI_119':'IBC', 
#     'ROI_120':'IBC',
#     'ROI_121':'DCIS',
#     'ROI_122':'DCIS', 
#     'ROI_123':'DCIS',
#     'ROI_124':'DCIS',
#     'ROI_125':'IBC', 
#     'ROI_126':'DCIS',
#     'ROI_127':'DCIS',
#     'ROI_128':'Normal'  
# }

# #good_peptides = data.columns[4:]
# #good_peptides = [1767.9236, 1326.6325, 1681.8167, 1421.7859, 1386.6873,1248.5677, 1609.7969, 1212.6219, 1205.5797, 1588.7813, 1458.7006 ] 
# good_peptides = [758.4519, 797.4264, 976.4517, 999.5946, 1060.5269, 1068.5069, 1082.6317, 1084.5018, 1089.4847, 1098.5062,	
#                  1098.5902, 1102.564, 1115.544, 1125.5283, 1128.528, 1128.532, 1142.5073, 1154.5073, 1166.5324, 1175.5473, 
#                  1205.5797, 1212.6219, 1226.6124, 1229.5909, 1248.5677, 1253.6008, 1257.5818, 1273.6747, 1280.6845,	1283.5862, 
#                  1286.66, 1291.6641, 1299.5811, 1326.6325, 1349.6444, 1351.6965, 1365.6474, 1379.7165, 1386.6873, 1399.6448,
#                  1406.705, 1406.7063, 1421.7859, 1458.7006, 1472.6434, 1480.7503, 1552.8078, 1588.7813, 1593.7867, 1595.7772, 
#                  1599.7358, 1609.7969, 1681.8167, 1692.8188, 1692.83, 1767.9236]

# def get_roi_stats():
#     '''
#     Calculates spatial centroid for each roi and the mean for every peptide in each roi. 
#     Return is used in generate_spatial_heatmap and get_hierarchical_clusters.  

#     Returns:
#         roi_stats (np array) contains roi name, spatial centroid (x & y), class, and means for every peptide
#         roi_stats is used as an argument in generate_spatial_heatmap and get_hierarchical_clusters
# '''
#     print('Calculating spatial centroids for each ROI...')
#     # Calculate spatial centroid of each roi
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





# import numpy as np 
# from matplotlib import pyplot as plt
# from pyslingshot import Slingshot
# from anndata import AnnData

# load = True
# num_cells = 1000
# num_dims_reduced = 2
# num_branches = 1

# K=10    # cluster labels
# filename = r'C:\Users\AngelLab\Downloads\fakedata-1branch.npy'
# start_node = 4

# if load:
#     data = np.load(filename, allow_pickle=True).item()
#     cluster_labels = data['cluster_labels']
#     data = data['data']

# #plt.scatter(data[:,0], data[:,1], c=cluster_labels)
# #plt.show()

# num_genes = 500
# ad = AnnData(np.zeros((num_cells, num_genes)))
# ad.obsm['X_umap'] = data
# ad.obs['celltype'] = cluster_labels
# ad

# fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,10))
# custom_xlim = (-12,12)
# custom_ylim = (-12,12)

# slingshot = Slingshot(ad, celltype_key='celltype', obsm_key='X_umap', start_node=start_node, is_debugging='verbose')
# slingshot.fit(num_epochs=1, debug_axes=axes)

# fig,axes = plt.subplots(ncols=2, figsize=(12,4))
# axes[0].set_title('Cluster')
# axes[1].set_title('Pseudotime')
# slingshot.plotter.curves(axes[0],slingshot.curves)
# slingshot.plotter.clusters(axes[0], labels=np.arange(slingshot.num_clusters), s=4, alpha=0.5)
# slingshot.plotter.clusters(axes[1], color_mode='pseudotime', s=5 )
# plt.show()



import json

def get_roi_labels_from_json(path):
    file = json.load(open(path))
    
    roi_labels = {}
    for entry in file['roi_labels']:
        roi_labels.update(entry)

    return roi_labels

PATH = r'C:\Users\AngelLab\Documents\GitHub\spatial-proteomics-analyzer\examples\example_roi_labels.json'
get_roi_labels_from_json(PATH)







'''
References for random forest model:
https://scikit-learn.org/stable/modules/permutation_importance.html#permutation-importance
https://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html
'''
