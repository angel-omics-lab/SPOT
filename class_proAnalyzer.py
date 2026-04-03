'''
This script contains class modules for spatial proteomics analysis 
'''
import pandas as pd 
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests 
from matplotlib import pyplot as plt
import os 
import math
import seaborn as sns
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage 
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import scale
import anndata as ad
import scanpy as sc
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx 
import time 

'''
TODO
-update roi_labels to be read from a json file at data_path rather than a passed argument 
-optional: have spatial roi map be overlayed the H&E picture, this would require (or make this conditional on its existence) user to have a .png or .jpeg in data folder
-have user pass classes (or just read unique classes from roi_labels.json)--would need to update a couple functions where 'DCIS'/'IBC' are hardcoded (get_roi_map and diff_expression_test notably)
-update plot formatting so all peptides are rounded to 3 decimals
-maybe: change get_random_forest_model_ranking to train random forest model on all peptides then rank them by their feature importance rather than one at a time (could also do both?)
-tag results folder name with date and time so it's unique? 
-add try except in allPipeline for non-critical/independent methods (ex. generate_boxplots)
-axes of dendrogram/fores barplot/heatmaps labeled False, should just be blank
-Dr. Angel feedback: add numbers to boxplots 
-edge case to fix: if all peptides are filtered out, then we get an error in generate boxplots, need to put an escape if all peptides are removed/nonsignificant
-separate box plots so they are saved individually rather than in array 
-forest barplot not in Arial for some reason
-change roi map so labels are larger and outside of dots
-remove legend from forest barplot 
-figure out why UMAP is taking so long
-remove 'class' title from pca/umap plots
-add points into box plots so can visualize potential clusters
'''

class SpatialOmicsAnalyzer:
    def __init__(self, data_path, roi_labels):
        self.data_path = data_path

        self.roi_labels = roi_labels #json.load(open('data.json', 'roi_labels'))       # Contains sheet name of each ROI and its label     
        self.good_peptides = None
        self.ann_obj = None
        try: 
            self.data = pd.read_excel(data_path, sheet_name=None)
        except FileNotFoundError as e:
            print('Issue loading spreadsheet, check file path')
            print(e)
        # Univeral plot formatting
        plt.rcParams.update({
            'font.family': 'Arial',
            'font.size': 12 
        })
    
    def check_sheet_format(self):
        '''
        Cleans each ROI sheet before any analysis:
            - Drops columns with non-numeric values
            - Drops fully empty columns and rows
            - Fills NaN intensities with 0 
            - Ensures the metadata (first 4 columns) are present
        
        '''
        print('Checking spreadsheet formatting...')
        for region in list(self.roi_labels.keys()):
            if region not in self.data:
                print(f' WARNING: {region} in roi_labels is not found in spreadsheet.')

            df = self.data[region]
            # Drop empty rows and columns
            df = df.dropna(how='all').reset_index(drop=True)
            df = df.loc[:, df.notna().any()]

            # Keep only numeric peptide column names
            metadata_cols = df.columns[:4].tolist()
            peptide_cols = df.columns[4:]

            numeric_peptide_cols = [
                c for c in peptide_cols 
                if isinstance(c, (int,float)) and not isinstance(c,list)
            ]
            dropped = set(peptide_cols) - set(numeric_peptide_cols)
            if dropped:
                print(f'{region}: dropping {len(dropped)} non-numeric columns')
            
            self.data[region] =  df[metadata_cols + numeric_peptide_cols]
        print('Spreadsheet formatting check complete.')


        
    def filter_rois(self):
        '''
        Load ROIs, filter out the ones with >25% 0s, then normalize the intensities via log transform
        
        Returns: 
            None, but roi_labels (dict) is updated from filtering. And intensities are normalized. 
        '''
        print('Filtering ROIs...')
        to_remove = []
        for region, label in self.roi_labels.items():
            data = self.data[region]
            intensities = data.iloc[:, 4:]      # Skip first 4 columns 
            print(f'Processing {region}...')
            intensities = intensities.fillna(0)        # Change all NaN values to 0 
            zeros = (intensities==0).astype(int).sum().sum()      # Number of 0 intensities
            total = intensities.shape[0]*data.shape[1]     # Total number of intensities
            
            if zeros > (0.25*total):       # Remove ROI from list if >25% of peptides are 0 
                to_remove.append(region)
                print(f'Removed {region} with {label} label from list: >25% zero intensities.')
            else: 
                print(f'{region} accepted')
        for region in to_remove:
            del self.roi_labels[region] 
        print('ROI filtering complete. Filtered list: ', self.roi_labels.keys())
        # return self.roi_labels


    
    def get_roi_map(self):
        '''
        Creates a spatial dot plot of the ROIs where each dot is the location of the ROIs spatial centroid. 
        Dots are colored according to their class. 
        NOTE potential ~upgrade~ is to have the dots overlay a picture of the tissue, This would require cv2 AND for user to pass a png/jpeg
        
        ''' 
        print('Generating spatial dot plot for regions...')
        # Calculate spatial centroid of each roi
        
        # Concatenate all ROI sheets, tagging each row with its ROI label
        combined = pd.concat(
            [df.assign(roi=roi) for roi, df in self.data.items() if roi in self.roi_labels],
            ignore_index=True
        )
            # Groupby computes centroid (vector operation)
        agg_dict = {'x': 'mean', 'y': 'mean'}
        roi_stats = combined.groupby('roi').agg(agg_dict).reset_index()     # roi_stats columns: ['roi', 'x', 'y', peptide_1, peptide_2, ...]
        roi_stats['class'] = roi_stats['roi'].map(self.roi_labels)
        roi_stats['roi'] = roi_stats['roi'].str.removeprefix('ROI_')    # remove prefix so label can fit inside dot
                
            # Generate colored scatter plot
        sns.scatterplot(data=roi_stats, x='x', y='y', 
                        hue='class', palette={'DCIS':'dodgerblue', 'IBC':'orange', 'Normal':'green'},      
                        s=250
        )
        for x, y, roi in zip(roi_stats['x'], roi_stats['y'], roi_stats['roi']):
            plt.text(x, y, str(roi),
                    ha='right', va='bottom', 
                    fontsize=13, fontweight='bold', color='black')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.xticks([])
        plt.yticks([])
        plt.xlabel(None)
        plt.ylabel(None)
        plt.title('Spatial ROI Plot')
        plt.tight_layout()

        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/roi_map.png'))
        plt.close()
        sns.reset_orig()

        print('Spatial dot plot generation successful. Saved to: ', os.path.dirname(self.data_path))
        
    

    def peptide_sparsity_filter(self): 
        '''
        Removes peptides only sparsely expressed across ROIs. Those with <20% expression in >10% of ROIs are removed. 

        Returns: 
            None, but updates good_peptides(list): list of all peptides that are well expressed across ROIs. 
        '''
        
        print('Filtering peptides...')
        # peptide_list = self.data[region].columns[4:]
        # Get union of all peptide columns across ROIs (handles ROIs w different peptides)       
        all_peptides = set()    
        for region in self.roi_labels:
            all_peptides.update(self.data[region].columns[4:])        
        
        zero_counts_per_peptide = pd.Series(0, index=list(all_peptides))
        
        for region in self.roi_labels:
            intensities = self.data[region].iloc[:, 4:]     # Grabs peptide columns for the roi
            zero_pct = (intensities==0).mean()      # Calculate percentage of 0 intensities for each peptide
            # Identify peptides as T/F based on zero_pct, and aligns Series to full peptide index (allows for peptides existing in some ROIs but not others)
            flagged = (zero_pct > 0.2).reindex(zero_counts_per_peptide.index, fill_value=0)     
            zero_counts_per_peptide = zero_counts_per_peptide.add(flagged.astype(int))      # Tells how many ROIs found each peptide to be too sparse, used by threshold

        # Keep peptides that are present in >= 10% ROIs 
        threshold = len(self.roi_labels)*0.10       
        self.good_peptides = zero_counts_per_peptide[zero_counts_per_peptide <= threshold].index.tolist()

        print('Peptide filtering complete. Filtered list: ', self.good_peptides)
        # return [float(p) for p in self.good_peptides]
    

    def normalize_intensities(self):
        '''
        Natural log transforms peptide intensities. Zeros are replaced with NaN before log transform to avoid -inf values. This function also rounds peptides to 3 decimals. 
        '''
        print('Normalizing intensities...')
        for region in self.roi_labels:
            data = self.data[region]
            cols = [p for p in self.good_peptides if p in data.columns]
            data[cols] = data[cols].replace(0,np.nan)
            data[cols] = np.log(data[cols])
            data[cols] = data[cols].replace(np.nan, 0)

        print('Normalization complete.')


    def diff_expression_test(self): 
        '''
        Identifies peptides that have differential expression between classes (in this instance DCIS v IBC) via 
        Kruskal-Wallis test with Benjamini-Hochberg FDR validation.
        
        Returns:
            good_peptides(list) updated so it is now only those differentially expressed.
        '''
        results_list = []
        # Identify significant peptide 
        print('Identifying differentially expressed peptides...')
        for peptide in self.good_peptides:
            group1 = [self.data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'DCIS']     # Mean intensities of peptide in DCIS (class 1) rois        
            group2 = [self.data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'IBC']      # Mean Intensities of peptide in IBC (class 2) rois
            # Run KW
            H_statistic, p_value = stats.kruskal(group1, group2)
            # Add to data frame 
            results_list.append({
            'peptide': peptide, 
            'pvalue_kw': p_value
            })
        results_df = pd.DataFrame(results_list) 
        # Run BH FDR, add q val to df 
        reject, q_values, _, _ = multipletests(results_df['pvalue_kw'], alpha=0.05, method='fdr_bh')
        results_df['qvalue_bh'] = q_values
        # Save peptides whose q value is <= 0.05
        self.good_peptides = results_df[results_df['qvalue_bh'] <= 0.05]['peptide'].tolist()
        if not self.good_peptides: 
            raise ValueError('No peptides were identified to be significant with 95 percent confidence. Pipeline stopped.') 
        print('Differential analysis complete. Significant peptides: ', self.good_peptides)
        print(f'{len(reject)} peptides removed: non-significant')
        
    

    def generate_boxplots(self): 
        '''
        Takes the good_peptides list and generates a grid of boxplots where each boxplot compares mean intensities between 2 classes for each peptide. 
        In this case, the two classes we choose are IBC and DCIS, but this can be changed. 
        '''
        print('Generating box plots for significant peptides...')

        for peptide in self.good_peptides:
            plt.figure()
            plot_data = (
                [{'Intensity': val, 'Class':'DCIS'} for val in [self.data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'DCIS' ]]
                +
                [{'Intensity': val, 'Class':'IBC'} for val in [self.data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'IBC' ]]
            )
            plot_df = pd.DataFrame(plot_data)

            ax = sns.boxplot(
                data=plot_df, 
                x='Class', y='Intensity', 
                hue='Class', 
                palette={'DCIS':'dodgerblue', 'IBC':'orange'}
            )

            ax.set_title(peptide)
            ax.set_ylabel('ln(Median Intensity)')
            ax.set_xlabel('')
            ax.legend_ = None
            plt.tight_layout()

            plt.savefig(os.path.join(os.path.dirname(self.data_path), f'results/boxplot_{peptide}.png'))
            plt.close()

        print('Box plots generation successful!')
        sns.reset_orig()


    def get_roi_stats(self):
        '''
        Calculates spatial centroid for each roi and the mean for every peptide in each roi. 
        Return is used in generate_spatial_heatmap and get_hierarchical_clusters.  

        Returns:
            roi_stats (np array) contains roi name, spatial centroid (x & y), class, and means for every peptide
                roi_stats is used as an argument in generate_spatial_heatmap and get_hierarchical_clusters
        '''
        print('Calculating spatial centroids for each ROI...')
        # Calculate spatial centroid of each roi
            # Concatenate all ROI sheets, tagging each row with its ROI label
        combined = pd.concat(
            [df.assign(roi=roi) for roi, df in self.data.items() if roi in self.roi_labels],
            ignore_index=True
        )
            # Groupby computes centroid AND all peptide means in one pass (vector operation)
        agg_dict = {'x': 'mean', 'y': 'mean', **{p: 'median' for p in self.good_peptides}}
        roi_stats = combined.groupby('roi').agg(agg_dict).reset_index()     # roi_stats columns: ['roi', 'x', 'y', peptide_1, peptide_2, ...]
        roi_stats['class'] = roi_stats['roi'].map(self.roi_labels)

        print('ROI calculations successful!')
        return roi_stats


    def generate_spatial_heatmap(self, roi_stats):
        '''
        Generates a spatial heatmap (on the original image) for each peptide based on the centroid of each ROI and mean intensity of each peptide.

        Returns:
            heatmap_X (png) where X is the peak, makes one for every good_peptide 
            roi_stats (np array) 
        '''
        print('Generating heatmaps for each peptide...')
        # Generate heatmaps for each peptide as a scatterplot where color indicates mean intensity
        for peptide in self.good_peptides:
            heatmap = plt.figure()
            plt.scatter(roi_stats['x'], roi_stats['y'], 
                        c=roi_stats[peptide], cmap='jet', s=70)
            cb = plt.colorbar()
            plt.xlabel(None)
            plt.xticks([])
            plt.ylabel(None)
            plt.yticks([])
            plt.title(peptide, font='Arial', fontsize=14)
            plt.tight_layout()
            plt.savefig(os.path.join(os.path.dirname(self.data_path), f'results/heatmap_{peptide}.png'))
            cb.remove()
            plt.close()
            
        print('Spatial heatmap generation successful.')
        return roi_stats        



    def make_hierarchical_clusters(self, roi_stats):
        '''
        Generates a dendrogram showing hierarchical clusters based on Ward's (+ Euclidean) method. 
        i.e. ROIs with similar peptide profiles will be clustered together
        
        Returns:
            dendrogram (png)
        '''
        print('Generating hierarchical clusters...')
        # Extract peptide intensities per roi 
        feature_matrix = roi_stats[self.good_peptides].values
        # Create dendrogram
        print('Constructing dendrogram for visualization...')
        linkage_data = linkage(feature_matrix, method='ward', metric='euclidean')
        dend = dendrogram(linkage_data, labels=roi_stats['roi'].values)
 
        # Formatting
        plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['mediumorchid', 'lightcoral', 'lightblue'])
        color_map = {'DCIS': 'dodgerblue', 'IBC': 'orange', 'Normal': 'green'}
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
        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/dendrogram.png'))
        plt.close()
        plt.rcParams['axes.prop_cycle'] = plt.rcParamsDefault['axes.prop_cycle']

        print('Dendrogram generation successful.')



    def get_random_forest_ranking(self, roi_stats):
        '''
        Runs a random forest classification model for all peptides then outputs a bar plot of peptides with their 'discriminating' score. 
        
        Returns
        random_forest_bar_graph (png)
        '''
        print('Generating random forest models per peptide...')
        oob_list = []
        # Convert peptides to string so yticks interpreted as categorical
        good_peptides_str = [str(p) for p in self.good_peptides]
        roi_stats_str = roi_stats
        roi_stats_str.columns = [str(c) for c in roi_stats_str.columns]
        for peptide in good_peptides_str : 
            X = roi_stats_str[[peptide]]    
            y = roi_stats_str['class']
            # Randomly split the dataset so 10% is reserved for model validation 
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)
            
            forest_model = RandomForestClassifier(n_estimators=100, oob_score=True, random_state=42)
            forest_model.fit(X_train, y_train)
            oob_list.append(forest_model.oob_score_)
            
        # Plot each peptide oob on bar plot
        print('Constructing oob bar plot...')
        # Sort them so plot can have them listed as ascending
        sorted_pairs = sorted(zip(oob_list, good_peptides_str), key=lambda x:x[0])
        sorted_scores, sorted_peptides = zip(*sorted_pairs)
        
        plt.barh(sorted_peptides, sorted_scores, color='mediumorchid')
        plt.title('Peptide discrimination score in random forest models')
        plt.xlabel('OOB score')
        plt.xlim(0,1)
        plt.xticks()
        plt.yticks(sorted_peptides)
        plt.tight_layout()
        plt.legend_ = None

        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/forest_barplot.png'))
        plt.close()
        print('Bar plot generation successful.')


    def create_anndata_object(self):
        '''
        Creates an AnnData object for pixel-level analysis. Concatenates all accepted ROI sheets, filters to good_peptides, z-score scales intensities, 
        and attached per-pixel metadata (ROI, class, x-coord, y-coord). This object is used as the input for PCA and UMAP analysis. 
        
        Returns: 
            ann_obj (AnnData): X (n_pixels, n_peptides), obs metadata, and var peptide index.
        
        '''
        print('Constructing AnnData object...')
        # Concat ROI sheets into a single flat df
        combined = pd.concat([df.assign(ROI=roi) for roi, df in self.data.items() if roi in self.roi_labels], 
                             ignore_index=True)
        
        # Extract and scale peptide intensity matrix 
        combined_intensities = combined[self.good_peptides].values.astype(float)
        combined_intensities = scale(combined_intensities)      # z-score scaling so peptides contribute equally to PCA
        
        # Build per-pixel metadata (obs)
        pixel_metadata = pd.DataFrame({
            'sample': combined['ROI'].values, 
            'class': combined['ROI'].map(self.roi_labels).values, 
            'x': combined['x'].values.astype(float),
            'y': combined['y'].values.astype(float)
        })
        
        # Build per-peptide metadata (var)
        peptide_metadata = pd.DataFrame(
            {'peptide':self.good_peptides}, 
            index=[str(p) for p in self.good_peptides]
        )
        
        self.ann_obj = ad.AnnData(
            X = combined_intensities,   # (n_pixels, n_peptides)
            obs = pixel_metadata,       # one row per pixel
            var = peptide_metadata      # one row per peptide
        )
        
        print('AnnData object created successfully. Shape:', self.ann_obj.shape)
    
    
    def run_pixel_dim_reduction(self):
        '''
        Runs PCA and UMAP on a pixel-level AnnData object then saves the resulting plots, colored by class label. 
        We use scanpy rather than sklearn because it saves computational time when identifying num of PCs to retain. 
        
        Args:
            ann_obj (AnnData): output of create_anndata_object()
        
        Returns:
            ann_obj (AnnData): updated with PCA and UMAP embeddings
            Not returned, but saves both pca and umap plots
        '''
        sc.settings.figdir = os.path.join(os.path.dirname(self.data_path), "results")
        # Register class colors so scanpy can use them directly 
        class_colors = {'Normal':'green', 'DCIS': 'dodgerblue', 'IBC': 'orange'}
        self.ann_obj.obs['class'] = pd.Categorical(self.ann_obj.obs['class'])
        self.ann_obj.uns['class_colors'] = [
            class_colors[c] for c in self.ann_obj.obs['class'].cat.categories
        ]

        # Calculate number of components to retain for PCA with Explained Variance Threshold heuristic
        #print('Calculating optimal number of PCA components...')
        # max_comps = min(self.ann_obj.n_obs, self.ann_obj.n_vars) - 1
        # sc.pp.pca(self.ann_obj, n_comps=max_comps)       # Run pca with max comps
        # cumsum = np.cumsum(self.ann_obj.uns['pca']['variance_ratio'])
        
        # n_comps = int(np.argmax(cumsum>=0.99) + 1)      # Retain comp
        # print(f'Retaining {n_comps} principal components (explain >= 99% of variance)')
        n_comps=5
        print('Running PCA...') 
        sc.pp.pca(self.ann_obj, n_comps=n_comps)
        
        print('Running UMAP analysis; this may take a while...')
        sc.pp.neighbors(self.ann_obj, use_rep='X_pca', n_pcs=n_comps)
        sc.tl.umap(self.ann_obj)
        
        print('Generating PCA and UMAP figures...')
        sc.pl.pca(self.ann_obj, color='class', 
                  show=False, save='pca.png')
        print('PCA plot generated successfully.')
        sc.pl.umap(self.ann_obj, color='class', 
                  show=False, save='umap.png')
        print('UMAP plot generated successfully.')

    
    def compute_mst(self):
        '''
        Generates a minimum spanning tree (MST) connecting ROI centroids in PCA space. 
        It is used as the 'backbone' of TSCAN and is overlaid on the UMAP generated from it. 
        
        Args:
            ann_obj (AnnData) : output of run_pixel_analysis(), is original ann_obj with added X_pca and X_umap in obs 
        
        Returns: 
            ann_obj (AnnData) : updated with MST stored in ann_obj.uns['mst']
                uns['mst']['graph'] -- networkx graph of MST
                uns['mst']['centroids_pca'] -- df, ROI centroids in PCA space
        '''
        if self.ann_obj is None:
            raise RuntimeError('ann_obj is empty. Run create_anndata_object() and run_pixel_analysis() first')
        
        print('Computing MST...')
        # Compute roi centroids in pca space
        n_comps = self.ann_obj.obsm['X_pca'].shape[1]
        pca_df = pd.DataFrame(
            self.ann_obj.obsm['X_pca'], 
            index = self.ann_obj.obs_names, 
            columns=[f'PC{i+1}' for i in range(n_comps)]
        )
        pca_df['sample']=self.ann_obj.obs['sample'].values
        centroids_pca = pca_df.groupby('sample').mean()
        
        # Build mst on centroids
        dist_matrix = cdist(centroids_pca.values, centroids_pca.values, metric='euclidean')     # Gets Eucl. distance between every ROI centroid pair in PCA space
        mst_sparse = minimum_spanning_tree(dist_matrix)         # Actual MST algorithm -- finding edge subset that minimizes total distance, outputs sparse matrix identifying those edges
        
        graph = nx.from_scipy_sparse_array(mst_sparse)      # Converts sparse matrix into graph object 
        roi_names = list(centroids_pca.index)
        mst_graph = nx.relabel_nodes(graph, {i: name for i, name in enumerate(roi_names)})      # Swaps default node names for actual ROIs
        
        print(f'MST built with {mst_graph.number_of_nodes()} nodes and {mst_graph.number_of_edges()} edges')
        
        # Add MST data to ann_obj
        self.ann_obj.uns['mst'] = {
            'graph': mst_graph, 
            'centroids' : centroids_pca
        }
            



    def run_pseudotime_slingshot(self):
        '''
        Runs Slingshot pseudotime reconstruction. Outputs UMAP with Slingshot principal curves and ROI centroids overlaid. 
        '''
        from pyslingshot import Slingshot
        print('Running pseudotime reconstruction. This may take a while...')
        # Step 1: Clustering
        sc.pp.neighbors(self.ann_obj, use_rep='X_umap', n_neighbors=15)     # Overwrite the neighbor graph to be based in UMAP space, making leiden clusters more coherent
        sc.tl.leiden(self.ann_obj, resolution=0.02)     # resolution up - more clusters V down = fewer clusters, ideally want # clusters to = # disease stages
        print(f"Leiden found {self.ann_obj.obs['leiden'].nunique()} clusters")

        # Step 2: Set root (identify cluster that Normal cells belong to)
        normal_mask = self.ann_obj.obs['class']=='Normal'
        normal_clusters = self.ann_obj.obs.loc[normal_mask, 'leiden']
        root_cluster = int(normal_clusters.mode().iloc[0])

        #Step 3: Slingshot fitting
        slingshot = Slingshot(
            self.ann_obj, 
            celltype_key='leiden', 
            obsm_key='X_umap', 
            start_node=root_cluster
        )
        slingshot.fit()
        self.ann_obj.obs['pseudotime'] = slingshot.unified_pseudotime
        print('Pseudotime fitting successful. Now proceeding to plotting...')

        # Step 4: Plot
        umap_coords = self.ann_obj.obsm['X_umap']
        self.ann_obj.obs['umap1'] = umap_coords[:,0]
        self.ann_obj.obs['umap2'] = umap_coords[:,1]

        class_colors = {'Normal':'green','DCIS':'dodgerblue', 'IBC':'orange'}
        point_colors = self.ann_obj.obs['class'].map(class_colors).values

        roi_centroids = self.ann_obj.obs.groupby('sample')[['umap1', 'umap2']].mean()
        roi_class = self.ann_obj.obs.groupby('sample')['class'].first()
        centroid_colors = roi_class.map(class_colors).values

        fig,ax = plt.subplots()
        # Plot UMAP as scatterplot
        scatter = ax.scatter(
            self.ann_obj.obs['umap1'], 
            self.ann_obj.obs['umap2'], 
            c=point_colors, 
            s=5
        )

        # Overlay centroids
        ax.scatter(
            roi_centroids['umap1'],
            roi_centroids['umap2'],
            c=centroid_colors,
            s=80, 
            edgecolors='black',
            linewidths=0.8,
            zorder=5
        )
        # Label centroids
        for roi, row in roi_centroids.iterrows():
            ax.annotate(
                roi, 
                xy=(row['umap1'], row['umap2']),
                fontsize=7,
                fontweight='bold',
                ha='center',
                va='bottom',
                xytext=(0,5),
                textcoords='offset points'
            )
        # Overlay slingshot curves
        for curve in slingshot.curves:
            points = curve.points
            ax.plot(
                points[:,0],
                points[:,1],
                color='black',
                linewidth = 1.5
            )
        
        ax.set_xlabel('UMAP1'),
        ax.set_ylabel('UMAP2')
        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/umap_pseudotime.png'))
        plt.close()
    
        
    
    
##### Entire pipeline ##### 
    def allPipeline(self):
        start = time.time()
        os.makedirs((os.path.join(os.path.dirname(self.data_path), 'results')), exist_ok=True)
        try : 
            self.check_sheet_format()
            self.filter_rois()
            self.get_roi_map()
            self.peptide_sparsity_filter()
            self.normalize_intensities()
            self.diff_expression_test()
            self.generate_boxplots()
            roi_stats = self.get_roi_stats()
            self.generate_spatial_heatmap(roi_stats)
            self.make_hierarchical_clusters(roi_stats)
            self.get_random_forest_ranking(roi_stats)
            self.create_anndata_object()
            self.run_pixel_dim_reduction()
            # self.compute_mst()
            self.run_pseudotime_slingshot()
        except:
            print('Pipeline broke before finishing.')
        finally: 
            end = time.time()
            duration = (end-start) // 60    # in minutes
            print(f'Total duration: {duration:.2f} minutes')
