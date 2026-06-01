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
import json
from pcurvepy2 import PrincipalCurve


class SpatialOmicsToolkit:
    def __init__(self, data_path, json_path):
        self.data_path = data_path 
        self.good_peptides = None
        self.adata = None
        self.roi_stats = None
        try: 
            self.data = pd.read_excel(data_path, sheet_name=None)
        except FileNotFoundError as e:
            print('Issue loading spreadsheet, check file path.')
            print(e)

        try: 
            file = json.load(open(json_path))
            self.roi_labels = {}
            for entry in file['roi_labels']:
                self.roi_labels.update(entry)
<<<<<<< Updated upstream
=======
            for entry in file['classes']:
                self.roi_class_colors.update(entry)
            self.roi_classes = list(self.roi_class_colors.keys())
            print('ROI classes: ', self.roi_classes)
            print('ROI class colors: ', self.roi_class_colors)

            self.normal_class = 'Normal' if 'Normal' in self.roi_classes else None
            self.analysis_classes = [c for c in self.roi_classes if c != self.normal_class]

>>>>>>> Stashed changes
        except Exception as e:
            print('Issue loading and/or reading json file, check file path.')
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
            metadata_cols = df.columns[:2].tolist()
            peptide_cols = df.columns[2:]

            numeric_peptide_cols = [
                c for c in peptide_cols 
                if isinstance(c, (int,float)) and not isinstance(c,list)
            ]
            dropped = set(peptide_cols) - set(numeric_peptide_cols)
            if dropped:
                print(f'{region}: dropping {len(dropped)} non-numeric columns')
            
            df = df[metadata_cols + numeric_peptide_cols]
            
            rename_map = {c: round(c,3) for c in numeric_peptide_cols}
            df = df.rename(columns=rename_map)
            
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
            intensities = data.iloc[:, 2:]      # Skip x,y columns 
            print(f'Processing {region}...')
<<<<<<< Updated upstream
            intensities = intensities.fillna(0)        # Change all NaN values to 0 
            # zeros = (intensities==0).astype(int).sum().sum()      # Number of 0 intensities
            # total = intensities.shape[0]*data.shape[1]     # Total number of intensities
            
            # if zeros > (0.25*total):       # Remove ROI from list if >25% of peptides are 0 
            #     to_remove.append(region)
            #     print(f'Removed {region} with {label} label from list: >25% zero intensities.')
            # else: 
            #     print(f'{region} accepted')
            NA_prop = (intensities == 0).mean()
            if (NA_prop > 0.25).mean() > 0.25:
=======
            
            zero_rows = (intensities == 0).all(axis=1)
            zero_row_prop = zero_rows.mean() 

            if (zero_row_prop > 0.25).mean() > 0.25:
>>>>>>> Stashed changes
                to_remove.append(region)
                print(f'Removed {region} with {label} label from list: >25% zero intensities.')
            else: 
                print(f'{region} accepted')
        for region in to_remove:
            del self.roi_labels[region] 
        print('ROI filtering complete. Filtered list: ', self.roi_labels.keys())


    
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
            all_peptides.update(self.data[region].columns[2:])        
        
        zero_counts_per_peptide = pd.Series(0, index=list(all_peptides))
        
        for region in self.roi_labels:
            intensities = self.data[region].iloc[:, 2:]     # Grabs peptide columns for the roi
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
            # data[cols] = data[cols].replace(0,np.nan)
            # data[cols] = np.log(data[cols])
            # data[cols] = data[cols].replace(np.nan, 0)
            data[cols] = np.log1p(data[cols])
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
<<<<<<< Updated upstream
        for peptide in self.good_peptides:
            group1 = [self.data[region][peptide].median() for region in self.roi_labels if self.roi_labels[region] == 'DCIS']     # Mean intensities of peptide in DCIS (class 1) rois        
            group2 = [self.data[region][peptide].median() for region in self.roi_labels if self.roi_labels[region] == 'IBC']      # Mean Intensities of peptide in IBC (class 2) rois
=======
        for peptide in self.good_peptides:                                                                           
            group1 = [self.data[region][peptide].median() 
                      for region in self.roi_labels 
                      if self.roi_labels[region] == self.analysis_classes[0]
                      and peptide in self.data[region].columns]
            group2 = [self.data[region][peptide].median() 
                      for region in self.roi_labels 
                      if self.roi_labels[region] == self.analysis_classes[1]
                      and peptide in self.data[region].columns]             
>>>>>>> Stashed changes
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
        self.good_peptides = results_df[results_df['pvalue_kw'] <= 0.05]['peptide'].tolist()
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
            fig, ax = plt.subplots(figsize=(2.5, 4))  # create axes explicitly
            plt.rcParams['font.serif'] = ['Arial']
<<<<<<< Updated upstream
            plot_data = (
                [{'Intensity': val, 'Class':'DCIS'} for val in [self.data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'DCIS' ]]
                +
                [{'Intensity': val, 'Class':'IBC'} for val in [self.data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'IBC' ]]
=======

            medians_class1 = [self.data[region][peptide].median() 
                              for region in self.roi_labels 
                              if self.roi_labels[region] == self.analysis_classes[0] 
                              and peptide in self.data[region].columns]
            medians_class2 = [self.data[region][peptide].median() 
                              for region in self.roi_labels 
                              if self.roi_labels[region] == self.analysis_classes[1] 
                              and peptide in self.data[region].columns] 

            plot_data = (                                    
                [{'Intensity': val, 'Class':self.analysis_classes[0]} for val in medians_class1]
                +
                [{'Intensity': val, 'Class':self.analysis_classes[1]} for val in medians_class2]
>>>>>>> Stashed changes
            )
            plot_df = pd.DataFrame(plot_data)

            sns.boxplot(
                data=plot_df, 
                x='Class', y='Intensity', 
                hue='Class', 
<<<<<<< Updated upstream
                palette={'DCIS':'dodgerblue', 'IBC':'orange'},
=======
                palette=self.roi_class_colors,    
>>>>>>> Stashed changes
                ax=ax
            )

            sns.swarmplot(
                data=plot_df,
                x='Class', y='Intensity',
                color='black',
                ax=ax
            )

            ax.set_title(peptide)
            ax.set_ylabel('ln(Median Intensity)')
            ax.set_xlabel('')
            ax.legend_ = None
            fig.tight_layout()

            fig.savefig(os.path.join(os.path.dirname(self.data_path), f'results/boxplot_{peptide}.png'))
            plt.close(fig)

        print('Box plots generation successful!')
        sns.reset_orig()


    def get_roi_stats(self):
        '''
        Calculates spatial centroid for each roi and the mean for every peptide in each roi. 
        Peptide medians are computed only over ROIs that actually contain that peptide. 
        
        self.roi_stats (np array) is created and contains roi name, spatial centroid (x,y), class, and medians. 
        It is used in generate_spatial_heatmap and get_hierarchical_clusters.  

        '''
        print('Calculating spatial centroids for each ROI...')
        
        roi_rows = []
        for roi, df in self.data.items():
            if roi not in self.roi_labels:
                continue
            
            row = {'roi': roi}
            
            # Spatial centroid 
            row['x'] = df['x'].mean()
            row['y'] = df['y'].mean()
            
            # Peptide medians
            for p in self.good_peptides:
                row[p] = df[p].median() if p in df.columns else np.nan
                
            roi_rows.append(row)
            
            
        self.roi_stats = pd.DataFrame(roi_rows)
        self.roi_stats['class'] = self.roi_stats['roi'].map(self.roi_labels)

        print('ROI calculations successful!')
<<<<<<< Updated upstream
=======
        print('Saving medians to datafram...')                
        class1_rows = self.roi_stats[self.roi_stats['class'] == self.analysis_classes[0]][self.good_peptides]
        class2_rows  = self.roi_stats[self.roi_stats['class'] == self.analysis_classes[1]][self.good_peptides]

        plot_data_df = pd.DataFrame({
            'median_class1_allrois': class1_rows.median(),      # NOTE
            'median_class2_allrois':  class2_rows.median(),       # NOTE
        }) 

        self.output_excel = pd.concat([self.output_excel, plot_data_df], axis=1)
>>>>>>> Stashed changes


    def generate_spatial_heatmap(self):
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
            plt.rcParams['font.serif'] = ['Arial']
            plt.scatter(self.roi_stats['x'], self.roi_stats['y'], 
                        c=scale(self.roi_stats[peptide]), cmap='jet', s=200)
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



    def make_hierarchical_clusters(self):
        '''
        Generates a dendrogram showing hierarchical clusters based on Ward's (+ Euclidean) method. 
        i.e. ROIs with similar peptide profiles will be clustered together
        
        Returns:
            dendrogram (png)
        '''
        print('Generating hierarchical clusters...')
        # Extract peptide intensities per roi 
        feature_matrix = self.roi_stats[self.good_peptides].values
        
        # Check for NaN/inf values
        if not np.isfinite(feature_matrix).all():
            n_bad=(~np.isfinite(feature_matrix)).sum()
            print(f'WARNING : Feature matrix contains {n_bad} non-finite values. Filling with column medians.')
            df_fm = pd.DataFrame(feature_matrix, columns=self.good_peptides)
            feature_matrix = df_fm.fillna(df_fm.median()).values
        
        # Create dendrogram
        print('Constructing dendrogram for visualization...')
        plt.rcParams['lines.linewidth'] = 2.5
        plt.rcParams['font.serif'] = ['Arial']
        linkage_data = linkage(feature_matrix, method='ward', metric='euclidean')
        dend = dendrogram(linkage_data, labels=self.roi_stats['roi'].values)
 
        # Formatting
        plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['mediumorchid', 'lightcoral', 'lightblue'])
        color_map = {'DCIS': 'dodgerblue', 'IBC': 'orange', 'Normal': 'green'}
        leaf_order = dend['ivl']  # ROIs in dendrogram order
        label_colors = [color_map[self.roi_stats.set_index('roi').loc[roi, 'class']] for roi in leaf_order]
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



    def get_random_forest_ranking(self):
        '''
        Runs a random forest classification model for all peptides then outputs a bar plot of peptides with their 'discriminating' score. 
        
        No returns, but generates a bar plot with every significant peptide and its feature importance (bar plot). Saved as a png
        '''
        print('Running random forest classification...')
        from sklearn.inspection import permutation_importance

        good_peptides_str = [str(p) for p in self.good_peptides]
        roi_stats_str = self.roi_stats.copy()
        roi_stats_str.columns = [str(c) for c in roi_stats_str.columns]
        roi_stats_str = roi_stats_str[roi_stats_str['class'] != self.normal_class]

        X = roi_stats_str[good_peptides_str]
        y = roi_stats_str['class']
        
        print('Fitting forest model...')
        model = RandomForestClassifier(n_estimators=1000, random_state=42, oob_score= True, bootstrap=True)
        model.fit(X, y)
        oob_score_pct = round(model.oob_score_*100 , 2)                                        
        print(f'Model accuracy (OOB score): {oob_score_pct}%')

        # Feature importance based on mean decrease in impurity (MDI)
        importances = model.feature_importances_
        variable_importances = pd.Series(importances, index=good_peptides_str).sort_values()
        # perm_imp = permutation_importance(model, X, y, n_repeats=30, scoring='accuracy', random_state=42)
        # variable_importances = pd.Series(perm_imp.importances_mean, index=good_peptides_str).sort_values()


        print('Plotting...')
        fig, ax = plt.subplots(figsize=(8, max(4, len(self.good_peptides) * 0.5 + 1)))
        plt.rcParams['font.serif'] = ['Arial']
        variable_importances.plot.barh(
            ax=ax,
            color='mediumorchid',
            #xerr=pd.Series(perm_imp.importances_std, index=good_peptides_str).reindex(variable_importances.index), 
            #ecolor='gray'
        )

        ax.set_title(f'Ranking of peptide importance in random forest classification model\n Model accuracy: {oob_score_pct}%')
        ax.set_xlabel('Feature importance')
        ax.set_ylabel('Peptide m/z')
        plt.tight_layout()
        plt.show()
        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/barplot.png'))
        plt.close()

    def create_anndata_object(self):
        '''
        Creates an AnnData object for pixel-level analysis. Concatenates all accepted ROI sheets, filters to good_peptides, z-score scales intensities, 
        and attached per-pixel metadata (ROI, class, x-coord, y-coord). This object is used as the input for PCA and UMAP analysis. 
        
        Returns: 
            adata (AnnData): X (n_pixels, n_peptides), obs metadata, and var peptide index.
        
        '''
        print('Constructing AnnData object...')
        # Concat ROI sheets into a single flat df
        combined = pd.concat([df.assign(ROI=roi) for roi, df in self.data.items() if roi in self.roi_labels], 
                             ignore_index=True)
        
        # Extract peptide intensity matrix 
        combined_intensities = combined[self.good_peptides].values.astype(float)
        
        # Impute NaNs with per-peptide medians
        if np.isnan(combined_intensities).any():
            col_medians = np.nanmedian(combined_intensities, axis=0)
            nan_mask = np.isnan(combined_intensities)
            combined_intensities[nan_mask] = np.take(col_medians, np.where(nan_mask)[1])
            print(f'Imputed NaNs in {nan_mask.sum()} pixel-peptide entries with per-peptide medians.')
        
        # Scale peptide intensity matrix 
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
        
        self.adata = ad.AnnData(
            X = combined_intensities,   # (n_pixels, n_peptides)
            obs = pixel_metadata,       # one row per pixel
            var = peptide_metadata      # one row per peptide
        )
        
        print('AnnData object created successfully. Shape:', self.adata.shape)
    
    
    def run_pixel_dim_reduction(self):
        '''
        Runs PCA and UMAP on a pixel-level AnnData object then saves the resulting plots, colored by class label. 
        We use scanpy rather than sklearn because it saves computational time when identifying num of PCs to retain. 
        
        Args:
            adata (AnnData): output of create_anndata_object()
        
        Returns:
            adata (AnnData): updated with PCA and UMAP embeddings
            Not returned, but saves both pca and umap plots
        '''
        sc.settings.figdir = os.path.join(os.path.dirname(self.data_path), "results")
        # Register class colors so scanpy can use them directly 
        class_colors = {'Normal':'green', 'DCIS': 'dodgerblue', 'IBC': 'orange'}
        self.adata.obs['class'] = pd.Categorical(self.adata.obs['class'])
        self.adata.uns['class_colors'] = [
            class_colors[c] for c in self.adata.obs['class'].cat.categories
        ]

        #Calculate number of components to retain for PCA with Explained Variance Threshold heuristic
        print('Calculating optimal number of PCA components...')
        max_comps = min(self.adata.n_obs, self.adata.n_vars) - 1
        sc.pp.pca(self.adata, n_comps=max_comps)       # Run pca with max comps
        cumsum = np.cumsum(self.adata.uns['pca']['variance_ratio'])
        
<<<<<<< Updated upstream
        n_comps = int(np.argmax(cumsum>=0.95) + 1)      # Retain comp
        print(f'Retaining {n_comps} principal components (explain >= 99% of variance)')

        print('Running PCA...') 
        sc.pp.pca(self.adata, n_comps=3)
=======
        n_comps_ret = int(np.argmax(cumsum>=0.95) + 1)      # Retain comp
        print(f'Retaining {n_comps_ret} principal components (explain >= 95% of variance)')

        print('Running PCA...') 
        sc.pp.pca(self.adata, n_comps=n_comps_ret)
>>>>>>> Stashed changes
        
        print('Running UMAP analysis; this may take a while...')
        sc.pp.neighbors(self.adata, use_rep='X_pca', n_pcs=3, n_neighbors=50)
        sc.tl.umap(self.adata, min_dist=0.75, n_components=2)
        
        print('Generating PCA and UMAP figures...')
        sc.pl.pca(self.adata, color='class', 
                  show=False, save='.png')
        print('PCA plot generated successfully.')
        sc.pl.umap(self.adata, color='class',
                  show=False, save='.png')
        print('UMAP plot generated successfully.')
        
        # # Calculate principal curve to be added to UMAP projection 
        # pc = PrincipalCurve()
        # pc.fit(self.adata.obsm['X_umap'])
        # curve_points= pc.points

        # Plot ROI centroid UMAP projection 
        print('Generated ROI centroid UMAP projection...')
        umap_coords = pd.DataFrame(
          self.adata.obsm['X_umap'], 
          columns=['UMAP1', 'UMAP2'],
          index=self.adata.obs.index
        )
        umap_coords['sample'] = self.adata.obs['sample'].values
        umap_coords['class'] = self.adata.obs['class'].values
        
        # Average UMAP coords across all pixes in each roi
        roi_umap = umap_coords.groupby('sample')[['UMAP1', 'UMAP2']].mean().reset_index()
        roi_umap['class'] = roi_umap['sample'].map(self.roi_labels)
        roi_umap['label'] = roi_umap['sample'].str.removeprefix('ROI_')
        
        
        fig,ax = plt.subplots(figsize=(8,6))
        
        # Plot  pixels
        ax.scatter(umap_coords['UMAP1'], umap_coords['UMAP2'],
            c=[class_colors[c] for c in umap_coords['class']],
            s=0.6, alpha=0.15, rasterized=True           
        )
        # Overlay ROI centroids
        for _, row in roi_umap.iterrows():
            ax.scatter(row['UMAP1'], row['UMAP2'],
                color=class_colors[row['class']],
                s=180, edgecolors='black', linewidths=0.5, zorder=5)
            ax.text(row['UMAP1'], row['UMAP2'], row['label'],
                    fontsize=9, fontweight='bold', ha='left', va='bottom', zorder=6)
            
        # # Overlay principal curve
        # ax.plot(
        #     curve_points[:, 0], 
        #     curve_points[:, 1],
        #     color='black', linewidth=2, zorder=4
        # )
        
        ax.set_title('UMAP with ROI Centroids')
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.set_xticks([])
        ax.set_yticks([])
        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/umap_roi_centroids.png'), dpi=150)
        plt.close()
        print('UMAP with ROI centroids generation successful.')


<<<<<<< Updated upstream
=======
        pca_coords = pd.DataFrame(
            self.adata.obsm['X_pca'][:, :3],
            columns=['PC1', 'PC2', 'PC3'],
            index=self.adata.obs.index
        )

        pca_coords['class'] = self.adata.obs['class'].values
        pca_coords['sample'] = self.adata.obs['sample'].values

        classes= pca_coords['class'].unique()

        fig = go.Figure()

        # Plot pixels per class
        for cls in classes:
            mask = pca_coords['class'] == cls
            subset = pca_coords[mask]
            fig.add_trace(go.Scatter3d(
                x=subset['PC1'],
                y=subset['PC2'],
                z=subset['PC3'],
                mode='markers',
                marker=dict(size=2,color=self.roi_class_colors[cls], opacity=0.15)
            ))
        
        # Save interactive html
        fig.write_html(os.path.join(os.path.dirname(self.data_path), 'results/pca_3d.html')) 

    # def get_3d_umap(self):


    def get_roi_level_mst(self):
        '''
        Builds a pseudotime MST on ROI spatial centroids. 
        Nodes = ROIs (from roi_stats), edges weighted by Euclidean distance between centroids
        Pseudotime is computed as shortest path distance from Normal root in the MST

        Returns:
            None, but saves mst_pseudotime.png and stores pseusotime in self.roi_stats['mst_pseudotime']
        '''
        print('Building ROI-level MST...')

        # Step 1: Build distance matrix on spatial centroids
        coords = self.roi_stats[self.good_peptides].values
        dist_matrix = cdist(coords, coords, metric='Euclidean')

        # Step 2: Compute MST
        mst_sparse = minimum_spanning_tree(dist_matrix)
        mst_array = mst_sparse.toarray()

        # Step 3: Build Networkx graph from MST
        G = nx.Graph()
        roi_names = self.roi_stats['roi'].values
        for i, name in enumerate(roi_names):
            G.add_node(i, roi=name, cls=self.roi_stats['class'].iloc[i])
        
        n = len(roi_names)
        for i in range(n):
            for j in range(n):
                w = mst_array[i, j]
                if w > 0: 
                    G.add_edge(i, j, weight=w)

        # Store
        self.mst_graph = G
        self.mst_coords = coords    # median intensities, shape (n,2) 
        
        # Step 4: Select root node
        if self.normal_class and (self.roi_stats['class'] == self.normal_class).any(): 
            normal_mask = self.roi_stats['class'] == 'Normal'
            normal_indices = self.roi_stats[normal_mask].index.tolist()
            normal_coords = coords[normal_mask]
            normal_centroid = normal_coords.mean(axis=0)
            dists_to_centroid = np.linalg.norm(normal_coords - normal_centroid, axis=1)
            local_root = np.argmin(dists_to_centroid)
            root_node = normal_indices[local_root]
            print('Root node identified as point closest to Normal centroid')
            
        else:
            global_centroid = coords.mean(axis=0)
            dists_to_global = np.linalg.norm(coords - global_centroid, axis=1)
            root_node = int(np.argmin(dists_to_global))
            print('Root node identified as point closest to global centroid')

        # Step 5: Pseudotime = shortest path length from root (in MST edge weights)
        path_lengths = nx.single_source_dijkstra_path_length(G, root_node, weight='weight')
        max_len = max(path_lengths.values()) or 1       # avoid dividing by zero
        self.roi_stats['mst_pseudotime'] = [path_lengths.get(i, np.nan) / max_len for i in range(n)]


        # Step 6: Plot
        class_colors = self.roi_class_colors
        pos = {i: (coords[i, 0], coords[i, 1]) for i in range(n)}
        labels = {i: roi_names[i].removeprefix('ROI_') for i in range(n)}
        node_colors = [class_colors.get(G.nodes[i]['cls'], 'black') for i in G.nodes]

        fig,ax = plt.subplots(figsize=(8,6))

        # Draw MST edges manually so we control the axes object
        for i, j in G.edges():
            ax.plot([coords[i, 0], coords[j, 0]],
                    [coords[i, 1], coords[j, 1]],
                    color='black', linewidth=1, zorder=1)

        # Draw ROI centroid nodes
        for i in range(n):
            ax.scatter(coords[i,0], coords[i,1],
                       color=node_colors[i], s=300, edgecolors='black', linewidths=0.5, zorder=5)
            ax.text(coords[i,0], coords[i,1], labels[i], 
                    fontsize=8, fontweight='bold', ha='left', va='bottom', zorder=6)

        # Legend
        for cls, col in class_colors.items():
            ax.scatter([], [], color=col, label=cls, s=80)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        ax.set_title('ROI-level MST')
        ax.set_xticks([])
        ax.set_yticks([])
        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/mst_pseudotime.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()
        print('ROI-level MST pseudotime complete. Saved to results/mst_pseudotime.png')


    def project_mst_onto_umap(self):
        '''
        Projects the ROI-level MST onto UMAP space. 
        The MST topology/connectivity is reconstructed identically, but nodes are at their mean UMAP coordinates instead.

        Returns:
            None, but saves mst_umap_projection.png to results/
        '''
        print('Projecting ROI-level MST onto UMAP space...')

        # Check existance of required columns
        if not hasattr(self, 'mst_graph') or not hasattr(self, 'mst_coords'):
            raise RuntimeError('MST construction data not found, Run get_roi_level_mst() first.')
        if 'UMAP1' not in self.roi_stats.columns: 
            raise ValueError(f'roi_stats is missing columns UMAP1, UMAP2. Make sure run_pixel_dim_reduction() is run first.')
        
        # Create plotting assets
        G = self.mst_graph
        umap_coords = self.roi_stats[['UMAP1', 'UMAP2']].values
        roi_names = self.roi_stats['roi'].values
        n = len(roi_names)

        class_colors = self.roi_class_colors
        node_colors = [class_colors.get(G.nodes[i]['cls'], 'black') for i in G.nodes]
        labels = {i: roi_names[i].removeprefix('ROI_') for i in range(n)}

        # Plot
        fig, ax = plt.subplots(figsize=(8,6))

        for i, j in G.edges():
            ax.plot([umap_coords[i, 0], umap_coords[j, 0]],
                    [umap_coords[i, 1], umap_coords[j, 1]],
                    color='black', linewidth=1, zorder=1)

        for i in range(n):
            ax.scatter(umap_coords[i, 0], umap_coords[i, 1],
                    color=node_colors[i], s=100,
                    edgecolors='black', linewidths=0.5, zorder=5)
            ax.text(umap_coords[i, 0], umap_coords[i, 1], labels[i],
                    fontsize=8, fontweight='bold', ha='left', va='bottom', zorder=6)
            
        ax.scatter(self.adata.obsm['X_umap'][:,0], self.adata.obsm['X_umap'][:,1], 
                   c=[class_colors[c] for c in self.adata.obs['class']], 
                   s=0.6, alpha=0.15, rasterized=True
                )
            
        # Legend 
        for cls, col in class_colors.items():
            ax.scatter([], [], color=col, label=cls, s=80)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
>>>>>>> Stashed changes

    
    def run_diffusion_pseudotime(self):
        print('Running diffusion pseudotime...')
        # Identify Normal root
        print('Identifying normal root...')
        normal_idx = self.adata.obs[self.adata.obs['class']=='Normal'].index[0]
        self.adata.uns['iroot'] = self.adata.obs_names.get_loc(normal_idx)

<<<<<<< Updated upstream
        # Calculate then store diffusion pseudotime
        print('Calculating pseudotime...')
        sc.pp.neighbors(self.adata, n_pcs=10) 
        sc.tl.diffmap(self.adata)
        sc.tl.dpt(self.adata)

        # UMAP + DPT 
        X = self.adata.obsm['X_umap']
        pseudotime = 1- (self.adata.obs['dpt_pseudotime'].values)
        
        fig, ax = plt.subplots(figsize=(8,5))
        plot = plt.scatter(X[:,0], X[:,1], 
                         c=pseudotime, cmap='magma', s=4, alpha=0.4, rasterized=True
                         )
        plt.colorbar(plot, ax=ax, shrink=0.8)
        ax.set_title('Diffusion Pseudotime')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_xticks([])
        ax.set_yticks([])

        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/dpt_pseudotime.png'))
        plt.close()
        print('Diffusion pseudotime reconstruction successful!')
=======
        ax.set_title('ROI-level MST (UMAP Projection)')
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.set_xticks([])
        ax.set_yticks([])
        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/mst_umap_projection.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()
        print('MST UMAP projection complete. Saved to results/mst_umap_projection.png')
    
    
    def run_pseudotime_phate(self):
            import phate
            print('Running PHATE Pseudotime...')
            plt.rcParams['font.serif'] = ['Arial']
            phate_operator = phate.PHATE(knn=10, t=45, verbose=True)
            self.adata.obsm['X_phate'] = phate_operator.fit_transform(self.adata.X)
            sc.pl.embedding(self.adata, basis='phate', color='class', save='_class.png')    # Color by class
            self.adata.obs['phate1'] = self.adata.obsm['X_phate'][:, 0]
            sc.pl.embedding(self.adata, basis='phate', color='phate1', cmap='plasma', save='_pseudotime.png')   # Color by phate trajectory

            print('PHATE generated successfully.')
>>>>>>> Stashed changes


    def run_diffusion_pseudotime_w_curve(self):
        print('Running diffusion pseudotime...')
        sc.settings.verbosity = 3
        from sklearn.neighbors import NearestNeighbors
        # Identify Normal root
        print('Identifying normal root...')
        normal_mask = self.adata.obs['class'] == 'Normal'
        normal_coords = self.adata.obsm['X_pca'][normal_mask]
        normal_centroid = normal_coords.mean(axis=0, keepdims=True)
        # Find the Normal pixel closest to the centroid
        nbrs = NearestNeighbors(n_neighbors=1).fit(normal_coords)
        _, local_idx = nbrs.kneighbors(normal_centroid)
        # Map back to global adata index
        normal_global_indices = np.where(normal_mask)[0]
        self.adata.uns['iroot'] = normal_global_indices[local_idx[0][0]]


        # Calculate then store diffusion pseudotime
        print('Calculating pseudotime...')
        sc.pp.neighbors(self.adata, use_rep='X_pca', n_pcs=3) 
        sc.tl.diffmap(self.adata)
        sc.tl.dpt(self.adata)

        # Fit principal curve in diffusion map space (DC1, DC2) 
        dc_coords = self.adata.obsm['X_diffmap'][:, 1:3]
        pc = PrincipalCurve(k=2)
        pc.fit(dc_coords)

        
        pseudotime_interp, point_inter, order = pc.unpack_params()  


        # Project curve onto UMAP space
        # For each point on pc (in diff space), find nearest neighbor, then look up that neighbor's UMAP coordinate
        neighbors = NearestNeighbors(n_neighbors=1).fit(dc_coords)
        _, indices = neighbors.kneighbors(pc.points)
        curve_umap = self.adata.obsm['X_umap'][indices.flatten()]
        

        # Plot UMAP w DPT colors and DPT principal curve overlaid 
        X = self.adata.obsm['X_umap']
        dpt_pseudotime = self.adata.obs['dpt_pseudotime'].values
        
        fig, ax = plt.subplots(figsize=(8,5))
        plot = plt.scatter(X[:,0], X[:,1], 
                         c=dpt_pseudotime, cmap='jet', s=0.6, alpha=0.4, rasterized=True
                         )
        plt.colorbar(plot, ax=ax, shrink=0.8, alpha=0.4)
        ax.set_title('Diffusion Pseudotime')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_xticks([])
        ax.set_yticks([])

        ax.plot(
            curve_umap[:, 0],
            curve_umap[:, 1],
            color='black', linewidth=2, zorder=4
        )

        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'results/dpt_pseudotime.png'))
        plt.close()
        print('Diffusion pseudotime reconstruction successful!')


    
##### Entire pipeline ##### 
    def allAnalysis(self):
        start = time.time()
        os.makedirs((os.path.join(os.path.dirname(self.data_path), 'results')), exist_ok=True)
        try : 
            self.check_sheet_format()
            self.filter_rois()
            self.get_roi_map()
            self.peptide_sparsity_filter()
            self.normalize_intensities()
            self.diff_expression_test()
            # self.generate_boxplots()
            self.get_roi_stats()
            # self.generate_spatial_heatmap()
            # self.make_hierarchical_clusters()
            # self.get_random_forest_ranking()
            self.create_anndata_object()
            self.run_pixel_dim_reduction()
            # self.compute_mst()
            #self.run_pseudotime_slingshot()
            self.run_diffusion_pseudotime_w_curve()
            #self.run_pseudotime_phate()
        except Exception as e:
            import traceback
            print('Pipeline broke before finishing.')
            traceback.print_exc()
        finally: 
            end = time.time()
            duration = (end-start) // 60    # in minutes
            print(f'Total duration: {duration:.2f} minutes')
