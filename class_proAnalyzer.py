'''
This script contains class modules for spatial proteomics analysis 
'''
import pandas as pd 
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests 
from matplotlib import pyplot as plt
import os 
import cv2
import math
import seaborn as sns

''' Example of roi_label format
roi_labels = {
    'roi_1':'DCIS', 
    'roi_2':'IBC', 
    'roi_3':'Normal', 
    .
    .
    . 
}'''

class SpatialProteomicsAnalyzer:
    def __init__(self, data_path, roi_labels):
        self.data_path = data_path
        self.roi_labels = roi_labels        # Contains sheet name of each ROI and its label
        # self.full_data = None
        self.good_peptides = None
        
    def load_and_preprocess(self):
        '''
        Load ROIs, filter out the ones with >25% 0s, then normalize the intensities via log transform
        
        Returns: 
            None, but roi_labels (dict) is updated from filtering. And intensities are normalized. 
        '''
        print('Filtering ROIs...')
        for region, label in self.roi_labels.items():
            data = pd.read_excel(self.data_path, sheet_name=region)
            intensities = data.iloc[:, 4:]      # Skip first 4 columns 
            print(f'Processing {region}...')
            intensities = intensities.fillna(0)        # Change all NaN values to 0 
            zeros = (intensities==0).astype(int).sum().sum()      # Number of 0 intensities
            total = intensities.shape[0]*data.shape[1]     # Total number of intensities
            
            if zeros > (0.25*total):       # Remove ROI from list if >25% of peptides are 0 
                del self.roi_labels[region]
                print(f'Removed {region} with {label} label from list: >25% zero intensities.')
            else: 
                print(f'{region} accepted')
                print(f'Normalizing intensities in {region}')
                data.iloc[:, 4:] = np.log(data.iloc[:, 4:])  # Natural log intensities
                ### NOTE save normalized intensities somewhere?? 
        print('ROI filtering complete. Filtered list: ', self.roi_labels.keys())
        return self.roi_labels


    
    def get_roi_map(self): 
        data = pd.read_excel(self.data_path, sheet_name=None)
        print('Generating roi map...')
        combined = pd.concat(
            [df.assign(roi=roi) for roi, df in data.items() if roi in self.roi_labels],
            ignore_index=True
        )
        agg_dict = {'x': 'mean', 'y': 'mean', **{r: 'mean' for r in self.roi_labels}}
        roi_stats = combined.groupby('Class').agg(agg_dict).reset_index()     # roi_stats columns: ['roi', 'class', 'x', 'y']
        
    

    def compute_peptide_sparsity(self): 
        '''
        Calculate which peptides are expressed across all ROIs, then remove those that arent. 

        Returns: 
            good_peptides(list): list of all peptides that are expressed across ROIs. 
        '''
        
        print('Filtering peptides...')
        peptide_list = [758.4519, 797.4264, 976.4517]       ## NOTE need to update this
        roi_data = {region: pd.read_excel(self.data_path, sheet_name=region).iloc[:, 4] for region in self.roi_labels}       # Read all ROIs into a dict
        all_data = pd.concat(roi_data, names=['ROI', 'sample'])     # Stack all ROI data with ROI labels
        zero_counts_per_peptide = pd.Series(dtype=int) 

        for region in self.roi_labels:
            intensities = pd.read_excel(self.data_path, sheet_name=region).iloc[:, 4:]
            zero_pct = (intensities==0).mean()      # Percentage of 0 intensities for each peptide
            zero_counts_per_peptide = zero_counts_per_peptide.add((zero_pct > 0.20).astype(int), fill_value=0)      # Label peptide for region with 0 if > 20% of values are 0 

        # Keep peptides that are present in >= 10% ROIs 
        threshold = len(self.roi_labels)*0.10       
        self.good_peptides = zero_counts_per_peptide[zero_counts_per_peptide <= threshold].index.tolist()

        print('Peptide filtering complete. Filtered list: ', self.good_peptides)
        return [float(p) for p in self.good_peptides]
    
    


    def diff_expression_test(self): 
        '''
        Identifies peptides that have differential expression between classes (in this instance DCIS v IBC) via 
        Kruskal-Wallis test with Benjamini-Hochberg FDR validation. Generates a boxplot of 
        
        Returns:
            good_peptides(list) updated so it is now only those differentially expressed.
        '''
        data = pd.read_excel(self.data_path, sheet_name=None)
        results_list = []
        # Identify significant peptide 
        print('Identifying differentially expressed peptides...')
        for peptide in self.good_peptides:
            group1 = [data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'DCIS']     # Mean intensities of peptide in DCIS (class 1) rois
            group2 = [data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'IBC']      # Mean Intensities of peptide in IBC (class 2) rois
            # Run KW
            H_statistic, p_value = stats.kruskal(group1, group2)
            # Add to data frame 
            results_list.append({
            'peptide': peptide, 
            'pvalue_ols': p_value
            })
        results_df = pd.DataFrame(results_list) 
        # Run BH FDR, add q val to df 
        reject, q_values, _, _ = multipletests(results_df['pvalue_ols'], alpha=0.05, method='fdr_bh')
        results_df['qvalue_bh'] = q_values
        # Save peptides whose q value is <= 0.05
        self.good_peptides = self.good_peptides = results_df[results_df['qvalue_bh'] <= 0.05]['peptide'].tolist()
        print('Differential analysis complete. Significant peptides: ', self.good_peptides)
        print(f'{len(reject)} peptides removed: non-significant')
        
    

    def generate_boxplots(self): 
        '''
        Takes the good_peptides list and generates a grid of boxplots where each boxplot compares mean intensities between 2 classes for each peptide. 
        In this case, the two classes we choose are IBC and DCIS, but this can be changed. 
        '''
        data=pd.read_excel(self.data_path, sheet_name=None)
        print('Generating box plots for significant peptides...')
        n_peptides = len(self.good_peptides)
        n_cols = 5
        n_rows = math.ceil(n_peptides/n_cols)

        fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols*4, n_rows*4))        # Grid of box plots
        axs = axs.flatten()     # Flatten so indexing is easier
        # Give overall plot name, axis, legend 
        for i, peptide in enumerate(self.good_peptides):
            plot_data = (
                [{'Intensity': val, 'Class':'DCIS'} for val in [data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'DCIS' ]]
                +
                [{'Intensity': val, 'Class':'IBC'} for val in [data[region][peptide].mean() for region in self.roi_labels if self.roi_labels[region] == 'IBC' ]]
            )
            plot_df = pd.DataFrame(plot_data)

            sns.boxplot(data=plot_df, x='Class', y='', 
                        hue='Class', palette={'DCIS':'blue', 'IBC':'orange'},
                         ax=axs[i], legend=False
            )
            axs[i].set_title(peptide)

        # Hide unused subplots
        for j in range(i+1, len(axs)): axs[j].set_visible(False)
            
        # Shared legend
        labels = [plt.Rectangle((0,0),1,1, color=c) for c in ['blue', 'orange']]
        fig.legend(labels, ['DCIS', 'IBC'], loc='upper right')

        # Shared title, yaxislabel
        fig.suptitle('Differentially Expressed Peptides (DCIS v IBC)')
        fig.supylabel('ln(Mean Intensity)')
        fig.tight_layout()

        plt.savefig(os.path.join(os.path.dirname(self.data_path), 'sig_peptide_boxplots.png'))
        plt.show()

        return self.good_peptides
    



    def plot_spatial_heatmap(self):
        '''
        Generates a spatial heatmap (on the original image) for each peptide based on the centroid of each ROI and mean intensity of each peptide.

        Returns:
            heatmap_X (png) where X is the peak, makes one for every good_peptide 
            roi_stats (np array) 
        '''
        data = pd.read_excel(self.data_path, sheet_name=None)
        print('Generating heatmaps for each peptide...')
        # Calculate spatial centroid of each roi
            # Concatenate all ROI sheets, tagging each row with its ROI label
        combined = pd.concat(
            [df.assign(roi=roi) for roi, df in data.items() if roi in self.roi_labels],
            ignore_index=True
        )
            # Groupby computes centroid AND all peptide means in one pass (vector operation)
        agg_dict = {'x': 'mean', 'y': 'mean', **{p: 'mean' for p in self.good_peptides}}
        roi_stats = combined.groupby('roi').agg(agg_dict).reset_index()     # roi_stats columns: ['roi', 'x', 'y', peptide_1, peptide_2, ...]
        
        # Generate heatmaps for each peptide
        for peptide in self.good_peptides:
            heatmap = plt.figure()
            plt.scatter(roi_stats['x'], roi_stats['y'], 
                        c=roi_stats[peptide], cmap='RdBu_r', s=50)
            plt.colorbar()
            plt.xlabel(False)
            plt.ylabel(False)
            plt.title(peptide)
            plt.savefig(os.path.join(os.path.dirname(self.data_path), f'heatmap_{peptide}.png'))
        print('Spatial heatmap generation successful. Saved to: ', os.path.dirname(self.data_path))
        return roi_stats        


    # def make_hierarchical_clusters(self, roi_stats):
    #     '''
        

    #     '''

            
        