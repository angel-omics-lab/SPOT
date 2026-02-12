'''
This script contains class modules for spatial proteomics analysis 
'''
import pandas as pd 
import numpy as np

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

    
        
    
    def compute_peptide_sparsity(self): 
        '''
        Calculate which peptides are expressed across all ROIs, then remove those that arent. 
        
        Parameters: 
            threshold(float)[0-1]: Percentage of ROIs where specified peptide is absent in. Those below it are not added to 'good' list

        
        Returns: 
            good_peptides(list): list of all peptides that are expressed across ROIs. 
        '''
        
        print('Filtering peptides...')
        peptide_list = [758.4519, 797.4264, 976.4517]
        roi_data = {region: pd.read_excel(self.data_path, sheet_name=region).iloc[:, 4] for region in self.roi_labels}       # Read all ROIs into a dict
        all_data = pd.concat(roi_data, names=['ROI', 'sample'])     # Stack all ROI data with ROI labels

        zero_pct_per_roi = (all_data == 0).groupby(level='ROI').mean()      # Calculate zero percentage per peptide per ROI 

        rois_with_low_zeros = (zero_pct_per_roi < 0.2).sum(axis=0)     # Count ROIs where each peptide has <20% zeros

        # Count peptides that are abundant in > 10% ROIs
        threshold = len(self.roi_labels)*0.1 
        self.good_peptides = rois_with_low_zeros[rois_with_low_zeros>=threshold].index

        print('Peptide filtering complete. Filtered list: ', self.good_peptides)
        return self.good_peptides
    
    
    
            
        