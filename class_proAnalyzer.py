'''
This script contains class modules for spatial proteomics analysis 
'''
import pandas as pd 

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
        Load ROIs, filter out the ones with >80% 0s, then normalize the intensities via log transform
        
        Returns: 
            None, but roi_labels (list) is updated from filtering. And intensities are normalized. 
        '''
        print('Filtering ROIs...')
        for region, label in self.roi_labels:
            data = pd.read_excel(self.data_path, sheet_name=region)     ### NOTE may need to change dep on how sheets are formatted
            print(f'Processing {region}...')
            na = data.isnull().sum()        # Number of NaN intensities     ### NOTE maybe add or zero to na count??
            
            if na > (0.25*data.shape[1]):       # Remove ROI from list if >25% of peptides are null 
                self.roi_labels.remove(region)
                print(f'Removed {region} with {label} label from list: >25% NaN intensities.')
                break       ### NOTE check to make sure this does what I think it does
            else: print(f'{region} accepted')
        print('ROI filtering complete. Filtered list: ', self.roi_labels)

            ### NOTE add normalizing here
    
    
    
    def compute_sparsity(self, threshold=0.2): 
        '''
        Calculate which peptides are expressed across all ROIs, then remove those that arent. 
        
        Parameters: 
            threshold(float)[0-1]: Percentage of ROIs where specified peptide is absent in. Those below it are not added to 'good' list

        
        Returns: 
            good_peptides(1d-array): list of all peptides that are expressed across ROIs. 
        '''
        
        print('Filtering peptides...')
        for peptide in peptide_list and region in self.roi_labels:
            intensities = []
            for region in self.roi_labels:
                # calculate intensity of peptide in region
                intensities.append(value)
            if intensities.isnull().sum(0) > (threshold*intensities.shape[0]):    # Only add peptide to good list if it appears in > threshold% of ROIs
                print(f'Peptide {peptide} did not meet threshold for ROI presence')
                break
            else: 
                self.good_peptides.append(peptide)
                print(f'Peptide {peptide} accepted')
            
        print('Peptide filtering complete. Filtered list: ', self.good_peptides)
    
    
    
            
        