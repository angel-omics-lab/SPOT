import pandas as pd 
import numpy as np 

DATA_PATH = r'C:\Users\AngelLab\Documents\GitHub\spatial-proteomics-analyzer\fake_data_for_testing.xlsx'

roi_labels = {
    'ROI_101':'DCIS', 
    'ROI_102':'IBC', 
    'ROI_103':'Normal', 
}

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
print('Filtering peptides...')
peptide_list = [758.4519, 797.4264, 976.4517]
zero_counts_per_peptide = pd.Series(dtype=int)      # Will track how many ROIs each peptide 'failed' in 

for region in roi_labels: 
    intensities = pd.read_excel(DATA_PATH, sheet_name=region).iloc[:, 4:]
    zero_pct = (intensities==0).mean()      # Percentage of 0 intensities for each peptide 
    zero_counts_per_peptide = zero_counts_per_peptide.add((zero_pct > 0.20).astype(int), fill_value=0)      # Label peptide for region with 0 if > 20% of values are 0

# Keep peptides above 10% threshold
threshold = len(roi_labels)*0.1
good_peptides = zero_counts_per_peptide[zero_counts_per_peptide <= threshold].index.tolist()


print('Peptide filtering complete. Filtered list: ', good_peptides)

