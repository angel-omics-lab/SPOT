import pandas as pd 
import numpy as np 
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests 

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


##### differential expression test #####
data = pd.read_excel(DATA_PATH, sheet_name=None)
good_peptides = [758.4519, 797.4264, 976.4517]
results_list = []

print('Identifying differentially expressed peptides...')
for peptide in good_peptides:
    group1 = [data[region][peptide].mean() for region in roi_labels if roi_labels[region] == 'DCIS']     # Mean intensities of peptide in DCIS (class 1) rois
    group2 = [data[region][peptide].mean() for region in roi_labels if roi_labels[region] == 'IBC']      # Mean Intensities of peptide in IBC (class 2) rois
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
good_peptides = good_peptides = results_df[results_df['qvalue_bh'] <= 0.05]['peptide'].tolist()
print('Differential analysis complete. Significant peptides: ', good_peptides)