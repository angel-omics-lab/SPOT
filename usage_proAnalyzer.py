
import json
from class_proAnalyzer import SpatialOmicsAnalyzer

DATA_PATH = r'F:\Bryn\DCIS_Souvik\dcis_data_full.xlsx'
#DATA_PATH = r'C:\Users\AngelLab\Desktop\pipeline-data_PMA\DCIS065_A5_HE_total_modified.xlsx'
#JSON_PATH = r'C:\Users\AngelLab\Documents\GitHub\spatial-proteomics-analyzer\roi_labels.xlsx'

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
    'ROI_111':'IBC', 
    'ROI_112':'IBC',
    'ROI_113':'IBC', 
    'ROI_114':'IBC',
    'ROI_115':'IBC', 
    'ROI_116':'IBC',
    'ROI_117':'IBC', 
    'ROI_118':'IBC',
    'ROI_119':'IBC', 
    'ROI_120':'IBC',
    'ROI_121':'DCIS',
    'ROI_122':'DCIS', 
    'ROI_123':'DCIS',
    'ROI_124':'DCIS',
    'ROI_125':'IBC', 
    'ROI_126':'DCIS',
    'ROI_127':'DCIS',
    'ROI_128':'Normal'  
}

SPOT = SpatialOmicsAnalyzer(data_path=DATA_PATH, roi_labels=roi_labels)

SPOT.allPipeline()