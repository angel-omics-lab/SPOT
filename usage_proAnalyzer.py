
import json
from class_dataAnalysis import SpatialOmicsToolkit

DATA_PATH = r'F:\Bryn\DCIS_Souvik\dcis_data_full.xlsx'
JSON_PATH = r'C:\Users\AngelLab\Documents\GitHub\spatial-proteomics-analyzer\examples\example_roi_labels.json'

spot = SpatialOmicsToolkit(data_path=DATA_PATH, json_path=JSON_PATH)

spot.allAnalysis()