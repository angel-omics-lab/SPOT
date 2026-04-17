
from SPOT.class_dataAnalysis import SpatialOmicsToolkit

DATA_PATH = r'/Users/bryngerding/Desktop/MUSC/dcis_data_full.xlsx'
JSON_PATH = r'/Users/bryngerding/Desktop/MUSC/example_roi_labels.json'

spot = SpatialOmicsToolkit(data_path=DATA_PATH, json_path=JSON_PATH)

spot.allAnalysis()