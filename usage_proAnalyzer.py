
from SPOT.class_dataAnalysis import SpatialOmicsToolkit
from SPOT.class_dataPrep import ImzmlConverter

# DATA_PATH = r'/Users/bryngerding/Desktop/MUSC/dcis_data_full.xlsx'
# JSON_PATH = r'/Users/bryngerding/Desktop/MUSC/example_roi_labels.json'
# IMZML_FOLDER = r'/Users/bryngerding/Desktop/MUSC/imzmls'

# DATA_PATH = r'C:\Users\AngelLab\Desktop\SPOT_stuff\example_data.xlsx'
# JSON_PATH = r'C:\Users\AngelLab\Desktop\SPOT_stuff\example_roi_labels.json'

IMZML_FOLDER = r'C:\Users\AngelLab\Desktop\SPOT_LuHan_data'
DATA_PATH=r'C:\Users\AngelLab\Desktop\SPOT_LuHan_data\data_combined.xlsx'
JSON_PATH=r'C:\Users\AngelLab\Desktop\SPOT_LuHan_data\roi_labels.json'

spot = SpatialOmicsToolkit(data_path=DATA_PATH, json_path=JSON_PATH)
spot.allAnalysis() 

# dp = ImzmlConverter(input_path= IMZML_FOLDER)
# dp.generateWorksheet()
