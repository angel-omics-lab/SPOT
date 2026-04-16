'''
This script is an example of how to use the data analysis component of SPOT starting from a worksheet and JSON file of ROIs.
'''
from SPOT.class_dataAnalysis import SpatialOmicsToolkit          # later: from class_dataAnalysis import SpatialOmicsToolkit


WORKSHEET_PATH = r''      # Paste the path to the worksheet
JSON_PATH = r''           # Paste the path to the json file

# NOTE all results will be saved to the same directory as the passed worksheet path within a 'results' folder the pipeline will generate. 

spot = SpatialOmicsToolkit(data_path = WORKSHEET_PATH, roi_labels = JSON_PATH)

spot.allAnalysis()
