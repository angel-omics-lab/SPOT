'''
This script is an example of how to use the data preparation component of SPOT to convert a set of imzML files into a singular worksheet.
'''

from SPOT.class_dataPrep import ImzmlConverter

FOLDER_PATH = r'F:\Bryn\imzML_extraction\heavy dcis\region_imzmls'       # Paste the path to the folder containing all relevant imzML files

converter = ImzmlConverter(input_path = FOLDER_PATH)

converter.generateWorksheet()

# NOTE before you run analysis, you may need to change the names of the individual sheets so they match what is in your JSON file. 
# Currently, the above function is written so the sheets are named after the imzML files they came from. 
