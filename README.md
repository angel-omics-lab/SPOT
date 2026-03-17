# README

SPOT (Spatial Omics Toolkit) is an end-to-end spatial omics pipeline 

## Installation (in your terminal; will only do this once)

```python
# Create and activate a virtual environment
python -m venv venv_spot
source venv_spot/bin/activate   # for Mac/Linux
venv_spot\Scripts\activate      # for Windows

# Install the SPOT package
pip install https://github.com/bgerd02/spatial-proteomics-analyzer.git

```

## Usage 

# 1. Copy usage_proAnalyzer.py and change the DATA_PATH variable to the path to your spreadsheet
# 2. Create a .json file with a class listed for every ROI. Make sure the names of the ROIs match the names of the sheets in the spreadsheet. You may copy the example json file and replace it with your information. This json file should be in the same directory as the spreadsheet. 
# 3. With the virtual environment activated, run the usage script: 
```python
python usage_proAnalyzer.py
```

## Issues
If you are having problems running SPOT, please add an issue to this repository and/or contact gerdingb@musc.edu
