<img width="1052" height="333" alt="Image" src="https://github.com/user-attachments/assets/0c105274-b155-44af-934f-71df72cf7b29" />

SPOT (Spatial Omics Toolkit) is an end-to-end spatial omics pipeline. open-sourced, pre-existing statistical methods, like random forest models and principal component analysis (PCA), to provide researchers with a straightforward and comprehensive pipeline for both class based and spatial omics analysis.

## Configuring
To make sure all the necessary packages are installed without conflict, it is advised to run SPOT in a virtual environment.
This means that every time you want to run SPOT, it needs to be done with the virtual environment activated. 
```python
# Create and activate a virtual environment
python -m venv venv_spot
source venv_spot/bin/activate   # for Mac/Linux
venv_spot\Scripts\activate      # for Windows

# Once it is activated, its name should appear at the start of the command prompt. 

# Install the SPOT package 
pip install git+https://github.com/angel-omics-lab/SPOT.git
```

## Usage
1. Prepare the input data. If you have already have an excel file with ROIs separated into unique sheets within the file, you can skip this step and proceed to Step 2. 
    a. Given imzML files for each region of interest, they can be converted to a pipeline-specific worksheet with the methods in class_dataPrep. You can see an example script using these methods in example_dataPrep. You can directly download a copy of the file and change the path variables or incorporate it a script of your own.
    b. Once your worksheet has been created, create a JSON file that lists the ROIs and their disease classification. The ROI name in the JSON file and the individual sheet names contain ROI-specific data need to be identical. You can see examples of both a properly formatted JSON file and worksheet in the example folder. 
2. Create a script to run the pipeline. You can see an example script using these methods in example_dataAnalysis. You can directly download a copy of the file and change the path variables or incorporate it a script of your own.
3. Run SPOT analysis. Make sure to run the script while in the virtual environment you created earlier. 
```python
# Within the activated virtual environment, enter the following command:
python your/path/to/script.py    # NOTE: on a Mac, you may need to have python3 as the prefix instead of python 
```

## Issues
If you are having problems running SPOT, please add an Issue to this repository and/or contact gerdingb@musc.edu
