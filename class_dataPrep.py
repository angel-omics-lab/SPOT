'''
This script contains class modules for data preparation in SPOT 
'''

import os
import pandas as pd 
import glob
from pyimzml.ImzMLParser import ImzMLParser

class ImzmlConverter:
    def __init__(self, input_path, output_path):
        self.input_path = input_path
        self.output_path = output_path
        self.imzml_dict = {}

    def check_existance(self):
        for file in os.listdir(self.input_path): 
            if file.endswith('.imzML'):
                self.imzml_dict[os.path.basename(file)] = os.path.join(self.input_path, file)
        print(f'{len(self.imzml_dict.keys())} imzML files found in folder.')
    

    def convert_imzml_to_df(self, imzml_path):
        parser = ImzMLParser(imzml_path)
        pixel_level_spectra = []

        for idx, (x,y,z) in enumerate(parser.coordinates):
            mzs, intensities = parser.getspectrum(idx)      # Unpack mzs and intensities so can zip loop over them
            for mz, intensity in zip(mzs, intensities):
                pixel_level_spectra.append({'x':x, 'y':y, 'mz':mz, 'intensity':intensity}) 
        
        df = pd.DataFrame(pixel_level_spectra)
        df = pd.pivot_table(df, 
            index=['x', 'y'], columns=['mz'], values='intensity', fill_value=0, aggfunc='first')
        
        return df 


    def combine_dfs_to_excel(self):
        writer = pd.ExcelWriter(os.path.join(self.output_path,'data_combined.xlsx'), engine='xlsxwriter')
        print('Converting all imzML files to DataFrames...')
        for name, path in self.imzml_dict.items():
            df = self.convert_imzml_to_df(path)
            sheet_name = os.path.splitext(name)[0][:31]     # remove file extension and cap at 31 characters so compatible with excel
            df.to_excel(writer, sheet_name=sheet_name, index=False)
        writer.close()
        print('Successfully saved imzMLs to excel sheet.')
    

    def generateWorksheet(self):
        self.check_existance()
        self.combine_dfs_to_excel()
