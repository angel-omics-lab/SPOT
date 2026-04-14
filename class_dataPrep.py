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
    


    def aggregate_columns(path):
        print('Filtering spreadsheet to aggregate duplicate columns...')
        spreadsheet = pd.read_excel(path, sheet_name=None)

        # Iteratively go through each sheet in spreadsheet
        for sheet_name, df  in spreadsheet.items():
            print(f'Processing sheet: {sheet_name}')
            output_cols = []

            # Convert peak col names to float, then create groups of peaks based on their rounded (1) values
            mz_cols = df.columns[2:]
            mz_groups = df[mz_cols].columns.groupby(pd.Series(mz_cols, dtype=float).round(1))
                
            # Iteratively go through each group and decide which ones to combine
            for peak, col_names in mz_groups.items():
                # Check for nonzeros in each row
                group = df[col_names]
                row_nonzero_count = (group != 0).sum(axis=1)    # Series w a nonzero count per row 
                # Aggregate columns if there is less than 10% overlap 
                conflict_rate = (row_nonzero_count > 1).mean()
                if conflict_rate < 0.10:
                    combined = group.max(axis=1)
                    combined_name = np.median(col_names.astype(float))
                    output_cols.append(combined)
                else:
                    output_cols.append(group)

            df2 = pd.concat(output_cols, axis=1)
            df = df2
            print(f'Column aggregation for {sheet_name} complete. Combined {len(mz_cols)} peaks into {df2.shape[1]} peaks')
            

        print('Worksheet filtering complete. Data is ready to to passed to SPOT analysis.')


    def generateWorksheet(self):
        self.check_existance()
        self.combine_dfs_to_excel()
        self.aggregate_columns()