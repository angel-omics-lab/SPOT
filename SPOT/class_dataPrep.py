'''
This script contains class modules for data preparation in SPOT 
'''

import os
import pandas as pd 
import numpy as np
import glob
from pyimzml.ImzMLParser import ImzMLParser
import xml.etree.ElementTree as ET

class ImzmlConverter:
    def __init__(self, input_path):
        self.input_path = input_path
        self.imzml_dict = {}

    def check_existance(self):
        for file in os.listdir(self.input_path): 
            if file.endswith('.imzML'):
                self.imzml_dict[os.path.basename(file)] = os.path.join(self.input_path, file)
        print(f'{len(self.imzml_dict.keys())} imzML files found in folder.')
    

    def get_3d_positions(self, imzml_path):
        tree=ET.parse(imzml_path)
        root = tree.getroot()
        ns = {'ms': 'http://psi.hupo.org/ms/mzml'}  # namespace from file's xmlns
        positions = []
        for scan in root.findall('.//ms:scan', ns):
            x = scan.find('ms:userParam[@name="3DPositionX"]', ns)
            y = scan.find('ms:userParam[@name="3DPositionY"]', ns)
            if x is not None and y is not None: 
                positions.append((float(x.get('value')), float(y.get('value'))))
        return positions


    def convert_imzml_to_df(self, imzml_path):
        parser = ImzMLParser(imzml_path)
        pixel_level_spectra = []

        positions = self.get_3d_positions(imzml_path)
        
        print(f"Scans from XML: {len(positions)}, Spectra from parser: {len(parser.coordinates)}")      # Sanity check

        for idx, (x,y) in enumerate(positions):
            mzs, intensities = parser.getspectrum(idx)      # Unpack mzs and intensities so can zip loop over them
            for mz, intensity in zip(mzs, intensities):
                pixel_level_spectra.append({'x':x, 'y':y, 'mz':mz, 'intensity':intensity}) 
        
        df = pd.DataFrame(pixel_level_spectra)
        df = pd.pivot_table(df, 
            index=['x', 'y'], columns=['mz'], values='intensity', fill_value=0, aggfunc='first')
        
        return df 


    def combine_dfs_to_excel(self):
        writer = pd.ExcelWriter(os.path.join(self.input_path,'data_combined.xlsx'), engine='xlsxwriter')
        print('Converting all imzML files to DataFrames...')
        for name, path in self.imzml_dict.items():
            df = self.convert_imzml_to_df(path)
            df = df.reset_index()
            sheet_name = os.path.splitext(name)[0][:31]     # remove file extension and cap at 31 characters so compatible with excel
            df.to_excel(writer, sheet_name=sheet_name, index=False)
        writer.close()
        print('Successfully saved imzMLs to excel sheet.')
    

    def aggregate_columns(self):
        print('Filtering spreadsheet to aggregate duplicate columns...')
        path = os.path.join(self.input_path,'data_combined.xlsx')
        spreadsheet = pd.read_excel(path, sheet_name=None)
        
        writer = pd.ExcelWriter(path, engine='xlsxwriter')

        # Iteratively go through each sheet in spreadsheet
        for sheet_name, df  in spreadsheet.items():
            print(f'Processing sheet: {sheet_name}')
            output_cols = []
            
            # Preserve x, y columns explicitly
            xy_cols = df.columns[:2]
            mz_cols = df.columns[2:]

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
                if conflict_rate >= 0.10 : 
                    # Pairwise compatibility check
                    ungrouped = list(col_names)
                    while ungrouped:
                        safe_group = [ungrouped.pop(0)]
                        i = 0 
                        while i < len(ungrouped):
                            candidate = ungrouped[i]
                            test_group = df[safe_group + [candidate]]
                            row_nonzero_count = (test_group != 0).sum(axis=1)
                            if (row_nonzero_count > 1).mean() < 0.10:
                                safe_group.append(candidate)
                                ungrouped.pop(i)
                            else:
                                i += 1
                        combined = df[safe_group].max(axis=1)
                        combined.name = np.median(np.array(safe_group, dtype=float))
                        output_cols.append(combined)
                
                else :      # conflict rate < 0.10
                    combined = group.max(axis=1)
                    combined.name = np.median(col_names.astype(float))
                    output_cols.append(combined)
                
            # Re-attach x, y columns at the front
            df2 = pd.concat([df[xy_cols].reset_index(drop=True),
                        pd.concat(output_cols, axis=1).reset_index(drop=True)], axis=1)

            df2.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f'Column aggregation for {sheet_name} complete. Combined {len(mz_cols)} peaks into {df2.shape[1]-2} peaks')
            
        writer.close()  
        print('Worksheet filtering complete. Data is ready to be passed to SPOT analysis.')


    def generateWorksheet(self):
        self.check_existance()
        self.combine_dfs_to_excel()
        self.aggregate_columns()