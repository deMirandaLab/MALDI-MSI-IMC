import os
import pandas as pd
import numpy as np
import re


script_dir = os.path.dirname(os.path.abspath(__file__))
print(script_dir)

# Calculate the percentages per cell per image
MsiFolder = os.path.join(script_dir, "Input/MSI_pixel_data/")
cell_intensity_folder = os.path.join(script_dir, "Input/Intensityfiles/")
PhenotypeFolder = os.path.join(script_dir, "Input/Export_categories_24/")

MsiFiles = [f for f in os.listdir(MsiFolder) if f.endswith(".csv")]
PhenotypeFiles = [f for f in os.listdir(PhenotypeFolder) if f.endswith(".csv")]
cell_intensity_files = [f for f in os.listdir(cell_intensity_folder) if f.endswith(".csv")]


All_cells = pd.DataFrame()

for i in range(len(PhenotypeFiles)):
    phenotypes = pd.read_csv(os.path.join(PhenotypeFolder, PhenotypeFiles[i]), header=None)
    metrics = pd.read_csv(os.path.join(cell_intensity_folder, cell_intensity_files[i]))
    MSI_pixel = pd.read_csv(os.path.join(MsiFolder, MsiFiles[i]))

    cell_id_list = range(1, len(phenotypes) + 1)
    #print(len(cell_id_list))
    MSI_all = pd.DataFrame()

    MSI_perc = MSI_pixel.copy()
    MSI_perc.iloc[:, 4:115] = MSI_perc.iloc[:, 4:115].div(MSI_perc.iloc[:, 4:115].sum(axis=1), axis=0)

    for v in cell_id_list:
        cell_n = v
        cell_end = f" {cell_n}]"
        cell_middle = f" {cell_n},"
        cell_start = f"[{cell_n},"
    
        regex_pattern = f"{re.escape(cell_end)}|{re.escape(cell_middle)}|{re.escape(cell_start)}"
        #print(regex_pattern)
        
        MSI_perc_cell = MSI_perc[MSI_perc['Cell_idx'].str.contains(regex_pattern, regex=True)]
        #print(MSI_perc_cell)
        if not MSI_perc_cell.empty:
            cell_count = MSI_perc_cell['Cell_idx'].str.count(str(cell_n))
            MSI_cell_comb = MSI_perc_cell.iloc[:, 4:116].mul(cell_count, axis=0).sum() / cell_count.sum()

            MSI_cell_comb['Cell_idx'] = cell_n
            MSI_cell_comb = MSI_cell_comb.to_frame().transpose().head()
            #print(MSI_cell_comb)
            MSI_all = pd.concat([MSI_all, MSI_cell_comb])


    phenotypes.columns = ["phenotype"]
    phenotype_intensity = pd.concat([phenotypes, metrics], axis=1)
    phenotype_intensity['Cell_idx'] = range(1, len(phenotype_intensity) + 1)
    print(MSI_all)
    IMC_MSI_Merge = pd.merge(phenotype_intensity, MSI_all, on="Cell_idx")

    sample = os.path.splitext(PhenotypeFiles[i])[0]
    IMC_MSI_Merge['Sample'] = sample
    All_cells = pd.concat([All_cells, IMC_MSI_Merge])


All_cells.to_csv(os.path.join(script_dir, "All_samples_IMC_MSI_merge_per_cell_perc.csv"))