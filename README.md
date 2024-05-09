# MALDI-MSI-IMC
Code to integrate MALDI-MSI and Imaging mass cytometry data obtained on a single section.

The here provided code is part of the combined MALDI-MSI IMC methodology proposed in: 

"Integration of Mass Cytometry and Mass Spectrometry Imaging for Spatially Resolved Single Cell Metabolic Profiling" Joana B Nunes, Marieke E Ijsselsteijn, Tamim Abdelaal, Rick Ursem, Manon van der Ploeg, Martin Giera, Bart Everts, Ahmed Mahfouz, Bram Heijs, Noel FCC de Miranda

The methodology combines mass spectrometry imaging-based metabolomics and imaging mass cytometry-based immunophenotyping on the same single tissue section to reveal metabolic heterogeneity within tissues and its association with specific cell populations like cancer cells or immune cells. 

## Use
Two scripts are used for the integration of MALDI-MSI images with Imaging mass cytometry images. The following input is required: 

- MALDI MSI TIFF files for each identified metabolite
- IMC images & cell segmentation masks

The MALDI_MSI_IMC_coregistration.py code is used to coregister the data after which the MALDI_MSI_IMC_integration.py code can be applied to generate metabolite abundances per cell identified by IMC. 


## Example Data
Example data to get started with integrating MALDI-MSI and IMC data is available via https://figshare.com/s/c58ddc70fa8dc0602842




## Cite 
When using the methodology and/or code for any publications, please cite: 

"Integration of Mass Cytometry and Mass Spectrometry Imaging for Spatially Resolved Single Cell Metabolic Profiling" Joana B Nunes, Marieke E Ijsselsteijn, Tamim Abdelaal, Rick Ursem, Manon van der Ploeg, Martin Giera, Bart Everts, Ahmed Mahfouz, Bram Heijs, Noel FCC de Miranda


