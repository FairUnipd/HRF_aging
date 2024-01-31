# HRF_aging
Repository including sample code and additional material to reproduce the results outlined in the manuscript "Revealing the spatial pattern of brain hemodynamic sensitivity to healthy aging through sparse DCM"

Data: This folder contains essential data required for running the scripts.

Scripts: Within this folder, you'll find sample code aimed at generating figures (prefixed with "Figure..") and executing the modeling pipeline. 
It is recommended to follow this order when running the scripts:

1. extract_HRF_features.m
2. main_multicollinearity_multivariate_selection.m
3. save_data_for_Figure4.m

Functions: Here, you'll find functions utilized by the scripts to streamline processes.

Schaefer_62Parcels_7Networks_order_FSLMNI152_2mm.nii.gz: Clustered functional atlas whose labels are included in Table 1-1 (MNI 152 FSL standard space).
