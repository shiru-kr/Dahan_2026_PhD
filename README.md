# An Investigation into Brain Structure-Function Connectivity Patterns in Bipolar Disorder

This repository contains the code, scripts, and analysis pipelines developed as part of my PhD thesis:

**An Investigation into Brain Structure-Function Connectivity Patterns in Bipolar Disorder**

The project investigates how structural and functional brain connectivity interact in bipolar disorder using multimodal MRI, network neuroscience, and imaging genetics approaches. The repository is organised according to the individual studies presented in the thesis, allowing each analysis pipeline to be reproduced independently.

This repository is intended to promote transparency and reproducibility by providing the analysis scripts used throughout the thesis.

The neuroimaging and genetic datasets analysed in this thesis are subject to data-sharing agreements and ethical approvals. As a result, the raw data cannot be publicly distributed through this repository.

Researchers with appropriate access permissions should obtain the data directly from the original data providers before running these analyses.

This work was completed as part of my PhD research at the University of Galway.

## Repository Structure

### Folder 1 - Study 1
**Treatment Resistance in Bipolar Disorder is Associated with Intrinsic Network Power**

Contains scripts and code for Study 1 sensitivity analyses. The study performed Independent Component Analysis (ICA) on resting-state functional MRI to investigate intrinsic brain network activity and its relationship with treatment response in bipolar disorder. Preprocessing was performed in SPM and ICA was performed with the GIFT toolbox (GUI). The python scripts in the folder were used to assess the stability of the power spectra results and for plotting.

Dahan, S., Casburn, M., Fernandes, A., Quirke, J., Corley, E., Broin, P. Ó., ... & Sweet, J. (2026). Treatment Resistance in Bipolar Disorder is Associated with Intrinsic Network Power. Journal of Affective Disorders Reports, 101107.

### Folder 2 - Study 2
**Altered Structure-Function Coupling and its Genetic Underpinnings in Bipolar Disorder**

Contains the analysis pipeline for regional structure-function coupling. Python scripts are provided to calculate regional structure-function coupling from structural and functional connectivity matrices, statistical analysis code, plotting code and genome-wide association analyses (GWAS) scripts.

Dahan, S., Ressin, L., Baltramonaityte, V., Corley, E., Laighneach, A., Walton, E., ... & Cannon, D. M. (2025). Altered Frontoparietal, Temporal and Sensorimotor Structure-Function Coupling and Its Genetic Underpinnings in Bipolar Disorder. bioRxiv, 2025-12.

### Folder 3 - Study 3
**Reduced Cross-Modal Modularity and Flexibility of Network Organisation in Bipolar Disorder**

Contains FSL scripts to preprocess and rescale T1-weighted MRI in FSL, MATLAB and SPM scripts for obtaining grey matter similarity networks, and Python code for multiplex network analysis integrating white matter connectivity, grey matter similarity, and functional connectivity to investigate multiplex community organisation.
