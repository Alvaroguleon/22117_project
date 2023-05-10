# 22117_project
This repository contains the necessary code to process the output from MutateX and ClinVar information into the principal component biplots. 

The contents are structured as follows:
<pre> 
├── 00_datasets
│   ├── processed
│   │   ├── clinvar.txt
│   │   ├── data.csv
│   │   ├── metadata.csv
│   │   └── pca_data.csv
│   └── raw
│       ├── E_Sequence.netsurfp.txt
│       ├── F_Sequence.netsurfp.txt
│       ├── clinvar_result_hba1.txt
│       ├── clinvar_result_hba2.txt
│       ├── clinvar_result_hbb.txt
│       ├── energies_alpha.csv
│       └── energies_beta.csv
├── 01_scripts
│   ├── 00_run.R
│   ├── 01_processing_clinvar.R
│   ├── 02_processing_energies.R
│   ├── 03_PCA.R
│   └── 04_Extra_plots_exploratory.R
├── 02_figures
│   ├── expo_full.png
│   ├── expo_zoom.png
│   ├── res_clin_full.png
│   ├── res_clin_zoom.png
│   ├── res_mut_full.png
│   ├── res_mut_zoom.png
│   ├── res_wt_full.png
│   ├── res_wt_zoom.png
│   ├── residues.png
│   └── residues_pca.png
├── LICENSE
└── README.md

</pre>
