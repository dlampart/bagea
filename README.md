Implements BAGEA (Bayesian Annotation Guided eQTL Analysis), which allows to fit eQTL data to directed and undirected genomic annotations.

An installation procedure depending on the devtools package can be found at

	`bagea_install.R`

Importantly, for the installation to succeed  an environment variable BAGEA_PATH has to be set and point to a folder where auxilliary data will be installed (consult the scripts on details).
By default the installation, will just download allready preprocessed files. But the essential data can also be downloaded from the original source and preprossed (see documentation of the bagea functions `install_external_data` and `install_1KG_data`)
The full installation downloads and preprocesses 1KG data, so it can take a while.
To run BAGEA efficiently, we recommend to run R using on optimized BLAS library (openBLAS, MKL, osX veclib, etc). This will make some matrix operations much more efficient. How to do this is system dependent.

An a full example can be found at.

	`bagea_tutorial.R`

 Stepping throuh this example will:
1)  Download GTEx and preprocess summary statistics, 
2)  Download and preprocess directed annotations (from Reshef et al. 2018) 
3)  Prepare a bagea input object from those input data
4)  Run BAGEA
5)  calculate MSE_dir from the results.

To get additional information of bagea functions please consult their help page e.g.:

`help(run_bagea)`


