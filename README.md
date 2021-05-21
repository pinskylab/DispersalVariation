# DispersalVariation
Project measuring variability in larval dispersal for clownfish Amphiprion clarkii, published here: https://doi.org/10.1111/mec.15732 

This repository is organized into 3 folders:
1. code
  -This contains the R code to produce the results in the manuscript “Quantifying dispersal variability among nearshore marine populations”. My original coding was done in Jupyter Notebooks (IR Kernel). These notebooks are included in the folder "WorkingJupyterNotebooks", and can be viewed as HTML using the web site https://nbviewer.jupyter.org/. The .R files were produced by downloading these notebooks as .R files. 
2. data
  -This contains the data to run the analysis using the code. The genomic data here is a filtered to use in Colony2. It also contains a backup of the SQL database with all data for the clownfish field and wetlab work in the Pinsky lab.
3. genomics
  -This is the version of https://github.com/pinskylab/genomics as it was at the time of this project. It contains all molecular wetlab metadata, raw genomic data, bioinformatics code, and protocols.
  - The final VCF created from the "genomics" repository filtering pipeline is [FIL_6.recode.vcf](https://github.com/pinskylab/genomics/blob/master/data/FIL_6.recode.vcf.zip)
  - The final genepop file used to make Colony2 input files is [seq33-03_norecap.gen](https://github.com/pinskylab/DispersalVariation/blob/master/data/seq33-03_norecap.gen)
  
  Please contact Katrina at kat.catalano@rutgers.edu with any questions.
  
  ## Workflow
  1. Protocols for sample collection in the field can be found within this separate repository: https://github.com/pinskylab/field
  2. Protocols for sequencing libraries made from tissue samples can be found in the "genomics" directory within this repository, or wiithin this separate repository: https://github.com/pinskylab/genomics
  3. Parentage analysis
  - [Make input files](https://github.com/pinskylab/DispersalVariation/blob/master/code/colony_prep2012_2018.R)
  - [Add metadata and summarise parentage results](https://github.com/pinskylab/DispersalVariation/blob/master/code/process_colony_dyad_2.0.R)

  4. Dispersal analysis
  - [Fit kernels](https://github.com/pinskylab/DispersalVariation/blob/master/code/kernel_fit_file_prep.R)
  - [Plot kernels](https://github.com/pinskylab/DispersalVariation/blob/master/code/kernel_plotting.R)
  - [Plot mean dispersal distance, meadian, and kurtosis](https://github.com/pinskylab/DispersalVariation/blob/master/code/kernel_plotting.R)
  - [Test for temporal correlations](https://github.com/pinskylab/DispersalVariation/blob/master/code/dispersal_route_correlation.R)
  - [Test significance of differences with a null model](https://github.com/pinskylab/DispersalVariation/blob/master/code/null_dispersal.R)
  - [Test for changes in direction](https://github.com/pinskylab/DispersalVariation/blob/master/code/dispersal_directionality.R)
