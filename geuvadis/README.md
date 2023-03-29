### Overview

This directory contains scripts for replicating the simulation and GEUVADIS results in the isoQTL paper.

The scripts in the 'scripts/' directory are all the scripts required to run this experiment. The cython file 'cis_pass.pyx' will need to be compiled first, as described in the top-level README for this repo.

However, there are data requirements, 'data/genotypes_yri.vcf.gz' and 'data/qqnorm_all.bed.gz' that cannot be provided as they contain individual genotype and phenotype information from the GEUVADIS study. As can be seen in the submit scripts listed below, the simulations and real data analysis both require the the genotypes stored in 'data/genotypes_yri.vcf.gz', the simulation scripts require the phenotypes stored in 'data/qqnorm_all.bed.gz' while the real data analysis requires chromosome-specific bed files, and the real data analysis also requires a covariates file stored in 'data/covariates.txt'.

The GEUVADIS genotype data can be obtained online from https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/ and the expression data can be obtained online from https://www.internationalgenome.org/data-portal/data-collection/geuvadis. It was preprocessed as described in the Supplemental Material of the manuscript -- see parse_data.R.

### Simulations

The scripts 'submit_mvn_array.sh' and 'submit_noncis_array.sh' are the scripts used to generate the simulation results used in the paper and were run on a computing cluster. These scripts require that QTLtools is installed and in the system path.

The plots were generated using plot_geuvadis_sim_res.ipynb.

### Real data analysis

The script 'submit_geuvadis_real.sh' runs the analysis.

The plot showing the distribution of eGenes found by each method was generated with plot_metagene.ipynb. This requires as input results files for different methods, the genotype VCF file, and a GTF file, in our case 'Homo_sapiens.GRCh37.87.gtf.gz', which can be downloaded from http://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/. 

