# Overview

This directory contains the scripts used to produce the GTEx results in the isoQTL manuscript.

In addition to the scripts in the top-level 'scripts/' directory, this directory also needs the scripts in the 'scripts/' directory in this folder.

GTEx results for QTLtools and p-value aggregation methods run on QTLtools results were generated using 'submit_preprocessing.sh' to preprocess (filter and transform) the data as described in the manuscript and then 'submit_gtex_thyroid.sh' to generate the results.

For all other methods (covered by the 'isoqtl.py' and 'cis_pass.pyx' scripts), this was computationall infeasible due to some genes with many isoforms. Consequently, the methods were run on each gene individually, using the analogous two scripts, 'submit_preprocessing_gene_level.sh' and 'submit_per_gene_gtex.sh'.

These scripts require GTEx genotype and phenotype data to run, as well as covariates and population identifiers. GTEx genotypes must be obtained by request as described in https://gtexportal.org/home/protectedDataAccess. The GTEx expression data, covariates, and population information can be obtained from https://gtexportal.org/home/datasets.


# Enrichment plots

Code for generating the functional enrichment plot is located in the enrichment_plots/ subdirectory. 

First 'submit_ftest_nominal.sh' and 'submit_qtltools_nominal.sh' were used to generate all nominally-significant eQTLs. 

Then 'prepare_torus_input.py' and 'prepare_torus_input_isoqtl.py' were used to generate files suitable for input to torus. 

Torus was then run to generate functional enrichments, as described in the torus wiki (https://github.com/xqwen/torus/tree/master/examples/GEUVADIS). 

The script plot_enrichments.ipynb was used to generate the figure seen in the text.
 
