# IsoQTL

This is the repository containing software implementation and information for replicating the results from the manuscript "Accounting for Isoform Expression in eQTL Mapping Substantially Increases Power". 

For more information, or if you use the software, please cite our paper:

LaPierre N, Pimentel H. Accounting for Isoform Expression in eQTL Mapping Substantially Increases Power. bioRxiv. 2023:2023-06.



### Installation and Overview

Most methods discussed in the paper (apart from QTLtools) are implemented in the 'isoqtl.py' and 'cis_pass.pyx' scripts in the 'scripts/' directory. Compilation of the cython file is required but is not very complicated and should take less than a minute. Instructions:

```
git clone https://github.com/nlapier2/isoQTL.git
cd isoQTL/scripts/
python setup.py build_ext --inplace
```

We also require the following dependencies: bcftools, statsmodels python package

The QTLtools results in the paper were obtained by running the QTLtools software package. See: https://qtltools.github.io/qtltools/

The results running p-value aggregation methods on QTLtools results were obtained using combine_qtltools_iso_pvals.py.



### Replication folders

The 'geuvadis' folder has scripts and basic instructions for replicating the simulation and GEUVADIS real data results. The 'gtex' folder has scripts and basic instructions for replicating the GTEx results. Please be aware that these require access to GEUVADIS and GTEx data, as explained in those respective directories.

The code used to plot the main text figures apart from Figure 1 is also provided in the 'geuvadis' and 'gtex' directories.



### Example run

```
python scripts/isoqtl.py --vcf dosages.vcf.gz --pheno iso_exp.bed.gz --output results_isoqtl.tsv
```

If multiple genes are being tested for association, we highly recommend running an FDR control method on these results. We recommend the Q value, which is available in the R "qvalue" package. We have included a helper script to run this. For example:

```
Rscript scripts/perform_qvalue_fdr.R results_isoqtl.tsv 0.1 IsoQTL fdr_results_isoqtl.tsv
```

The arguments are the input file, the FDR control level, the method that generated the results (use 'IsoQTL' if running isoqtl.py', and the output file.

By default, the script will run all methods it is capable of running. The '--methods' flag can be used to run only a single method, e.g. for F-test only, append '--methods ftest'.

For more details on the input file formats, output file formats, and arguments for these scripts, see below.



### Input file formats

IsoQTL takes as input a VCF file for the genotypes and a BED file for the phenotypes (isoform expression levels).

The input VCF is expected to be in dosage ("DS") format. If you have dosage format in addition to genotypes, you can use bcftools query to extract the dosage field. See: https://samtools.github.io/bcftools/bcftools-man.html#query

If you have genotypes but not dosage, you can convert to dosage format using the bcftools dosage plugin. This requires a little workaround to make a properly-formatted VCF file. For example:

```
bcftools view genos.vcf.gz --header-only > dosages.vcf
bcftools +dosage genos.vcf.gz | tail -n +2 >> dosages.vcf
```

We also require that VCFs are compressed and indexed with bgzip and tabix, so that we can extract information from them quickly. This can be done as follows:

```
bgzip -c dosages.vcf > dosages.vcf.gz
tabix -p vcf dosages.vcf.gz
```

The BED file should be a tab-delimited file with each row corresponding to an isoform and with six header columns (chromosome, start position, end position, isoform name, gene name, and strand) plus one additional column for each sample containing the isoform's expression level for that sample. The format is identical to that specified by QTLtools, so we refer to their documentation for more information: https://qtltools.github.io/qtltools/


The rest of the options are as follows (viewable in terminal by running "python isoqtl.py -h"):

```
usage: isoqtl.py [-h] --vcf VCF --pheno PHENO [--bcftools BCFTOOLS]
                 [--covariates COVARIATES]
                 [--methods {wilks,fisher,min,cauchy,ftest} [{wilks,fisher,min,cauchy,ftest} ...]]
                 [--nominal NOMINAL] [--output OUTPUT] [--permute PERMUTE]
                 [--window WINDOW]

isoQTL main script.

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF             VCF or BCF file with genotypes. Required.
  --pheno PHENO         BED file or tsv file with transcript expression
                        levels.
  --bcftools BCFTOOLS   Path to bcftools executable ("bcftools" by default).
  --covariates COVARIATES
                        tsv file listing covariates to adjust phenotypes for.
  --methods {wilks,fisher,min,cauchy,ftest} [{wilks,fisher,min,cauchy,ftest} ...]
                        Specify method to obtain gene-level p-value.
  --nominal NOMINAL     Print genes with a nominal assocation p-value below
                        this cutoff.
  --output OUTPUT       Base name for results.
  --permute PERMUTE     Number of permutations to do for a permutation pass.
                        Set to 0 to do nominal pass only.
  --window WINDOW       Size of window in bp around start position of
                        phenotypes.
```



### Output format

IsoQTL's output format should look like this:

```
#Gene   SNP     Statistic       Nominal P-value Dir. Perm. P-value      Beta Perm. P-value
ENSG00000123       rs987      2.59822 0.0110439       0.19802 0.176881
ENSG00000456      rs654      35.6631 8.82281e-08     0.00990099      2.47607e-05
```

The SNP field contains the top associated SNP with the Gene in the first column. The Statistic and Nominal P-value are computed using the hypothesis test described in our paper. The Direct Permutation P-value is computed by a standard permutation test, while the Beta Permutation P-value is used to approximate highly significant permutation p-values, as described in our paper. We highly recommend using this last column as the p-value input to FDR control.

The output from our FDR control helper script should look like:

```
ENSG00000732 rs598 49.0966 2.18176e-11 0.00990099 3.01129e-08 3.01129e-07
ENSG00000188 rs745 -6.67684 2.37684e-09 0.00990099 5.81673e-06 4.23034909090909e-05
```

Essentially, it is the same format, except without the header row, and with a new last column, which contains the Q value. All genes in this output are eGenes, as they have a Q value lower than the specified threshold. Again, for reporting p-values we recommend the Beta Permutation P-value, which is the second-to-last column here.

