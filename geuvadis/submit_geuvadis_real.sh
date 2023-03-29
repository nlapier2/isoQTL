#!/bin/sh
#$ -S /bin/bash  # run in a bash shell
#$ -N job-geuvadis-isoqtl-array_nlapier2  # this is the (N)ame of your job
#$ -cwd  # this (c)hanges the (w)orking (d)irectory to the directory with this script
#$ -o stdout-geuvadis-isoqtl-array.out  # this is the file that standard (o)utput will be written to
#$ -l h_data=16G,highp,h_rt=72:00:00  
#$ -t 1-22:1  # run an array job, with job numbers ranging from 1 to 22 in increments of 1

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile
R_LIBS_USER=/u/home/n/nlapier2/R/x86_64-pc-linux-gnu-library/4.0/  # setting R library path

i=${SGE_TASK_ID}
bedfile=data/beds/iso_level_norm/qqnorm_${SGE_TASK_ID}.bed.gz
genebedfile=data/beds/gene_level_norm/qqnorm_${SGE_TASK_ID}.bed.gz
covarfile=data/covariates.txt
vcffile=data/genotypes_yri.vcf.gz

outdir=2results/geuvadis_chr${SGE_TASK_ID}
[ ! -d ${outdir} ] && mkdir ${outdir}  # make outdir if it doesn't exist
touch ${outdir}/start_timestamp.txt

python scripts/isoqtl.py --vcf ${vcffile} --pheno ${bedfile} --covariates ${covarfile} --nominal 1.0 --window 50000 --output ${outdir}/results_chr_${i} --permute 100 --methods wilks fisher min cauchy ftest 

QTLtools cis --vcf ${vcffile} --bed ${genebedfile} --cov ${covarfile} --window 50000 --nominal 1.0 --out ${outdir}/results_qtltools_gene_chr_${i}.tsv --permute 100
QTLtools cis --vcf ${vcffile} --bed ${bedfile} --cov ${covarfile} --window 50000 --nominal 1.0 --out ${outdir}/results_qtltools_iso_chr_${i}.tsv --permute 100
QTLtools cis --vcf ${vcffile} --bed ${bedfile} --cov ${covarfile} --window 50000 --nominal 1.0 --grp-best --out ${outdir}/results_qtltools_iso_grpbest_chr_${i}.tsv --permute 100
QTLtools cis --vcf ${vcffile} --bed ${bedfile} --cov ${covarfile} --window 50000 --nominal 1.0 --grp-pca1 --out ${outdir}/results_qtltools_iso_grppca1_chr_${i}.tsv --permute 100
QTLtools cis --vcf ${vcffile} --bed ${bedfile} --cov ${covarfile} --window 50000 --nominal 1.0 --grp-mean --out ${outdir}/results_qtltools_iso_grpmean_chr_${i}.tsv --permute 100

python scripts/combine_qtltools_iso_pvals.py --qtltools ${outdir}/results_qtltools_iso_chr_${i}.tsv --tx2gene data/transcripts_to_genes.txt --method fisher --output ${outdir}/results_fisher_on_qtltools_iso_chr_${i}.tsv
python scripts/combine_qtltools_iso_pvals.py --qtltools ${outdir}/results_qtltools_iso_chr_${i}.tsv --tx2gene data/transcripts_to_genes.txt --method min --output ${outdir}/results_min_on_qtltools_iso_chr_${i}.tsv
python scripts/combine_qtltools_iso_pvals.py --qtltools ${outdir}/results_qtltools_iso_chr_${i}.tsv --tx2gene data/transcripts_to_genes.txt --method cauchy --output ${outdir}/results_cauchy_on_qtltools_iso_chr_${i}.tsv

Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_chr_${i}.wilks.tsv 0.1 IsoQTL ${outdir}/fdr_results_wilks_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_chr_${i}.ftest.tsv 0.1 IsoQTL ${outdir}/fdr_results_ftest_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_chr_${i}.fisher.tsv 0.1 fisher_perm ${outdir}/fdr_results_fisher_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_chr_${i}.min.tsv 0.1 fisher_perm ${outdir}/fdr_results_min_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_chr_${i}.cauchy.tsv 0.1 fisher_perm ${outdir}/fdr_results_cauchy_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_qtltools_gene_chr_${i}.tsv 0.1 QTLtools ${outdir}/fdr_results_qtltools_gene_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_qtltools_iso_grpbest_chr_${i}.tsv 0.1 QTLtools_grp ${outdir}/fdr_results_qtltools_grpbest_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_qtltools_iso_grppca1_chr_${i}.tsv 0.1 QTLtools_grp ${outdir}/fdr_results_qtltools_grppca1_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_qtltools_iso_grpmean_chr_${i}.tsv 0.1 QTLtools_grp_mean ${outdir}/fdr_results_qtltools_grpmean_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_fisher_on_qtltools_iso_chr_${i}.tsv 0.1 combined_qtltools ${outdir}/fdr_results_fisher_on_qtltools_iso_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_min_on_qtltools_iso_chr_${i}.tsv 0.1 combined_qtltools ${outdir}/fdr_results_min_on_qtltools_iso_chr_${i}.tsv
Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_cauchy_on_qtltools_iso_chr_${i}.tsv 0.1 combined_qtltools ${outdir}/fdr_results_cauchy_on_qtltools_iso_chr_${i}.tsv
#
