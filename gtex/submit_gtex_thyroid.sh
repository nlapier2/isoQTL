#!/bin/sh
#$ -S /bin/bash  # run in a bash shell
#$ -N job-gtex_nlapier2  # this is the (N)ame of your job
#$ -cwd  # this (c)hanges the (w)orking (d)irectory to the directory with this script
#$ -o stdout-gtex.out  # this is the file that standard (o)utput will be written to
#$ -l h_data=16G,highp,h_rt=23:59:00  
#$ -t 1-22:1  # run an array job, with job numbers ranging from 1 to 22 in increments of 1

source ~/.bash_profile  # load your account settings stored in your bash profile

i=${SGE_TASK_ID}
vcffile=data/thyroid/gtex_thyroid_genos.vcf.gz
bedfile=data/thyroid/beds/gtex_thyroid_chrom_${i}.bed.gz
covarfile=data/thyroid/gtex_thyroid.covariates.txt
outdir=results/thyroid/gtex_chrom_${i}
[ ! -d ${outdir} ] && mkdir ${outdir}  # make outdir if it doesn't exist

# run wilks, ftest, fisher, min, cauchy
#python scripts/isoqtl.py --vcf ${vcffile} --pheno ${bedfile} --covariates ${covarfile} --window 1000000 --permute 100 --output ${outdir}/gtex_res_chrom_${i}

# qtltools runs
QTLtools cis --vcf ${vcffile} --bed ${bedfile} --cov ${covarfile} --window 1000000 --out ${outdir}/gtex_res_qtltools_iso_chrom_${i}.tsv --permute 100
QTLtools cis --vcf ${vcffile} --bed ${bedfile} --cov ${covarfile} --window 1000000 --grp-best --out ${outdir}/gtex_res_qtltools_grpbest_chrom_${i}.tsv --permute 100
QTLtools cis --vcf ${vcffile} --bed ${bedfile} --cov ${covarfile} --window 1000000 --grp-pca1 --out ${outdir}/gtex_res_qtltools_grppca1_chrom_${i}.tsv --permute 100
QTLtools cis --vcf ${vcffile} --bed ${bedfile} --cov ${covarfile} --window 1000000 --grp-mean --out ${outdir}/gtex_res_qtltools_grpmean_chrom_${i}.tsv --permute 100

# fisher/min/cauchy on qtltools isoform level results
python scripts/combine_qtltools_iso_pvals.py --qtltools ${outdir}/gtex_res_qtltools_iso_chrom_${i}.tsv --bed ${bedfile} --method fisher --output ${outdir}/gtex_res_fisher_on_qtltools_iso_chrom_${i}.tsv
python scripts/combine_qtltools_iso_pvals.py --qtltools ${outdir}/gtex_res_qtltools_iso_chrom_${i}.tsv --bed ${bedfile} --method min --output ${outdir}/gtex_res_min_on_qtltools_iso_chrom_${i}.tsv
python scripts/combine_qtltools_iso_pvals.py --qtltools ${outdir}/gtex_res_qtltools_iso_chrom_${i}.tsv --bed ${bedfile} --method cauchy --output ${outdir}/gtex_res_cauchy_on_qtltools_iso_chrom_${i}.tsv

#
