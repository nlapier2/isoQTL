#!/bin/sh
#$ -S /bin/bash  # run in a bash shell
#$ -N job-gtex_nlapier2  # this is the (N)ame of your job
#$ -cwd  # this (c)hanges the (w)orking (d)irectory to the directory with this script
#$ -o stdout-gtex.out  # this is the file that standard (o)utput will be written to
#$ -l h_data=16G,highp,h_rt=167:59:00  
#$ -t 1-22:1  # run an array job, with job numbers ranging from 1 to 22 in increments of 1

source ~/.bash_profile  # load your account settings stored in your bash profile

i=${SGE_TASK_ID}
vcffile=data/thyroid/gtex_thyroid_genos.vcf.gz
bedfile=data/thyroid/beds/gtex_thyroid_chrom_${i}.bed.gz
covarfile=data/thyroid/gtex_thyroid.covariates.txt
outdir=results/thyroid/gtex_chrom_${i}
[ ! -d ${outdir} ] && mkdir ${outdir}  # make outdir if it doesn't exist

python scripts/ftest_allsnps_nominal.py --vcf ${vcffile} --pheno ${bedfile} --covariates ${covarfile} --window 1000000 --nominal 0.05 --output ${outdir}/nominal_gtex_res_chrom_${i}

#
