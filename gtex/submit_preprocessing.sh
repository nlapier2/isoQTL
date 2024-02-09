#!/bin/sh
#$ -S /bin/bash  # run in a bash shell
#$ -N job-preprocess_nlapier2  # this is the (N)ame of your job
#$ -cwd  # this (c)hanges the (w)orking (d)irectory to the directory with this script
#$ -o stdout-preprocess.out  # this is the file that standard (o)utput will be written to
#$ -M nathanl2012@gmail.com  # e(M)ail address to which hoffman2 should send notifications about job status
#$ -m abe  # e(m)ail you at the address from -M when job is (a)borted, (b)egins, or (e)nds
#$ -l h_data=16G,highp,h_rt=4:00:00 

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile

# tissue-specific names
basename=data/thyroid/gtex_thyroid
basename_beds=data/thyroid/beds/gtex_thyroid
gene_gct=data/gct/gene_reads_2017-06-05_v8_thyroid.gct.gz

# fixed names
gct_iso_reads=data/gct/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz
gct_iso_tpm=data/gct/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
covar=data/covariates_from_web/Thyroid.v8.covariates.txt
eur=data/population/EUR.list
bfile=data/bfile/GTEx_MAF_0.01
gtf=data/gtf/Homo_sapiens.GRCh38.88.chr.gtf.gz

# preprocessing scripts
python scripts/subset_files.py --gene_gct ${gene_gct} --iso_gct_tpm ${gct_iso_tpm} --iso_gct_reads ${gct_iso_reads} --covariates ${covar} --eur ${eur} --output ${basename}
python scripts/create_vcf.py --bfile ${bfile} --ids_file ${basename}.ids.txt --samples_file ${basename}.samples.txt --output ${basename}
python scripts/create_beds.py --gct_tpm ${basename}.tpm.gct.gz --gct_reads ${basename}.reads.gct.gz --gtf ${gtf} --output ${basename_beds}
#
