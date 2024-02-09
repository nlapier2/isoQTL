#!/bin/sh
#$ -S /bin/bash  # run in a bash shell
#$ -N job-sim-isoqtl-array_nlapier2  # this is the (N)ame of your job
#$ -cwd  # this (c)hanges the (w)orking (d)irectory to the directory with this script
#$ -o stdout-sim-isoqtl-array.out  # this is the file that standard (o)utput will be written to
#$ -l h_data=16G,highp,h_rt=72:00:00  
#$ -t 1-25:1  # run an array job, with job numbers ranging from 1 to 25 in increments of 1

. /u/local/Modules/default/init/modules.sh  # allows you to load modules
source ~/.bash_profile  # load your account settings stored in your bash profile
R_LIBS_USER=/u/home/n/nlapier2/R/x86_64-pc-linux-gnu-library/4.0/  # setting R library path

h2cis=(0.01 0.05 0.1 0.15 0.2)
negpct=(0.1 0.3 0.5 0.7 0.9)

index=$((${SGE_TASK_ID}-1))
h2c_index=$((${index}%${#h2cis[@]}))
neg_index=$((${index}/${#h2cis[@]}))
this_h2c=${h2cis[$h2c_index]}
this_neg=${negpct[$neg_index]}
outdir=results/mvn_genomewide_numiso_4_ftest/h2cis_${this_h2c}_negpct_${this_neg}
[ ! -d ${outdir} ] && mkdir ${outdir}  # make outdir if it doesn't exist

for i in {1..1} #10}
do
    python scripts/mvn_simulation.py --vcf data/genotypes_yri.vcf.gz --pheno data/qqnorm_all.bed.gz --dropout_rate 0.5 --dropout_iso 0.0 --h2cis ${this_h2c} --h2noncis 0.0 --neg_pct ${this_neg} --max_corr 0.49 --min_corr 0.09 --max_corr_env 0.49 --min_corr_env 0.09 --output ${outdir}/sim_phenos_${i} --window 50000 --num_iso 4 # --num_causal 3 --num_iso 2

    python scripts/isoqtl.py --vcf data/genotypes_yri.vcf.gz --pheno ${outdir}/sim_phenos_${i}.iso.bed.gz --nominal 0.5 --window 50000 --output ${outdir}/results_sim_${i} --permute 100 --methods wilks fisher min cauchy ftest 

    QTLtools cis --vcf data/genotypes_yri.vcf.gz --bed ${outdir}/sim_phenos_${i}.gene.bed.gz --window 50000 --out ${outdir}/results_qtltools_gene_perm_sim_${i}.tsv --permute 100
    QTLtools cis --vcf data/genotypes_yri.vcf.gz --bed ${outdir}/sim_phenos_${i}.iso.bed.gz --window 50000 --out ${outdir}/results_qtltools_iso_perm_sim_${i}.tsv --permute 100
    QTLtools cis --vcf data/genotypes_yri.vcf.gz --bed ${outdir}/sim_phenos_${i}.iso.bed.gz --window 50000 --grp-best --out ${outdir}/results_qtltools_iso_grpbest_perm_sim_${i}.tsv --permute 100
    QTLtools cis --vcf data/genotypes_yri.vcf.gz --bed ${outdir}/sim_phenos_${i}.iso.bed.gz --window 50000 --grp-pca1 --out ${outdir}/results_qtltools_iso_grppca1_perm_sim_${i}.tsv --permute 100
    QTLtools cis --vcf data/genotypes_yri.vcf.gz --bed ${outdir}/sim_phenos_${i}.iso.bed.gz --window 50000 --grp-mean --out ${outdir}/results_qtltools_iso_grpmean_perm_sim_${i}.tsv --permute 100

    python scripts/combine_qtltools_iso_pvals.py --qtltools ${outdir}/results_qtltools_iso_perm_sim_${i}.tsv --tx2gene data/transcripts_to_genes.txt --method fisher --output ${outdir}/results_fisher_on_qtltools_iso_sim_${i}.tsv
    python scripts/combine_qtltools_iso_pvals.py --qtltools ${outdir}/results_qtltools_iso_perm_sim_${i}.tsv --tx2gene data/transcripts_to_genes.txt --method min --output ${outdir}/results_min_on_qtltools_iso_sim_${i}.tsv
    python scripts/combine_qtltools_iso_pvals.py --qtltools ${outdir}/results_qtltools_iso_perm_sim_${i}.tsv --tx2gene data/transcripts_to_genes.txt --method cauchy --output ${outdir}/results_cauchy_on_qtltools_iso_sim_${i}.tsv

    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_sim_${i}.wilks.tsv 0.1 IsoQTL ${outdir}/fdr_results_wilks_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_sim_${i}.ftest.tsv 0.1 IsoQTL ${outdir}/fdr_results_ftest_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_sim_${i}.fisher.tsv 0.1 fisher_perm ${outdir}/fdr_results_fisher_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_sim_${i}.min.tsv 0.1 fisher_perm ${outdir}/fdr_results_min_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_sim_${i}.cauchy.tsv 0.1 fisher_perm ${outdir}/fdr_results_cauchy_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_qtltools_gene_perm_sim_${i}.tsv 0.1 QTLtools ${outdir}/fdr_results_qtltools_gene_perm_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_qtltools_iso_grpbest_perm_sim_${i}.tsv 0.1 QTLtools_grp ${outdir}/fdr_results_qtltools_grpbest_perm_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_qtltools_iso_grppca1_perm_sim_${i}.tsv 0.1 QTLtools_grp ${outdir}/fdr_results_qtltools_grppca1_perm_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_qtltools_iso_grpmean_perm_sim_${i}.tsv 0.1 QTLtools_grp_mean ${outdir}/fdr_results_qtltools_grpmean_perm_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_fisher_on_qtltools_iso_sim_${i}.tsv 0.1 combined_qtltools ${outdir}/fdr_results_fisher_on_qtltools_iso_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_min_on_qtltools_iso_sim_${i}.tsv 0.1 combined_qtltools ${outdir}/fdr_results_min_on_qtltools_iso_sim_${i}.tsv
    Rscript scripts/perform_qvalue_fdr.R ${outdir}/results_cauchy_on_qtltools_iso_sim_${i}.tsv 0.1 combined_qtltools ${outdir}/fdr_results_cauchy_on_qtltools_iso_sim_${i}.tsv

    python scripts/eval_sim_mvn.py --gene_info ${outdir}/sim_phenos_${i}.causal.txt --wilks ${outdir}/fdr_results_wilks_sim_${i}.tsv --ftest ${outdir}/fdr_results_ftest_sim_${i}.tsv --qtltools_gene ${outdir}/fdr_results_qtltools_gene_perm_sim_${i}.tsv --qtltools_grpbest ${outdir}/fdr_results_qtltools_grpbest_perm_sim_${i}.tsv --qtltools_grppca1 ${outdir}/fdr_results_qtltools_grppca1_perm_sim_${i}.tsv --qtltools_grpmean ${outdir}/fdr_results_qtltools_grpmean_perm_sim_${i}.tsv --fisher_perm ${outdir}/fdr_results_fisher_sim_${i}.tsv --cauchy_perm ${outdir}/fdr_results_cauchy_sim_${i}.tsv --min_perm ${outdir}/fdr_results_min_sim_${i}.tsv --fisher_qtltools_iso ${outdir}/fdr_results_fisher_on_qtltools_iso_sim_${i}.tsv --min_qtltools_iso ${outdir}/fdr_results_min_on_qtltools_iso_sim_${i}.tsv --cauchy_qtltools_iso ${outdir}/fdr_results_cauchy_on_qtltools_iso_sim_${i}.tsv --threshold 0.1 --use_perm > ${outdir}/eval_perm_${i}.txt
done

python scripts/tabulate_results.py ${outdir}/ > ${outdir}/tabulated_results.txt
#
