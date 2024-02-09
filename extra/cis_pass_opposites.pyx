import os
import subprocess
import sys
from io import StringIO
import pandas as pd
import statsmodels.formula.api as smf
from scipy.special import betaln as betaln_func
from scipy.special import ncfdtr
from scipy.stats import linregress, beta, chi2, norm, t, f, combine_pvalues
from scipy.optimize import minimize
import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t


cdef (double, double, double) get_beta_stderr(x, double var_x, double mean_x, y, double mean_y, double var_y, int dof):
    cdef int n = len(x)
    cdef double cov_xy = 0.0
    cdef int i
    for i in range(n):
        cov_xy += (x[i] - mean_x) * (y[i] - mean_y)
    cov_xy /= n
    cdef double r = cov_xy / np.sqrt(var_x * var_y)
    cdef double stderr = np.sqrt((1 - r**2) * var_y / var_x / dof)
    cdef double beta = cov_xy / var_x
    cdef double intercept = mean_y - beta * mean_x
    return beta, stderr, intercept


cdef double get_pval(double beta, double stderr, int dof):
    cdef double stat = beta / stderr
    if stat <= 0:
        return 2 * t.cdf(stat, dof)
    else:
        return 2 * t.sf(stat, dof)


def stacked_regress(snp_genos, tx2expr, txlist, genos_var, genos_mean, dof):
    # perform multiple ("stacked") linear regressions
    betas, stderrs, pvals, residuals = [], [], [], []
    for tx in txlist:
        #regression_results = linregress(snp_genos, tx2expr[tx])
        #beta, stderr = regression_results[0], regression_results[4]
        tx_arr = np.array(tx2expr[tx])
        tx_mean, tx_var = np.mean(tx_arr), np.var(tx_arr)
        beta, stderr, intercept = get_beta_stderr(np.array(snp_genos), genos_var, genos_mean, tx_arr, tx_mean, tx_var, dof)
        pval = get_pval(beta, stderr, dof)
        resid = tx2expr[tx] - (snp_genos * beta) - intercept
        betas.append(beta)
        stderrs.append(stderr)
        pvals.append(pval)
        residuals.append(resid)
    return np.array(betas), np.array(stderrs), np.array(pvals), np.array(residuals)


def get_regr_res(tx2expr, snp2geno, txlist, nominal_thresh, permuted_tx2expr, n_perms, methods):
    cdef int dof = len(snp2geno.iloc[:, 0]) - 2

    for snp in snp2geno:
        genos = np.array(snp2geno[snp]).astype(np.double)
        if len(txlist) == 1:
            betas, inter, r, all_pvals, stderrs = linregress(tx2expr[txlist[0]], genos)
        else:
            genos_var = np.var(genos)
            genos_mean = np.mean(genos)
            if genos_var == 0.0:
                continue
            betas, stderrs, all_pvals, residuals = stacked_regress(snp2geno[snp], tx2expr, txlist, genos_var, genos_mean, dof)
        regr_res = [betas, stderrs, all_pvals]
    return regr_res


def preprocess_snp2geno(tx2expr, snp2geno):
    # remove genotypes for non-phenotyped individuals and remove "nan" entries
    snp2geno = (snp2geno.T[list(tx2expr.index)]).T  # remove genotypes for non-phenotyped individuals
    new_snp2geno = pd.DataFrame()
    for i in snp2geno:
        if snp2geno[i].isnull().values.any():
            continue  # exclude SNPs with null genotypes
        new_snp2geno[i] = pd.to_numeric(snp2geno[i])
    return new_snp2geno


def get_gene_window(txlist, tx2info, window):
    chrom = ''  # chromosome this gene is on
    first_start, last_end = 10**20, 0  # first start post & last end pos among all transcripts
    for tx in txlist:
        chrom, start, end, gid, strand = tx2info[tx]
        if start < first_start:
            first_start = start
        if end > last_end:
            last_end = end
    # center window around TSS
    # window_start, window_end = first_start - window, last_end + window
    window_start, window_end = first_start - window, first_start + window
    return int(window_start), int(window_end), chrom


def get_cis_snps(vcf, bcftools, txlist, tx2info, tx2expr, meta_lines, window):
    window_start, window_end, chrom = get_gene_window(txlist, tx2info, window)
    window_str = str(chrom) + ':' + str(window_start) + '-' + str(window_end)
    bcf_proc = subprocess.Popen([bcftools, 'view', vcf, '-r', window_str], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bcf_out = StringIO(bcf_proc.communicate()[0].decode('utf-8'))
    gene_window_snps = pd.read_csv(bcf_out, sep='\t', header=meta_lines, index_col=2).T
    # read in SNPs in the window and clear out null SNPs and non-phenotyped individuals
    snp2info, snp2geno = gene_window_snps.iloc[:8], gene_window_snps.iloc[8:]
    snp2geno = preprocess_snp2geno(tx2expr, snp2geno)
    return snp2geno


def nominal_pass(prev, vcf, tx2info, tx2expr, gene_to_tx, meta_lines, window, bcftools, n_perms, nominal_thresh, methods, qtltools):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    gene_results = {}
    for gene in gene_to_tx:
        txlist = gene_to_tx[gene]
        # get window around gene and subset the vcf for that window using bcftools
        snp2geno = get_cis_snps(vcf, bcftools, txlist, tx2info, tx2expr, meta_lines, window)
        if qtltools:
            gene = gene.split('.')[0]
        if gene not in prev:  # not a gene whose nominal results we want to print
            continue
        target_snp = prev[gene].item()
        gene_snp_pairname = gene + '_' + target_snp
        snp2geno = pd.DataFrame({target_snp: snp2geno[target_snp]})

        regr_res = get_regr_res(tx2expr, snp2geno, txlist, nominal_thresh, {}, n_perms, methods)
        gene_results[gene] = regr_res
    fnull.close()
    return gene_results
