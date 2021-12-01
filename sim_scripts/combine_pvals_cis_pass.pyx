import os
import subprocess
import sys
import pandas as pd
from io import StringIO

from scipy.special import betaln as betaln_func
from scipy.special import ncfdtr
from scipy.stats import linregress, beta, chi2, norm, combine_pvalues
from scipy.optimize import minimize

import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t


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


cdef (double, double) get_beta(x, double var_x, double mean_x, y, double mean_y):
    cdef int n = len(x)
    cdef double cov_xy = 0.0
    cdef int i
    for i in range(n):
        cov_xy += (x[i] - mean_x)*(y[i] - mean_y)
    cov_xy /= n
    cdef double beta = cov_xy / var_x
    cdef double intercept = mean_y - beta * mean_x
    return beta, intercept


def get_betas_stderrs(snp_genos, tx2expr, txlist, genos_var, genos_mean):
    betas, residuals = [], []
    for tx in txlist:
        #regression_results = linregress(snp_genos, tx2expr[tx])
        #beta, stderr = regression_results[0], regression_results[4]
        tx_arr = np.array(tx2expr[tx])
        beta, intercept = get_beta(np.array(snp_genos), genos_var, genos_mean, tx_arr, np.mean(tx_arr))
        resid = tx2expr[tx] - (snp_genos * beta) - intercept
        betas.append(beta)
        residuals.append(resid)
    return np.array(betas), np.array(residuals)


def cauchy_acat(pval_vec):
    weights = [1.0 / len(pval_vec) for i in pval_vec]  # equal weights for now
    stat = sum([weights[i] * np.tan((0.5 - pval_vec[i]) * np.pi) for i in range(len(pval_vec))])
    pval = 0.5 - (np.arctan(stat / sum(weights))) / np.pi
    return stat, pval


def find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, permuted_tx2expr, n_perms, combine_method):
    cdef double peak_stat = 0.0
    cdef double peak_pval = 1.0
    cdef double pval = 0.0
    cdef double stat = 0.0
    cdef np.ndarray[DTYPE_t, ndim=1] betas
    cdef np.ndarray[DTYPE_t, ndim=1] stderrs
    cdef np.ndarray[DTYPE_t, ndim=2] residuals
    cdef np.ndarray[DTYPE_t, ndim=1] genos
    cdef np.ndarray[DTYPE_t, ndim=1] perm_peak_pvals = np.ones(len(permuted_tx2expr)) * 0.99999999
    cdef int i
    peak_snp = 'none'
    # peak_snp, peak_stat, peak_pval = 'rsid', 0.0, 1.0
    # perm_peak_pvals = [1.0 for i in range(len(permuted_tx2expr))]
    for snp in snp2geno:
        genos = np.array(snp2geno[snp]).astype(np.double)
        #print(snp)
        all_pvals = [linregress(tx2expr[tx], genos)[3] for tx in txlist]
        if combine_method == 'cauchy':
            stat, pval = cauchy_acat(all_pvals)
        elif combine_method == 'min':
            stat, pval = 0.0, min(all_pvals)
        else:
            stat, pval = combine_pvalues(all_pvals, method=combine_method, weights=None)
        if pval < peak_pval:
            peak_pval = pval
            peak_stat = stat
            peak_snp = snp
        if pval < nominal_thresh:
            for i in range(len(permuted_tx2expr)):
                all_pvals = [linregress(permuted_tx2expr[i][tx], genos)[3] for tx in txlist]
                if combine_method == 'cauchy':
                    stat, pval = cauchy_acat(all_pvals)
                elif combine_method == 'min':
                    stat, pval = 0.0, min(all_pvals)
                else:
                    stat, pval = combine_pvalues(all_pvals, method=combine_method, weights=None)
                if pval < perm_peak_pvals[i]:
                    perm_peak_pvals[i] = pval
    return peak_snp, peak_stat, peak_pval, perm_peak_pvals


cdef (double, double) set_initial_beta_params(pval_vec):
    pval_mean = np.mean(pval_vec)
    pval_var = np.var(pval_vec)
    beta_shape1 = pval_mean * (pval_mean * (1.0 - pval_mean) / pval_var - 1.0)
    beta_shape2 = beta_shape1 * (1.0 / pval_mean - 1.0)
    return beta_shape1, beta_shape2


# def beta_likelihood_function(beta_shape1, beta_shape2, pval_vec):
cdef double beta_likelihood_function(beta_params, other_params):
    beta_shape1, beta_shape2 = beta_params
    likelihood_sum1, likelihood_sum2, n_var = other_params
    beta_func_res = betaln_func(beta_shape1, beta_shape2)
    return -1.0 * ((beta_shape1 - 1) * likelihood_sum1 + (beta_shape2 - 1) * likelihood_sum2 - n_var * beta_func_res)


cdef double compute_beta_perm_pvals(nominal_pval, perm_peak_pvals):
    beta_shape1, beta_shape2 = set_initial_beta_params(perm_peak_pvals)
    beta_params = np.array([beta_shape1, beta_shape2])
    likelihood_sum1 = sum([np.log(i) for i in perm_peak_pvals])
    likelihood_sum2 = sum([np.log(1 - i) for i in perm_peak_pvals])
    n_var = len(perm_peak_pvals)
    other_params = [likelihood_sum1, likelihood_sum2, n_var]
    res = minimize(beta_likelihood_function, x0=beta_params, args=other_params, method='nelder-mead')
    beta_shape1, beta_shape2 = res.x
    adjusted_pval = 1 - beta.sf(nominal_pval, beta_shape1, beta_shape2)
    return adjusted_pval


def permute_transcripts(tx2expr, txlist, n_perms):
    tx_perms = []
    for i in range(n_perms):
        permuted_dict = {}
        for tx in txlist:
            permuted_dict[tx] = np.copy(tx2expr[tx])
            np.random.shuffle(permuted_dict[tx])
        tx_perms.append(permuted_dict)
    return np.transpose(np.array(tx_perms))


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


def nominal_pass(vcf, tx2info, tx2expr, gene_to_tx, meta_lines, window, bcftools, n_perms, nominal_thresh, combine_method):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    gene_results, perm_pvals = {}, {}
    for gene in gene_to_tx:
        txlist = gene_to_tx[gene]
        # get window around gene and subset the vcf for that window using bcftools
        snp2geno = get_cis_snps(vcf, bcftools, txlist, tx2info, tx2expr, meta_lines, window)
        
        # find and store the peak SNP for this gene and its association statistic and p-value
        if n_perms > 0:
            permuted_tx2expr = permute_transcripts(tx2expr, txlist, n_perms)
            snp, fstat, pval, perm_peak_pvals = find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, permuted_tx2expr, n_perms, combine_method)
            if pval < nominal_thresh:
                dir_perm_pval = (1.0 + sum([pval >= perm_i for perm_i in perm_peak_pvals])) / float(n_perms + 1)
                #if np.var(perm_peak_pvals) == 0.0:
                if all(perm_peak_pvals > 0.999) and pval < 0.999:
                    beta_perm_pval = sys.float_info.min # dir_perm_pval
                elif dir_perm_pval == 1.0:
                    beta_perm_pval = 1.0
                else:
                    beta_perm_pval = compute_beta_perm_pvals(pval, perm_peak_pvals)

                perm_pvals[gene] = [dir_perm_pval, beta_perm_pval]
        else:
            snp, fstat, pval, perm_peak_pvals = find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, {}, n_perms, combine_method)
        gene_results[gene] = [snp, fstat, pval]
    fnull.close()
    return gene_results, perm_pvals
