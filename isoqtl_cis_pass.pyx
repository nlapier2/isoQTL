import os
import subprocess
import sys
import pandas as pd
from io import StringIO

from scipy.special import betaln as betaln_func
from scipy.special import ncfdtr
from scipy.stats import linregress, beta, chi2
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
    # center window around first_start and last_end
    window_start, window_end = first_start - window, last_end + window
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


cdef (double, double, double) precomp_wilks_bartlett(tx2expr, txlist, sample_size, num_iso):
    resid_intercept = []
    #Y_pred = []
    #print(len(txlist))
    for tx in txlist:
        intercept = np.mean(tx2expr[tx])
        resid = tx2expr[tx] - intercept
        resid_intercept.append(resid)
        #Y_pred.append(intercept)
        #print(tx2expr[tx]) ; print(intercept) ; print(resid)
    #resid_intercept = np.array(resid_intercept)
    resid_intercept = np.transpose(np.array(resid_intercept))
    #print(np.matmul(np.transpose(Y_pred), resid_intercept)) ; sys.exit()
    #print(resid_intercept.shape)
    #print(is_invertible(np.matmul(np.transpose(resid_intercept), resid_intercept)))
    cdef double det_resid_intercept = np.linalg.det(np.matmul(np.transpose(resid_intercept), resid_intercept) / sample_size)
    cdef double multiplier = -(sample_size - 2 - 0.5 * num_iso)
    cdef double chi2_dof = num_iso
    #print(det_resid_intercept, multiplier, chi2_dof) ; import sys ; sys.exit()
    return det_resid_intercept, multiplier, chi2_dof


cdef (double, double) wilks_bartlett(np.ndarray[DTYPE_t, ndim=2] resid_full, double det_resid_intercept, double multiplier, double chi2_dof, double sample_size):
    cdef double det_resid_full = np.linalg.det(np.matmul(np.transpose(resid_full), resid_full) / sample_size)
    #print(multiplier, det_resid_full, det_resid_intercept) ; sys.exit()
    cdef double chisq_stat = multiplier * np.log(det_resid_full / det_resid_intercept)
    #cdef double chisq_stat = multiplier * np.log(det_resid_intercept / det_resid_full)
    cdef double pval = chi2.sf(chisq_stat, chi2_dof)
    #print(multiplier, det_resid_full, det_resid_intercept)
    #print(chisq_stat, chi2_dof, pval)
    #print()
    return chisq_stat, pval


def is_invertible(a):
    return a.shape[0] == a.shape[1] and np.linalg.matrix_rank(a) == a.shape[0]


# implement test proposed by Hashimoto & Ohtani, Econometrics Letters (1990)
cdef (double, double) hashimoto_ohtani_t2_stat(np.ndarray[DTYPE_t, ndim=1] betas, np.ndarray[DTYPE_t, ndim=2] residuals, np.ndarray[DTYPE_t, ndim=1] genos, np.ndarray[DTYPE_t, ndim=2] z, double genos_dot_inv):
    cdef int i
    cdef double p = len(betas) * 1.0  # number of isoforms
    cdef int T = len(genos)  # sample size
    cdef double K = 1.0  # num regressors
    cdef int r = int((T-K) / K)
    cdef double dof = r - p + 1

    cdef np.ndarray[DTYPE_t, ndim=2] Sigma = np.matmul(residuals, np.transpose(residuals)) / (T - K)
    cdef np.ndarray[DTYPE_t, ndim=2] resid_star = np.matmul(residuals, z)
    cdef np.ndarray[DTYPE_t, ndim=2] S = np.matmul(resid_star, np.transpose(resid_star)) * genos_dot_inv / r
    
    sig_div_genos = Sigma / np.dot(genos, genos)
    if not is_invertible(S):
        if max(betas) < 10**-10:  # singular matrix due to zero association
            return 0.0, 1.0
        else:  # singular matrix due to extreme association
            return sys.float_info.max, sys.float_info.min
    if not is_invertible(sig_div_genos):
        if max(betas) < 10**-10:  # singular matrix due to zero association
            return 0.0, 1.0
        else:  # singular matrix due to extreme association
            return sys.float_info.max, sys.float_info.min

    cdef double hotelling_t2_stat = (np.matmul(np.matmul(np.transpose(betas), np.linalg.inv(S)), betas) / r) * ((r - p + 1) / p)
    ncp = np.matmul(np.matmul(betas, np.linalg.inv(sig_div_genos)), betas)
    cdef double pval = 1 - ncfdtr(p, dof, ncp, hotelling_t2_stat)
    # print(p, dof, ncp, hotelling_t2_stat, pval)
    return hotelling_t2_stat, pval


def find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, permuted_tx2expr, statistic):
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
        genos_var = np.var(genos)
        genos_mean = np.mean(genos)
        betas, residuals = get_betas_stderrs(snp2geno[snp], tx2expr, txlist, genos_var, genos_mean)
        genos = np.array(snp2geno[snp]).astype(np.double)
        if statistic == 'hashimoto':
            genos_dot_inv = 1 / np.dot(genos, genos)
            N = np.identity(len(genos)) - np.outer(genos * genos_dot_inv, np.transpose(genos))
            eigvals, eigvecs = np.linalg.eigh(N)
            z = np.delete(eigvecs, 0, axis=1)
            #print(betas) ; print(residuals) ; print(genos) ; print()
            stat, pval = hashimoto_ohtani_t2_stat(betas, residuals, genos, z, genos_dot_inv)
        else:
            det_resid_intercept, multiplier, chi2_dof = precomp_wilks_bartlett(tx2expr, txlist, len(genos), len(txlist))
            stat, pval = wilks_bartlett(np.transpose(residuals), det_resid_intercept, multiplier, chi2_dof, len(genos))
        if pval < peak_pval:
            peak_snp = snp
            peak_stat = stat
            peak_pval = pval
        if pval < nominal_thresh:
            for i in range(len(permuted_tx2expr)):
                betas, residuals = get_betas_stderrs(genos, permuted_tx2expr[i], txlist, genos_var, genos_mean)
                if statistic == 'hashimoto':
                    stat, pval = hashimoto_ohtani_t2_stat(betas, residuals, genos, z, genos_dot_inv)
                else:
                    stat, pval = wilks_bartlett(np.transpose(residuals), det_resid_intercept, multiplier, chi2_dof, len(genos))
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


def nominal_pass(vcf, tx2info, tx2expr, gene_to_tx, meta_lines, window, bcftools, n_perms, nominal_thresh, statistic):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    gene_results, perm_pvals = {}, {}
    for gene in gene_to_tx:
        txlist = gene_to_tx[gene]
        # get window around gene and subset the vcf for that window using bcftools
        window_start, window_end, chrom = get_gene_window(txlist, tx2info, window)
        window_str = str(chrom) + ':' + str(window_start) + '-' + str(window_end)
        bcf_proc = subprocess.Popen([bcftools, 'view', vcf, '-r', window_str],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bcf_out = StringIO(bcf_proc.communicate()[0].decode('utf-8'))
        gene_window_snps = pd.read_csv(bcf_out, sep='\t', header=meta_lines+4, index_col=2).T
        # read in SNPs in the window and clear out null SNPs and non-phenotyped individuals
        snp2info, snp2geno = gene_window_snps.iloc[:8], gene_window_snps.iloc[8:]
        snp2geno = preprocess_snp2geno(tx2expr, snp2geno)
        # find and store the peak SNP for this gene and its association statistic and p-value
        if n_perms > 0:
            permuted_tx2expr = permute_transcripts(tx2expr, txlist, n_perms)
            snp, fstat, pval, perm_peak_pvals = find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, permuted_tx2expr, statistic)
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
            snp, fstat, pval, perm_peak_pvals = find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, {}, statistic)
        gene_results[gene] = [snp, fstat, pval]
    fnull.close()
    return gene_results, perm_pvals
