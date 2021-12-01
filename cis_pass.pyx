import os
import subprocess
import sys
import pandas as pd
from io import StringIO

from scipy.special import betaln as betaln_func
from scipy.special import ncfdtr
from scipy.stats import linregress, beta, chi2, norm
from scipy.optimize import minimize

import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t


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


def steiger_test(tx2expr, txlist, snp2geno, snp, permuted_tx2expr, n_perms, steiger_thresh):
    # get corrs of each tx with the snp
    tx_mat = [np.array(tx2expr[tx]) for tx in txlist]
    snp_corrs = [np.corrcoef(tx_mat[i], snp2geno[snp])[0][1] for i in range(len(tx_mat))]
    tx_corr_mat = np.corrcoef(tx_mat)
    # compute phi, c, and z_1* (eq. 3, 10, 12 in steiger)
    for k in range(len(tx_mat)):
        for h in range(len(tx_mat)):
            if h <= k:  # matrix is symmetrical, so only compute half
                continue
            r_kh = tx_corr_mat[k][h]
            r_jk, r_jh = snp_corrs[k], snp_corrs[h]
            # compute fisher-transformed z-statistics (eq. 8 in steiger)
            z_jk = -0.5 * np.log((1 + r_jk) / (1 - r_jk))
            z_jh = -0.5 * np.log((1 + r_jh) / (1 - r_jh))
            # compute phi (eq. 3 in steiger)
            phi = r_kh * (1 - r_jk**2 - r_jh**2) - 0.5 * (r_jk * r_jh) * (1 - r_jk**2 - r_jh**2 - r_kh**2)
            # compute s and z_1* (eq 10 and 12 in steiger)
            s_jkh = phi / (1 - r_jk**2) * (1 - r_jh**2)
            z1 = np.sqrt(len(tx_mat[0]) - 3) * (z_jk - z_jh) * -np.sqrt((2 - 2 * s_jkh))
            # compute pval and if significant than proceed to multivariate test
            if z1 >= 0:
                pval = 2 * norm.sf(z1)
            else:
                pval = 2 * norm.cdf(z1)
            #print(z1, pval)
            if pval < steiger_thresh:
                return tx2expr, txlist, permuted_tx2expr
    # if no sig diff corrs, combine transcripts and proceed to univariate test
    this_txlist = [txlist[0]]
    this_tx2expr = {this_txlist[0]: sum(tx_mat)}
    if n_perms > 0:
        this_perm = permute_transcripts(tx2expr, txlist, n_perms)
    else:
        this_perm = {}
    return this_tx2expr, this_txlist, this_perm


def find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, permuted_tx2expr, n_perms, steiger_thresh):
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
        if len(txlist) > 1:
            this_tx2expr, this_txlist, this_perm = steiger_test(tx2expr, txlist, snp2geno, snp, permuted_tx2expr, n_perms, steiger_thresh)
        else:
            this_tx2expr, this_txlist, this_perm = tx2expr, txlist, permuted_tx2expr
        #print('')
        if len(this_txlist) == 1:
            beta, inter, r, pval, stderr = linregress(this_tx2expr[this_txlist[0]], genos)
            stat = beta / stderr
        else:
            genos_var = np.var(genos)
            genos_mean = np.mean(genos)
            betas, residuals = get_betas_stderrs(snp2geno[snp], tx2expr, txlist, genos_var, genos_mean)
            genos = np.array(snp2geno[snp]).astype(np.double)
            det_resid_intercept, multiplier, chi2_dof = precomp_wilks_bartlett(tx2expr, txlist, len(genos), len(txlist))
            stat, pval = wilks_bartlett(np.transpose(residuals), det_resid_intercept, multiplier, chi2_dof, len(genos))
        if pval < peak_pval:
            peak_snp = snp
            peak_stat = stat
            peak_pval = pval
        if pval < nominal_thresh:
            for i in range(len(permuted_tx2expr)):
                if len(this_txlist) == 1:
                    beta, inter, r, pval, stderr = linregress(this_perm[i][this_txlist[0]], genos)
                else:
                    betas, residuals = get_betas_stderrs(genos, permuted_tx2expr[i], txlist, genos_var, genos_mean)
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


def nominal_pass(vcf, tx2info, tx2expr, gene_to_tx, meta_lines, window, bcftools, n_perms, nominal_thresh, steiger_thresh):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    gene_results, perm_pvals = {}, {}
    for gene in gene_to_tx:
        txlist = gene_to_tx[gene]
        # get window around gene and subset the vcf for that window using bcftools
        snp2geno = get_cis_snps(vcf, bcftools, txlist, tx2info, tx2expr, meta_lines, window)
        
        # find and store the peak SNP for this gene and its association statistic and p-value
        if n_perms > 0:
            permuted_tx2expr = permute_transcripts(tx2expr, txlist, n_perms)
            snp, fstat, pval, perm_peak_pvals = find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, permuted_tx2expr, n_perms, steiger_thresh)
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
            snp, fstat, pval, perm_peak_pvals = find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, {}, n_perms, steiger_thresh)
        gene_results[gene] = [snp, fstat, pval]
    fnull.close()
    return gene_results, perm_pvals
