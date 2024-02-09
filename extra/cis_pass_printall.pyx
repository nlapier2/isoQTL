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


cdef (double, double, double) precomp_wilks_bartlett(tx2expr, txlist, sample_size, num_iso):
    resid_intercept = []
    for tx in txlist:
        intercept = np.mean(tx2expr[tx])
        resid = tx2expr[tx] - intercept
        resid_intercept.append(resid)
    if np.linalg.cond(resid_intercept) >= 1 / sys.float_info.epsilon:  # if singular matrix
        for ind in range(len(resid_intercept)):
            noise = np.random.normal(0, 0.1, len(resid_intercept[ind]))
            resid_intercept[ind] += noise
    resid_intercept = np.transpose(np.array(resid_intercept))
    cdef double det_resid_intercept = np.linalg.det(np.matmul(np.transpose(resid_intercept), resid_intercept) / sample_size)
    cdef double multiplier = -(sample_size - 2 - 0.5 * num_iso)
    cdef double chi2_dof = num_iso
    return det_resid_intercept, multiplier, chi2_dof


cdef (double, double) wilks_bartlett(np.ndarray[DTYPE_t, ndim=2] resid_full, double det_resid_intercept, double multiplier, double chi2_dof, double sample_size):
    cdef double det_resid_full = np.linalg.det(np.matmul(np.transpose(resid_full), resid_full) / sample_size)
    cdef double chisq_stat = multiplier * np.log(det_resid_full / det_resid_intercept)
    cdef double pval = chi2.sf(chisq_stat, chi2_dof)
    return chisq_stat, pval


cdef (double, double) cauchy_acat(pval_vec):
    weights = [1.0 / len(pval_vec) for i in pval_vec]  # equal weights for now
    stat = sum([weights[i] * np.tan((0.5 - pval_vec[i]) * np.pi) for i in range(len(pval_vec))])
    pval = 0.5 - (np.arctan(stat / sum(weights))) / np.pi
    return stat, pval


cdef (double, double) multi_regress_f_test(tx2expr, txlist, snp2geno, snp):
    # perform reverse (multiple) regression of snp on isoforms and get f statistic
    this_tx2expr = tx2expr
    this_tx2expr[snp] = snp2geno[snp]  # merge snp into tx2expr df so we can run statsmodels regression
    formula = snp + ' ~ ' + ' + '.join(['Q("' + i + '")' for i in txlist])  # formula for statsmodels regression
    reg_res = smf.ols(formula, data=this_tx2expr).fit()
    return reg_res.fvalue, reg_res.f_pvalue


def precomp_xtx_inv(tx2expr, txlist):
    mat = [tx2expr[tx] for tx in txlist]
    mat.append(np.ones(len(mat[0])))
    x = np.array(mat).T  # add intercept term
    # x = np.array([tx2expr[tx] for tx in txlist])
    xtx = np.matmul(x.T, x)
    if np.linalg.cond(xtx) >= 1 / sys.float_info.epsilon:  # if singular matrix
        for ind in range(len(xtx)):
            noise = np.random.normal(0, 0.1, len(xtx[ind]))
            xtx[ind] += noise
    xtx_inv = np.linalg.inv(xtx)
    return x, xtx_inv


cdef (double, double) multi_regress_f_test_given_xtx(y, x, xtx_inv):
    # compute OLS effect estimates for multiple regression
    dof1, dof2 = len(x[0]) - 1, len(y) - len(x[0])
    xty = np.matmul(x.T, y)
    betahat = np.matmul(xtx_inv, xty)
    # compute standard errors of OLS estimates
    yhat = np.matmul(x, betahat)
    resid = y - yhat
    resid_var = sum([i**2 for i in resid]) / dof2
    mat = resid_var * xtx_inv
    betahat_stderrs = np.array([mat[i][i] ** 2 for i in range(len(mat))])
    # compute F statistic
    rss_full = sum([i**2 for i in resid])
    resid_intercept = y - np.mean(y)
    rss_intercept = sum([i**2 for i in resid_intercept])
    numerator = (rss_intercept - rss_full) / dof1
    denominator = rss_full / dof2
    stat = numerator / denominator
    pval = f.sf(stat, dof1, dof2)
    return stat, pval


def find_peak_snp(tx2expr, snp2geno, txlist, nominal_thresh, permuted_tx2expr, n_perms, methods):
    # initialize returned variables indicating peak associations
    peak_snps, peak_stats, peak_pvals, perm_peak_pvals = {}, {}, {}, {}
    for m in methods:
        peak_snps[m] = 'none'
        peak_stats[m] = 0.0
        peak_pvals[m] = 1.0
        perm_peak_pvals[m] = np.ones(len(permuted_tx2expr)) * 0.99999999
    if len(snp2geno) == 0:
        return peak_snps, peak_stats, peak_pvals, perm_peak_pvals

    # initialize variables for computations
    cdef double pval = 0.0
    cdef double stat = 0.0
    cdef np.ndarray[DTYPE_t, ndim=1] betas
    cdef np.ndarray[DTYPE_t, ndim=1] stderrs
    cdef np.ndarray[DTYPE_t, ndim=2] residuals
    cdef np.ndarray[DTYPE_t, ndim=1] all_pvals
    cdef np.ndarray[DTYPE_t, ndim=1] genos
    cdef int i
    cdef int n_genos = len(snp2geno.iloc[:, 0])
    cdef int dof = n_genos - 2
    det_resid_intercept, multiplier, chi2_dof = precomp_wilks_bartlett(tx2expr, txlist, n_genos, len(txlist))
    if 'ftest' in methods:
        x, xtx_inv = precomp_xtx_inv(tx2expr, txlist)
        permuted_xtx_inv = [precomp_xtx_inv(permuted_tx2expr[i], txlist) for i in range(len(permuted_tx2expr))]
    else:
        xtx_inv, permuted_xtx_inv = 0, 0

    # go through each SNP to find strongest p-value for each method
    for snp in snp2geno:
        genos = np.array(snp2geno[snp]).astype(np.double)

        # compute nominal pass
        if len(txlist) == 1:
            beta, inter, r, pval, stderr = linregress(tx2expr[txlist[0]], genos)
            stat = beta / stderr
        else:
            genos_var = np.var(genos)
            genos_mean = np.mean(genos)
            if genos_var == 0.0:
                continue
            betas, stderrs, all_pvals, residuals = stacked_regress(snp2geno[snp], tx2expr, txlist, genos_var, genos_mean, dof)
        for m in methods:
            if len(txlist) != 1:
                if m == 'wilks':
                    # det_resid_intercept, multiplier, chi2_dof = precomp_wilks_bartlett(tx2expr, txlist, len(genos), len(txlist))
                    stat, pval = wilks_bartlett(np.transpose(residuals), det_resid_intercept, multiplier, chi2_dof, len(genos))
                elif m == 'ftest':
                    #stat, pval = multi_regress_f_test(tx2expr, txlist, snp2geno, snp)
                    stat, pval = multi_regress_f_test_given_xtx(snp2geno[snp], x, xtx_inv)
                elif m == 'min':
                    stat, pval = 0.0, min(all_pvals)
                elif m == 'cauchy':
                    stat, pval = cauchy_acat(all_pvals)
                else:  # m == 'fisher'
                    stat, pval = combine_pvalues(all_pvals, method=m, weights=None)
            if pval < peak_pvals[m]:
                peak_pvals[m] = pval
                peak_stats[m] = stat
                peak_snps[m] = snp

        # compute permutation pass
        if min([peak_pvals[key] for key in peak_pvals]) < nominal_thresh:
            for i in range(len(permuted_tx2expr)):
                if len(txlist) == 1:
                    beta, inter, r, pval, stderr = linregress(permuted_tx2expr[i][txlist[0]], genos)
                else:
                    betas, stderrs, all_pvals, residuals = stacked_regress(genos, permuted_tx2expr[i], txlist, genos_var, genos_mean, dof)
                for m in methods:
                    if len(txlist) != 1:
                        if m == 'wilks':
                            stat, pval = wilks_bartlett(np.transpose(residuals), det_resid_intercept, multiplier, chi2_dof, len(genos))
                        elif m == 'ftest':
                            #stat, pval = multi_regress_f_test(permuted_tx2expr[i], txlist, snp2geno, snp)
                            perm_x, perm_xtx_inv = permuted_xtx_inv[i]
                            stat, pval = multi_regress_f_test_given_xtx(snp2geno[snp], perm_x, perm_xtx_inv)
                        elif m == 'min':
                            stat, pval = 0.0, min(all_pvals)
                        elif m == 'cauchy':
                            stat, pval = cauchy_acat(all_pvals)
                        else:  # m == 'fisher'
                            stat, pval = combine_pvalues(all_pvals, method=m, weights=None)
                    if pval < perm_peak_pvals[m][i]:
                        perm_peak_pvals[m][i] = pval
    return peak_snps, peak_stats, peak_pvals, perm_peak_pvals


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


def nominal_pass(vcf, tx2info, tx2expr, gene_to_tx, meta_lines, window, bcftools, n_perms, nominal_thresh, methods):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    gene_results, perm_pvals = {}, {}
    for m in methods:
        gene_results[m] = {}
        perm_pvals[m] = {}

    for gene in gene_to_tx:
        txlist = gene_to_tx[gene]
        # get window around gene and subset the vcf for that window using bcftools
        snp2geno = get_cis_snps(vcf, bcftools, txlist, tx2info, tx2expr, meta_lines, window)

        # here we run find_peak_snp for each snp separately so that each
        #   snp-gene result will be printed
        for snp in snp2geno:
            tmp_snp2geno = pd.DataFrame({snp: snp2geno[snp]})
            gene_snp_pairname = gene + '_' + snp

            # find and store the peak SNP for this gene and its association statistic and p-value
            if n_perms > 0:
                permuted_tx2expr = permute_transcripts(tx2expr, txlist, n_perms)

                snp, stat, pval, perm_peak_pvals = find_peak_snp(tx2expr, tmp_snp2geno, txlist, nominal_thresh, permuted_tx2expr, n_perms, methods)

                for m in methods:
                    if pval[m] < nominal_thresh:
                        dir_perm_pval = (1.0 + sum([pval[m] >= perm_i for perm_i in perm_peak_pvals[m]])) / float(n_perms + 1)
                        #if np.var(perm_peak_pvals) == 0.0:
                        if all(perm_peak_pvals[m] > 0.999) and pval[m] < 0.999:
                            beta_perm_pval = sys.float_info.min # dir_perm_pval
                        elif dir_perm_pval == 1.0:
                            beta_perm_pval = 1.0
                        else:
                            beta_perm_pval = compute_beta_perm_pvals(pval[m], perm_peak_pvals[m])
                        perm_pvals[m][gene_snp_pairname] = [dir_perm_pval, beta_perm_pval]
            else:
                snp, stat, pval, perm_peak_pvals = find_peak_snp(tx2expr, tmp_snp2geno, txlist, nominal_thresh, {}, n_perms, methods)
            for m in methods:
                gene_results[m][gene_snp_pairname] = [snp[m], stat[m], pval[m]]
    fnull.close()
    return gene_results, perm_pvals
