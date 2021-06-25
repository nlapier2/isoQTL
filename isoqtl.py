import argparse
import gzip
import sys
import os
import subprocess
import pandas as pd
import numpy as np
from scipy.special import betaln as betaln_func
from scipy.special import ncfdtr
from scipy.stats import linregress, beta
from scipy.optimize import minimize
# import statsmodels.api as sm
import statsmodels.formula.api as smf


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='isoQTL main script.')
    parser.add_argument('--vcf', required=True, help='VCF or BCF file with genotypes. Required.')
    parser.add_argument('--pheno', required=True, help='BED file or tsv file with transcript expression levels.')
    parser.add_argument('--bcftools', default='bcftools', help='Path to bcftools executable ("bcftools" by default).')
    parser.add_argument('--covariates', default='NONE', help='tsv file listing covariates to adjust phenotypes for.')
    parser.add_argument('--nominal', default=0.05, type=float,
                        help='Print transcripts with a nominal assocation p-value below this cutoff.')
    parser.add_argument('--output', default='isoqtl_results.tsv', help='Where to output results to.')
    parser.add_argument('--permute', default=1000, type=int,
                        help='Number of permutations to do for a permutation pass. Set to 0 to do nominal pass only.')
    parser.add_argument('--window', default=50000, type=int,
                        help='Size of window in bp around start position of phenotypes.')
    args = parser.parse_args()
    return args


def check_pheno_ids_in_covars(tx2expr, covar_df):
    # Ensure that all people in tx2expr are in covariates file
    id_list = list(tx2expr.index)
    for iid in id_list:
        if iid not in covar_df:
            sys.exit('Error: ID ' + iid + ' not found in covariates file.')


def regress_out_covars(covariates_fname, tx2expr):
    # add covariates into expression dataframe while maintaining variable names
    covar_df = pd.read_csv(covariates_fname, sep='\t', index_col=0)
    check_pheno_ids_in_covars(tx2expr, covar_df)
    individuals = list(tx2expr.index)
    tx_names = list(tx2expr.columns)
    covar_names = list(covar_df.index)
    all_var_names = tx_names + covar_names
    df_with_covars = tx2expr.T.merge(covar_df, how='outer')[individuals].T
    df_with_covars.columns = all_var_names
    # coerce pandas datatype from object to numeric and standardize, then regress out covariates
    covars_formula = ' + '.join(covar_names)
    for i in tx_names:
        df_with_covars[i] = pd.to_numeric(df_with_covars[i])  # coerce to numeric
        mean = np.mean(df_with_covars[i])
        stddev = np.std(df_with_covars[i])
        df_with_covars[i] = (df_with_covars[i] - mean) / stddev  # standardize
        formula = 'Q("' + i + '") ~ ' + covars_formula
        df_with_covars[i] = smf.ols(formula, data=df_with_covars).fit().resid  # regress out covars
    return df_with_covars[tx_names]


def get_gene_to_tx(tx2info):
    gene_to_tx = {}
    for tx in tx2info:
        gene = tx2info[tx][3]
        if gene not in gene_to_tx:
            gene_to_tx[gene] = []
        gene_to_tx[gene].append(tx)
    return gene_to_tx


def read_pheno_file(pheno_fname, covariates):
    pheno_df = pd.read_csv(pheno_fname, sep='\t', index_col=3).T  # read in phenotype file
    # split into info about the tx and expression levels into separate dataframes
    tx2info, tx2expr = pheno_df.iloc[:5], pheno_df.iloc[5:]
    if covariates != 'NONE':
        tx2expr = regress_out_covars(covariates, tx2expr)  # regress out covariates from phenotypes
    else:  # still need to convert to numeric
        tx2expr = tx2expr.apply(pd.to_numeric, errors='coerce')
    gene_to_tx = get_gene_to_tx(tx2info)
    return tx2info, tx2expr, gene_to_tx


def check_vcf(vcf, tx2expr):
    # check how many metadata lines there are and ensure that all phenotyped IDs are also genotyped
    compressed = vcf.endswith('.gz')
    if compressed:
        infile = gzip.open(vcf, 'r')
    else:
        infile = open(vcf, 'r')
    meta_lines = 0
    for line in infile:
        if compressed:
            line = line.decode('ascii')  # convert binary to string
        if not line.startswith('#CHROM'):
            meta_lines += 1
        else:  # line starts with #CHROM -- in other words, this is the header line
            splits = line.strip().split('\t')[9:]  # get a list of the people in this vcf
            id_dict = {iid: True for iid in splits}  # convert to dict for indexing speed
            # ensure that all phenotyped individuals are also genotyped
            pheno_iid_list = list(tx2expr.index)
            for iid in pheno_iid_list:
                if iid not in id_dict:
                    sys.exit('Error: ID ' + iid + ' not found in genotypes file.')
            break  # here we are only looking for the header, not reading the whole file, so we are done
    infile.close()
    return meta_lines


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


def get_betas_stderrs(snp_genos, tx2expr, txlist):
    betas, stderrs, residuals = [], [], []
    for tx in txlist:
        # regression_results = linregress(tx2expr[tx], snp_genos)
        regression_results = linregress(snp_genos, tx2expr[tx])
        beta, stderr = regression_results[0], regression_results[4]
        resid = tx2expr[tx] - (snp_genos * beta)
        betas.append(beta)
        stderrs.append(stderr)
        residuals.append(resid)
    return np.array(betas), np.array(stderrs), np.array(residuals)


# implement test proposed by Hashimoto & Ohtani, Econometrics Letters (1990)
def hashimoto_ohtani_t2_stat(betas, residuals, genos):
    p = M = len(betas)  # number of isoforms
    # R = np.identity(M)  # linear restrictions on betas
    T = len(genos)  # sample size
    K = 1  # num regressors
    r = int((T-K) / K)
    dof = r - p + 1
    X = genos  # (genos - np.mean(genos)) / np.std(genos)
    Q = np.sqrt(1 / np.dot(X, X))
    # Q = np.array([1 / np.dot(X, X)])
    # Q = np.invert(np.matmul(np.transpose(X), X))
    Sigma = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            Sigma[i][j] = np.dot(residuals[i], residuals[j]) / (T - K)

    N = np.identity(T) - np.outer(X * (1 / np.dot(X, X)), np.transpose(X))
    eigvals, eigvecs = np.linalg.eigh(N)
    eigvals_1 = [i for i in range(len(eigvals)) if 0.999 < eigvals[i] < 1.001]
    z = np.array([v[eigvals_1] for v in eigvecs])
    resid_star = np.array([np.matmul(np.transpose(z), residuals[i]) for i in range(len(residuals))])

    S = np.zeros((M, M))
    for j in range(r):  # T):
        vec_e_j = np.array([e[j] for e in resid_star])
        # kron_prod = np.kron(np.identity(M), np.transpose(Q))
        kron_prod = np.identity(M) * Q
        eta_j = np.matmul(kron_prod, np.transpose(vec_e_j))
        S += np.outer(eta_j, np.transpose(eta_j))
    S /= r

    hotelling_t2_stat = (np.matmul(np.matmul(np.transpose(betas), np.linalg.inv(S)), betas) / r) * ((r - p + 1) / p)
    ncp = np.matmul(np.matmul(betas, np.linalg.inv(Sigma / np.dot(X, X))), betas)
    pval = 1 - ncfdtr(p, dof, ncp, hotelling_t2_stat)
    # print(p, dof, ncp, hotelling_t2_stat, pval)
    return hotelling_t2_stat, pval


def find_peak_snp(tx2expr, snp2geno, txlist, permuted_tx2expr):
    peak_snp, peak_stat, peak_pval = 'rsid', 0, 1
    perm_peak_pvals = [1.0 for i in range(len(permuted_tx2expr))]
    for snp in snp2geno:
        betas, stderrs, residuals = get_betas_stderrs(snp2geno[snp], tx2expr, txlist)
        stat, pval = hashimoto_ohtani_t2_stat(betas, residuals, np.array(snp2geno[snp]))
        if pval < peak_pval:
            peak_snp = snp
            peak_stat = stat
            peak_pval = pval
        for i in range(len(permuted_tx2expr)):
            betas, stderrs, residuals = get_betas_stderrs(snp2geno[snp], permuted_tx2expr[i], txlist)
            stat, pval = hashimoto_ohtani_t2_stat(betas, residuals, np.array(snp2geno[snp]))
            if pval < perm_peak_pvals[i]:
                perm_peak_pvals[i] = pval
    return peak_snp, peak_stat, peak_pval, perm_peak_pvals


def set_initial_beta_params(pval_vec):
    pval_mean = np.mean(pval_vec)
    pval_var = np.var(pval_vec)
    beta_shape1 = pval_mean * (pval_mean * (1.0 - pval_mean) / pval_var - 1.0)
    beta_shape2 = beta_shape1 * (1.0 / pval_mean - 1.0)
    return beta_shape1, beta_shape2


# def beta_likelihood_function(beta_shape1, beta_shape2, pval_vec):
def beta_likelihood_function(beta_params, other_params):
    beta_shape1, beta_shape2 = beta_params
    likelihood_sum1, likelihood_sum2, n_var = other_params
    beta_func_res = betaln_func(beta_shape1, beta_shape2)
    return -1.0 * ((beta_shape1 - 1) * likelihood_sum1 + (beta_shape2 - 1) * likelihood_sum2 - n_var * beta_func_res)


def compute_beta_perm_pvals(nominal_pval, perm_peak_pvals):
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


def nominal_pass(vcf, tx2info, tx2expr, gene_to_tx, meta_lines, window, bcftools, n_perms):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    gene_results, perm_pvals = {}, {}
    for gene in gene_to_tx:
        txlist = gene_to_tx[gene]
        # get window around gene and subset the vcf for that window using bcftools
        window_start, window_end, chrom = get_gene_window(txlist, tx2info, window)
        window_str = str(chrom) + ':' + str(window_start) + '-' + str(window_end)
        subprocess.Popen([bcftools, 'view', vcf, '-r', window_str, '-o', 'TEMP_isoqtl'], stderr=fnull).wait()
        # read in SNPs in the window and clear out null SNPs and non-phenotyped individuals
        gene_window_snps = pd.read_csv('TEMP_isoqtl', sep='\t', header=meta_lines+4, index_col=2).T
        subprocess.Popen(['rm', 'TEMP_isoqtl']).wait()
        snp2info, snp2geno = gene_window_snps.iloc[:8], gene_window_snps.iloc[8:]
        snp2geno = preprocess_snp2geno(tx2expr, snp2geno)
        # find and store the peak SNP for this gene and its association statistic and p-value
        if n_perms > 0:
            permuted_tx2expr = permute_transcripts(tx2expr, txlist, n_perms)
            snp, fstat, pval, perm_peak_pvals = find_peak_snp(tx2expr, snp2geno, txlist, permuted_tx2expr)
            dir_perm_pval = (1.0 + sum([pval >= perm_i for perm_i in perm_peak_pvals])) / float(n_perms + 1)
            if dir_perm_pval < 1.0:
                beta_perm_pval = compute_beta_perm_pvals(pval, perm_peak_pvals)
            else:
                beta_perm_pval = 1.0
            perm_pvals[gene] = [dir_perm_pval, beta_perm_pval]
        else:
            snp, fstat, pval, perm_peak_pvals = find_peak_snp(tx2expr, snp2geno, txlist, {})
        gene_results[gene] = [snp, fstat, pval]
    fnull.close()
    return gene_results, perm_pvals


def write_results(results, perm_results, out_fname, nominal_cutoff):
    with(open(out_fname, 'w')) as outfile:
        if len(perm_results) > 0:
            outfile.write('#Gene\tSNP\tF-Statistic\tP-value\tDir. Perm. P-value\tBeta Perm. P-value\n')
        else:
            outfile.write('#Gene\tSNP\tF-Statistic\tP-value\n')
        for gene in results:
            snp, fstat, pval = results[gene]
            if pval > nominal_cutoff:
                continue
            fields = [gene, snp]
            fields.extend(['{:.6}'.format(i) for i in [float(fstat), float(pval)]])
            if len(perm_results) > 0:
                dir_perm_pval, beta_perm_pval = perm_results[gene]
                fields.append('{:.6}'.format(dir_perm_pval))
                fields.append('{:.6}'.format(beta_perm_pval))
            outfile.write('\t'.join(fields) + '\n')


if __name__ == "__main__":
    args = parseargs()
    # print('Reading in phenotypes file...')
    tx2info, tx2expr, gene_to_tx = read_pheno_file(args.pheno, args.covariates)
    meta_lines = check_vcf(args.vcf, tx2expr)
    # print('Performing nominal pass...')
    nominal_results, perm_results = nominal_pass(
        args.vcf, tx2info, tx2expr, gene_to_tx, meta_lines, args.window, args.bcftools, args.permute)
    write_results(nominal_results, perm_results, args.output, args.nominal)
