import argparse
import gzip
import sys
import os
import subprocess
import pandas as pd
import numpy as np
from scipy.stats import f, linregress  # , pearsonr
from scipy.special import ncfdtr
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
    parser.add_argument('--window', default=100000, type=int,
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


def hotelling_test(betas, covar_mat, sample_size):
    # compute the parameters and inputs for the hotelling, then compute the stat
    num_params = len(betas)
    dof = sample_size - num_params - 1  # sample_size * params - params
    # covar_mat = np.identity(num_params)  # np.outer(stderrs, stderrs)
    # hotelling_stat = dof * np.matmul(np.matmul(betas, np.linalg.inv(covar_mat)), betas)
    hotelling_stat = np.matmul(np.matmul(betas, np.linalg.inv(covar_mat)), betas)
    # convert hotelling stat to f stat for scipy-supported hypothesis testing
    f_stat = (dof - num_params + 1) / (dof * num_params) * hotelling_stat
    f_dfn = num_params
    f_dfd = dof - num_params + 1
    pval = f.sf(f_stat, f_dfn, f_dfd)
    return f_stat, pval


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
    '''
    # convert hotelling stat to f stat for scipy-supported hypothesis testing
    f_stat = (dof - M + 1) / (dof * M) * hotelling_t2_stat
    f_dfn = M
    f_dfd = dof - M + 1
    pval = f.sf(f_stat, f_dfn, f_dfd)
    return f_stat, pval
    '''
    return hotelling_t2_stat, pval


'''
def get_covar_mat(tx2expr, txlist):
    covar_mat = np.zeros((len(txlist), len(txlist)))
    for i in range(len(txlist)):
        for j in range(i+1):
            tx1 = tx2expr[txlist[i]]
            tx2 = tx2expr[txlist[j]]
            corr = pearsonr(tx1, tx2)[0]
            covar_mat[i][j] = corr
            covar_mat[j][i] = corr
    return covar_mat
'''


'''
def get_covar_mat(stderrs, genos):
    covar_mat = np.outer(stderrs, stderrs)
    # xtx = np.dot(genos, genos)
    sample_size = len(genos)
    adjust_dof = (sample_size - 2) / (sample_size - 3)
    for i in range(len(covar_mat)):
        for j in range(len(covar_mat)):
            if i != j:
                covar_mat[i][j] = covar_mat[i][j] * 0  # adjust_dof  # * xtx
    return covar_mat
'''


def get_covar_mat(betas, residuals, genos):
    n_vars = len(betas)
    sample_size = len(genos)
    xtx = np.dot(genos, genos)
    xtx_inv = 1.0 / xtx
    covar_mat = np.zeros((n_vars, n_vars))
    for i in range(len(covar_mat)):
        for j in range(len(covar_mat)):
            # print(residuals[i].shape)
            covar_mat[i][j] = (1 / sample_size) * np.dot(residuals[i], residuals[j]) * xtx_inv
    return covar_mat


def find_peak_snp(tx2expr, snp2geno, txlist):
    peak_snp, peak_stat, peak_pval = 'rsid', 0, 1
    # covar_mat = get_covar_mat(tx2expr, txlist)
    for snp in snp2geno:
        sample_size = len(snp2geno[snp])
        betas, stderrs, residuals = get_betas_stderrs(snp2geno[snp], tx2expr, txlist)
        # covar_mat = get_covar_mat(stderrs, snp2geno[snp])
        # covar_mat = get_covar_mat(betas, residuals, snp2geno[snp])
        # stat, pval = hotelling_test(betas, stderrs, sample_size)
        # stat, pval = hotelling_test(betas, covar_mat, sample_size)
        stat, pval = hashimoto_ohtani_t2_stat(betas, residuals, np.array(snp2geno[snp]))
        if pval < peak_pval:
            peak_snp = snp
            peak_stat = stat
            peak_pval = pval
    return peak_snp, peak_stat, peak_pval


def nominal_pass(vcf, tx2info, tx2expr, gene_to_tx, meta_lines, window, bcftools):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    gene_results = {}
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
        snp, fstat, pval = find_peak_snp(tx2expr, snp2geno, txlist)
        gene_results[gene] = [snp, fstat, pval]
    fnull.close()
    return gene_results


def write_results(results, out_fname, nominal_cutoff):
    with(open(out_fname, 'w')) as outfile:
        outfile.write('#Gene\tSNP\tF-Statistic\tP-value\n')
        for gene in results:
            snp, fstat, pval = results[gene]
            if pval > nominal_cutoff:
                continue
            fields = [str(i) for i in [gene, snp, fstat, pval]]
            outfile.write('\t'.join(fields) + '\n')


if __name__ == "__main__":
    args = parseargs()
    # print('Reading in phenotypes file...')
    tx2info, tx2expr, gene_to_tx = read_pheno_file(args.pheno, args.covariates)
    meta_lines = check_vcf(args.vcf, tx2expr)
    # print('Performing nominal pass...')
    nominal_results = nominal_pass(args.vcf, tx2info, tx2expr, gene_to_tx, meta_lines, args.window, args.bcftools)
    write_results(nominal_results, args.output, args.nominal)
