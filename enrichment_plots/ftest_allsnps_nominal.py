import argparse
import gzip
import sys
import os
import subprocess
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
#from cis_pass import nominal_pass

from io import StringIO
from scipy.special import betaln as betaln_func
from scipy.special import ncfdtr
from scipy.stats import linregress, beta, chi2, norm, t, f, combine_pvalues
from scipy.optimize import minimize


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Run F-test nominal pass.')
    parser.add_argument('--vcf', required=True, help='VCF or BCF file with genotypes. Required.')
    parser.add_argument('--pheno', required=True, help='BED file or tsv file with transcript expression levels.')
    parser.add_argument('--bcftools', default='bcftools', help='Path to bcftools executable ("bcftools" by default).')
    parser.add_argument('--covariates', default='NONE', help='tsv file listing covariates to adjust phenotypes for.')
    parser.add_argument('--nominal', default=0.5, type=float,
                        help='Print genes with a nominal assocation p-value below this cutoff.')
    parser.add_argument('--output', default='eqtl_results', help='Base name for results.')
    parser.add_argument('--window', default=50000, type=int,
                        help='Size of window in bp around start position of phenotypes.')
    args = parser.parse_args()
    return args


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


def multi_regress_f_test_given_xtx(y, x, xtx_inv):
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


def nominal_pass(vcf, tx2info, tx2expr, gene_to_tx, meta_lines, window, bcftools, nominal_thresh):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    gene_results = {}

    for gene in gene_to_tx:
        txlist = gene_to_tx[gene]
        gene_results[gene] = []
        # get window around gene and subset the vcf for that window using bcftools
        snp2geno = get_cis_snps(vcf, bcftools, txlist, tx2info, tx2expr, meta_lines, window)
        n_genos = len(snp2geno.iloc[:, 0])
        dof = n_genos - 2
        x, xtx_inv = precomp_xtx_inv(tx2expr, txlist)

        for snp in snp2geno:
            genos = np.array(snp2geno[snp]).astype(np.double)
            if len(txlist) == 1:
                beta, inter, r, pval, stderr = linregress(tx2expr[txlist[0]], genos)
                stat = beta / stderr
            else:
                genos_var = np.var(genos)
                if genos_var == 0.0:
                    continue
                stat, pval = multi_regress_f_test_given_xtx(snp2geno[snp], x, xtx_inv)
            if pval < nominal_thresh:
                gene_results[gene].append([snp, stat, pval])
    fnull.close()
    return gene_results



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


def check_vcf(vcf, bcftools, tx2expr):
    # check how many metadata lines there are and ensure that all phenotyped IDs are also genotyped
    bcf_proc = subprocess.Popen([bcftools, 'view', vcf, '--header-only'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    meta_lines = bcf_proc.communicate()[0].decode('utf-8').count('\n') - 1
    compressed = vcf.endswith('.gz')
    if compressed:
        infile = gzip.open(vcf, 'r')
    else:
        infile = open(vcf, 'r')
    for line in infile:
        if compressed:
            line = line.decode('ascii')  # convert binary to string
        if not line.startswith('#CHROM'):
            continue
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


def write_results(results, out_fname, nominal_cutoff):
    with(open(out_fname, 'w')) as outfile:
        outfile.write('#Gene\tSNP\tStatistic\tNominal P-value\n')
        for gene in results:
            for vec in results[gene]:
                snp, stat, pval = vec
                fields = [gene, snp]
                fields.extend(['{:.6}'.format(i) for i in [float(stat), float(pval)]])
                outfile.write('\t'.join(fields) + '\n')


if __name__ == "__main__":
    args = parseargs()
    # print('Reading in phenotypes file...')
    tx2info, tx2expr, gene_to_tx = read_pheno_file(args.pheno, args.covariates)
    meta_lines = check_vcf(args.vcf, args.bcftools, tx2expr)
    # print('Performing nominal pass...')
    nominal_results = nominal_pass(
        args.vcf, tx2info, tx2expr, gene_to_tx, meta_lines, args.window, args.bcftools, args.nominal)
    outname = args.output + '.ftest.tsv'
    write_results(nominal_results, outname, args.nominal)
