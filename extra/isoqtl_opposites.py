# given a list of SNP-gene pairs from a previous run, prints out the
#   individual isoform-SNP association results for those SNP-gene pairs

import argparse
import gzip
import sys
import os
import pickle
import subprocess
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from cis_pass_opposites import nominal_pass


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='isoQTL main script.')
    parser.add_argument('--vcf', required=True, help='VCF or BCF file with genotypes. Required.')
    parser.add_argument('--pheno', required=True, help='BED file or tsv file with transcript expression levels.')
    parser.add_argument('--results', required=True, help='Results from a previous run. Required.')
    parser.add_argument('--qtltools', action='store_true', help='Use if results are qtltools.')
    parser.add_argument('--bcftools', default='bcftools', help='Path to bcftools executable ("bcftools" by default).')
    parser.add_argument('--covariates', default='NONE', help='tsv file listing covariates to adjust phenotypes for.')
    parser.add_argument('--methods', nargs='+', default=['wilks', 'fisher', 'min', 'cauchy', 'ftest'],
        choices=['wilks', 'fisher', 'min', 'cauchy', 'ftest'], help='Specify method to obtain gene-level p-value.')
    parser.add_argument('--nominal', default=0.5, type=float,
                        help='Print genes with a nominal assocation p-value below this cutoff.')
    parser.add_argument('--output', default='eqtl_results', help='Base name for results.')
    parser.add_argument('--permute', default=100, type=int,
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


if __name__ == "__main__":
    args = parseargs()
    tx2info, tx2expr, gene_to_tx = read_pheno_file(args.pheno, args.covariates)
    meta_lines = check_vcf(args.vcf, args.bcftools, tx2expr)
    if args.qtltools:
        prev = pd.read_csv(args.results, sep=' ', header=None, index_col=0, usecols=[0,7]).transpose()
    else:
        prev = pd.read_csv(args.results, sep=' ', header=None, index_col=0, usecols=[0,1]).transpose()
    gene_results = nominal_pass(
        prev, args.vcf, tx2info, tx2expr, gene_to_tx, meta_lines, args.window, args.bcftools, args.permute, args.nominal, args.methods, args.qtltools)
    with open(args.output, 'wb') as handle:
      pickle.dump(gene_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
