import argparse
import gzip
import sys
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, linregress
# import statsmodels.api as sm
import statsmodels.formula.api as smf


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='isoQTL main script.')
    parser.add_argument('--vcf', required=True, help='VCF or BCF file with genotypes. Required.')
    parser.add_argument('--pheno', required=True, help='BED file or tsv file with transcript expression levels.')
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


def read_pheno_file(pheno_fname, covariates):
    pheno_df = pd.read_csv(pheno_fname, sep='\t', index_col=3).T  # read in phenotype file
    # split into info about the tx and expression levels into separate dataframes
    tx2info, tx2expr = pheno_df.iloc[:5], pheno_df.iloc[5:]
    if covariates != 'NONE':
        tx2expr = regress_out_covars(covariates, tx2expr)  # regress out covariates from phenotypes
    else:  # still need to convert to numeric
        tx2expr = tx2expr.apply(pd.to_numeric, errors='coerce')
    return tx2info, tx2expr


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


def create_tx_windows(tx2info, window):
    # made for efficient mapping of a SNP location to transcripts in that location
    half_window = window / 2.0
    location_to_tx = {}  # maps chromosome to megabase to coordinates of transcript windows overlapping with that Mb
    for tx in tx2info:
        chrom, start, end, gid, strand = tx2info[tx]
        start, end = int(start), int(end)
        window_start = max(start - half_window, 0)
        window_end = end + half_window
        window_start_mb = int(window_start / 1000000)
        window_end_mb = int(window_end / 1000000)
        if chrom not in location_to_tx:
            location_to_tx[chrom] = {}
        if window_start_mb not in location_to_tx[chrom]:
            location_to_tx[chrom][window_start_mb] = []
        location_to_tx[chrom][window_start_mb].append([tx, window_start, window_end])
        if window_end_mb != window_start_mb:
            if window_end_mb not in location_to_tx[chrom]:
                location_to_tx[chrom][window_end_mb] = []
            location_to_tx[chrom][window_end_mb].append([tx, window_start, window_end])
    return location_to_tx


def preprocess_snp2geno(snp2geno):
    # remove genotypes for non-phenotyped individuals and remove "nan" entries
    snp2geno = (snp2geno.T[list(tx2expr.index)]).T  # remove genotypes for non-phenotyped individuals
    new_snp2geno = pd.DataFrame()
    for i in snp2geno:
        if snp2geno[i].isnull().values.any():
            continue  # exclude SNPs with null genotypes
        new_snp2geno[i] = pd.to_numeric(snp2geno[i])
    return new_snp2geno


def find_peak_snps(snp2info, snp2geno, location_to_tx, tx_to_peak_snp):
    # find transcript windows that each snp is in and see if this snp is the most correlated snp (so far) with the tx
    for snp in snp2geno:
        chrom, pos = snp2info[snp]['#CHROM'], snp2info[snp]['POS']
        pos = int(pos)
        pos_mb = int(pos / 1000000)
        if chrom in location_to_tx:
            if pos_mb in location_to_tx[chrom]:
                # iterate through transcripts and find which this snp is in the window for
                tx_window_list = location_to_tx[chrom][pos_mb]
                for tx, window_start, window_end in tx_window_list:
                    if window_start <= pos <= window_end:  # snp is within this transcript's window
                        # print(tx2expr[tx])
                        # print(snp2geno[snp])
                        r_squared = pearsonr(tx2expr[tx], snp2geno[snp])[0] ** 2
                        # update peak snp if this tx doesn't have a peak or this r_sq is stronger than previous high
                        if tx not in tx_to_peak_snp or r_squared > tx_to_peak_snp[tx][1]:
                            tx_to_peak_snp[tx] = [snp, r_squared, list(snp2geno[snp])]
    return tx_to_peak_snp


def regress_peak_snps(tx_to_peak_snp):
    # now perform full linear regression of transcripts against peak snps
    for tx in tx_to_peak_snp:
        snp, peak_r2, snp_genos = tx_to_peak_snp[tx]
        regression_results = linregress(tx2expr[tx], snp_genos)
        beta, stderr, pval = regression_results[0], regression_results[4], regression_results[3]
        tx_to_peak_snp[tx] = [snp, peak_r2, beta, stderr, pval]
    return tx_to_peak_snp


def nominal_pass(vcf, tx2expr, location_to_tx, meta_lines):
    # read thru chunks of vcf, finding which SNPs are in which transcripts' windows and keeping an updated dict mapping
    #    transcripts to their peak SNPs, then perform linear regression of transcripts against peak SNP genotypes
    tx_to_peak_snp = {}
    lines_read = 0
    chunksize = 10 ** 4
    with pd.read_csv(vcf, chunksize=chunksize, sep='\t', header=meta_lines, index_col=2) as reader:
        for chunk in reader:
            chunk = chunk.T
            snp2info, snp2geno = chunk.iloc[:8], chunk.iloc[8:]
            # print(snp2geno.head())
            # print(tx2expr.head())
            snp2geno = preprocess_snp2geno(snp2geno)
            tx_to_peak_snp = find_peak_snps(snp2info, snp2geno, location_to_tx, tx_to_peak_snp)
            lines_read += chunksize
            print('Lines processed from vcf: ' + str(lines_read))
    tx_to_peak_snp = regress_peak_snps(tx_to_peak_snp)
    return tx_to_peak_snp


def write_results(results, out_fname, nominal_cutoff):
    with(open(out_fname, 'w')) as outfile:
        outfile.write('#Transcript\tSNP\tR^2\tBeta\tStdErr\tP-value\n')
        for tx in results:
            snp, peak_r2, beta, stderr, pval = results[tx]
            if pval > nominal_cutoff:
                continue
            fields = [str(i) for i in [tx, snp, peak_r2, beta, stderr, pval]]
            outfile.write('\t'.join(fields) + '\n')


if __name__ == "__main__":
    args = parseargs()
    print('Reading in phenotypes file...')
    tx2info, tx2expr = read_pheno_file(args.pheno, args.covariates)
    meta_lines = check_vcf(args.vcf, tx2expr)
    print('Computing cis windows around isoforms...')
    location_to_tx = create_tx_windows(tx2info, args.window)
    print('Performing nominal pass...')
    nominal_results = nominal_pass(args.vcf, tx2expr, location_to_tx, meta_lines)
    write_results(nominal_results, args.output, args.nominal)
