import argparse
import gzip
import os
import random
import subprocess
import sys
import pandas as pd
import numpy as np
from scipy.stats import matrix_normal
from io import StringIO


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='isoQTL main script.')
    parser.add_argument('--vcf', required=True, help='VCF or BCF file with genotypes. Required.')
    parser.add_argument('--pheno', required=True, help='BED file or tsv file with transcript expression levels.')
    parser.add_argument('--bcftools', default='bcftools', help='Path to bcftools executable ("bcftools" by default).')
    parser.add_argument('--dropout_rate', default=0.25, type=float,
                        help='Percent of genes to randomly set to zero heritability.')
    parser.add_argument('--dropout_iso', default=0.0, type=float,
                        help='Percent of isoforms to randomly set to zero heritability.')
    parser.add_argument('--h2cis', default=0.1, type=float, help='Heritability of cis SNPs on gene expression.')
    parser.add_argument('--h2noncis', default=0.0, type=float,
                        help='Heritability of non-cis effects on gene expression.')
    parser.add_argument('--max_corr', default=0.64, type=float, help='Max r^2 betweeen isoforms.')
    parser.add_argument('--max_corr_env', default=0.64, type=float, help='Max r^2 between environments.')
    parser.add_argument('--min_corr', default=0.04, type=float, help='Min r^2 between isoforms.')
    parser.add_argument('--min_corr_env', default=0.04, type=float, help='Min r^2 between environments.')
    parser.add_argument('--neg_pct', default=0.5, type=float, help='Percent of isoforms negatively correlated.')
    parser.add_argument('--num_causal', default=-1, type=int,
        help='Specify a fixed number of causal SNPs per gene.')
    parser.add_argument('--num_iso', default=-1, type=int,
                        help='Use to include genes only with a certain number of isoforms.')
    parser.add_argument('--output', default='simulated_phenos', help='Output file base name.')
    parser.add_argument('--window', default=50000, type=int,
                        help='Size of window in bp around start position of phenotypes.')
    args = parser.parse_args()
    return args


def get_gene_to_tx(tx2info):
    gene_to_tx = {}
    for tx in tx2info:
        gene = tx2info[tx][3]
        if gene not in gene_to_tx:
            gene_to_tx[gene] = []
        gene_to_tx[gene].append(tx)
    return gene_to_tx


def read_pheno_file(pheno_fname):
    pheno_df = pd.read_csv(pheno_fname, sep='\t', index_col=3).T  # read in phenotype file
    # split into info about the tx and expression levels into separate dataframes
    tx2info, tx2expr = pheno_df.iloc[:5], pheno_df.iloc[5:]
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


def preprocess_snp2geno(tx2expr, snp2geno):
    # remove genotypes for non-phenotyped individuals and remove "nan" entries
    snp2geno = (snp2geno.T[list(tx2expr.index)]).T  # remove genotypes for non-phenotyped individuals
    new_snp2geno = pd.DataFrame()
    for i in snp2geno:
        if snp2geno[i].isnull().values.any():
            continue  # exclude SNPs with null genotypes
        snp_i = pd.to_numeric(snp2geno[i])
        snp_i = (snp_i - snp_i.mean()) / snp_i.std()  # standardize genotype
        new_snp2geno[i] = snp_i
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
    return int(window_start), int(window_end), int(last_end), chrom
 

def create_covar_mat(num_iso, args_min_corr, args_max_corr, args_neg_pct):
    # simulate MVN covariance matrix, given min and max covars
    mean_vec, covar_mat = np.zeros(num_iso), np.zeros((num_iso, num_iso))
    sqrt_min, sqrt_max = np.sqrt(args_min_corr), np.sqrt(args_max_corr)
    for i in range(num_iso):
        for j in range(num_iso):
            if i == j:
                covar_mat[i][j] = 1.0
            elif j > i:
                choice = np.random.uniform(sqrt_min, sqrt_max)
                if np.random.uniform(0.0, 1.0) <= args_neg_pct:
                    choice *= -1.0  # random sign
                covar_mat[i][j] = choice
            else:
                covar_mat[i][j] = covar_mat[j][i]
    return mean_vec, covar_mat


def sim_iso_effect(num_iso, h2, min_corr, max_corr, neg_pct):
    # simulate correlated effect sizes on all isoforms
    min_eff = h2 / 2
    if num_iso == 1:
        iso_effect = np.random.normal(0, np.sqrt(h2))
        while abs(iso_effect) < min_eff:
            iso_effect = np.random.normal(0, np.sqrt(h2))
        iso_effect = [iso_effect]
    else:
        if h2 == 0.0:
            iso_effect = np.zeros(num_iso)
        else:
            mean_vec, covar_mat = create_covar_mat(num_iso, min_corr, max_corr, neg_pct)
            while min(np.linalg.eigh(covar_mat)[0]) < 0:  # ensure covar_mat is positive semidefinite
                mean_vec, covar_mat = create_covar_mat(num_iso, min_corr, max_corr, neg_pct)
            iso_effect = np.random.multivariate_normal(mean_vec, covar_mat) * np.sqrt(h2)
            while min([abs(i) for i in iso_effect]) < min_eff:
                iso_effect = np.random.multivariate_normal(mean_vec, covar_mat) * np.sqrt(h2)
    return iso_effect
    
    
def sim_noncis_matrix(num_iso, num_people, h2, min_corr, max_corr, neg_pct):
    if h2 == 0.0:
        return np.zeros((num_iso, num_people))
    mean_vec, iso_covar_mat = create_covar_mat(num_iso, min_corr, max_corr, neg_pct)
    mean_vec, ppl_covar_mat = create_covar_mat(num_people, 0.99, 0.99, neg_pct)
    iso_covar_mat = np.dot(iso_covar_mat, iso_covar_mat.T)
    iso_covar_mat = iso_covar_mat / np.max(iso_covar_mat)
    ppl_covar_mat = np.dot(ppl_covar_mat, ppl_covar_mat.T)
    ppl_covar_mat = ppl_covar_mat / np.max(ppl_covar_mat)
    # iso_covar_mat, ppl_covar_mat = np.sqrt(h2) * iso_covar_mat, np.sqrt(h2) * ppl_covar_mat  # scale by h2
    # we don't use mean_vec here; scipy matrix normal default mean is zero matrix
    noncis_effect_mat = matrix_normal.rvs(rowcov=iso_covar_mat, colcov=ppl_covar_mat)
    return noncis_effect_mat * np.sqrt(h2)


def simulate_phenotypes(args, txlist, causal_snp_df, snp_h2_eff):
    sim_tx2expr, sim_gene2expr = {}, np.zeros(len(causal_snp_df))
    all_snp_h2 = [snp_h2_eff[snp] for snp in snp_h2_eff]
    num_iso, total_h2 = len(txlist), sum(all_snp_h2)
    num_people = len(causal_snp_df)
    noncis_effect_mat = sim_noncis_matrix(num_iso, num_people, args.h2noncis, args.min_corr_env, args.max_corr_env, args.neg_pct)
    # first sim genetic component: the effects of the SNPs on genes & isoforms
    for snp in causal_snp_df:
        snp_h2 = snp_h2_eff[snp]
        cis_effect = sim_iso_effect(num_iso, snp_h2, args.min_corr, args.max_corr, args.neg_pct)
        # noncis_effect = sim_iso_effect(num_iso, args.h2noncis, args.min_corr_env, args.max_corr_env)
        # num_people = len(causal_snp_df[snp])
        # noncis_effect_mat = sim_noncis_matrix(num_iso, num_people, args.h2noncis, args.min_corr_env, args.max_corr_env)
        
        # genetic component of expression = SNP * genetic effect for each tx
        all_iso_dropped = True
        for i in range(len(txlist)):
            tx = txlist[i]
            if random.random() < args.dropout_iso and (not all_iso_dropped or i + 1 != len(txlist)):
                # randomly drop out isoforms, but keep at least one to avoid having zero effect on all isoforms
                cis_component = causal_snp_df[snp] * 0.0
            else:
                cis_component = causal_snp_df[snp] * cis_effect[i]
                all_iso_dropped = False
            # noncis_component = noncis_effect[i] * np.ones(len(causal_snp_df[snp]))
            # noncis_component = noncis_effect_mat[i]
            if tx not in sim_tx2expr:
                sim_tx2expr[tx] = cis_component #+ noncis_component
            else:
                sim_tx2expr[tx] += (cis_component)# + noncis_component)
    # add in noncis effects
    for i in range(len(txlist)):
        tx = txlist[i]
        if tx in sim_tx2expr:
            sim_tx2expr[tx] += noncis_effect_mat[i]

    # now add in the noise and shared components
    for tx in sim_tx2expr:
        cis_plus_noncis = sim_tx2expr[tx]
        h2noise = 1.0 - total_h2 - args.h2noncis
        noise_component = np.random.normal(0, np.sqrt(h2noise), size=len(cis_plus_noncis))  # draw noise for each person
        sim_expr = cis_plus_noncis + noise_component
        # std_sim_expr = (sim_expr - np.mean(sim_expr)) / np.std(sim_expr)
        sim_tx2expr[tx] = sim_expr  # final iso expression is sum of components
        sim_gene2expr += sim_expr  # gene expression is sum of iso expressions
    return pd.DataFrame(sim_tx2expr), sim_gene2expr
    
    
def make_gene_info(gene, tx2info, txlist, window_start, window_end, gene_end, window):
    tx_info = tx2info[txlist[0]]
    gene_chr, gene_strand = tx_info[0], tx_info[-1]
    # gene_start, gene_end = window_start + window, window_end - window
    gene_start = window_start + window
    gene_info = [gene_chr, gene_start, gene_end, gene, gene, gene_strand]
    return gene_info


def sim_causal_snps(snp2geno, num_causal, h2, dropout_rate):
    if num_causal < 1:  # if num_causal not specified, pick randomly from 1-3
        r = random.random()
        if r < 0.3333:
            num_causal = 1
        elif r < 0.6666:
            num_causal = 2
        else:
            num_causal = 3
    # randomly pick num_causal snps in cis to be causal
    snp_set = list(snp2geno.keys())
    causal_snp_names = np.random.choice(snp_set, size=num_causal, replace=False)
    causal_snp_df = snp2geno[causal_snp_names]
    # set relative amount of variance explained by each SNP randomly
    relative_effects = np.array([random.random() for i in range(num_causal)])
    relative_effects /= sum(relative_effects)  # make sure they sum to 1
    # randomly set effects on some genes to 0
    dropout = False
    if random.random() < dropout_rate:
        dropout = True
        relative_effects *= 0.0
    # now scale relative effects to heritability parameters
    h2_rel_eff = h2 * relative_effects
    snp_h2_eff = {causal_snp_names[i]: h2_rel_eff[i] for i in range(num_causal)}
    return causal_snp_names, causal_snp_df, snp_h2_eff, dropout


def get_cis_snps(args, gene, gene2info, txlist, tx2info, tx2expr, meta_lines):
    window_start, window_end, gene_end, chrom = get_gene_window(txlist, tx2info, args.window)
    gene2info[gene] = make_gene_info(gene, tx2info, txlist, window_start, window_end, gene_end, args.window)
    window_str = str(chrom) + ':' + str(window_start) + '-' + str(window_end)
    bcf_proc = subprocess.Popen([args.bcftools, 'view', args.vcf, '-r', window_str], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bcf_out = StringIO(bcf_proc.communicate()[0].decode('utf-8'))
    gene_window_snps = pd.read_csv(bcf_out, sep='\t', header=meta_lines, index_col=2).T
    # read in SNPs in the window and clear out null SNPs and non-phenotyped individuals
    snp2info, snp2geno = gene_window_snps.iloc[:8], gene_window_snps.iloc[8:]
    snp2geno = preprocess_snp2geno(tx2expr, snp2geno)
    return snp2geno, gene2info


def simulation_pass(args, tx2info, tx2expr, gene_to_tx, meta_lines):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp = {}, {}, {}, {}
    for gene in gene_to_tx:
        txlist = gene_to_tx[gene]
        if args.num_iso != -1 and args.num_iso != len(txlist):
            continue
        # get window around gene and subset the vcf for that window using bcftools
        snp2geno, gene2info = get_cis_snps(args, gene, gene2info, txlist, tx2info, tx2expr, meta_lines)
        if len(list(snp2geno.keys())) < 5:
            del gene2info[gene]
            continue  # ignore genes with insufficient number of cis SNPs
        # randomly select causal snps
        causal_snp_names, causal_snp_df, snp_h2_eff, dropout = sim_causal_snps(snp2geno, args.num_causal, args.h2cis, args.dropout_rate)
        # simultate phenotypes
        if dropout:
            gene2causalsnp[gene] = [causal_snp_names, '0.0']
        else:
            gene2causalsnp[gene] = [causal_snp_names, str(args.h2cis)]
        sim_phenos, sim_gene_pheno = simulate_phenotypes(args, txlist, causal_snp_df, snp_h2_eff)
        for tx in sim_phenos:
            sim_tx2expr[tx] = sim_phenos[tx]
        sim_gene2expr[gene] = sim_gene_pheno
    fnull.close()
    sim_tx2expr = pd.DataFrame(sim_tx2expr).astype(str)
    sim_gene2expr = pd.DataFrame(sim_gene2expr).astype(str)
    return sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp


def write_results(tx2info, sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp, output):
    with(open(output + '.iso.bed', 'w')) as outfile:
        info_fields = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand']
        outfile.write('\t'.join(info_fields + list(sim_tx2expr.index)) + '\n')
        for tx in tx2info:
            if tx not in sim_tx2expr:
                continue
            tx_info_fields = [str(i) for i in tx2info[tx]]
            tx_info_fields = tx_info_fields[:3] + [tx] + tx_info_fields[3:]  # add tx name
            outfile.write('\t'.join(tx_info_fields) + '\t')
            outfile.write('\t'.join(list(sim_tx2expr[tx])) + '\n')
    with(open(output + '.gene.bed', 'w')) as outfile:
        info_fields = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand']
        outfile.write('\t'.join(info_fields + list(sim_tx2expr.index)) + '\n')
        for gene in gene2info:
            gene_info_fields = [str(i) for i in gene2info[gene]]
            outfile.write('\t'.join(gene_info_fields) + '\t')
            outfile.write('\t'.join(list(sim_gene2expr[gene])) + '\n')
    with(open(output + '.causal.txt', 'w')) as causal_outfile:
        causal_outfile.write('#Gene\tCausal SNP\th2\n')
        for gene in gene2causalsnp:
            causal_snps, h2 = gene2causalsnp[gene]
            causal_snps = ','.join(causal_snps)
            causal_outfile.write('\t'.join([gene, causal_snps, str(h2)]) + '\n')


def compress_and_index(output):
    # compress bed files with bgzip, then index with tabix
    iso_bed_name = output + '.iso.bed'
    gene_bed_name = output + '.gene.bed'
    iso_gz_outname = iso_bed_name + '.gz'
    gene_gz_outname = gene_bed_name + '.gz'
    with(open(iso_gz_outname, 'w')) as gz_outf:
        subprocess.Popen(['bgzip', '-c', iso_bed_name], stdout=gz_outf).wait()
    with(open(gene_gz_outname, 'w')) as gz_outf:
        subprocess.Popen(['bgzip', '-c', gene_bed_name], stdout=gz_outf).wait()
    subprocess.Popen(['tabix', '-p', 'bed', iso_gz_outname]).wait()
    subprocess.Popen(['tabix', '-p', 'bed', gene_gz_outname]).wait()
    subprocess.Popen(['rm', iso_bed_name, gene_bed_name]).wait()


if __name__ == "__main__":
    args = parseargs()
    corr_all = [args.max_corr, args.max_corr_env, args.min_corr, args.min_corr_env]
    if min(corr_all) < 0.0 or max(corr_all) > 1.0:
        sys.exit('Error: correlation arguments must be between 0.0 and 1.0, inclusive.')
    h2all = [args.h2cis, args.h2noncis]
    if min(h2all) < 0.0 or sum(h2all) > 1.0:
        sys.exit('Error: h2 args must be between 0.0 and 1.0 and sum to 1 or less.')
    
    tx2info, tx2expr, gene_to_tx = read_pheno_file(args.pheno)
    meta_lines = check_vcf(args.vcf, args.bcftools, tx2expr)
    sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp = simulation_pass(
        args, tx2info, tx2expr, gene_to_tx, meta_lines)
    write_results(tx2info, sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp, args.output)
    compress_and_index(args.output)
