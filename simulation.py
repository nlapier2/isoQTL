import argparse
import gzip
import os
import random
import subprocess
import sys
import pandas as pd
import numpy as np
from io import StringIO


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='isoQTL main script.')
    parser.add_argument('--vcf', required=True, help='VCF or BCF file with genotypes. Required.')
    parser.add_argument('--pheno', required=True, help='BED file or tsv file with transcript expression levels.')
    parser.add_argument('--bcftools', default='bcftools', help='Path to bcftools executable ("bcftools" by default).')
    parser.add_argument('--dropout_rate', default=0.25, type=float, help='Percent of genes to randomly set to zero heritability.')
    parser.add_argument('--h2g', default=0.1, type=float, help='Gene-level heritability.')
    parser.add_argument('--h2i', default=0.05, type=float, help='Isoform-level heritability.')
    parser.add_argument('--h2shared', default=0.25, type=float, help='Heritability shared across the gene, but not due to QTL.')
    parser.add_argument('--num_causal', default=-1, type=int,
        help='Specify a fixed number of causal SNPs per gene.')
    parser.add_argument('--num_iso', default=-1, type=int, help='Use to include genes only with a certain number of isoforms.')
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
    # center window around first_start and last_end
    window_start, window_end = first_start - window, last_end + window
    return int(window_start), int(window_end), chrom


def simulate_phenotypes(tx2expr, txlist, causal_snp_df, snp_h2g_eff, snp_h2i_eff, args_h2shared):
    sim_tx2expr, sim_gene2expr = {}, np.zeros(len(causal_snp_df))
    total_h2g, total_h2i = 0.0, 0.0
    # first sim genetic component: the effects of the SNPs on genes & isoforms
    for snp in causal_snp_df:
        snp_h2g = snp_h2g_eff[snp]
        snp_h2i = snp_h2i_eff[snp]
        total_h2g, total_h2i = total_h2g + snp_h2g, total_h2i + snp_h2i
        min_eff = min(snp_h2g, snp_h2i) / 2
        gene_eff = np.random.normal(0, np.sqrt(snp_h2g))
        for tx in txlist:
            iso_eff = np.random.normal(0, np.sqrt(snp_h2i))
            while abs(gene_eff + iso_eff) < min_eff:
                # don't let effects cancel each other out
                iso_eff = np.random.normal(0, np.sqrt(snp_h2i))
            genetic_component = causal_snp_df[snp] * (gene_eff + iso_eff)
            if tx not in sim_tx2expr:
                sim_tx2expr[tx] = genetic_component
            else:
                sim_tx2expr[tx] += genetic_component
    # now add in the noise and shared components
    for tx in sim_tx2expr:
        genetic_component = sim_tx2expr[tx]
        if total_h2g + total_h2i == 0.0:
            h2noise = 1.0 - args_h2shared
            this_h2shared = args_h2shared
        else:  # pick noise based on realized genetic variance to preserve h2
            h2noise = np.var(genetic_component) * ((1.0 - args_h2shared) / (total_h2g + total_h2i) - 1.0)
            this_h2shared = np.var(genetic_component) * (args_h2shared / (total_h2g + total_h2i))
        shared_component = np.random.normal(0, np.sqrt(this_h2shared)) * np.ones(len(genetic_component))  # draw once, then shared
        noise_component = np.random.normal(0, np.sqrt(h2noise), size = len(genetic_component))  # draw separately for each person
        sim_expr = genetic_component + noise_component + shared_component
        # std_sim_expr = (sim_expr - np.mean(sim_expr)) / np.std(sim_expr)
        sim_tx2expr[tx] = sim_expr  # final iso expression is sum of components
        sim_gene2expr += sim_expr  # gene expression is sum of iso expressions
    return pd.DataFrame(sim_tx2expr), sim_gene2expr
    
    
def make_gene_info(gene, tx2info, txlist, window_start, window_end, window):
    tx_info = tx2info[txlist[0]]
    gene_chr, gene_strand = tx_info[0], tx_info[-1]
    gene_start, gene_end = window_start + window, window_end - window
    gene_info = [gene_chr, gene_start, gene_end, gene, gene, gene_strand]
    return gene_info


def sim_causal_snps(snp2geno, num_causal, h2g, h2i, dropout_rate):
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
    h2g_rel_eff = h2g * relative_effects
    h2i_rel_eff = h2i * relative_effects
    snp_h2g_eff = {causal_snp_names[i]: h2g_rel_eff[i] for i in range(num_causal)}
    snp_h2i_eff = {causal_snp_names[i]: h2i_rel_eff[i] for i in range(num_causal)}
    return causal_snp_names, causal_snp_df, snp_h2g_eff, snp_h2i_eff, dropout


def get_cis_snps(args, gene, gene2info, txlist, tx2info, tx2expr, meta_lines):
    window_start, window_end, chrom = get_gene_window(txlist, tx2info, args.window)
    gene2info[gene] = make_gene_info(gene, tx2info, txlist, window_start, window_end, args.window)
    window_str = str(chrom) + ':' + str(window_start) + '-' + str(window_end)
    bcf_proc = subprocess.Popen([args.bcftools, 'view', args.vcf, '-r', window_str], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bcf_out = StringIO(bcf_proc.communicate()[0].decode('utf-8'))
    gene_window_snps = pd.read_csv(bcf_out, sep='\t', header=meta_lines+4, index_col=2).T
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
        # randomly select causal snps
        causal_snp_names, causal_snp_df, snp_h2g_eff, snp_h2i_eff, dropout = sim_causal_snps(snp2geno, args.num_causal, args.h2g, args.h2i, args.dropout_rate)
        # simultate phenotypes
        if dropout:
            gene2causalsnp[gene] = [causal_snp_names, '0.0', '0.0']
        else:
            gene2causalsnp[gene] = [causal_snp_names, str(args.h2g), str(args.h2i)]
        sim_phenos, sim_gene_pheno = simulate_phenotypes(tx2expr, txlist, causal_snp_df, snp_h2g_eff, snp_h2i_eff, args.h2shared)
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
        causal_outfile.write('#Gene\tCausal SNP\th2g\th2i\n')
        for gene in gene2causalsnp:
            causal_snps, h2g, h2i = gene2causalsnp[gene]
            causal_snps = ','.join(causal_snps)
            h2g, h2i = str(h2g), str(h2i)
            causal_outfile.write('\t'.join([gene, causal_snps, h2g, h2i]) + '\n')


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
    h2all = [args.h2g + args.h2i + args.h2shared]
    if min(h2all) < 0.0 or sum(h2all) > 1.0:
        sys.exit('Error: h2 parameters must be positive and sum to less than 1.0.')
    tx2info, tx2expr, gene_to_tx = read_pheno_file(args.pheno)
    meta_lines = check_vcf(args.vcf, tx2expr)
    sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp = simulation_pass(
        args, tx2info, tx2expr, gene_to_tx, meta_lines)
    write_results(tx2info, sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp, args.output)
    compress_and_index(args.output)
