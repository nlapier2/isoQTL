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
    parser.add_argument('--h2g', default=0.05, type=float, help='Gene-level heritability.')
    parser.add_argument('--h2i', default=0.02, type=float, help='Isoform-level heritability.')
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


def simulate_phenotypes(tx2expr, txlist, causal_snp, h2g, h2i):
    sim_tx2expr, sim_gene2expr = {}, np.zeros(len(causal_snp))
    h2e = 1.0 - h2g - h2i  # residual variance not explained by gene or iso effects
    gene_eff = np.random.normal(0, h2g)
    for tx in txlist:
        iso_eff = np.random.normal(0, h2i)
        genetic_component = causal_snp * (gene_eff + iso_eff)
        noise_component = np.random.normal(0, h2e, size = len(genetic_component))
        sim_tx2expr[tx] = genetic_component + noise_component
        sim_gene2expr += sim_tx2expr[tx]
    return pd.DataFrame(sim_tx2expr), sim_gene2expr
    
    
def make_gene_info(gene, tx2info, txlist, window_start, window_end, window):
    tx_info = tx2info[txlist[0]]
    gene_chr, gene_strand = tx_info[0], tx_info[-1]
    gene_start, gene_end = window_start + window, window_end - window
    gene_info = [gene_chr, gene_start, gene_end, gene, gene, gene_strand]
    return gene_info


def simulation_pass(vcf, tx2info, tx2expr, gene_to_tx, meta_lines, dropout_rate, h2g, h2i, window, bcftools):
    fnull = open(os.devnull, 'w')  # used to suppress some annoying bcftools warnings
    sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp = {}, {}, {}, {}
    for gene in gene_to_tx:
        txlist = gene_to_tx[gene]
        # get window around gene and subset the vcf for that window using bcftools
        window_start, window_end, chrom = get_gene_window(txlist, tx2info, window)
        gene2info[gene] = make_gene_info(gene, tx2info, txlist, window_start, window_end, window)
        window_str = str(chrom) + ':' + str(window_start) + '-' + str(window_end)
        bcf_proc = subprocess.Popen([bcftools, 'view', vcf, '-r', window_str],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bcf_out = StringIO(bcf_proc.communicate()[0].decode('utf-8'))
        gene_window_snps = pd.read_csv(bcf_out, sep='\t', header=meta_lines+4, index_col=2).T
        # read in SNPs in the window and clear out null SNPs and non-phenotyped individuals
        snp2info, snp2geno = gene_window_snps.iloc[:8], gene_window_snps.iloc[8:]
        snp2geno = preprocess_snp2geno(tx2expr, snp2geno)
        
        # randomly select a causal snp
        snp_set = snp2geno.keys()
        causal_snp = snp2geno[snp_set[int(random.random() * len(snp_set))]]
        #print(causal_snp) ; print(tx2expr[txlist[0]]) ; sys.exit()
        # simultate phenotypes
        if random.random() < dropout_rate:
            # randomly select some genes to have no causal effect
            gene2causalsnp[gene] = [causal_snp.name, '0.0', '0.0']
            sim_phenos, sim_gene_pheno = simulate_phenotypes(tx2expr, txlist, causal_snp, 0.0, 0.0)
        else:
            gene2causalsnp[gene] = [causal_snp.name, str(h2g), str(h2i)]
            sim_phenos, sim_gene_pheno = simulate_phenotypes(tx2expr, txlist, causal_snp, h2g, h2i)
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
            causal_snp, h2g, h2i = gene2causalsnp[gene]
            h2g, h2i = str(h2g), str(h2i)
            causal_outfile.write('\t'.join([gene, causal_snp, h2g, h2i]) + '\n')


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
    tx2info, tx2expr, gene_to_tx = read_pheno_file(args.pheno)
    meta_lines = check_vcf(args.vcf, tx2expr)
    sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp = simulation_pass(
        args.vcf, tx2info, tx2expr, gene_to_tx, meta_lines, args.dropout_rate, args.h2g, args.h2i, args.window, args.bcftools)
    write_results(tx2info, sim_tx2expr, sim_gene2expr, gene2info, gene2causalsnp, args.output)
    compress_and_index(args.output)
