import argparse
import numpy as np
from scipy.stats import combine_pvalues


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Script for combining iso-level p-values from QTLtools.')
    parser.add_argument('--qtltools', required=True, help='QTLtools permutation results file. Required.')
    parser.add_argument('--tx2gene', required=True, help='Tab-separated mapping transcripts (1st col) to their genes (2nd col).')
    parser.add_argument('--method', default='fisher',
        choices=['fisher', 'min', 'tippett', 'pearson', 'stouffer', 'cauchy'],
        help='Specify fisher or other method for combining pvals.')
    parser.add_argument('--output', default='combined_pvals_qtltools.txt',
        help='Output file name.')
    args = parser.parse_args()
    return args


def read_tx2gene(tx2gene_fname):
    tx2gene = {}
    with(open(tx2gene_fname, 'r')) as infile:
        for line in infile:
            splits = line.strip().split('\t')
            tx, gene = splits[0], splits[1]
            tx2gene[tx] = gene
    return tx2gene


def read_qtltools_res(qtltools_fname, tx2gene):
    gene2pvals = {}
    with(open(qtltools_fname, 'r')) as infile:
        for line in infile:
            splits = line.strip().split(' ')
            tx, pval = splits[0], splits[-1]
            if pval == 'NA':
                continue
            pval = float(pval)
            gene = tx2gene[tx]
            if gene not in gene2pvals:
                gene2pvals[gene] = [pval]
            else:
                gene2pvals[gene].append(pval)
    return gene2pvals


def cauchy_acat(pval_vec):
    weights = [1.0 / len(pval_vec) for i in pval_vec]  # equal weights for now
    stat = sum([weights[i] * np.tan((0.5 - pval_vec[i]) * np.pi) for i in range(len(pval_vec))])
    pval = 0.5 - (np.arctan(stat / sum(weights))) / np.pi
    return pval


def combine_and_write(gene2pvals, method, out_fname):
    with(open(out_fname, 'w')) as outfile:
        for gene in gene2pvals:
            if method == 'min':
                combined_pval = min(gene2pvals[gene])
            elif method == 'cauchy':
                combined_pval = cauchy_acat(gene2pvals[gene])
            else:
                combined_pval = combine_pvalues(gene2pvals[gene], method=method, weights=None)[1]
            outfile.write(gene + '\t' + str(combined_pval) + '\n')


if __name__ == "__main__":
    args = parseargs()
    tx2gene = read_tx2gene(args.tx2gene)
    gene2pvals = read_qtltools_res(args.qtltools, tx2gene)
    combine_and_write(gene2pvals, args.method, args.output)
