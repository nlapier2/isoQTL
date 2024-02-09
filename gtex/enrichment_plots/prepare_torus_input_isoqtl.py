# given IsoQTL results file and annotations file, prepare torus input file
import argparse
import gzip
from scipy.stats import t


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Prepare torus input file.')
    parser.add_argument('--results', required=True, help='Assocation results file. Required.')
    parser.add_argument('--bim', required=True, help='Plink bim file for SNP info. Required.')
    parser.add_argument('--annotations', required=True, help='Annotations file. Required.')
    parser.add_argument('--output', default='torus_', help='Base name for output files.')
    parser.add_argument('--threshold', type=float, default=0.00000005, 
        help='P-value threshold for including SNPs.')
    parser.add_argument('--perm', action='store_true', help='Use if permutation pass was done.')
    parser.add_argument('--qvalue', action='store_true', help='Use if qvalue is included.')
    parser.add_argument('--gene_list', default='NONE', help='Read genes to retain SNPs from.')
    pargs = parser.parse_args()
    return pargs


# read bim file to map SNP location to rsid
def read_bim_file(bimfile):
    loc2name = {}
    with(open(bimfile, 'r')) as infile:
        for line in infile:
            splits = line.strip().split('\t')
            chrom, rsid, pos = splits[0], splits[1], splits[3]
            loc = 'chr' + chrom + '_' + pos
            loc2name[loc] = rsid
    return loc2name


# from annotation file, derive map of SNP location to full SNP name
def read_annotation_file(annotations, output, loc2name):
    name_map = {}
    with(gzip.open(annotations, 'r')) as infile:
        infile.readline()  # skip header
        for line in infile:
            snp_name = line.decode().split('\t')[0]
            snp_chr, snp_pos = snp_name.split('_')[:2]
            loc = snp_chr + '_' + snp_pos
            if loc not in loc2name:
                continue
            rsid = loc2name[loc]
            name_map[rsid] = snp_name
    return name_map


# read list of genes that passed qvalue fdr control
def read_gene_list(gene_list_file):
    gene_list = {}
    if gene_list_file == 'NONE':
        return gene_list
    with(open(gene_list_file, 'r')) as infile:
        for line in infile:
            gene = line.split(' ')[0]
            gene_list[gene] = True
    return gene_list


# convert IsoQTL format to matrix eQTL format, which torus can read
def make_assoc_file_matrixeqtl(resfile, output, thresh, perm, qvalue, name_map, gene_list):
    with(open(resfile, 'r')) as infile:
        infile.readline()  # skip header
        with(gzip.open(output + 'assoc.txt.gz', 'w')) as outfile:
            outfile.write('SNP\tgene\tbeta\tt-stat\tp-value\n'.encode())  # write header
            for line in infile:
                if qvalue:
                    gene, rsid, stat, nom, dirperm, pval, qval = line.strip().split(' ')
                elif perm:
                    gene, rsid, stat, nom, dirperm, pval = line.strip().split('\t')
                else:
                    gene, rsid, stat, pval = line.strip().split('\t')
                if len(gene_list) > 0 and gene not in gene_list:
                    continue
                if float(pval) > thresh:
                    continue
                stat = str(t.ppf(float(pval)/2, 492))
                beta = stat
                snp_name = name_map[rsid]
                outline = '\t'.join([snp_name, gene, beta, stat, pval]) + '\n'
                outfile.write(outline.encode())


if __name__ == '__main__':
    args = parseargs()
    loc_to_name = read_bim_file(args.bim)
    name_dict = read_annotation_file(args.annotations, args.output, loc_to_name)
    list_of_genes = read_gene_list(args.gene_list)
    make_assoc_file_matrixeqtl(args.results, args.output, args.threshold, 
        args.perm, args.qvalue, name_dict, list_of_genes)
