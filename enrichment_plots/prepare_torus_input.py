# given qtltools nominal results file and annotations file, prepare torus input files
import argparse
import gzip
import subprocess
from scipy.stats import t


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Prepare torus input files.')
    parser.add_argument('--qtltools', required=True, help='QTLtools results file. Required.')
    parser.add_argument('--annotations', required=True, help='Annotations file. Required.')
    parser.add_argument('--output', default='torus_', help='Base name for output files.')
    parser.add_argument('--threshold', type=float, default=0.00000005, 
        help='P-value threshold for including SNPs.')
    parser.add_argument('--gene_list', default='NONE', help='Read genes to retain SNPs from.')
    parser.add_argument('--map', action='store_true', help='Use to generate map files.')
    pargs = parser.parse_args()
    return pargs


# from annotation file, derive map of SNP location to full SNP name
def read_annotation_file(annotations, output):
    loc2name = {}
    with(gzip.open(annotations, 'r')) as infile:
        infile.readline()  # skip header
        for line in infile:
            snp_name = line.decode().split('\t')[0]
            snp_chr, snp_loc = snp_name.split('_')[:2]
            loc2name[snp_chr + '_' + snp_loc] = snp_name
    return loc2name


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


# convert qtltools format to matrix eQTL format, which torus can read
# also convert snp names and store gene location info for gene map file
def make_assoc_file_matrixeqtl(qtltools, output, threshold, loc2name, gene_list):
    gene2loc = {}
    pcol = 11
    betacol = 13
    snpcol = 7
    with(open(qtltools, 'r')) as infile:
        with(gzip.open(output + 'assoc.txt.gz', 'w')) as outfile:
            outfile.write('SNP\tgene\tbeta\tt-stat\tp-value\n'.encode())  # write header
            for line in infile:
                splits = line.strip().split(' ')
                gene, chr, start, end = splits[:4]
                if len(gene_list) > 0 and gene not in gene_list:
                    continue
                gene2loc[gene] = [chr, start, end]
                snp, snpchr, snploc = splits[snpcol:snpcol+3]
                snp_name = loc2name['chr' + snpchr + '_' + snploc]
                splits[snpcol] = snp_name
                pval = float(splits[pcol])
                if pval > threshold:
                    continue
                stat = t.ppf(pval/2, 492)
                fields = [snp_name, gene, splits[betacol], str(stat), splits[pcol]]
                outline = '\t'.join(fields) + '\n'
                outfile.write(outline.encode())
    return gene2loc


# write snp map for torus
def write_snp_map(output, loc2name):
    with(gzip.open(output + 'snp_map.txt.gz', 'w')) as outfile:
        for loc in loc2name:
            name = loc2name[loc]
            chr, pos = loc.split('_')
            chr = chr[3:]  # strip off 'chr'
            outline = name + '\t' + chr + '\t' + pos + '\n'
            outfile.write(outline.encode())
            
            
# write gene map for torus
def write_gene_map(output, gene2loc):
    with(gzip.open(output + 'gene_map.txt.gz', 'w')) as outfile:
        for gene in gene2loc:
            chr, start, end = gene2loc[gene]
            outline = gene + '\t' + chr + '\t' + start + '\t' + end + '\n'
            outfile.write(outline.encode())


if __name__ == '__main__':
    args = parseargs()
    loc_to_snp_name = read_annotation_file(args.annotations, args.output)
    #subprocess.Popen(['cp', args.annotations, args.output + 'annotations.txt.gz']).wait()
    list_of_genes = read_gene_list(args.gene_list)
    gene2location = make_assoc_file_matrixeqtl(args.qtltools, args.output, 
                        args.threshold, loc_to_snp_name, list_of_genes)
    if args.map:
        write_snp_map(args.output, loc_to_snp_name)
        write_gene_map(args.output, gene2location)
