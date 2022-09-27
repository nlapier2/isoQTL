# convert GCT files to BED phenotype files. use GTF for gene/isoform annotations.
import argparse
import gzip
import numpy as np
import pandas as pd
import qtl.norm
import subprocess


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Convert GCT/TSV file to BED files.')
    parser.add_argument('--gct_tpm', required=True, help='GCT TPM file name. Required.')
    parser.add_argument('--gct_reads', required=True, help='GCT read counts file name. Required.')
    parser.add_argument('--gtf', required=True, help='GTF file with gene information. Required.')
    parser.add_argument('--output', required=True, help='Base name for output BED files.')
    parser.add_argument('--tpm_thresh', type=float, default=0.1, help='Minimum TPM threshold.')
    parser.add_argument('--read_thresh', type=float, default=6.0, 
        help='Minimum read count threshold.')
    parser.add_argument('--frac_thresh', type=float, default=0.2, 
        help='Fraction of samples that must meet tpm and read count thresholds.')
    pargs = parser.parse_args()
    return pargs


# read gtf file which tells us chromosome & location for an isoform
def read_gtf(gtf_fname):
    tx2info = {}
    with(gzip.open(gtf_fname, 'r')) as infile:
        for line in infile:
            line = line.decode()  # convert binary to string
            if line.startswith('#'):
                continue
            splits = line.split('\t')
            line_type = splits[2]
            if line_type != 'transcript':
                continue
            chrom, start, end, strand = splits[0], splits[3], splits[4], splits[6]
            start, end = int(start), int(end)
            tx = splits[8].split('"')[5]
            tx2info[tx] = [chrom, start, end, strand]
    return tx2info


def prepare_expression(counts_df, tpm_df, sample_frac_threshold=0.2,
                       count_threshold=6, tpm_threshold=0.1, mode='tmm'):
    """
    Original author of this function: Francois Aguet
    Adapted from: https://github.com/broadinstitute/gtex-pipeline
    Genes are filtered using the following expression thresholds:
      TPM >= tpm_threshold in >= sample_frac_threshold * samples
      read counts >= count_threshold in sample_frac_threshold * samples
    The filtered counts matrix is then normalized using:
      TMM (mode='tmm'; default) or
      quantile normalization (mode='qn')
    """

    # expression thresholds
    ns = tpm_df.shape[1]
    mask = (
        (np.sum(tpm_df >= tpm_threshold, axis=1) >= sample_frac_threshold * ns) &
        (np.sum(counts_df >= count_threshold, axis=1) >= sample_frac_threshold * ns)
    ).values

    # apply normalization
    if mode.lower() == 'tmm':
        tmm_counts_df = qtl.norm.edger_cpm(counts_df, normalized_lib_sizes=True)
        norm_df = qtl.norm.inverse_normal_transform(tmm_counts_df[mask])
    elif mode.lower() == 'qn':
        qn_df = qtl.norm.normalize_quantiles(tpm_df.loc[mask])
        norm_df = qtl.norm.inverse_normal_transform(qn_df)
    else:
        raise ValueError(f'Unsupported mode {mode}')

    return norm_df


# read gct files and perform normalization and filtering
def process_gct(gct_tpm, gct_reads, frac_thresh, read_thresh, tpm_thresh):
    tpm_df = pd.read_csv(gct_tpm, sep='\t', skiprows=2, index_col=0)
    reads_df = pd.read_csv(gct_reads, sep='\t', skiprows=2, index_col=0)
    gene_ids = tpm_df['gene_id']
    tpm_df.drop('gene_id', axis=1, inplace=True)
    reads_df.drop('gene_id', axis=1, inplace=True)
    norm_df = prepare_expression(reads_df, tpm_df, sample_frac_threshold=frac_thresh,
                                    count_threshold=read_thresh, tpm_threshold=tpm_thresh)
    return norm_df, gene_ids


# Add annotation information, ordered in BED format
def annotate_df(norm_df, gene_ids, annotations):
    iso_list = norm_df.index
    iso_list = [iso.split('.')[0] for iso in iso_list]
    chroms = [annotations[iso][0] for iso in iso_list if iso in annotations]
    starts = [annotations[iso][1] for iso in iso_list if iso in annotations]
    ends = [annotations[iso][2] for iso in iso_list if iso in annotations]
    strands = [annotations[iso][3] for iso in iso_list if iso in annotations]
    norm_df.insert(0, '#Chr', chroms)
    norm_df.insert(1, 'start', starts)
    norm_df.insert(2, 'end', ends)
    norm_df.insert(3, 'pid', norm_df.index, allow_duplicates=True)
    norm_df.insert(4, 'gid', gene_ids)
    norm_df.insert(5, 'strand', strands)
    return norm_df


# sort, write, compress, and index bed file for each chromosome
def write_beds(out_basename, bed_df):
    for i in range(22):
        target_chrom = str(i+1)
        chrom_df = bed_df[bed_df['#Chr'] == target_chrom]
        chrom_df = chrom_df.sort_values(by=['start'])
        outname = out_basename + '_chrom_' + str(i+1) + '.bed'
        chrom_df.to_csv(outname, sep='\t', index=False)
        subprocess.Popen(['bgzip', outname]).wait()
        subprocess.Popen(['tabix', '-p', 'bed', outname + '.gz']).wait()
    for i in ['X', 'Y', 'MT']:
        chrom_df = bed_df[bed_df['#Chr'] == i]
        chrom_df = chrom_df.sort_values(by=['start'])
        outname = out_basename + '_chrom_' + i + '.bed'
        chrom_df.to_csv(outname, sep='\t', index=False)
        subprocess.Popen(['bgzip', outname]).wait()
        subprocess.Popen(['tabix', '-p', 'bed', outname + '.gz']).wait()


if __name__ == '__main__':
    args = parseargs()
    gtf_annotations = read_gtf(args.gtf)
    normalized_df, gene_id_col = process_gct(args.gct_tpm, args.gct_reads, args.frac_thresh,
                                                args.read_thresh, args.tpm_thresh)
    annotated_df = annotate_df(normalized_df, gene_id_col, gtf_annotations)
    write_beds(args.output, annotated_df)
