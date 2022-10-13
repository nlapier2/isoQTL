# intersect samples from gct, covariates file, and EUR file, and subset files accordingly
import argparse
import gzip


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(
        description='Subset gct and covariates files to only contain shared samples.')
    parser.add_argument('--gene_gct_tpm', required=True, 
        help='Gene-level TPM gct file for a specific tissue. Required.')
    parser.add_argument('--gene_gct_reads', required=True, 
        help='Gene-level read counts gct file for a specific tissue. Required.')
    parser.add_argument('--covariates', required=True, 
        help='Covariates file for a tissue. Required.')
    parser.add_argument('--eur', required=True, 
        help='File with names of samples from Europeans. Required.')
    parser.add_argument('--output', default='subset', help='Output base name.')
    pargs = parser.parse_args()
    return pargs


def read_gct_samples(gct_gene_fname):
    gct_samples, id2sample = {}, {}
    with(gzip.open(gct_gene_fname, 'r')) as infile:
        infile.readline() ; infile.readline()  # skip metadata lines
        samples = infile.readline().decode().strip().split('\t')[2:]
        for samp in samples:
            gct_samples[samp] = True
            idval = '-'.join(samp.split('-')[:2])  # shorten sample name to donor ID
            id2sample[idval] = samp
    return gct_samples, id2sample


def read_covariates_file(covar_fname):
    covar_ids = {}
    with(open(covar_fname, 'r')) as infile:
        idvals = infile.readline().strip().split('\t')[1:]
        for i in idvals:
            covar_ids[i] = True
    return covar_ids


def read_eur_file(eur_fname):
    eur_ids = {}
    with(open(eur_fname, 'r')) as infile:
        for line in infile:
            name = line.strip().split()[0]
            name = name.replace('.', '-')
            eur_ids[name] = True
    return eur_ids


def intersect(gct_samples, covar_ids, eur_ids):
    intersect_samples, intersect_ids = {}, {}
    for samp in gct_samples:
        idval = '-'.join(samp.split('-')[:2])  # shorten sample name to donor ID
        if idval in covar_ids and idval in eur_ids:
            intersect_samples[samp] = True
            intersect_ids[idval] = True
    return intersect_samples, intersect_ids


def write_samples_and_ids(basename, intersect_samples, intersect_ids):
    with(open(basename + '.ids.txt', 'w')) as outfile:
        for idval in intersect_ids:
            outfile.write(idval + '\n')
    with(open(basename + '.samples.txt', 'w')) as outfile:
        for samp in intersect_samples:
            outfile.write(samp + '\n')


def subset_gct(basename, gct_fname, intersect_samples, tpm_or_reads='tpm'):
    cols_to_keep = [1]
    with(gzip.open(gct_fname, 'r')) as infile:
        with(gzip.open(basename + '.' + tpm_or_reads + '.gct.gz', 'w')) as outfile:
            # write meta lines, then figure out which columns in the header to keep
            outfile.write(infile.readline())
            outfile.write(infile.readline())
            header = infile.readline().decode().strip().split('\t')
            cols_to_keep += [i for i in range(3, len(header)) if header[i] in intersect_samples]
            outline = '\t'.join([header[i] for i in cols_to_keep]) + '\n'
            outfile.write(outline.encode())
            
            # main pass through gct
            for line in infile:
                splits = line.decode().strip().split('\t')
                outline = '\t'.join([splits[i] for i in cols_to_keep]) + '\n'
                outfile.write(outline.encode())


def subset_covariates(basename, covar_file, intersect_samples, id2sample):
    cols_to_keep = [0]
    with(open(covar_file, 'r')) as infile:
        with(open(basename + '.covariates.txt', 'w')) as outfile:
            # figure out which columns in the header to keep
            header = infile.readline().strip().split('\t')
            header = [header[0]] + [id2sample[i] for i in header[1:]]  # convert donor ID to sample
            cols_to_keep += [i for i in range(1, len(header)) if header[i] in intersect_samples]
            outline = '\t'.join([header[i] for i in cols_to_keep]) + '\n'
            outfile.write(outline)
            
            # main pass through gct
            for line in infile:
                splits = line.strip().split('\t')
                outline = '\t'.join([splits[i] for i in cols_to_keep]) + '\n'
                outfile.write(outline)


if __name__ == '__main__':
    args = parseargs()
    gct_samp, id2samp = read_gct_samples(args.gene_gct_tpm)
    covar_id = read_covariates_file(args.covariates)
    eur_id = read_eur_file(args.eur)
    intersect_samp, intersect_id = intersect(gct_samp, covar_id, eur_id)
    write_samples_and_ids(args.output, intersect_samp, intersect_id)
    subset_gct(args.output, args.gene_gct_tpm, intersect_samp, tpm_or_reads='tpm')
    subset_gct(args.output, args.gene_gct_reads, intersect_samp, tpm_or_reads='reads')
    subset_covariates(args.output, args.covariates, intersect_samp, id2samp)
