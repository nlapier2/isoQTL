# convert GTEx PLINK files to VCF, subset, convert to dosage format, and rename columns
import argparse
import subprocess


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Convert GTEx plink files to VCFs.')
    parser.add_argument('--bfile', required=True, help='Plink bfile name. Required.')
    parser.add_argument('--ids_file', required=True, 
        help='File with donor ID names to keep. Required.')
    parser.add_argument('--samples_file', required=True, 
        help='File with sample names to keep. Required.')
    parser.add_argument('--output', default='gtex_genos', help='Output base name.')
    pargs = parser.parse_args()
    return pargs


def get_id2sample(samples_file):
    id2sample = {}
    with(open(samples_file, 'r')) as infile:
        for line in infile:
            sample = line.strip()
            idval = '-'.join(sample.split('-')[:2])  # shorten sample name to donor ID
            id2sample[idval] = sample
    return id2sample


def reformat_vcf(in_vcf, outname, id2sample):
    with(open(in_vcf, 'r')) as infile:
        with(open(outname, 'w')) as outfile:
            for line in infile:
                if line.startswith('##') and not line.startswith('##FORMAT'):
                    outfile.write(line)
                elif line.startswith('##FORMAT'):  # replace GT format with DS format
                    outfile.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Imputed genotype dosage">\n')
                elif line.startswith('#'):
                    splits = line.strip().split('\t')
                    for i in range(9, len(splits)):
                        iid = splits[i].split('_')[0]
                        sample = id2sample[iid]
                        splits[i] = sample
                    outfile.write('\t'.join(splits) + '\n')
                else:
                    splits = line.strip().split('\t')
                    splits[8] = 'DS'  # converting from GT to DS (genotype to dosage)
                    for i in range(9, len(splits)):
                        geno = int(splits[i][0]) + int(splits[i][2])  # e.g. 0/1 --> 1, etc
                        dosage = float(geno)
                        splits[i] = str(dosage)
                    outfile.write('\t'.join(splits) + '\n')


if __name__ == '__main__':
    args = parseargs()
    unf_name = args.output + '_unformatted'
    formatted_name = args.output + '_genos.vcf'
    subprocess.Popen(['plink', '--bfile', args.bfile, '--recode', 'vcf', 
        '--keep-fam', args.ids_file, '--out', unf_name]).wait()
    id2samp = get_id2sample(args.samples_file)
    reformat_vcf(unf_name + '.vcf', formatted_name, id2samp)
    subprocess.Popen(['bgzip', formatted_name]).wait()
    subprocess.Popen(['tabix', '-p', 'vcf', formatted_name + '.gz']).wait()
    subprocess.Popen(['rm', unf_name + '.vcf']).wait()
