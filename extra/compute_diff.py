# subset an isoform-aware fdr results file with only the genes not detected by qtltools

import argparse
import pandas as pd

def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Compute isoform-aware-exclusive genes.')
    parser.add_argument('--results', required=True, 
        help='Results from a previous isoform-aware run. Required.')
    parser.add_argument('--qtltools_results', required=True, 
        help='Results from a previous qtltools run. Required.')
    parser.add_argument('--output', default='eqtl_results', help='Output file name.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parseargs()
    res_ia = pd.read_csv(args.results, sep=' ', header=None)
    res_qt = pd.read_csv(args.qtltools_results, sep=' ', header=None)
    tmp = [i.split('.')[0] for i in res_ia[0]]
    uniq_rows = [i for i in range(len(tmp)) if tmp[i] not in list(res_qt[0])]
    subset = res_ia.iloc[uniq_rows]
    subset.to_csv(args.output, sep=' ', header=False, index=False)
