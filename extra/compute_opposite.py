# compute genes whose isoform nominal associations have opposite signs

import pickle
import argparse


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Find genes with opposite-sign isoforms.')
    parser.add_argument('--tld', required=True, 
        help='Top level directory name containing pickled isoform regression results. Required.')
    parser.add_argument('--ignore_singleton', action='store_true', 
        help='Use to ignore single-isoform genes.')
    parser.add_argument('--req_signif', action='store_true', 
        help='Use to require that opposite-sign isoforms have significant p-values.')
    args = parser.parse_args()
    return args
  
  
def get_combined_dict(tld, method_name):
    all_res = {}
    for i in range(1, 23):
        dirname = tld + 'geuvadis_chr' + str(i) + '/'
        fname = dirname + 'results_' + method_name + '_chr_' + str(i)
        with(open(fname, 'rb')) as infile:
            chr_res = pickle.load(infile)
        for j in chr_res:
            all_res[j] = chr_res[j]
    return all_res


def get_opposite(res, req_signif = False, ignore_singleton = False):
    total = 0
    opposite = 0
    for gene in res:
        if not hasattr(res[gene][0], "__len__"):  # integer not array, so single-isoform gene
            if not ignore_singleton:  # optionally ignore 1-isoform genes
                total += 1
            continue
        betas, pvals = res[gene][0], res[gene][2]
        total += 1
        # search for isoforms with different signs, optionally requiring significance
        pos = False
        neg = False
        for i in range(len(betas)):
            if float(betas[i]) > 0 and (not req_signif or float(pvals[i]) < 0.05):
                pos = True
            elif float(betas[i]) < 0 and (not req_signif or float(pvals[i]) < 0.05):
                neg = True
        if pos and neg:
            opposite += 1
    return total, opposite


if __name__ == '__main__':
    args = parseargs()
    if not args.tld.endswith('/'):
        args.tld += '/'
    all_res_ftest = get_combined_dict(args.tld, 'ftest')
    all_res_diff = get_combined_dict(args.tld, 'diff_ftest_qtlgene')
    #all_res_qtltools = get_combined_dict(args.tld, 'qtltools_gene')
    ftest_tot, ftest_opp = get_opposite(all_res_ftest, 
        req_signif = args.req_signif, ignore_singleton = args.ignore_singleton)
    diff_tot, diff_opp = get_opposite(all_res_diff, 
        req_signif = args.req_signif, ignore_singleton = args.ignore_singleton)
    #qtl_tot, qtl_opp = get_opposite(all_res_qtltools, 
    #    req_signif = args.req_signif, ignore_singleton = args.ignore_singleton)
    print('F-test opposite / total genes: ' +  str(ftest_opp) + ' / ' + str(ftest_tot))
    # print('QTLtools-sum opposite / total genes: ' +  str(qtl_opp) + ' / ' + str(qtl_tot))
    print('F-test-only opposite / total genes: ' +  str(diff_opp) + ' / ' + str(diff_tot))
