import argparse


def parseargs():  # handle user arguments
    parser = argparse.ArgumentParser(description='isoQTL main script.')
    parser.add_argument('--gene_info', required=True, help='File with gene info from simulation script. Required.')
    parser.add_argument('--wilks', help='Wilks-Bartlett results file.')
    parser.add_argument('--ftest', help='F-test results file.')
    parser.add_argument('--qtltools_gene', help='Path to qtltools gene-level results.')
    parser.add_argument('--qtltools_grpbest', help='Path to qtltools grpbest gene-level results.')
    parser.add_argument('--qtltools_grppca1', help='Path to qtltools grppca1 gene-level results.')
    parser.add_argument('--qtltools_grpmean', help='Path to qtltools grpmean gene-level results.')
    parser.add_argument('--fisher_perm', help='Fisher permutation results file.')
    parser.add_argument('--min_perm', help='Min permutation results file.')
    parser.add_argument('--cauchy_perm', help='Cauchy permutation results file.')
    parser.add_argument('--fisher_qtltools_iso', help='Fisher on QTLtools results file.')
    parser.add_argument('--min_qtltools_iso', help='Min-pval on QTLtools results file.')
    parser.add_argument('--cauchy_qtltools_iso', help='Cauchy on QTLtools results file.')
    parser.add_argument('--threshold', default=0.05, type=float,
                        help='P-value threshold for calling an eGene. Default: 0.05.')
    parser.add_argument('--use_perm', action='store_true',
                        help='Use permuted p-value instead of nominal p-value, if available.')
    args = parser.parse_args()
    return args


def read_gene_info(gene_info):
    egenes, null_genes = {}, {}
    with(open(gene_info, 'r')) as infile:
        infile.readline()  # skip header line
        for line in infile:
            splits = line.strip().split('\t')
            gene, h2cis = splits[0], float(splits[2])
            if h2cis == 0.0:
                null_genes[gene] = True
            else:
                egenes[gene] = True
    return egenes, null_genes


def read_isoqtl(isoqtl_res, use_perm, threshold):
    isoqtl_egenes = {}
    with(open(isoqtl_res, 'r')) as infile:
        infile.readline()  # skip header
        for line in infile:
            if use_perm:
                gene = line.strip().split(' ')[0]
                isoqtl_egenes[gene] = True
            else:
                splits = line.strip().split('\t')
                gene = splits[0]
                pval = float(splits[3])
                if pval < threshold:
                    isoqtl_egenes[gene] = True
    return isoqtl_egenes


def read_qtltools_gene(qtltools_gene_res, use_perm, threshold):
    qtltools_gene_egenes = {}
    with(open(qtltools_gene_res, 'r')) as infile:
        for line in infile:
            splits = line.strip().split(' ')
            gene = splits[0]
            if use_perm:
                qtltools_gene_egenes[gene] = True
                # pval = float(splits[-1])
            else:
                pval = float(splits[11])
                if pval < threshold:
                    qtltools_gene_egenes[gene] = True
    return qtltools_gene_egenes


def read_combined_pvals(combined_gene_res, use_perm, threshold):
    combined_egenes = {}
    with(open(combined_gene_res, 'r')) as infile:
        infile.readline()  # skip header
        for line in infile:
            if use_perm:
                gene = line.strip().split(' ')[0]
                combined_egenes[gene] = True
            else:
                splits = line.strip().split('\t')
                gene = splits[0]
                pval = float(splits[1])
                if pval < threshold:
                    combined_egenes[gene] = True
    return combined_egenes


def calc_tp_fp_tn_fn(true_egenes, null_genes, pred_egenes):
    tp, fp, tn, fn = 0.0, 0.0, 0.0, 0.0
    for gene in pred_egenes:
        if gene in true_egenes:
            tp += 1.0
        else:
            fp += 1.0
    fn = len(true_egenes) - tp
    tn = len(null_genes) - fp
    return tp, fp, tn, fn


def calc_metrics_and_print(method_name, tp, fp, tn, fn):
    if tp + fp == 0.0:
        precision = 'N/A'
    else:
        precision = tp / (tp + fp)
    if tp + fn == 0.0:
        recall = 'N/A'
    else:
        recall = tp / (tp + fn)
    if precision == 'N/A' or recall == 'N/A' or precision + recall == 0.0:
        f1_score = 'N/A'
    else:
        f1_score = 2 * precision * recall / (precision + recall)

    print('Results for ' + str(method_name) + ':\n')
    print('Precision: ' + str(precision))
    print('Recall: ' + str(recall))
    print('F1 Score: ' + str(f1_score) + '\n')
    str_vals = ', '.join([str(i) for i in [tp, fp, tn, fn]])
    print('TP FP TN FN: ' + str_vals + '\n')


if __name__ == '__main__':
    args = parseargs()
    true_egenes, null_genes = read_gene_info(args.gene_info)

    if args.wilks is not None:
        wilks_egenes = read_isoqtl(args.wilks, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, wilks_egenes)
        calc_metrics_and_print('Wilks-Bartlett', tp, fp, tn, fn)

    if args.qtltools_gene is not None:
        qtltools_gene_egenes = read_qtltools_gene(args.qtltools_gene, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, qtltools_gene_egenes)
        calc_metrics_and_print('QTLtools_gene', tp, fp, tn, fn)

    if args.qtltools_grpbest is not None:
        qtltools_grpbest_egenes = read_qtltools_gene(args.qtltools_grpbest, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, qtltools_grpbest_egenes)
        calc_metrics_and_print('QTLtools_grpbest', tp, fp, tn, fn)

    if args.qtltools_grppca1 is not None:
        qtltools_grppca1_egenes = read_qtltools_gene(args.qtltools_grppca1, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, qtltools_grppca1_egenes)
        calc_metrics_and_print('QTLtools_grppca1', tp, fp, tn, fn)

    if args.qtltools_grpmean is not None:
        qtltools_grpmean_egenes = read_qtltools_gene(args.qtltools_grpmean, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, qtltools_grpmean_egenes)
        calc_metrics_and_print('QTLtools_grpmean', tp, fp, tn, fn)

    if args.ftest is not None:
        ftest_egenes = read_isoqtl(args.ftest, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, ftest_egenes)
        calc_metrics_and_print('F-test', tp, fp, tn, fn)

    if args.fisher_perm is not None:
        fisher_perm_egenes = read_isoqtl(args.fisher_perm, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, fisher_perm_egenes)
        calc_metrics_and_print('fisher_perm', tp, fp, tn, fn)

    if args.min_perm is not None:
        min_perm_egenes = read_isoqtl(args.min_perm, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, min_perm_egenes)
        calc_metrics_and_print('min_perm', tp, fp, tn, fn)

    if args.cauchy_perm is not None:
        cauchy_perm_egenes = read_isoqtl(args.cauchy_perm, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, cauchy_perm_egenes)
        calc_metrics_and_print('cauchy_perm', tp, fp, tn, fn)

    if args.fisher_qtltools_iso is not None:
        fisher_qtltools_iso_egenes = read_combined_pvals(args.fisher_qtltools_iso, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, fisher_qtltools_iso_egenes)
        calc_metrics_and_print('fisher_qtltools_iso', tp, fp, tn, fn)

    if args.min_qtltools_iso is not None:
        min_qtltools_iso_egenes = read_combined_pvals(args.min_qtltools_iso, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, min_qtltools_iso_egenes)
        calc_metrics_and_print('min_qtltools_iso', tp, fp, tn, fn)

    if args.cauchy_qtltools_iso is not None:
        cauchy_qtltools_iso_egenes = read_combined_pvals(args.cauchy_qtltools_iso, args.use_perm, args.threshold)
        tp, fp, tn, fn = calc_tp_fp_tn_fn(true_egenes, null_genes, cauchy_qtltools_iso_egenes)
        calc_metrics_and_print('cauchy_qtltools_iso', tp, fp, tn, fn)
