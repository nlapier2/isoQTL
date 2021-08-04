# Usage: "python tabulate_results.py my_results_dir/
import glob
import sys
    

def parse_results_from_file(fname):
    # parse one results file to get results for each metric for each method
    results = {}
    method = ''
    with(open(fname, 'r')) as infile:
        for line in infile:
            if line.startswith('Results'):
                method = line.split('Results for ')[1].split(':')[0]
                results[method] = {'Precision': 0.0, 'Recall': 0.0, 'F1 Score': 0.0}
            elif len(line) < 3:
                continue  # skip blank line
            else:
                metric, value = line.strip().split(': ')
                if value == 'N/A':
                    continue
                results[method][metric] += float(value)
    return results


def get_tabulated_results(all_fnames):
    # average over results for all files
    tabulated_results = {}
    for fname in all_fnames:
        this_res = parse_results_from_file(fname)
        for method in this_res:
            if method not in tabulated_results:
                tabulated_results[method] = {'Precision': 0.0, 'Recall': 0.0, 'F1 Score': 0.0}
            for metric in this_res[method]:  # first we sum the results
                tabulated_results[method][metric] += this_res[method][metric]
    # now we average the results
    for method in tabulated_results:
        for metric in tabulated_results[method]:
            tabulated_results[method][metric] /= len(all_fnames)
    return tabulated_results


res_dir = sys.argv[1]
if not res_dir.endswith('/'):
    res_dir += '/'
nom_fnames = glob.glob(res_dir + 'eval_nom_*.txt')
perm_fnames = glob.glob(res_dir + 'eval_perm_*.txt')
nom_results = get_tabulated_results(nom_fnames)
perm_results = get_tabulated_results(perm_fnames)
print('\nTabulated results for ' + res_dir)
for method in nom_results:
    print('\n' + method, 'nominal pass:')
    for metric in nom_results[method]:
        print(metric + ': ' + str(nom_results[method][metric]))
for method in perm_results:
    print('\n' + method, 'permutation pass:')
    for metric in perm_results[method]:
        print(metric + ': ' + str(perm_results[method][metric]))
print('\n\n')
