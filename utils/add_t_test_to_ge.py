
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from scipy.stats import ttest_ind

import numpy as np
import os
import infra
import constants
from numpy import log10


def calc_ttest(dataset="TNFa", gene_expression_file_name="ge.tsv"):
    constants.update_dirs(DATASET_NAME="TNFa")
    h_rows, h_cols, ge_dataset = infra.separate_headers(infra.load_gene_expression_profile_by_genes(gene_expression_file_name=gene_expression_file_name))

    pvals = []
    rows_to_delete = []
    pval_dict = {}

    for i, cur in enumerate(list(h_rows)):
        pval_dict[cur] = ttest_ind(ge_dataset[i,:3], ge_dataset[i,3:]).pvalue
        if np.isnan(pval_dict[cur]):
            print "case: {}, wt: {}".format(ge_dataset[i,:3], ge_dataset[i,3:])
            rows_to_delete.append(i)
        else:
            pvals.append(pval_dict[cur])
    ind = np.ones((len(h_rows),), bool)
    ind[rows_to_delete] = False
    h_rows = h_rows[ind]
    ge_dataset = ge_dataset[ind, :]
    # print pvals
    qvals = fdrcorrection0(pvals, alpha=0.05, method='indep', is_sorted=False)[1]
    qscores = []
    for i, cur in enumerate(h_rows):
        qscores.append(-log10(qvals[i]))

    output_h_cols = ["ensemblids"] +  list(h_cols) +[ "pval", "qval", "qscore"]


    output_matrix = np.c_[h_rows, ge_dataset, pvals,qvals,qscores]
    output_matrix = np.r_[np.reshape(output_h_cols, (1, len(output_h_cols))), output_matrix]

    lines = []
    for i, cur in enumerate(output_matrix):
        lines.append("\t".join(cur))

    file(os.path.join(constants.OUTPUT_DIR, "gescore.tsv"), "w+").write("\n".join(lines))


calc_ttest()



