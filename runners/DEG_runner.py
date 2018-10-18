#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import os
import pandas as pd

import numpy as np
# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()


import constants
import infra

from utils.r_runner import run_rscript


def run_DEG(method, conditions, data, genes, group):
    script = file(os.path.join(constants.REPO_DIR,"r","scripts","{}.r".format(method))).read()
    return run_rscript(script=script, data=data, genes=genes, conditions=conditions, group=group)


def prepare_input(gene_expression_file_name="ge.tsv"):
    ge_raw = infra.load_gene_expression_profile_by_genes(gene_expression_file_name=gene_expression_file_name)
    genes, conditions, data = infra.separate_headers(ge_raw)
    group = [1, 1, 1, 2, 2, 2]
    conditions = np.array(conditions)
    group = np.array(group)
    data = pd.DataFrame(data, index=genes, columns=conditions, dtype=np.int, )
    return conditions, data, genes, group



def main(method=constants.DEG_EDGER):
    conditions, data, genes, group = prepare_input()

    results = run_DEG(method, conditions, data, genes, group)
    if method==constants.DEG_EDGER:
        res = results["result"][["PValue", "FDR"]]
        res = res.rename(columns={"PValue": "pval", "FDR": "qval"})
    else:
        res = results["result"][["pvalue", "padj"]]
        res = res.rename(columns={"pvalue": "pval", "padj": "qval"})

    combined_DEG = pd.concat([data, res], axis=1, sort=False)
    combined_DEG = combined_DEG.sort_values("pval".format(method))
    combined_DEG.to_csv(os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method)), sep="\t", index_label="id")

if __name__ == "__main__":
    main()

