#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys
sys.path.insert(0, '../')

import os
import time
import shutil
# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()



import constants

from utils.r_runner import run_rscript
from utils.server import get_parameters
from utils.network import remove_subgraph_by_nodes
from utils.network import build_all_reports
from utils.server import get_score
import utils.server as server
import DEG_runner


import utils.go
import pandas as pd

from utils.network import output_modules

ALGO_NAME = "bionet"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)
NETWORK_NAME = "dip"


def run_bionet(deg_file_name, is_pval_score, network_file_name, fdr=0.05):
    script = file(os.path.join(ALGO_DIR, "{}.r".format(ALGO_NAME))).read()
    return run_rscript(script=script, output_vars = ["module_genes", "bg_genes"], network_file_name=network_file_name, deg_file_name=deg_file_name, fdr=fdr, is_pval_score=is_pval_score)


def init_specific_params(network_file_name=os.path.join(constants.NETWORKS_DIR, "{}.sif".format(NETWORK_NAME)), omitted_genes = [], ts=str(time.time())):
    return remove_subgraph_by_nodes(omitted_genes, network_file_name, ts=ts)


def get_module(network_file_name, score_file_name, is_pval_score, omitted_genes, ts=str(time.time()),fdr=0.05):
    network_file_name = init_specific_params(network_file_name=network_file_name, omitted_genes=omitted_genes, ts=ts)
    results = run_bionet(score_file_name, is_pval_score, network_file_name,fdr)
    module_genes = results["module_genes"]
    bg_genes = results["bg_genes"]
    file(os.path.join(constants.OUTPUT_DIR, "{}_module_genes_{}.txt".format(ALGO_NAME, ts)), "w+").write(
        "\n".join(module_genes))
    file(os.path.join(constants.OUTPUT_DIR, "{}_bg_genes_{}.txt".format(ALGO_NAME, ts)), "w+").write(
        "\n".join(bg_genes))

    sys.stdout.write("module gene size: {}. ratio: {}\n".format(module_genes.shape[0], module_genes.shape[0]/float(bg_genes.shape[0])))
    return list(module_genes), list(bg_genes)


def run_bionet_for_all_modules(fdr, network_file_name, score_file_name, is_pval_score):
    omitted_genes = []
    modules = []
    all_bg_genes = []
    for x in range(5):
        module_genes, bg_genes = get_module(network_file_name, score_file_name, is_pval_score, omitted_genes, str(x), fdr=fdr)
        omitted_genes += list(module_genes)
        modules.append(module_genes)
        all_bg_genes.append(bg_genes)
    return all_bg_genes, modules

def main(dataset_name=constants.DATASET_NAME, disease_name=None, expected_genes = None, fdr=0.05):
    global NETWORK_NAME
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    network_file_name, score_file_name, score_method, bg_genes = server.init_common_params(NETWORK_NAME)

    all_bg_genes, modules = run_bionet_for_all_modules(fdr, network_file_name, score_file_name, score_method!=constants.PREDEFINED_SCORE)

    output_base_dir = ""
    if constants.REPORTS:
        output_base_dir = build_all_reports(ALGO_NAME, modules, all_bg_genes, score_file_name, network_file_name, disease_name, expected_genes)

    output_file_name = os.path.join(constants.OUTPUT_DIR,
                                    "{}_client_output.txt".format(ALGO_NAME))
    output_modules(output_file_name, modules, score_file_name, output_base_dir)




if __name__ == "__main__":
    main()





