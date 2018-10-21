#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys
sys.path.insert(0, '../')


import os
import numpy as np
import pandas as pd
import subprocess

# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import infra

import DEG_runner
import utils.go

from utils.scripts import format_script
from utils.network import get_bg_genes

from utils.server import get_parameters


ALGO_NAME = "jactivemodules"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)
NETWORK_NAME = "dip"

if __name__ == "__main__":
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params

    score_method = constants.DEG_EDGER
    deg_file_name = os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(score_method))
    if not os.path.exists(deg_file_name):
        DEG_runner.main(method=score_method)

    network_file_name = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(NETWORK_NAME))
    bg_genes = get_bg_genes()

    for cur_search_method in ["sa"]:

        results_file_name = "{}/{}_{}_results.txt".format(constants.OUTPUT_DIR, ALGO_NAME, cur_search_method)
        format_script(os.path.join(constants.SH_DIR, "run_{}.sh".format(ALGO_NAME)), BASE_FOLDER=constants.BASE_PROFILE, DATASET_DIR=constants.DATASET_DIR,
                                                                          ALGO_DIR=ALGO_DIR, NETWORK_NAME=NETWORK_NAME, RANK_METHOD=score_method,
                                                                          IS_GREEDY=str(cur_search_method=="greedy"), OUTPUT_FILE=results_file_name)

        subprocess.Popen("bash {}/run_{}.sh".format(constants.SH_DIR, ALGO_NAME), shell=True,
                               stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()

        results = file(results_file_name.format(constants.OUTPUT_DIR, ALGO_NAME)).readlines()

        module_genes = [y for x in results for y in x.split()[:-1]]
        file(os.path.join(constants.OUTPUT_DIR,"{}_{}_module_genes.txt".format(ALGO_NAME, cur_search_method)), "w+").write("\n".join(module_genes))

        # utils.go.check_group_enrichment(module_genes, bg_genes)

        sys.stdout.write(os.path.join(constants.OUTPUT_DIR,"{}_{}_module_genes.txt".format(ALGO_NAME, cur_search_method)))









