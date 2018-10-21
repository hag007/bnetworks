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

import DEG_runner

# import utils.go


ALGO_NAME = "bionet"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)
NETWORK_NAME = "dip"

def run_bionet(deg_file_name, network_file_name):
    script = file(os.path.join(ALGO_DIR, "{}.r".format(ALGO_NAME))).read()
    return run_rscript(script=script, output_vars = ["module_genes", "bg_genes"], network_file_name=network_file_name, deg_file_name=deg_file_name)


def prepare_input(method=constants.DEG_EDGER):
    deg_file_name = str(os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method)))
    if not os.path.exists(deg_file_name):
        DEG_runner.main(method=method)

    return deg_file_name




if __name__ == "__main__":
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params

    deg_file_name = prepare_input(method=constants.DEG_EDGER)

    network_file_name = str(os.path.join(constants.NETWORKS_DIR, "{}.sif".format(NETWORK_NAME)))
    results = run_bionet(deg_file_name, network_file_name)
    module_genes = results["module_genes"]
    bg_genes = results["bg_genes"]

    file(os.path.join(constants.OUTPUT_DIR, "bionet_module_genes.txt"),"w+").write("\n".join(module_genes))
    file(os.path.join(constants.OUTPUT_DIR, "bionet_bg_genes.txt"), "w+").write("\n".join(bg_genes))

  #  utils.go.check_group_enrichment(module_genes,bg_genes)

    sys.stdout.write(os.path.join(constants.OUTPUT_DIR, "bionet_module_genes.txt"))





