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

from utils.r_runner import run_rscript
from utils.network import get_bg_genes

import infra

import utils.go

import DEG_runner

from utils.ensembl2entrez import ensembl2entrez_convertor
from utils.scripts import format_script
from utils.server import get_parameters

ALGO_NAME = "reactomefi"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)
NETWORK_NAME = "dip"



if __name__ == "__main__":
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params


    bg_genes = get_bg_genes()

    script_name = "run_{}.sh".format(ALGO_NAME)
    format_script(os.path.join(constants.SH_DIR, script_name), BASE_FOLDER=constants.BASE_PROFILE, DATASET_DIR=constants.DATASET_DIR, ALGO_DIR=ALGO_DIR, NETWORK_NAME=NETWORK_NAME)

    script_path = os.path.join(constants.SH_DIR, script_name)
    subprocess.Popen("bash {}/run_{}.sh".format(constants.SH_DIR, ALGO_NAME), shell=True,
                           stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()

    results = file(os.path.join(constants.OUTPUT_DIR, "reactomefi_modules.txt")).readlines()
    module_genes = [x.split()[0] for x in results[1:]]

    file(os.path.join(constants.OUTPUT_DIR,"{}_module_genes.txt".format(ALGO_NAME)), "w+").write("\n".join(module_genes))
    # utils.go.check_group_enrichment(module_genes, bg_genes)
    sys.stdout.write(os.path.join(constants.OUTPUT_DIR,"{}_module_genes.txt".format(ALGO_NAME)))







