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

ALGO_NAME = "netbox"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)


def prepare_input(method=constants.DEG_EDGER, network_name="dip"):
    deg_file_name = os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method))
    if not os.path.exists(deg_file_name):
        DEG_runner.main(method=method)


    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=deg_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)

    ind = np.where(h_cols=="qval")[0][0]
    ordered_ind = np.argsort(deg_data[:,ind])
    deg_data=deg_data[ordered_ind,:]
    h_rows=h_rows[ordered_ind]
    last_q_index = np.where(deg_data[:,np.where(h_cols=="qval")[0][0]]>0.05)[0][0]
    ge_list = os.path.join(ALGO_DIR, "ge_list.txt")
    file(os.path.join(ALGO_DIR, "ge_list.txt"), "w+").write("\n".join([x for x in h_rows[:last_q_index] if len(ensembl2entrez_convertor([x]))>0 ])) # ensembl2entrez_convertor([x])[0]

    conf_file = "conf.props"
    format_script(os.path.join(ALGO_DIR, conf_file), pval_threshold=0.05, sp_threshold=2, gene_file=ge_list)

    return conf_file




if __name__ == "__main__":
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params

    conf_file = prepare_input(method=constants.DEG_EDGER)

    bg_genes = get_bg_genes()

    script_name = "run_{}.sh".format(ALGO_NAME)
    format_script(os.path.join(constants.SH_DIR, script_name), BASE_FOLDER=constants.BASE_PROFILE, DATASET_DIR=constants.DATASET_DIR, NETBOX_DIR=ALGO_DIR)

    script_path = os.path.join(constants.SH_DIR, script_name)
    print subprocess.Popen("bash {}/run_{}.sh".format(constants.SH_DIR, ALGO_NAME), shell=True,
                           stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()

    results = file(os.path.join(ALGO_DIR, "modules.txt")).readlines()
    module_genes = [x.split()[0] for x in results[1:]]

    file(os.path.join(constants.OUTPUT_DIR,"netbox_module_genes.txt"), "w+").write("\n".join(module_genes))
    utils.go.check_group_enrichment(module_genes, bg_genes)
    sys.stdout.write(os.path.join(constants.OUTPUT_DIR,"netbox_module_genes.txt"))






