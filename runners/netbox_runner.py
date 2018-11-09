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
from utils.network import get_network_genes
from utils.network import build_all_reports
from utils.server import get_score
import infra

import utils.go

import DEG_runner

from utils.ensembl2entrez import ensembl2entrez_convertor
from utils.scripts import format_script
import utils.server as server
from utils.network import output_modules

ALGO_NAME = "netbox"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)

NETWORK_NAME = "dip"


def init_specific_params(score_file_name):

    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=score_file_name)
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

def extract_modules_and_bg(bg_genes):
    results = file(os.path.join(ALGO_DIR, "modules.txt")).readlines()
    modules = [[] for x in range(max([int(x.strip().split(" =")[1]) for x in results[1:]]) + 1)]
    for x in results[1:]:
        if int(x.strip().split(" =")[1]) != -1:
            modules[int(x.strip().split(" =")[1])].append(x.strip().split(" =")[0])
        else:
            modules.append([x.strip().split(" =")[0]])
    modules = filter(lambda x: len(x) > 0, modules)
    all_bg_genes = [bg_genes for x in modules]
    return modules, all_bg_genes


def main(dataset_name=constants.DATASET_NAME, disease_name=None, expected_genes = None):
    global NETWORK_NAME
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    network_file_name, score_file_name, score_method, bg_genes = server.init_common_params(NETWORK_NAME)

    init_specific_params(score_file_name)
    script_name = "run_{}.sh".format(ALGO_NAME)
    format_script(os.path.join(constants.SH_DIR, script_name), BASE_FOLDER=constants.BASE_PROFILE,
                  DATASET_DIR=constants.DATASET_DIR, NETBOX_DIR=ALGO_DIR)
    print subprocess.Popen("bash {}/run_{}.sh".format(constants.SH_DIR, ALGO_NAME), shell=True,
                           stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()

    modules, all_bg_genes = extract_modules_and_bg(bg_genes)
    output_base_dir = ""
    if constants.REPORTS:
        output_base_dir = build_all_reports(ALGO_NAME, modules, all_bg_genes, score_file_name, network_file_name, disease_name, expected_genes)
    output_file_name = os.path.join(constants.OUTPUT_DIR,
                                    "{}_client_output.txt".format(ALGO_NAME))
    output_modules(output_file_name, modules, score_file_name, output_base_dir )


if __name__ == "__main__":
    constants.update_dirs(DATASET_NAME_u="MCF7_2")
    main()






