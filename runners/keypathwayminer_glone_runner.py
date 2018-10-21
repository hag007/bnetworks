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

import shutil
import constants

from utils.scripts import format_script
from utils.network import get_bg_genes

from utils.server import get_parameters

import infra
import DEG_runner

import utils.go

ALGO_NAME = "keypathwayminer"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)

NETWORK_NAME = "dip"

def prepare_input(method=constants.DEG_EDGER):
    deg_file_name = os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method))
    if not os.path.exists(deg_file_name):
        DEG_runner.main(method=method)
    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=deg_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)
    ind = np.where(h_cols=="qval")[0][0]
    ordered_ind = np.argsort(deg_data[:,ind])
    deg_data=deg_data[ordered_ind,:]
    h_rows=h_rows[ordered_ind]
    sig_binary_col = deg_data[:,np.where(h_cols=="qval")[0][0]]<0.05
    sig_binary_output = np.c_[h_rows,  np.array(sig_binary_col, dtype=np.int)]
    file(os.path.join(constants.CACHE_DIR, "binary_score_{}.txt".format(method)), "w+").write("\n".join(["\t".join(x) for x in sig_binary_output]))

    bg_genes_file_name=os.path.join(constants.OUTPUT_DIR, "{}_bg_genes.txt".format(ALGO_NAME))
    file(os.path.join(constants.OUTPUT_DIR, bg_genes_file_name), "w+").write("\n".join(bg_genes))



def format_scripts(network_name="dip", algo_dir=ALGO_DIR,STRATEGY="ines"):
    format_script(os.path.join(constants.SH_DIR, "run_{}.sh".format(ALGO_NAME)), BASE_FOLDER=constants.BASE_PROFILE, DATASET_DIR=constants.DATASET_DIR, STRATEGY=STRATEGY)
    format_script(os.path.join(ALGO_DIR, "kpm.properties"), base_folder=constants.BASE_PROFILE, network_name=network_name, algo_dir=ALGO_DIR)
    format_script(os.path.join(ALGO_DIR, "datasets_file.txt"), base_folder=constants.BASE_PROFILE, dataset=constants.DATASET_NAME, score_method=score_method)


def extract_module_genes(STRATEGY):
    i = 1
    module_genes = []
    while os.path.exists(os.path.join(ALGO_DIR, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))):
        results = file(
            os.path.join(ALGO_DIR, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))).readlines()
        results = map(lambda x: x.strip(), results)
        module_genes = module_genes + results
        i += 1
    module_genes = list(set(module_genes))
    file(os.path.join(constants.OUTPUT_DIR, "{}_{}_module_genes.txt".format(ALGO_NAME, STRATEGY)), "w+").write("\n".join(module_genes))
    return module_genes, os.path.join(constants.OUTPUT_DIR, "{}_{}_module_genes.txt".format(ALGO_NAME,STRATEGY))


if __name__ == "__main__":
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params
    STRATEGY="GLONE"
    network_file_dir = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(NETWORK_NAME))
    bg_genes = get_bg_genes()

    if os.path.exists(os.path.join(ALGO_DIR, "results")):
        shutil.rmtree(os.path.join(ALGO_DIR, "results"))




    score_method=constants.DEG_EDGER
    prepare_input(method=score_method)

    format_scripts(network_name=NETWORK_NAME, algo_dir=ALGO_DIR,STRATEGY=STRATEGY)

    print subprocess.Popen("bash {}/run_{}.sh".format(constants.SH_DIR, ALGO_NAME), shell=True,
                           stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()

    module_genes, module_genes_file_name = extract_module_genes(STRATEGY)


    utils.go.check_group_enrichment(module_genes, bg_genes)

    sys.stdout.write(module_genes_file_name)






