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
import json

# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import infra

import utils.go

import DEG_runner
from utils.scripts import format_script
from utils.network import get_network_genes
from utils.network import build_all_reports
import utils.server as server
from utils.network import output_modules

ALGO_NAME = "jactivemodules"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)


def init_specific_params(search_method):
    results_file_name = "{}/{}_{}_results.txt".format(constants.OUTPUT_DIR, ALGO_NAME, search_method)
    return results_file_name


def extract_modules_and_bg(bg_genes, results_file_name, modules_genes_file_name):
    results = file(results_file_name.format(constants.OUTPUT_DIR, ALGO_NAME)).readlines()
    modules = [x.split()[:-1] for x in results]
    modules = [cur for cur in modules if len(cur) > 3]
    all_bg_genes = [bg_genes for x in modules]
    module_genes = [y for x in modules for y in x]
    file(modules_genes_file_name, "w+").write("\n".join(module_genes))
    print "extracted {} modules".format(len(modules))
    return all_bg_genes, modules


def main(dataset_name=constants.DATASET_NAME, disease_name=None, expected_genes = None, score_method=constants.DEG_EDGER, network_file_name="dip.sif"):
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    search_method = "greedy"
    network_file_name, score_file_name, score_method, bg_genes= server.init_common_params(network_file_name, score_method)

    results_file_name = init_specific_params(search_method)

    script_file_name=format_script(os.path.join(constants.SH_DIR, "run_{}.sh".format(ALGO_NAME)), BASE_FOLDER=constants.BASE_PROFILE,
                  DATASET_DIR=constants.DATASET_DIR,
                  ALGO_DIR=ALGO_DIR, NETWORK_NAME=network_file_name, SCORE_FILE_NAME=score_file_name,
                  IS_GREEDY=str(search_method == "greedy"), OUTPUT_FILE=results_file_name, NUM_OF_MODULES=10, OVERLAP_THRESHOLD=0)
    print "running :{}".format(script_file_name)
    subprocess.Popen("bash {}".format(script_file_name), shell=True,
                     stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()

    os.remove(script_file_name)
    modules_genes_file_name = os.path.join(constants.OUTPUT_DIR, "{}_{}_module_genes.txt".format(ALGO_NAME, search_method))
    all_bg_genes, modules = extract_modules_and_bg(bg_genes, results_file_name, modules_genes_file_name)

    output_base_dir = ""
    if constants.REPORTS:
        output_base_dir = build_all_reports(ALGO_NAME + "_" + search_method, dataset_name, modules, all_bg_genes, score_file_name, network_file_name, disease_name, expected_genes)
 
    output_file_name=os.path.join(constants.OUTPUT_DIR,
                  "{}_{}_client_output.txt".format(ALGO_NAME, search_method))
    output_modules(output_file_name, modules, score_file_name, output_base_dir)


if __name__ == "__main__":
    ds = ["GE_TNFa_2"] # NADAV" ,"GE_NADAV_2"]  # ["PASCAL_SUM_Breast_Cancer.G50", "PASCAL_SUM_Crohns_Disease.G50", "PASCAL_SUM_Schizophrenia.G50", "PASCAL_SUM_Triglycerides.G50", "PASCAL_SUM_Type_2_Diabetes.G50"]
    for cur in ds:
        constants.update_dirs(DATASET_NAME_u=cur)  # Type_2_Diabetes Crohns_Disease
        main(dataset_name=constants.DATASET_NAME, score_method=constants.DEG_EDGER)



