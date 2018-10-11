#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import os
import numpy as np
from numpy import log10
import pandas as pd
import json
import subprocess

# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

from r.r_runner import run_rscript

import DEG_runner

import infra

import utils.go

from utils.ensembl2gene_symbol import e2g_convertor
from utils.ensembl2entrez import ensembl2entrez_convertor


def prepare_input(method=constants.DEG_EDGER, network_name="dip"):
    deg_file_name = os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method))
    network_file_name = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(network_name))

    network_df = pd.read_csv(network_file_name, sep="\t")
    src = np.array(network_df["ID_interactor_A"])
    dst = np.array(network_df["ID_interactor_B"])

    vertices = set(np.append(src,dst))


    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=deg_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)
    total_deg_genes = h_rows
    ind = np.where(h_cols=="qval")[0][0]
    ordered_ind = np.argsort(deg_data[:,ind])
    deg_data=deg_data[ordered_ind,:]
    h_rows=h_rows[ordered_ind]
    last_q_index = np.where(deg_data[:,np.where(h_cols=="qval")[0][0]]>0.05)[0][0]

    bg_genes = list(vertices.intersection(set(total_deg_genes)))
    bg_genes_file_name=os.path.join(constants.OUTPUT_DIR, "jactivemodules_bg_genes.txt")
    file(os.path.join(constants.OUTPUT_DIR, "netbox_bg_genes.txt"), "w+").write("\n".join(bg_genes))
    return network_file_name, bg_genes


if __name__ == "__main__":
    score_method=constants.DEG_EDGER
    network_file_name, bg_genes = prepare_input(method=score_method)
    formatted_script = file("../sh/scripts/run_jactivemodules.sh.format").read().format(base_folder=constants.BASE_PROFILE, dataset="TNFa_2",
                                                                     network_name="dip", rank_method=score_method, is_greedy="true")
    file("../sh/scripts/run_jactivemodules.sh","w+").write(formatted_script)
    print subprocess.Popen("bash ../sh/scripts/run_jactivemodules.sh", shell=True,
                           stdout=subprocess.PIPE).stdout.read()  # cwd=dir_path

    results = file("/home/hag007/bnet/datasets/TNFa_2/output/jactivemodules_results.txt").readlines()

    module_genes = [y for x in results for y in x.split()[:-1]]
    file(os.path.join(constants.OUTPUT_DIR,"jactivemodules_sa_module_genes.txt"), "w+").write("\n".join(module_genes))

    utils.go.check_group_enrichment(module_genes, bg_genes)


    file(os.path.join(constants.OUTPUT_DIR,"jactivemodules_greedy_module_genes.txt"), "w+").write("\n".join(module_genes))

    formatted_script = file("../sh/scripts/run_jactivemodules.sh.format").read().format(base_folder=constants.BASE_PROFILE, dataset="TNFa_2",
                                                                     network_name="dip", rank_method=score_method, is_greedy="false")
    file("../sh/scripts/run_jactivemodules.sh","w+").write(formatted_script)
    print subprocess.Popen("bash ../sh/scripts/run_jactivemodules.sh", shell=True,
                           stdout=subprocess.PIPE).stdout.read()  # cwd=dir_path

    results = file("/home/hag007/bnet/datasets/TNFa_2/output/jactivemodules_results.txt").readlines()

    module_genes = [y for x in results for y in x.split()[:-1]]
    file(os.path.join(constants.OUTPUT_DIR,"jactivemodules_sa_module_genes.txt"), "w+").write("\n".join(module_genes))

    utils.go.check_group_enrichment(module_genes, bg_genes)






