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

import shutil
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
    sig_binary_col = deg_data[:,np.where(h_cols=="qval")[0][0]]<0.05
    sig_binary_output = np.c_[h_rows,  np.array(sig_binary_col, dtype=np.int)]
    file(os.path.join(constants.CACHE_DIR, "binary_score_{}.txt".format(method)), "w+").write("\n".join(["\t".join(x) for x in sig_binary_output]))
    bg_genes = list(vertices)# list(vertices.intersection(set(total_deg_genes)))
    bg_genes_file_name=os.path.join(constants.OUTPUT_DIR, "keypathwayminer_bg_genes.txt")
    file(os.path.join(constants.OUTPUT_DIR, bg_genes_file_name), "w+").write("\n".join(bg_genes))
    return network_file_name, bg_genes


if __name__ == "__main__":

    if os.path.exists("../repos/keypathwayminer/results"):
        results = shutil.rmtree("../repos/keypathwayminer/results")

    score_method=constants.DEG_EDGER
    network_file_name, bg_genes = prepare_input(method=score_method)
    formatted_script = file("../sh/scripts/run_keypathwayminer.sh.format").read().format(base_folder=constants.BASE_PROFILE, dataset="TNFa_2",
                                                                     network_name="dip", rank_method=score_method, is_greedy="true")
    file("../sh/scripts/run_keypathwayminer.sh","w+").write(formatted_script)

    formatted_script = file("../repos/keypathwayminer/kpm.properties.format").read().format(base_folder=constants.BASE_PROFILE, network_name="dip")
    file("../repos/keypathwayminer/kpm.properties","w+").write(formatted_script)

    formatted_script = file("../repos/keypathwayminer/datasets_file.txt.format").read().format(base_folder=constants.BASE_PROFILE, dataset="TNFa_2", score_method=score_method)
    file("../repos/keypathwayminer/datasets_file.txt","w+").write(formatted_script)

    print subprocess.Popen("bash ../../sh/scripts/run_keypathwayminer.sh", shell=True,
                           stdout=subprocess.PIPE, cwd="../repos/keypathwayminer").stdout.read()  # cwd=dir_path

    i=1
    module_genes=[]
    while os.path.exists("../repos/keypathwayminer/results/Pathway-{}-NODES-.txt".format("%02d" % (i,))):
        results = file("../repos/keypathwayminer/results/Pathway-{}-NODES-.txt".format("%02d" % (i,))).readlines()
        results = map(lambda x: x.strip(), results)
        module_genes = module_genes + results
        i+=1
    module_genes = list(set(module_genes))


    file(os.path.join(constants.OUTPUT_DIR,"keypathwayminer_module_genes.txt"), "w+").write("\n".join(module_genes))

    utils.go.check_group_enrichment(module_genes, bg_genes)






