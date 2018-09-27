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

def sif2hotnet2(network_name):
    network_file_name = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(network_name))

    network_df = pd.read_csv(network_file_name, sep="\t")
    src = np.array(network_df["ID_interactor_A"])
    dst = np.array(network_df["ID_interactor_B"])

    vertices = list(set(np.append(src,dst)))

    lns = ["{} {}".format(i+1, cur) for i, cur in enumerate(vertices)]
    file(os.path.join(constants.CACHE_DIR, "hotnet2_vertices.txt"), "w+").write("\n".join(lns))

    inds = map(lambda x: (vertices.index(x[0])+1, vertices.index(x[1])+1), zip(src,dst))


    lns = ["{} {} {}".format(cur[0], cur[1], 1) for cur in inds]
    file(os.path.join(constants.CACHE_DIR, "hotnet2_edges.txt"), "w+").write("\n".join(lns))

    print subprocess.Popen("bash ../sh/scripts/prepare_hotnet2.sh", shell=True, stdout=subprocess.PIPE).stdout.read() # cwd=dir_path

def run_hotnet2(deg_file_name, network_file_name):
    script = file("scripts/bionet.r").read()
    return run_rscript(script=script, output_vars = ["module_genes", "bg_genes"], network_file_name=network_file_name, deg_file_name=deg_file_name)


def prepare_input(method=constants.DEG_EDGER, network_name="dip"):
    deg_file_name = os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method))
    network_file_name = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(network_name))
    heat_file_name = os.path.join(constants.CACHE_DIR, "heatfile.txt")
    if not os.path.exists(os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method))):
        DEG_runner.main(method=method)

    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=deg_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)
    ind = np.where(h_cols=="qval")[0][0]

    lns = []
    for i, cur in enumerate(deg_data):
        lns.append(" ".join([str(h_rows[i]), str(-log10(cur[ind]))]))
    file(heat_file_name,"w+").write("\n".join(lns))

    sif2hotnet2(network_name=network_name)

    vertices = set(pd.read_csv(os.path.join(constants.CACHE_DIR, "hotnet2_vertices.txt"), sep=" ", index_col=False, header=None).ix[:,1])
    degs = set(pd.read_csv(deg_file_name, sep="\t").ix[:,0])
    bg_genes = list(vertices.intersection(degs))
    file(os.path.join(constants.OUTPUT_DIR, "hotnet2_bg_genes.txt"), "w+").write("\n".join(bg_genes))
    return heat_file_name, network_file_name, bg_genes


if __name__ == "__main__":
    heat_file_name, network_file_name, bg_genes = prepare_input(method=constants.DEG_EDGER)
    print subprocess.Popen("bash ../sh/scripts/run_hotnet2.sh", shell=True,
                           stdout=subprocess.PIPE).stdout.read()  # cwd=dir_path

    results = json.load(file(os.path.join(constants.OUTPUT_DIR,"results","consensus","subnetworks.json")))
    modules = [x["core"] for x in results["consensus"]]
    module_genes = reduce((lambda x, y: x + y), modules)
    file(os.path.join(constants.OUTPUT_DIR,"hotnet2_module_genes.txt"), "w+").write("\n".join(module_genes))

    utils.go.check_group_enrichment(module_genes, bg_genes)






