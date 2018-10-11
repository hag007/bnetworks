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
def sif2netbox(network_name):
    network_file_name = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(network_name))

    network_df = pd.read_csv(network_file_name, sep="\t")
    src = np.array(network_df["ID_interactor_A"])
    dst = np.array(network_df["ID_interactor_B"])

    vertices = list(set(np.append(src,dst)))

    lns = ["{} {}".format(i+1, cur) for i, cur in enumerate(vertices)]
    file(os.path.join(constants.CACHE_DIR, "netbox1_vertices.txt"), "w+").write("\n".join(lns))

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

    ge_list = os.path.join(constants.CACHE_DIR, "ge_list.txt")
    file(os.path.join(constants.CACHE_DIR, "ge_list.txt"), "w+").write("\n".join([x for x in h_rows[:last_q_index] if len(ensembl2entrez_convertor([x]))>0 ])) # ensembl2entrez_convertor([x])[0]

#    sif2netbox(network_name=network_name)
    conf_format_file = os.path.join("..", "repos", "netbox", "netbox1.props.format")
    conf_file = os.path.join(constants.CACHE_DIR, "netbox1.props")
    conf=file(conf_format_file).read()
    conf=conf.format(pval_threshold=0.05, sp_threshold=2, gene_file=ge_list)
    file(conf_file, "w+").write(conf)
    bg_genes = list(vertices.intersection(set(total_deg_genes)))
    bg_genes_file_name=os.path.join(constants.OUTPUT_DIR, "netbox_bg_genes.txt")
    file(os.path.join(constants.OUTPUT_DIR, "netbox_bg_genes.txt"), "w+").write("\n".join(bg_genes))
    return network_file_name, bg_genes  , conf_file


if __name__ == "__main__":
    network_file_name, bg_genes, conf_file = prepare_input(method=constants.DEG_EDGER)
    print subprocess.Popen("bash ../sh/scripts/run_netbox.sh", shell=True,
                           stdout=subprocess.PIPE).stdout.read()  # cwd=dir_path

    results = file("/home/hag007/repos/bnetworks/runners/modules.txt").readlines()

    module_genes = [x.split()[0] for x in results[1:]]
    file(os.path.join(constants.OUTPUT_DIR,"netbox_module_genes.txt"), "w+").write("\n".join(module_genes))

    utils.go.check_group_enrichment(module_genes, bg_genes)






