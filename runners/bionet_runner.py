#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import os

# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()


import constants

from r.r_runner import run_rscript

import DEG_runner

import utils.go

def run_bionet(deg_file_name, network_file_name):
    script = file("scripts/bionet.r").read()
    return run_rscript(script=script, output_vars = ["module_genes", "bg_genes"], network_file_name=network_file_name, deg_file_name=deg_file_name)


def prepare_input(method=constants.DEG_EDGER, network_name="dip"):
    deg_file_name = str(os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method)))
    network_file_name = str(os.path.join(constants.NETWORKS_DIR, "{}.sif".format(network_name)))
    if not os.path.exists(os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method))):
        DEG_runner.main(method=method)

    return deg_file_name, network_file_name


if __name__ == "__main__":
    deg_file_name, network_file_name = prepare_input(method=constants.DEG_DESEQ)

    results = run_bionet(deg_file_name, network_file_name)
    module_genes = results["module_genes"]
    bg_genes = results["bg_genes"]

    file(os.path.join(constants.OUTPUT_DIR, "bionet_module_genes.txt"),"w+").write("\n".join(module_genes))
    file(os.path.join(constants.OUTPUT_DIR, "bionet_bg_genes.txt"), "w+").write("\n".join(bg_genes))

    utils.go.check_group_enrichment(module_genes,bg_genes)





