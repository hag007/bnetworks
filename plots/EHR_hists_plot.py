import sys
sys.path.insert(0, '../')

import pandas as pd

import seaborn as sns
sns.set(color_codes=True)
import numpy as np

import constants

import argparse

import os

from matplotlib import pyplot as plt

from statsmodels.sandbox.stats.multicomp import fdrcorrection0

PVAL="emp_pval"
QVAL="qval"
METHOD='average'


def main(base_folder=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX"), file_format="emp_diff_modules_{}_{}_passed_oob.tsv", algo=None, dataset=None, ax=None):
    ax.set_facecolor('#ffffff')


    go_terms_emp_tables=pd.read_csv(os.path.join(base_folder, file_format.format(dataset,algo)),sep='\t', index_col=0).dropna()
    hg_hist = go_terms_emp_tables['hg_pval_max']

    go_terms_emp_tables=go_terms_emp_tables[go_terms_emp_tables["passed_oob_permutation_test"].str.contains('True', regex=False)==True]
    emp_hist = go_terms_emp_tables['hg_pval_max']


    ax=sns.distplot(hg_hist, norm_hist=False, kde=False, label="# HG enriched terms", bins=25, hist_kws=dict(alpha=0.5, range=(0,25)), ax=ax)
    ax=sns.distplot(emp_hist, norm_hist=False, kde=False, label="# EMP validated terms", bins=25, hist_kws=dict(alpha=0.5, range=(0,25)), ax=ax)
    ax.set_xlabel("enrichment score: -log10(pval)", fontsize=20)
    ax.set_ylabel("# of GO terms", fontsize=20)
    ax.legend(prop={"size": 20}, loc='upper right', facecolor='#ffffff')
    font={'size' : 17}
    # plt.rc('font' , **font)
    # subplot.legend()
    # subplot.set_title("algorithm: {}, dataset: {}\n"
    #           "EHR: {}".format(algo, dataset, round(len(emp_hist)/float(len(hg_hist)), 2)), fontdict={"size": 18})
    ax.set_title("Dataset: {}, Algorithm: {}\nEHR: {}".format(dataset, algo, round(len(emp_hist) / float(len(hg_hist)), 2)), fontdict={"size": 20})

    # plt.clf()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,Schizophrenia.G50") # Breast_Cancer.G50
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,netbox") # jactivemodules_sa,netbox,bionet,hotnet2,keypathwayminer_INES_GREEDY
    parser.add_argument('--pf', dest='pf', default=4)
    args = parser.parse_args()

    prefix = args.prefix
    datasets=args.datasets.split(",")
    algos = args.algos.split(",")

    pf=int(args.pf)
    ds_summary=pd.DataFrame()
    # figure, subplots = plt.subplots(7, 6, figsize=(35,30)) # figsize=(27,23)
    figure, subplots = plt.subplots(2, 2, figsize=(15, 15))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.05)
    for i, algo in enumerate(algos):
        for j, dataset in enumerate(datasets):
            main(algo=algo, dataset=dataset, ax=subplots[j][i])

    figure.text(0.01, 0.97, "A:", weight='bold', fontsize=22)
    figure.text(0.55, 0.97, "B:", weight='bold', fontsize=22)
    figure.text(0.01, 0.5, "C:", weight='bold', fontsize=22)
    figure.text(0.55, 0.5, "D:", weight='bold', fontsize=22)
    figure.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_12.png".format(",".join(algos), ",".join(datasets))))