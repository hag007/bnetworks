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


def main(base_folder="/home/hag007/Desktop/aggregate_report/oob", algo=None, datasets=None, subplot=None):

    go_terms_emp_tables={}
    go_terms_hg_tables = {}
    go_terms_emp = {}
    go_terms_hg = {}
    norm=0
    for cur in [x for x in os.listdir(base_folder) if not os.path.isdir(x) and x and x.endswith(".tsv")]:
        go_terms_emp_tables[cur[9:-15]]=pd.read_csv(os.path.join(base_folder, cur),sep='\t', index_col=0)
        norm=max(-np.log10(go_terms_emp_tables[cur[9:-15]][PVAL].min()),norm)

        go_terms_hg_tables[cur[9:-15]] = pd.read_csv(os.path.join(base_folder, cur), sep='\t', index_col=0)
        norm = max(-np.log10(go_terms_hg_tables[cur[9:-15]][PVAL].min()), norm)

    for k,v in go_terms_emp_tables.iteritems():
        go_terms_emp[k]= v[v["passed_oob_permutation_test"].str.contains('True', regex=False)==True]

    for k,v in go_terms_hg_tables.iteritems():
        go_terms_hg[k]= v.dropna()

    hg_hist = np.array([])
    for k, v in go_terms_hg.iteritems():
        if sum([d in k for d in datasets]) and algo in k:
            hg_hist = np.append(hg_hist, v['hg_pval_max'])
            # filtered_genes=v['hg_pval']
            # sorted_genes_hg = filtered_genes.sort_values( ascending=False)
            # hg_scores = np.append(sorted_genes_hg.values,
            #                               np.zeros(7435 - np.size(sorted_genes_hg.values)))
            # hg_pvals = [10 ** (-x) for x in hg_scores]
            # fdr_results = fdrcorrection0(hg_pvals, alpha=0.05, method='indep', is_sorted=False)
            # hg_fdr_adjusted_scores=hg_scores[fdr_results[0]]
            # hg_hist = np.append(hg_hist, hg_fdr_adjusted_scores)
            # # n_hg_true = len([cur for cur in fdr_results[0] if cur])

    emp_hist = np.array([])
    for k, v in go_terms_emp.iteritems():
        if sum([d in k for d in datasets]) and algo in k:
            emp_hist=np.append(emp_hist, v['hg_pval_max'])

    ax=sns.distplot(hg_hist, norm_hist=False, kde=False, label="# HG enriched terms", bins=25, hist_kws=dict(alpha=0.5, range=(0,25)), ax=subplot)
    ax=sns.distplot(emp_hist, norm_hist=False, kde=False, label="# EMP validated terms", bins=25, hist_kws=dict(alpha=0.5, range=(0,25)), ax=subplot)
    subplot.set_xlabel("enrichment score: -log10(pval)", fontsize=15)
    subplot.set_ylabel("counts", fontsize=15)
    subplot.legend(prop={"size": 15}, loc='upper right')
    font={'size' : 15}
    # plt.rc('font' , **font)
    # subplot.legend()
    subplot.set_title("algorithm: {}, dataset: {}\n"
              "EHR: {}".format(algo, ",".join(datasets), round(len(emp_hist)/float(len(hg_hist)), 2)), fontdict={"size": 15})
    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "TP_histogram_{}_{}.png".format(algo, ",".join(datasets))))
    # plt.clf()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,netbox,bionet,hotnet2,keypathwayminer_INES_GREEDY")
    parser.add_argument('--pf', dest='pf', default=4)
    args = parser.parse_args()

    prefix = args.prefix
    datasets=args.datasets.split(",")
    # datasets=["{}_{}".format(prefix,x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")

    # algos=["jactivemodules_greedy", "netbox"]
    # datasets=["TNFa_2", "ROR_1"]

    pf=int(args.pf)
    ds_summary=pd.DataFrame()
    figure, subplots = plt.subplots(7, 6, figsize=(35,30)) # figsize=(27,23)
    # figure, subplots = plt.subplots(2, 2, figsize=(15, 15))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.05)
    for i, algo in enumerate(algos):
        for j, dataset in enumerate(datasets):
            main(algo=algo, datasets=[dataset], subplot=subplots[j][i])
