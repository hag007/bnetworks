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

# from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

PVAL="emp_pval"
QVAL="qval"
METHOD='average'


def main(base_folder=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX"), file_format="emp_diff_modules_{}_{}_passed_oob.tsv", algo=None, dataset=None, ax=None):
    ax.set_facecolor('#ffffff')


    go_terms_emp_tables=pd.read_csv(os.path.join(base_folder, file_format.format(dataset,algo)),sep='\t', index_col=0).dropna()
    hg_hist = go_terms_emp_tables['hg_pval_max']

    if len(hg_hist):
        go_terms_emp_tables=go_terms_emp_tables[go_terms_emp_tables["passed_oob_permutation_test"].str.contains('True', regex=False)==True]
        emp_hist = go_terms_emp_tables['hg_pval_max']

    else:
        hg_hist = [0]
        emp_hist = [0]


    ax=sns.distplot(hg_hist, norm_hist=False, kde=False, label="# HG enriched terms", bins=25, hist_kws=dict(alpha=0.5, range=(0,25)), ax=ax)
    ax=sns.distplot(emp_hist, norm_hist=False, kde=False, label="# EMP validated terms", bins=25, hist_kws=dict(alpha=0.5, range=(0,25)), ax=ax)
    ax.set_xlabel("enrichment score: -log10(pval)", fontsize=20)
    ax.set_ylabel("# of GO terms", fontsize=20)
    # plt.legend()# prop={"size": 20}, loc='upper right', facecolor='#ffffff')
    font={'size' : 17}
    # plt.rc('font' , **font)
    # subplot.legend()
    # subplot.set_title("algorithm: {}, dataset: {}\n"
    #           "EHR: {}".format(algo, dataset, round(len(emp_hist)/float(len(hg_hist)), 2)), fontdict={"size": 18})
    ax.set_title("Dataset: {}, algorithm: {}".format(dataset, constants.ALGOS_ACRONYM[algo]), fontdict={"size": 18})

    x0, xmax = ax.get_xlim()
    y0, ymax = ax.get_ylim()
    ehr_score=round(len(emp_hist) / float(len(hg_hist)), 2)
    ax.text(xmax/2.0, ymax/2.0, "EHR: {}".format(ehr_score), fontdict={"size": 20})

    return ehr_score
    # plt.clf()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,Schizophrenia.G50") # "TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM" # TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM") # "Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50" # "TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="r") # "PASCAL_SUM" "GE"
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,netbox") # jactivemodules_sa,bionet,netbox,keypathwayminer_INES_GREEDY,hotnet2")
    parser.add_argument('--pf', dest='pf', default=4)
    args = parser.parse_args()

    prefix = args.prefix
    datasets=args.datasets.split(",")
    algos = args.algos.split(",")

    pf=int(args.pf)
    ds_summary=pd.DataFrame()
    # figure, subplots = plt.subplots(len(datasets), len(algos), figsize=(35,30)) # figsize=(27,23)
    figure, subplots = plt.subplots(2, 2, figsize=(15, 15))


    df_summary=pd.DataFrame(index=[constants.ALGOS_ACRONYM[algo] for algo in algos], columns=datasets)
    for i, algo in enumerate(algos):
        for j, dataset in enumerate(datasets):
            ehr_score=main(algo=algo, dataset=dataset, ax=subplots[j][i])
            df_summary.iloc[i, j] = ehr_score

    # figure.text(0.01, 0.97, "A:", weight='bold', fontsize=22)
    # figure.text(0.55, 0.97, "B:", weight='bold', fontsize=22)
    # figure.text(0.01, 0.5, "C:", weight='bold', fontsize=22)
    # figure.text(0.55, 0.5, "D:", weight='bold', fontsize=22)

    legends = ["# of HG enriched terms", "# of EV-terms"]
    # lgd=plt.legend(legends, loc = 'upper left', bbox_to_anchor = (-len(algos)*0.78, len(datasets)*1.528),  facecolor='#ffffff', ncol=2, prop={'size': 30})
    lgd = plt.figlegend(legends, loc='upper center',
                     facecolor='#ffffff', ncol=2, prop={'size': 14})
    figure.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    figure.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_12_{}.png".format(prefix)) , bbox_extra_artists=(lgd,), bbox_inches='tight')
    df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_12_{}_matrix.tsv").format(prefix), sep='\t', index_label="algos")