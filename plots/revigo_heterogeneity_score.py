import math
import random
import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.SemSim.SetSemSim import SetSemSim

import matplotlib
matplotlib.use("Agg")
import seaborn as sns
from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star

import argparse

import utils.go_hierarcies as go_hierarcies

import simplejson as json

import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def main(algos=None, cutoffs=[1.0, 2.0, 3.0, 4.0, 5.0], ax=None, axs_violin=None, title=""):

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')

    df_summary_agg = pd.DataFrame()
    for cutoff in cutoffs:
        df_summary=pd.read_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, "solution_richness_matrix_{}_{}.tsv".format(prefix, cutoff)),
            sep='\t', index_col=0).loc[algos]

        for a in set(df_summary.index).intersection(constants.ALGOS):
            for b in df_summary.columns:
                df_summary_agg=df_summary_agg.append({"algo": a, "dataset" : b, "cutoff": cutoff, "value" : df_summary.loc[a,b]}, ignore_index=True)


    for i_cutoff, cutoff in enumerate(cutoffs):
        my_order = df_summary_agg[df_summary_agg["cutoff"]==cutoff].groupby(by=["algo"])["value"].median().sort_values().index
        g=sns.violinplot(x="algo", y="value", data=df_summary_agg[df_summary_agg["cutoff"]==cutoff], ax=axs_violin[i_cutoff], order=my_order, palette={a: constants.COLORDICT[a] for a in my_order})
        g.set_xticklabels(g.get_xticklabels(), rotation=45)
    results = {}
    for cutoff in cutoffs:
        df_summary=pd.read_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, "solution_richness_matrix_{}_{}.tsv".format(prefix, cutoff)),
            sep='\t', index_col=0).loc[algos]

        for k, v in df_summary.iterrows():
            if k not in results:
                results[k] = []
            results[k].append(v.values)


    i=0
    df_average=pd.DataFrame()
    df_median=pd.DataFrame()
    df_std = pd.DataFrame()
    for k,v in sorted(list(results.iteritems()),key=lambda a: a[0]):
        if len(v)>0:
            print k, v
        y_average=[ np.average(cur_measurement) for cur_measurement in v]

        y_median=[ np.median(cur_measurement) for cur_measurement in v]
        y_stds= [np.std(cur_measurement) for cur_measurement in v]
        df_average = df_average.append(pd.Series(y_average,index=[1,2,3,4,5], name=k))
        df_median = df_median.append(pd.Series(y_median,index=[1,2,3,4,5], name=k))
        df_std = df_std.append(pd.Series(y_stds,index=[1,2,3,4,5], name=k))

        ax.plot(cutoffs,y_median,label=constants.ALGOS_ACRONYM[k], c=constants.COLORDICT[k])

        # for x, y, std in zip(cutoffs, ys, stds):
        #     ax.errorbar(x - 0.15 + 0.05 * i, y, yerr=std, linestyle='None', marker='^', c=colorlist[i])

        i += 1

    df_std.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"biological_richness_std_{}.tsv".format(prefix)), sep='\t')
    df_average.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"biological_richness_average_{}.tsv".format(prefix)), sep='\t')
    df_median.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"biological_richness_median_{}.tsv".format(prefix)), sep='\t')

    ax.set_xlabel("similarity cutoff", fontdict={"size": 22})
    ax.set_ylabel("median # of non-redundant terms", fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    ax.legend(fontsize=17, loc=(0,1.1), ncol=2, facecolor='#ffffff')



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,ERS_1,IEM,SHERA,SHEZH_1") # ",Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="dcem,jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td,keypathwayminer_INES_GREEDY,hotnet2") # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    parser.add_argument('--pf', dest='pf', default=3)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX"))
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    args = parser.parse_args()

    prefix = args.prefix
    file_format = args.file_format
    base_folder= args.base_folder
    sim_method = args.sim_method
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    cutoffs=[1.0,2.0,3.0,4.0,5.0]
    pf=int(args.pf)


    fig, axs = plt.subplots(1,2,figsize=(20,10))
    fig_violin, axs_violin = plt.subplots(2,len(cutoffs),figsize=(4*len(cutoffs)*2,10))
    prefix = "GE"
    datasets=["TNFa_2" ,"HC12","ROR_1","ERS_1","IEM","SHERA","SHEZH_1"]
    algos=["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY"]
    main(algos=algos, cutoffs=cutoffs, ax=axs[0], axs_violin=axs_violin[0], title="GE")

    prefix = "PASCAL_SUM"
    datasets=["Breast_Cancer.G50", "Crohns_Disease.G50", "Schizophrenia.G50", "Triglycerides.G50", "Type_2_Diabetes.G50"]
    algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox"]
    main(algos=algos, cutoffs=cutoffs, ax=axs[1], axs_violin=axs_violin[1], title="GWAS")

    fig.text(0.01, 0.97, "A:", weight='bold', fontsize=22)
    fig.text(0.5, 0.97, "B:", weight='bold', fontsize=22)
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"figure_13.png"))

    fig_violin.text(0.0, 0.97, "A:", weight='bold', fontsize=22)
    fig_violin.text(0.0, 0.5, "B:", weight='bold', fontsize=22)
    fig_violin.tight_layout()
    fig_violin.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_13_violin.png"))