import json
import matplotlib
from matplotlib_venn import venn2
matplotlib.use('Agg')

from pandas._libs.parsers import k
import sys
sys.path.insert(0, '../')

import argparse
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from utils.param_builder import build_gdc_params
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import shutil
from pandas.errors import EmptyDataError
from utils.permute_network import EdgeSwapGraph
import scipy
from scipy.optimize import least_squares
from runners.FDR_runner import run_FDR
import random
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import multiprocessing
from utils.daemon_multiprocessing import func_star
def calc_emp_pval(cur_rv, cur_dist):
    pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), cur_rv, side='left')


    return pos / float(np.size(cur_dist))



def main(algo_sample = None, dataset_sample = None,tsv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_{}_{}_md.tsv")):

    try:
        output_md = pd.read_csv(
            tsv_file_name.format(dataset_sample, algo_sample),
            sep='\t', index_col=0).dropna()
    except Exception, e:
        print e
        return 0, 0

    output_md = output_md.rename(columns={"filtered_pval": "hg_pval_max"})
    filtered_genes=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500]), ["GO name","hg_pval_max", "emp_pval_max", "passed_oob_permutation_test"]]


    print "total n_genes with pval:{}/{}".format(np.size(filtered_genes["hg_pval_max"].values), 7435)

    sorted_genes_hg = filtered_genes.sort_values(by=['hg_pval_max'], ascending=False)
    sig_genes_hg_pval=np.append(sorted_genes_hg["hg_pval_max"].values,np.zeros(7435-np.size(sorted_genes_hg["hg_pval_max"].values)))
    sig_genes_hg_pval = [10**(-x) for x in sig_genes_hg_pval]
    fdr_results = fdrcorrection0(sig_genes_hg_pval, alpha=0.05, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur])
    sig_hg_genes= sorted_genes_hg.iloc[:n_hg_true, :] if n_hg_true > 0 else 0
    if type(sig_hg_genes) == int:
        HG_CUTOFF = 0
        print "HG cutoff: {}, n=0".format(HG_CUTOFF)

    else:
        HG_CUTOFF = 10**(-sig_hg_genes.iloc[- 1]["hg_pval_max"])
        print "HG cutoff: {}, n={}".format(HG_CUTOFF, len(sig_hg_genes.index))


    sorted_genes_emp = filtered_genes.sort_values(by=['emp_pval_max'])
    sorted_genes_emp.loc[sorted_genes_emp['emp_pval_max']==0,'emp_pval_max']=1.0/5000
    sig_genes_emp_pval = sorted_genes_emp["emp_pval_max"].values
    fdr_results = fdrcorrection0(sig_genes_emp_pval, alpha=0.05, method='indep', is_sorted=False)
    n_emp_true = len([cur for cur in fdr_results[0] if cur])
    sig_emp_genes = sorted_genes_emp.iloc[:n_emp_true, :]
    EMP_CUTOFF = sig_emp_genes.iloc[- 1]["emp_pval_max"] if n_emp_true > 0 else 0
    print "EMP cutoff: {}, n={}".format(EMP_CUTOFF, len(sig_emp_genes.index))

    # genes_oob = filtered_genes.loc[filtered_genes["passed_oob_permutation_test"],:]
    # print "OOB genes: n={}".format(len(genes_oob.index))



    return sig_hg_genes, sig_emp_genes


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM") # Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50   TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM
    parser.add_argument('--prefix', dest='prefix', default="GE") # PASCAL_SUM   GE
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox,keypathwayminer_INES_GREEDY,hotnet2,my_netbox_td,dcem,dcem2,dcem3,dcem4") #
    omic_prefix=""
    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix

    terms_limit = 0
    df_rank_matrix=pd.DataFrame()
    df_matrix = pd.DataFrame()
    df_count_matrix = pd.DataFrame()
    df_all = pd.DataFrame()
    for cur_ds in datasets:
        df_ds=pd.DataFrame()
        for cur_alg in algos:
            sig_hg_genes, sig_emp_genes=main(cur_alg, cur_ds,tsv_file_name=os.path.join("/home/hag007/Desktop/aggregate{}_report/oob".format(omic_prefix), "emp_diff_modules_{}_{}_passed_oob.tsv")) # _gwas
            # venn2([set(sig_hg_genes.index), set(sig_emp_genes.index)], set_labels = ('HG', 'EMP'))
            # plt.title("EMP/HG ratio: {}".format(round(float(len(sig_emp_genes.index))/len(sig_hg_genes.index),3)))
            # plt.savefig("/home/hag007/Desktop/venn/venn_{}_{}.png".format(cur_ds, cur_alg))
            # plt.clf()
            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "algo"]=cur_alg
            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "dataset"]=cur_ds
            if type(sig_emp_genes)==int or type(sig_hg_genes)==int:
                df_ds.loc["{}_{}".format(cur_ds, cur_alg), "n_emp"] = 0
                df_ds.loc["{}_{}".format(cur_ds, cur_alg), "n_hg"] = 0
            else:
                df_ds.loc["{}_{}".format(cur_ds,cur_alg),"n_emp"]=np.sum(sig_emp_genes["passed_oob_permutation_test"].apply(lambda x: np.any(np.array(x[1:-1].split(', '),dtype=np.bool))).values)
                df_ds.loc["{}_{}".format(cur_ds,cur_alg), "n_hg"] = len(sig_hg_genes.index)
            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "ratio"] = round(float(df_ds.loc["{}_{}".format(cur_ds,cur_alg), "n_emp"])/max(df_ds.loc["{}_{}".format(cur_ds,cur_alg), "n_hg"],1),3)
            df_rank_matrix.loc[cur_alg, cur_ds]=df_ds.loc["{}_{}".format(cur_ds, cur_alg), "ratio"]
        df_ds["ratio_rank"]=df_ds["ratio"].rank(ascending=False)
        for cur_alg in algos:
            df_rank_matrix.loc[cur_alg, cur_ds] = df_ds.loc["{}_{}".format(cur_ds, cur_alg), "ratio_rank"]
            df_count_matrix.loc[cur_alg, cur_ds] = df_ds.loc["{}_{}".format(cur_ds, cur_alg), "n_emp"]
            df_matrix.loc[cur_alg, cur_ds] = df_ds.loc["{}_{}".format(cur_ds, cur_alg), "ratio"]
        df_all=pd.concat([df_all,df_ds])

    df_all.to_csv("/home/hag007/Desktop/aggregate{}_report/venn/summary.tsv".format(omic_prefix), sep='\t')
    df_rank_matrix.to_csv("/home/hag007/Desktop/aggregate{}_report/venn/rank_matrix.tsv".format(omic_prefix), sep='\t')
    df_matrix.to_csv("/home/hag007/Desktop/aggregate{}_report/venn/ratio_matrix.tsv".format(omic_prefix), sep='\t')
    df_count_matrix.to_csv("/home/hag007/Desktop/aggregate{}_report/venn/count_matrix.tsv".format(omic_prefix), sep='\t')