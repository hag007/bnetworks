import json
import matplotlib
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
    cur_dist=np.array(cur_dist, np.float32)
    emp_pvals=[]
    if type(cur_rv)==str:
        hg_pvals=np.array(cur_rv[1:-1].split(", "), dtype=np.float32)
    else:
        hg_pvals=np.array([cur_rv], dtype=np.float32)

    for hg_pval in hg_pvals:
        
        pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), hg_pval, side='left')
        emp_pval = pos / float(np.size(cur_dist))
        emp_pvals.append(emp_pval)

    return str(emp_pvals)


def main(algo_sample = None, dataset_sample = None, n_dist_samples = 300, n_total_samples = None, shared_list=None, n_start_i=None, limit = 10000,
         md_path=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_modules_{}_{}_md.tsv"), \
         dist_path=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX","emp_diff_modules_{}_{}.tsv")):


    output_md = pd.read_csv(md_path.format(dataset_sample, algo_sample), sep='\t', index_col=0).dropna()
    # constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_{}_{}_md.tsv"
    output_md = output_md.rename(columns={"filtered_pval": "hg_pval"})
    n_genes_pvals=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500]), "hg_pval"]
    

    if n_genes_pvals.shape[0]==0:
       n_modules=0
       n_terms=0
       n_genes_pvals=np.array([])
    else:
       n_modules=str(n_genes_pvals.iloc[0]).count(",")+1
       n_terms=np.size(n_genes_pvals)
       n_genes_pvals=n_genes_pvals.values
    print "start reduction.."
    n_genes_pvals = [np.power([10 for a in range(x.count(",")+1)],-np.array(x[1:-1].split(", ")).astype(np.float)) if type(x)==str else [10**(-x)] for x in n_genes_pvals]
    max_genes_pvals = reduce(lambda a, x : np.append(a,np.min(x)), n_genes_pvals , np.array([]))
    n_genes_pvals = reduce(lambda a, x : np.append(a,x), n_genes_pvals , np.array([]))
    n_total_pvals=n_modules*n_terms

    print "total n_genes with pval:{}/{}".format(np.size(max_genes_pvals), 7435)
    max_genes_pvals=np.append(max_genes_pvals,np.ones(7435-np.size(max_genes_pvals)))
    print max_genes_pvals 
    fdr_results = fdrcorrection0(max_genes_pvals, alpha=0.05, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur == True])
    HG_CUTOFF=np.sort(max_genes_pvals)[n_hg_true-1]
    print "HG cutoff: {}, (n={})".format(HG_CUTOFF, n_hg_true)

    output_md = output_md.loc[np.logical_and.reduce(
        [output_md["n_genes"].values > 5, output_md["n_genes"].values < 500,
         output_md["hg_pval_max"].values > np.log10(HG_CUTOFF)]), :]

    output = pd.read_csv(dist_path.format(dataset_sample, algo_sample),
        sep='\t', index_col=0).dropna()
    output = output.rename(columns={"filtered_pval": "hg_pval"})
    output = output.loc[output_md.index.values, :]
    counter = 0
    emp_dists = []
    emp_pvals = []

    n_total_samples=n_total_samples if n_total_samples is not None else len(output.iloc[0].loc["dist_n_samples"][1:-1].split(", "))
    if n_start_i is None:
        np.random.seed(int(random.random()*1000))        
        i_choice=np.random.choice(n_total_samples, n_dist_samples, replace=False)
        i_dist=i_choice[:n_dist_samples]
    else:
        i_dist=np.arange(n_start_i, n_start_i+n_fist_samples)



    for index, cur in output.iterrows():
        if counter == limit: break
        # print "cur iteration in real data: {}/{}".format(counter, len(output.index))
        pval = np.array([float(x) for x in cur["dist_n_samples"][1:-1].split(", ")])[i_dist]
            
        emp_pvals.append(calc_emp_pval(cur["hg_pval"], pval))
        output_md.loc[index,'emp_pval']=emp_pvals[-1]
        emp_dists.append(pval)
        counter += 1

    mask_ids=output.index.values    
    emp_pvals = [np.array([x]) if type(x)!=str else np.array(x[1:-1].split(", ")).astype(np.float32)
                       for x in emp_pvals]
 
    if len(emp_pvals)==0:
        n_modules=0
        n_terms=0
    else:
        n_modules=emp_pvals[0].shape[0]
        n_terms=len(emp_pvals)

    max_emp_pvals = reduce(lambda a, x: np.append(a, np.min(x)), emp_pvals, np.array([]))
    emp_pvals = reduce(lambda a, x: np.append(a, x), emp_pvals, np.array([]))
    dist_name = "emp"
    df_dists = pd.DataFrame(index=output.index)
    df_dists["emp"] = pd.Series(emp_dists, index=output.index[:limit])
    max_emp_pvals = np.sort([x if x != 0 else 1.0 / n_dist_samples for x in max_emp_pvals])
    print "max emp pvals len: {}".format(len(max_emp_pvals))
    # print "\n".join([str(a) for a in max_emp_pvals])
    file(os.path.join(constants.OUTPUT_GLOBAL_DIR , "consistent_pval.txt"),'w+').write("\n".join([str(a) for a in max_emp_pvals]))
    print "min vals", 1.0 / n_dist_samples, np.min(list(max_emp_pvals) + [1])
    zero_bool=[x<=0.004 for x in max_emp_pvals]
    print "max_genes_pvals: {}".format(max_genes_pvals.shape[0])
    file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "consistent.txt"), 'w+').write("\n".join(str(list(max_emp_pvals))[1:-1].split(", "))) 
    fdr_bh_results = fdrcorrection0(max_emp_pvals, alpha=0.05, method='indep', is_sorted=False)[0]
    # print len(max_emp_pvals), fdr_bh_results, np.sum(fdr_bh_results), max_emp_pvals
    n_emp_true=np.sum(fdr_bh_results)
    print "n_emp_true: {}".format(n_emp_true)

    if n_emp_true==0: 
       BH_TH=0
    else: 
       BH_TH=np.sort(max_emp_pvals)[n_emp_true-1]
     
    mask_terms=np.array([a<=BH_TH for a in emp_pvals])
    go_ids_result=np.array([])
    go_names_result=np.array([])
    n_emp_true=0
    if len(mask_terms) > 0:
       mask_terms=np.array(mask_terms).reshape(-1, n_modules)
       go_ids_result=output.index.values[mask_terms.any(axis=1)]
       go_names_result=output["GO name"].values[mask_terms.any(axis=1)]
       n_emp_true =np.sum(mask_terms)


    print "BH cutoff: {} # true terms passed BH cutoff: {}".format(BH_TH, n_emp_true)
    if shared_list is not None:
        shared_list.append((BH_TH, n_emp_true, HG_CUTOFF, n_hg_true, go_ids_result, go_names_result, mask_ids, mask_terms))
    return BH_TH, n_emp_true, HG_CUTOFF, n_hg_true, go_ids_result, go_names_result, mask_ids, mask_terms


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--n_iteration', dest='n_iteration', default=1)
    parser.add_argument('--n_total_samples', help="n_total_samples", dest='n_total_samples', default=5000)
    parser.add_argument('--n_dist_samples', help="n_dist_samples", dest='n_dist_samples', default=5000)
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=3)

    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix
    n_iteration=int(args.n_iteration)
    n_total_samples=int(args.n_total_samples)
    n_dist_samples = int(args.n_dist_samples)
    pf=int(args.pf)

    terms_limit = 0
    sig_terms_summary=pd.DataFrame()
    full_report=pd.DataFrame()
    p = multiprocessing.Pool(pf)

    for cur_ds in datasets:
        for cur_alg in algos:

            l_emp_cutoff=[]
            l_n_emp_true=[]
            l_hg_cutoff=[]
            l_n_hg_true=[]
            l_mask_terms=[]
            l_mask_ids=[]
            
            go_ids_results=[]
            go_names_results=[]
            go_names_intersection=[]
            ys1=[]
            ys2=[]
            print "{}_{}".format(cur_ds, cur_alg)

            shared_list = multiprocessing.Manager().list()
            params=[]
            print "about to start {} intersections".format(n_iteration)
            for cur in range(n_iteration):
                # params.append([main, [cur_alg, cur_ds, n_dist_samples, n_total_samples, shared_list, ]])
                main(cur_alg, cur_ds, n_dist_samples, n_total_samples, shared_list)
                                     

            # p.map(func_star, params)            
           
            for BH_TH, n_emp_true, HG_CUTOFF, n_hg_true, go_ids_result, go_names_result, mask_ids, mask_terms in list(shared_list):

                # print BH_TH, n_emp_true, HG_CUTOFF, n_hg_true, go_ids_result, go_names_result, mask_ids, mask_terms

                l_emp_cutoff.append(str(BH_TH))
                l_n_emp_true.append(str(n_emp_true))
                l_hg_cutoff.append(str(HG_CUTOFF))
                l_n_hg_true.append(str(n_hg_true))
                l_mask_terms.append(mask_terms)
                l_mask_ids.append(mask_ids)


                if len(go_names_result) >terms_limit:
                    go_ids_results.append(go_ids_result)
                    go_names_results.append(go_names_result)
                    go_names_intersection = list(set(go_names_results[0]).intersection(*go_names_results))
                   
                    print "total intersection size : {}".format(len(go_names_intersection))
                    ys1.append(len(go_names_intersection))
                    ys2.append(len(go_names_result))

            plt.clf()

            go_ids_intersection = list(set(go_ids_results[0]).intersection(*go_ids_results)) if len(go_ids_results) > 0 else []

            output_md = pd.read_csv(
                os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX",
                             "emp_diff_modules_{}_{}_md.tsv".format(cur_ds, cur_alg)),
                sep='\t', index_col=0)

            # print output_md.shape, l_mask_terms[0].shape 
            output_md.loc[l_mask_ids[0],"passed_oob_permutation_test"]="[False]"
            output_md.loc[l_mask_ids[0],"passed_oob_permutation_test"]=[str(list(a)) for a in l_mask_terms[0]]
            # output_md.loc[~np.isin(output_md.index.values, np.array(go_ids_intersection)), "passed_oob_permutation_test"] = False
            print "emp_pval:\n", output_md['emp_pval'].iloc[0], type(output_md['emp_pval'].iloc[0]) 
            if type(output_md['emp_pval'].iloc[0])!=str:
                output_md['emp_pval_max']=1.0
            else:
                print "output_md :{}".format(np.min(np.array(output_md['emp_pval'].iloc[0][1:-1].split(", "),dtype=np.float)))
                output_md['emp_pval_max']=output_md['emp_pval'].apply(lambda a: np.min(np.array(a[1:-1].split(", "),dtype=np.float)) if type(a)==str else "" )

            output_md.to_csv(
                os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX",
                             "emp_diff_modules_{}_{}_passed_oob.tsv".format(cur_ds, cur_alg)),
                sep='\t')

            sig_terms_summary.loc[cur_alg,cur_ds]=len(go_names_intersection)
            modules_summary=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"{}_{}".format(prefix,cur_ds),cur_alg,"modules_summary.tsv"), sep='\t')
            full_report=full_report.append({"index" : "{}_{}".format(cur_alg, cur_ds), "dataset" : cur_ds, "algo" : cur_alg, "n_mutual_sig_terms" : len(go_names_intersection), "mean_sig_terms" : np.mean(ys2) if len(ys2) > 0 else 0, "n_genes" : modules_summary["#_genes"].sum() , "n_modules" : modules_summary["#_genes"].count(), "modules_size_mean" : modules_summary["#_genes"].mean(), "module_size_std" : modules_summary["#_genes"].sum()}, ignore_index=True)

            fig,ax = plt.subplots()
            plt.plot(np.arange(len(ys1)),ys1, label='# intersected terms')
            plt.plot(np.arange(len(ys2)), ys2, label='# sig. terms in permutation')
            plt.title("min # terms above threshold ({}):  ({}/{}). ratio: {}/{}".format(terms_limit, len(ys1), n_iteration, ys1[-1] if len(ys1) > 0 else 0,
                                                                                min([len(x) for x in go_names_results]) if len(go_names_results) > 0 else 0 ))
            ax.set_xlabel("iteration index")
            ax.set_ylabel("# terms")

            plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "intersection_function_{}_{}.png".format(cur_ds, cur_alg)))

            fdr_test_summary=pd.DataFrame()
            fdr_test_summary["emp_cutoff"]=l_emp_cutoff
            fdr_test_summary["n_emp_true"]=l_n_emp_true
            fdr_test_summary["hg_cutoff"]=l_hg_cutoff
            fdr_test_summary["n_hg_true"]=l_n_hg_true

            fdr_test_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "fdr_test_summary_{}_{}.tsv".format(cur_ds, cur_alg)), sep='\t')

    sig_terms_summary.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "sig_terms_summary.tsv"),
        sep='\t')

    rank_col=[]
    for cur_ds_i, cur_ds in enumerate(datasets):
        rank_col = np.append(rank_col,full_report["n_mutual_sig_terms"].iloc[cur_ds_i:cur_ds_i+len(algos)].rank(ascending=1))

    full_report["mutual_sig_terms_rank"]=rank_col
    full_report.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "full_report.tsv"),
        sep='\t')

