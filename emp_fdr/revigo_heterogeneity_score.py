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

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star

import argparse

import utils.go_hierarcies as go_hierarcies

import simplejson as json

def main(datasets, algos, pf=10, base_folder='/home/hag007/Desktop/aggregate_report/oob', file_format="emp_diff_modules_{}_{}_passed_oob.tsv"):

    df = pd.DataFrame()
    for cur_ds in datasets:
        print "cur ds: {}".format(cur_ds) 
        constants.update_dirs(DATASET_NAME_u=cur_ds)
        algo_go_sim_score = []
        total_num_genes = []
        algos_signals = []
        for i_algo, cur_algo in enumerate(algos):
            print "current cur_algo: {}".format(cur_algo)
            try:
                emp_before_filter = pd.read_csv(
                    os.path.join(base_folder,
                                 file_format.format(cur_ds, cur_algo)), sep='\t', index_col=0)

                emp_after_filter = pd.read_csv(os.path.join(base_folder,
                             file_format.format(cur_ds, cur_algo)[:-4] + "_reduced.tsv"), sep='\t', index_col=0)
            except:
                total_num_genes.append(0)
                algos_signals.append(0)
                algo_go_sim_score.append(1)
                continue

            emp_results_fdr=emp_before_filter.dropna().loc[emp_before_filter.dropna()["passed_oob_permutation_test"].apply(lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values,:]
            df.loc[cur_algo, cur_ds]=float(emp_after_filter.shape[0])/max(emp_results_fdr.shape[0],1)
            print "{}/{}".format(emp_after_filter.shape[0], emp_results_fdr.shape[0])
    df.to_csv(os.path.join("/home/hag007/Desktop/aggregate_report", "ds_go_heterogeneity_matrix.tsv"), sep='\t')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM") # "Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50"
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox,hotnet2,keypathwayminer_INES_GREEDY") #"jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td"
    parser.add_argument('--pf', dest='pf', default=3)
    args = parser.parse_args()

    prefix = args.prefix
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    pf=int(args.pf)
    print "test" 
    ds_summary=pd.DataFrame()
    for cur_ds in datasets:
        print "current dataset: {}".format(cur_ds)
        main(datasets=datasets, algos=algos, pf=pf)
