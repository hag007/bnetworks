import sys

sys.path.insert(0, '../')

import argparse

import pandas as pd
import numpy as np
import shutil

from pandas.errors import EmptyDataError

import constants
import os

import utils.goids2gonames as goids2gonames

import multiprocessing
from utils.daemon_multiprocessing import func_star
from utils.randomize_data import get_permutation_name

from runners.datasets_multithread_runner import run_dataset


# def calc_dist1(algos, datasets, shared_list=None, is_max=True):
#     print "here!!!"
#     shared_list.append(pd.DataFrame())

def calc_dist(algos, datasets, shared_list=None, is_max=True):
    print "here!!"
    try:

        for cur_algo in algos:
            algos_filter = cur_algo

            df_go_pvals = pd.DataFrame()
            df_go_pvals.index.name = "GO id"
            for cur_ds in datasets:
                print "fetch permutation {}".format(cur_ds)
                n_modules = len(
                    pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, "modules_summary.tsv"),
                                sep='\t').index)
                go_results = [os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, cur_module) for cur_algo in
                              os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds))
                              if os.path.isdir(
                        os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) and cur_algo in algos_filter for
                              cur_module in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) if
                              "separated_modules" in cur_module and int(cur_module.split("_")[1]) < n_modules]

                for cur in go_results:
                    try:
                        df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']),
                                                axis=1)

                        df_go_pvals[df_go_pvals.isna()] = 1
                        if is_max:
                            df_go_pvals = df_go_pvals.min(axis=1).to_frame()

                    except EmptyDataError, e:
                        print e
                        pass
                if len(go_results) == 0:
                    df_go_pvals = pd.DataFrame(data=np.array([[1]]), index=["GO:0008150"])

            if not is_max:
                df_go_pvals[df_go_pvals.isna()] = 1

            if shared_list is not None:
                shared_list.append(df_go_pvals)
                print "done aggregate {} permutations".format(len(shared_list))
            return df_go_pvals
    except Exception, e:
        print Exception, e
        pass


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="Breast_Cancer2.G50")
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM")
    parser.add_argument('--algos', dest='algos', default="my_netbox_td")
    parser.add_argument('--network', dest='network', default="dip.sif")
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_end', default=1000)
    parser.add_argument('--network_permutations', dest='network_perm', default="false")
    parser.add_argument('--max_dist', help="takes max or all samples", dest='max_dist', default="true")
    parser.add_argument('--pval_measurement', help="how to calc pval. can be emp or beta", dest='pval_measurement',
                        default="emp")
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=1)
    parser.add_argument('--recalc_true_modules', dest='recalc_true_modules', default="false")

    args = parser.parse_args()

    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix
    network_file_name = args.network
    n_start = int(args.n_start)
    n_end = int(args.n_end)
    max_dist = args.max_dist.lower() == "true"
    network_perm = args.network_perm.lower() == "true"
    recalc_true_modules = args.recalc_true_modules.lower() == "true"
    pval_measurement = args.pval_measurement
    pf = args.pf

    for dataset in datasets:

        score_method = constants.PREDEFINED_SCORE
        if prefix == "GE":
            score_method = constants.DEG_EDGER
            if dataset.startswith("IE"):
                score_method = constants.DEG_T

        for algo in algos:
            try:

                df_all_terms = pd.DataFrame()

                ## empirical pvalues
                df_real_agg_pval = calc_dist([algo], ["{}_{}".format(prefix, dataset)], is_max=False)
                df_real_agg_list_pval = df_real_agg_pval.apply(lambda x: str(list(-np.log10(x))), axis=1).to_frame()
                df_real_agg_list_pval.index = df_real_agg_pval.index
                df_counter = 0

                df_results = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX",
                                                      "emp_diff_modules_{}_{}.tsv".format(dataset, algo)),
                                         sep='\t', index_col=0)

                df_results['hg_pval'] = df_real_agg_list_pval.iloc[:, 0]

                df_results.to_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_diff_modules_{}_{}.tsv".format(dataset, algo)),
                    sep='\t', index_label="GO id")

                if args.max_dist:
                    df_results.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX",
                                                   "emp_diff_modules_{}_{}.tsv".format(dataset, algo)), sep='\t',
                                      index_label="GO id")

                print "permutation shape: {}".format(df_all_terms.shape)
            except Exception, e:
                print Exception, e
                pass
