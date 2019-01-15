import sys
sys.path.insert(0, '../')
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from utils.go_pval_dist import create_random_ds
from multiprocessing import Process
import argparse
from pandas.errors import EmptyDataError




def plot_dist(pval ,algos_filter):



    fig, ax = plt.subplots(figsize=(12, 10))
    sns.distplot(pval, kde=False)
    plt.title("pval dist. -log10(pval) go terms. algo: {}".format("_".join(algos_filter)))
    real_n_term = np.size(pval)
    txt = "# terms={}".format(np.size(pval))
    fig.text(.5, .05, txt, ha='center')
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                             "go_pval_dist_{}_{}.png".format(prefix.lower(), "_".join(algos_filter))))
    plt.clf()

def calc_dist(algos, datasets,is_plot=False,empirical_th=None):
    for cur_algo in algos:
        algos_filter = cur_algo

        df_go = pd.DataFrame(columns=['qval', 'pval'])
        df_go_pvals = pd.DataFrame()
        df_go_pvals.index.name="GO id"
        for cur_ds in datasets:
            go_results = [os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, cur_module) for cur_algo in
                          os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds))
                          if os.path.isdir(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) and cur_algo in algos_filter for
                          cur_module in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) if
                          "separated_modules" in cur_module]

            for cur in go_results:
                try:
                    df_go = pd.concat((df_go, pd.read_csv(cur, sep='\t')))
                    df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']), axis=1)
                except EmptyDataError:
                    pass
        df_go_pvals[df_go_pvals.isna()]=1
        df_go = df_go[df_go['qval'] < 0.05]
        if empirical_th:
            df_go = df_go[df_go['pval'].apply(lambda x:-np.log10(x)) > empirical_th]

        # df_go = df_go[~df_go.index.duplicated(keep='first')]
        pval = -np.log10(df_go["pval"].values)
        if np.size(pval) == 0:
            pval = np.array([0])

        if is_plot:
            plot_dist(pval, algos_filter=algos+datasets)

        return pval, df_go, df_go_pvals


def empirical_dist_iteration(prefix, dataset, cur, algo):
    # from datasets_multithread_runner import run_dataset
    print "starting iteration: {}, {}, {}".format(prefix, dataset, cur)
    random_ds = create_random_ds(prefix, "{}_{}".format(prefix, dataset), cur, algo)
    permuted_network_file_name = "dip"  # _perm
    # if cur==0:
    #     permuted_network_file_name=create_permuted_network(network_file_name=network_file_name)
    run_dataset(random_ds, score_method=score_method,
                algos=[algo], network_file_name=permuted_network_file_name)
    cur_pval, df_terms, df_pval_terms = calc_dist([algo], [random_ds.format(prefix, dataset)])
    print "done iteration: {}, {}, {}".format(prefix, dataset, cur)
    return cur_pval, df_pval_terms


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SOC")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--network', dest='network', default="dip")
    parser.add_argument('--pf', help="parallelization_factor", dest='pf', default=10)
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n', default=100)
    # parser.add_argument('--max_dist', help="takes max or all samples", dest='max_dist', default="true")
    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix
    network_file_name = args.network
    parallelization_factor = args.pf
    n_start=args.n_start
    n_end=args.n_end
    # max_dist=args.max_dist.lower()=="true"


    summary = []
    for dataset in datasets:

        score_method = constants.PREDEFINED_SCORE
        if prefix == "GE":
            score_method = constants.DEG_EDGER
            if dataset.startswith("IE"):
                score_method = constants.DEG_T

        df_all_terms = pd.DataFrame()
        cur_real_ds= "{}_{}".format(prefix, dataset)

        for algo in algos:
            pval = np.array([])
            prcs = []
            for cur_parallel in range(n_start, n_end):
                prcs = []
                print "cur parallel: {}".format(cur_parallel)
                for cur in range(cur_parallel*parallelization_factor, cur_parallel*parallelization_factor+parallelization_factor):
                    prcs.append(Process(target=empirical_dist_iteration,
                            args=[prefix, dataset, cur, algo]))

                for cur in prcs:
                    cur.start()
                for cur in prcs:
                    cur.join()

