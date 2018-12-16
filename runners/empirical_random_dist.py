import json
from matplotlib import style
from pandas._libs.parsers import k

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
import utils.goids2gonames as goids2gonames
import shutil
from datasets_multithread_runner import run_dataset
from utils.go_pval_dist import create_random_ds
from utils.go_pval_dist import create_permuted_network
from multiprocessing import Process

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


def empirical_dist_iteration(prefix, dataset, cur):
    print "starting iteration: {}, {}, {}".format(prefix, dataset, cur)
    random_ds = create_random_ds(prefix, "{}_{}".format(prefix, dataset), cur)
    permuted_network_file_name = "dip"  # _perm
    # if cur==0:
    #     permuted_network_file_name=create_permuted_network(network_file_name=network_file_name)
    run_dataset(random_ds, score_method=score_method,
                algos=[algo], network_file_name=permuted_network_file_name)
    cur_pval, df_terms, df_pval_terms = calc_dist([algo], [random_ds.format(prefix, dataset)])
    print "done iteration: {}, {}, {}".format(prefix, dataset, cur)
    return cur_pval, df_pval_terms


if __name__ == "__main__":
    # prefix="GE"
    # datasets = ["{}_TNFa_2".format(prefix)]
    # algos = np.array(["jactivemodules_greedy"])
    #
    # calc_dist(algos, datasets)
    network_file_name="dip"
    prefix="GE"
    datasets = ["SOC"] # , "IEM", "IEN", "HC12", "MCF7_2", "TNFa_2"   alzheimers, schizophrenia
    algos = ["jactivemodules_greedy"] # "bionet"

    summary = []
    for dataset in datasets:


        score_method = constants.PREDEFINED_SCORE
        if prefix == "GE":
            score_method = constants.DEG_EDGER
            if dataset.startswith("IE"):
                score_method = constants.DEG_T

        total_n_terms=[]
        n_terms_filtered_in=[]
        n_terms_2_oom=[]
        summary = []
        diff_values=np.array([0])
        df_all_terms = pd.DataFrame()
        cur_real_ds= "{}_{}".format(prefix, dataset) # "{}_random_{}_{}".format(prefix, dataset, cur_real_from_rand)

        for algo in algos:
            pval = np.array([])
            random_ds = "{}_random_{}".format(prefix, dataset)
            prcs = []
            for cur in range(105, 110):
                # if cur==cur_real_from_rand: continue
                prcs.append(Process(target=empirical_dist_iteration,
                        args=[prefix, dataset, cur]))
            for cur in prcs:
                cur.start()
            for cur in prcs:
                cur.join()
