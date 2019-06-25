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
from utils.randomize_data import create_random_ds
from utils.randomize_data import permutation_output_exists
import argparse
from pandas.errors import EmptyDataError
from runners.datasets_multithread_runner import run_dataset
from utils.daemon_multiprocessing import MyPool, func_star
import shutil
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

def impute_dataset(dataset_name, network_file_name, imputation_fraction=0.5):
    if dataset_name.startswith("GE_"):
        data_file_path=os.path.join(constants.DATASETS_DIR, dataset_name,"data","ge.tsv")
    else:
        data_file_path = os.path.join(constants.DATASETS_DIR, dataset_name, "data", "score.tsv")

    df=pd.read_csv(data_file_path, index_col=0, sep='\t')

    network_genes = set(np.unique(pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_file_name), sep='\t').values))
    assay_genes=set(df.index.values)
    overlapped_genes=list(network_genes.intersection(assay_genes))
    genes_to_impute=np.random.choice(np.array(overlapped_genes), int(len(overlapped_genes)*imputation_fraction), replace=False)

    df = df.loc[genes_to_impute, :]
    df.to_csv(os.path.join(constants.DATASETS_DIR, dataset_name,"data","ge.tsv"), sep='\t')


def get_enriched_terms(algos, datasets, is_plot=False, empirical_th=None):
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

        fdr_results = fdrcorrection0(df_go['pval'].values, alpha=0.05, method='indep', is_sorted=False)[0]
        df_go = df_go[fdr_results]

        pval = -np.log10(df_go["pval"].values)
        if np.size(pval) == 0:
            pval = np.array([0])

        return pval, df_go, df_go_pvals


def recovery_iteration(prefix, dataset, cur, algo, network_file_name="dip.sif"):

    print "starting iteration: {}, {}, {}".format(prefix, dataset, cur)
    recovered_dataset_name="{}_{}_recovery_{}".format(prefix, dataset, cur)

    if not os.path.exists(os.path.join(constants.DATASETS_DIR, recovered_dataset_name)):

        shutil.copytree(os.path.join(constants.DATASETS_DIR, prefix+"_"+dataset),
                        os.path.join(constants.DATASETS_DIR, recovered_dataset_name))
        shutil.rmtree(os.path.join(constants.DATASETS_DIR, recovered_dataset_name,"output"))
        os.makedirs(os.path.join(constants.DATASETS_DIR, recovered_dataset_name, "output"))
        shutil.rmtree(os.path.join(constants.DATASETS_DIR, recovered_dataset_name,"cache"))
        os.makedirs(os.path.join(constants.DATASETS_DIR, recovered_dataset_name, "cache"))

        impute_dataset(recovered_dataset_name, network_file_name, imputation_fraction=0.5)

    permuted_network_file_name = network_file_name #   # _perm

    run_dataset(recovered_dataset_name, score_method=score_method,
                algos=[algo], network_file_name=permuted_network_file_name)
    cur_pval, df_terms, df_pval_terms = get_enriched_terms([algo], [recovered_dataset_name])

    df=pd.read_csv(os.path.join("/home/hag007/Desktop/aggregate_report/oob", "emp_diff_modules_{}_{}_passed_oob.tsv".format(dataset,algo)), sep='\t')
    full_sig_terms=df.loc[df["passed_oob_permutation_test"].dropna(axis=0).apply(
        lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :].sort_values(by=["hg_pval_max"], ascending=False)['GO id']

    recovery_fraction=len(set(full_sig_terms).intersection(df_terms["GO id"].values))/float(len(full_sig_terms))
    print "recovery_fraction: {}/{}={} ({}, {})".format(len(set(full_sig_terms).intersection(df_terms["GO id"].values)), float(len(full_sig_terms)), recovery_fraction, len(set(full_sig_terms)), len(set(df_terms["GO id"].values)))
    print "done iteration: {}, {}, {}".format(prefix, dataset, cur)
    return recovery_fraction


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="netbox,jactivemodules_greedy,my_netbox_td,jactivemodules_sa")
    parser.add_argument('--network', dest='network', default="dip.sif")
    parser.add_argument('--pf', help="parallelization_factor", dest='pf', default=1)
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=100)
    parser.add_argument('--override_permutations', help="takes max or all samples", dest='override_permutations', default="false")

    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix
    network_file_name = args.network
    parallelization_factor = int(args.pf)
    n_start=args.n_start
    n_end=args.n_end
    override_permutations=args.override_permutations.lower()=="true"


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
            recovery_fractions=[]
            prcs = []
            # p = MyPool(parallelization_factor)
            # params=[[recovery_iteration, [prefix, dataset, x, algo, network_file_name]] for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_output_exists(prefix, dataset, algo, x)]
            # p.map(func_star, params)

            for x in np.arange(n_start,n_end):
                recovery_fractions.append(recovery_iteration(prefix, dataset, x, algo, network_file_name))

            print "recovery_fractions for dataset {} and algo {}:\n{}".format(dataset,algo,recovery_fractions)
            open(os.path.join(constants.OUTPUT_GLOBAL_DIR, "recovery_results.txt"),'a+').write("recovery_fractions for dataset {} and algo {}:\n{}\n".format(dataset,algo,recovery_fractions))


