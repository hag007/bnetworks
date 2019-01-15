import numpy as np
from scipy.stats import friedmanchisquare
import pandas as pd
import os
import constants
import scikit_posthocs as sp
from scipy.stats import wilcoxon

if __name__ == "__main__":

    df=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "ds_go_rank_summary.tsv"), sep='\t', index_col=0)

    df.groupby(by="algo")["sig_terms", "sig_terms_rank", "variability", "variability_rank"].agg(['sum', 'mean', 'std'])\
    .to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "go_rank_sum.tsv"), sep='\t')



    dataset_ranks=[]
    datasets = df[df["algo"] == "bionet"]["dataset"].values
    for cur_dataset in datasets:
        dataset_ranks.append(df.loc[df["dataset"].values == cur_dataset, "sig_terms"].values)

    print dataset_ranks

    algo_ranks=[]
    algos = df[df["dataset"] == "GE_TNFa_2"]["algo"].values
    n_algos = np.size(algos)
    for cur_algo in algos:
        algo_ranks.append(df.loc[df["algo"].values==cur_algo, "sig_terms"].values)


    fried=friedmanchisquare(*algo_ranks)
    posthoc_conover=sp.posthoc_conover_friedman(dataset_ranks, p_adjust="fdr_bh")
    # posthoc_nemenyi = sp.posthoc_nemenyi_friedman(dataset_ranks) # , p_adjust="fdr_bh"
    # posthoc_miller = sp.posthoc_miller_friedman(dataset_ranks) # , p_adjust="fdr_bh"
    # posthoc_siegel = sp.posthoc_siegel_friedman(dataset_ranks, p_adjust="fdr_bh")
    # posthoc_siegel = sp.posthoc_siegel_friedman(dataset_ranks, p_adjust="fdr_bh")
    posthoc_dunn = sp.posthoc_dunn(algo_ranks, p_adjust="fdr_bh")
    manwhitney_posthoc=sp.posthoc_mannwhitney(algo_ranks)
    wilcoxon_posthoc = sp.posthoc_wilcoxon(algo_ranks)

    df_pw_signed_test=pd.DataFrame(columns=np.arange(n_algos), index=np.arange(n_algos))
    for i in np.arange(n_algos):
        for j in np.arange(n_algos):
            df_pw_signed_test.loc[i,j]=-1 if i==j else wilcoxon(algo_ranks[i],algo_ranks[j])[1]

    # print df[df["dataset"]=="GE_TNFa_2"]["algo"].values
    print fried
    # print posthoc_conover
    # print posthoc_nemenyi
    # print posthoc_miller
    # print posthoc_siegel
    # print posthoc_dunn
    print manwhitney_posthoc
    print wilcoxon_posthoc
    # print df_pw_signed_test


    # ranks = []
    # for cur_dataset in np.unique(df["dataset"].values):
    #     ranks.append(df.loc[df["dataset"].values == cur_dataset, "variability"].values)
    #
    # fried = friedmanchisquare(*ranks)
    # posthoc = sp.posthoc_conover_friedman(ranks, p_adjust="fdr_bh")
    # print df[df["dataset"] == "GE_TNFa_2"]["algo"].values
    # print fried
    # print posthoc



