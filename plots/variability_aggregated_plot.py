import matplotlib.pyplot as plt
import pandas as pd
from fastsemsim.SemSim import *
from rpy2.robjects import pandas2ri

pandas2ri.activate()
import constants

from matplotlib.lines import Line2D

import seaborn as sns


def variability_plot(algos, base_folder, homogeneity_file_format, homogeneity_std_file_format, heterogeneity_file_format, suffix, cutoffs, ax, axs_violin, title=""):

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')

    heterogeneities = pd.DataFrame()
    homogeneities = pd.DataFrame()
    homogeneities_std=pd.DataFrame()
    empties = pd.DataFrame()


    for i_cutoff, cutoff in enumerate(cutoffs):
        df_homogeneity = pd.read_csv(os.path.join(base_folder, homogeneity_file_format.format(suffix, cutoff)),
                                     sep='\t', index_col=0)
        df_homogeneity_std = pd.read_csv(os.path.join(base_folder, homogeneity_file_format.format(suffix, cutoff)),
                                     sep='\t', index_col=0)
        df_heterogeneity = pd.read_csv(os.path.join(base_folder, heterogeneity_file_format.format(suffix, cutoff)),
                                       sep='\t', index_col=0)

        df_homogeneity[df_homogeneity.isna()] = 0
        df_homogeneity_std[df_homogeneity_std.isna()] = 0
        df_heterogeneity[df_heterogeneity.isna()] = 0

        homogeneities = pd.concat([homogeneities, df_homogeneity.mean(axis=1)], axis=1)
        homogeneities_std = pd.concat([homogeneities_std, df_homogeneity_std.mean(axis=1)], axis=1)
        heterogeneities = pd.concat([heterogeneities, df_heterogeneity.mean(axis=1)], axis=1)

        df_summary_agg=pd.DataFrame()
        for algo in set(df_homogeneity.index).intersection(constants.ALGOS):
            for dataset in df_homogeneity.columns:
                df_summary_agg=df_summary_agg.append({"algo":algo, "dataset": dataset, "value": df_homogeneity.loc[algo,dataset]},ignore_index=True)

        df_summary_agg=df_summary_agg.dropna(axis=0)
        my_order = df_summary_agg.groupby(by=["algo"])["value"].median().sort_values().index

        g = sns.violinplot(x="algo", y="value", data=df_summary_agg,
                           ax=axs_violin[i_cutoff], order=my_order,
                           palette={a: constants.COLORDICT[a] for a in my_order})
        g.set_xticklabels(g.get_xticklabels(), rotation=45)

    homogeneities = homogeneities.loc[algos]
    homogeneities_std = homogeneities_std.loc[algos]

    # if algos is None:
    #     algos=list(homogeneities.index)


    i = 0
    homogeneities.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"homogeneity_summary_{}.tsv".format(title)), sep='\t')
    for cur in set(homogeneities.index).intersection(constants.ALGOS):
        ax.plot(cutoffs, homogeneities.loc[cur], c=constants.COLORDICT[cur])

        # for x, y, std in zip(cutoffs, homogeneities.iloc[cur_i], homogeneities_std.iloc[cur_i]):
        #
        #     ax.errorbar(x-0.15+0.05*i, y, yerr=std, linestyle='None', marker='^', c=colorlist[cur_i])
        i += 1

    ax.set_xlabel("edge (similarity) cutoff", fontsize=22)
    ax.set_ylabel("median intra-module homogeneity", fontsize=22)
    ax.set_title(title, fontsize=22)

    ax.legend(handles=constants.PATCHES, loc=(0,1.1), ncol=3, fontsize=22, facecolor='#ffffff')


if __name__=="__main__":

    base_folder=constants.OUTPUT_GLOBAL_DIR
    homogeneity_file_format = "homogeneity_avg_matrix_{}_{}.tsv"
    homogeneity_std_file_format = "homogeneity_std_matrix_{}_{}.tsv"
    heterogeneity_file_format = "heterogeneity_avg_matrix_{}_{}.tsv"
    cutoffs = [1.0, 2.0, 3.0, 4.0] #1.0, 2.0, 3.0, 4.0, 5.0  , 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]  #


    fig,axs=plt.subplots(1,2, figsize=(22,10))
    fig_violin, axs_violin = plt.subplots(2, len(cutoffs), figsize=(4 * len(cutoffs) * 2, 10))


    suffix = "GE"
    algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY"]
    variability_plot(algos, base_folder, homogeneity_file_format, homogeneity_std_file_format, heterogeneity_file_format, suffix, cutoffs, axs[0], axs_violin[0], title="GE")
    suffix = "PASCAL_SUM"
    algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox"]
    variability_plot(algos, base_folder, homogeneity_file_format, homogeneity_std_file_format, heterogeneity_file_format, suffix, cutoffs, axs[1], axs_violin[1], title="GWAS")

    fig.text(0.01,0.97, "A:", weight='bold',fontsize=22)
    fig.text(0.5, 0.97, "B:", weight='bold',fontsize=22)
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_15.png"))

    fig_violin.text(0.0, 0.97, "A:", weight='bold', fontsize=22)
    fig_violin.text(0.0, 0.5, "B:", weight='bold', fontsize=22)
    fig_violin.tight_layout()
    fig_violin.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_15_violin.png"))

