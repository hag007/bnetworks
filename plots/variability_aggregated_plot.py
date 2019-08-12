import matplotlib.pyplot as plt
import pandas as pd
from fastsemsim.SemSim import *
from rpy2.robjects import pandas2ri

pandas2ri.activate()
import constants

from matplotlib.lines import Line2D

import seaborn as sns




algos_acronym={"jactivemodules_greedy":"jAM_greedy",
               "jactivemodules_sa": "jAM_SA",
               "netbox": "netbox",
               "my_netbox_td": "my_netbox_td",
               "bionet": "bionet",
               "hotnet2": "hotnet2",
               "keypathwayminer_INES_GREEDY": "KPM",
               "dcem":"dcem"
}


def variability_plot(algos, base_folder, homogeneity_file_format, homogeneity_std_file_format, heterogeneity_file_format, suffix, cutoffs, ax, axs_violin):
    heterogeneities = pd.DataFrame()
    homogeneities = pd.DataFrame()
    homogeneities_std=pd.DataFrame()
    empties = pd.DataFrame()

    algos = sorted(list(algos))
    colorlist = [sns.color_palette("hls", n_colors=len(algos))[i] for i in
                 np.array(list(range(len(algos))))]

    patches = [Line2D([0], [0], marker='o', color='gray', label=algos_acronym[a], markersize=12,
                      markerfacecolor=c, alpha=0.7) for i, a, c in
               zip(list(range(len(algos))), algos, colorlist)]


    for i_cutoff, cutoff in enumerate(cutoffs):
        df_homogeneity = pd.read_csv(os.path.join(base_folder, homogeneity_file_format.format(suffix, cutoff)),
                                     sep='\t', index_col=0)
        df_homogeneity_std = pd.read_csv(os.path.join(base_folder, homogeneity_file_format.format(suffix, cutoff)),
                                     sep='\t', index_col=0)
        df_heterogeneity = pd.read_csv(os.path.join(base_folder, heterogeneity_file_format.format(suffix, cutoff)),
                                       sep='\t', index_col=0)

        homogeneities = pd.concat([homogeneities, df_homogeneity.mean(axis=1)], axis=1)
        homogeneities_std = pd.concat([homogeneities_std, df_homogeneity_std.mean(axis=1)], axis=1)
        heterogeneities = pd.concat([heterogeneities, df_heterogeneity.mean(axis=1)], axis=1)

        df_summary_agg=pd.DataFrame()
        for algo in df_homogeneity.index:
            for dataset in df_homogeneity.columns:
                df_summary_agg=df_summary_agg.append({"algo":algo, "dataset": dataset, "value": df_homogeneity.loc[algo,dataset]},ignore_index=True)

        df_summary_agg=df_summary_agg.dropna(axis=0)
        my_order = df_summary_agg.groupby(by=["algo"])["value"].mean().sort_values().index

        g = sns.violinplot(x="algo", y="value", data=df_summary_agg,
                           ax=axs_violin[i_cutoff], order=my_order,
                           palette={a: colorlist[algos.index(a)] for a in my_order})
        g.set_xticklabels(g.get_xticklabels(), rotation=45)

    homogeneities = homogeneities.loc[np.sort(homogeneities.index.values)]
    homogeneities_std = homogeneities_std.loc[np.sort(homogeneities.index.values)]

    if algos is None:
        algos=list(homogeneities.index)


    i = 0
    for cur_i in np.arange(homogeneities.shape[0]):
        ax.plot(cutoffs, homogeneities.iloc[cur_i], c=colorlist[algos.index(homogeneities.index[cur_i])])

        # for x, y, std in zip(cutoffs, homogeneities.iloc[cur_i], homogeneities_std.iloc[cur_i]):
        #
        #     ax.errorbar(x-0.15+0.05*i, y, yerr=std, linestyle='None', marker='^', c=colorlist[cur_i])
        i += 1

    ax.set_xlabel("edge (similarity) cutoff", fontsize=22)
    ax.set_ylabel("homogeneity", fontsize=22)
    ax.legend(handles=patches, loc='upper right')


if __name__=="__main__":

    base_folder=constants.OUTPUT_GLOBAL_DIR
    homogeneity_file_format = "homogeneity_avg_matrix_{}_{}.tsv"
    homogeneity_std_file_format = "homogeneity_std_matrix_{}_{}.tsv"
    heterogeneity_file_format = "heterogeneity_avg_matrix_{}_{}.tsv"
    cutoffs = [0.0,1.0, 2.0, 3.0, 4.0] # , 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]  #


    fig,axs=plt.subplots(1,2, figsize=(20,10))
    fig_violin, axs_violin = plt.subplots(2, len(cutoffs), figsize=(4 * len(cutoffs) * 2, 10))

    algos = ["dcem", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY",
             "hotnet2", "my_netbox_td"]
    suffix = "GE"
    variability_plot(algos, base_folder, homogeneity_file_format, homogeneity_std_file_format, heterogeneity_file_format, suffix, cutoffs, axs[0], axs_violin[0])
    suffix = "PASCAL_SUM"
    variability_plot(algos, base_folder, homogeneity_file_format, homogeneity_std_file_format, heterogeneity_file_format, suffix, cutoffs, axs[1], axs_violin[1])

    fig.text(0.01,0.97, "A:", weight='bold',fontsize=22)
    fig.text(0.5, 0.97, "B:", weight='bold',fontsize=22)
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_15.png"))

    fig_violin.text(0.0, 0.97, "A:", weight='bold', fontsize=22)
    fig_violin.text(0.0, 0.5, "B:", weight='bold', fontsize=22)
    fig_violin.tight_layout()
    fig_violin.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_15_violin.png"))

