import pandas as pd

from fastsemsim.SemSim import *

import matplotlib.pyplot as plt

from rpy2.robjects import pandas2ri
pandas2ri.activate()
from mpl_toolkits.mplot3d import Axes3D
import constants

from scipy.cluster import hierarchy
import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform


import utils.go
import utils.go_hierarcies
import math
import random
import matplotlib.cm as cm

from matplotlib.lines import Line2D
import matplotlib.colors as ml_colors

from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

import seaborn as sns



ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.01
SIM_TH= 0.4


algos_acronym={"jactivemodules_greedy":"jAM_greedy",
               "jactivemodules_sa": "jAM_SA",
               "netbox": "netbox",
               "keypathwayminer_INES_GREEDY": "KPM",
               "hotnet2": "hotnet2",
               "bionet": "bionet",
               "dcem": "dcem",
               "my_netbox_td": "netbox_td"
               }


def plot_scatter_with_size(df_size, df_ratio, grid_type, ax=None):
    df_size=df_size.loc[np.sort(df_size.index), np.sort(df_size.columns)]
    df_ratio=df_ratio.loc[np.sort(df_ratio.index), np.sort(df_ratio.columns)]
    algos = list(df_size.index)
    sns.set_palette("husl", len(algos))

    if ax is None:
        fig, ax = plt.subplots(figsize=(13, 13))

    ax.set_facecolor('#fffde3')
    ax.grid(color='gray')

    i_x = np.arange(len(df_size.columns))
    i_y = np.arange(len(df_size.index))
    size = np.array([df_size.iloc[b, a] for a in i_x for b in i_y])
    x = np.array([df_ratio.iloc[b, a] for a in i_x for b in i_y])
    y = np.array([df_size.iloc[b, a] for a in i_x for b in i_y])
    labels=np.array([sns.color_palette()[b] for a in i_x for b in i_y])
    sc = ax.scatter(x, y, s=200, c=labels, cmap='jet') # size/float(max(size))*2000+20

    # sorted_list = sorted([[x[i], y[i]] for i in range(len(x))], reverse=True)
    # pareto_front = [sorted_list[0]]
    # for pair in sorted_list[1:]:
    #     if True:
    #         if pair[1] >= pareto_front[-1][1]:
    #             pareto_front.append(pair)
    #     else:
    #         if pair[1] <= pareto_front[-1][1]:
    #             pareto_front.append(pair)
    #
    # pf_X = [pair[0] for pair in pareto_front]
    # pf_Y = [pair[1] for pair in pareto_front]
    # ax.plot(pf_X, pf_Y)


    colormap = cm.jet

    colorlist = [sns.color_palette()[i] for i in
                 np.array(list(range(len(algos))))]
    patches = [Line2D([0], [0], marker='o', markersize=12, color='gray', label=a,
                      markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
    ax.legend(handles=patches, loc='upper left', prop={'size': 22})
    ax.margins(0.03, 0.03)
    ax.set_xlabel("EHR", fontdict={"size":22})
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)
    ax.set_ylabel("term count", fontdict={"size":22})



def dot_grid(df_measurements, grid_type, y_label, filter_zeros=False, ax=None, algos=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    ax.set_facecolor('#fffde3')
    ax.grid(color='gray')

    if algos is None:
        algos=algos_acronym.keys()

    algos = sorted(list(algos))

    algo_acrs=[algos_acronym[a] for a in algos]

    colorlist = [sns.color_palette("hls", n_colors=len(algos))[i] for i in
                 np.array(list(range(len(algos))))]
    patches = [Line2D([0], [0], marker='o', color='gray', label=algos_acronym[a], markersize=12,
                      markerfacecolor=c, alpha=0.7) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]

    df_new = pd.DataFrame()
    for i, cur_row in df_measurements.iterrows():
        for cur_entry in cur_row:
            df_new=df_new.append({"algo": algos_acronym[i], y_label : cur_entry}, ignore_index=True)

    df_new=df_new.dropna(axis=0)
    my_order = df_new.groupby(by=["algo"])[y_label].mean().sort_values().index
    ax = sns.stripplot(x='algo', y=y_label, data=df_new, jitter=True, size=15, ax=ax, order=my_order, palette={a: colorlist[algo_acrs.index(a)] for a in my_order})
    ax.set_title("{}".format(y_label), fontdict={"size":22})

    patches += [Line2D([0], [0], marker='D', color='gray', markersize=0)]
    patches += [Line2D([0], [0], marker='D', color='gray', label='mean', markersize=12, markerfacecolor='red', alpha=0.7)]
    patches += [Line2D([0], [0], marker='P', color='gray', label='median', markersize=12, markerfacecolor='red', alpha=0.7)]
    i=0
    for index, row in df_measurements.iterrows():

        if algos_acronym[index] not in my_order:
            continue

        if filter_zeros:
            mn = row.values[row.values != 0].mean()
            mdn = np.median(row.values[row.values != 0])
        else:
            mn = row.values.mean()
            mdn = np.median(row.values)


        ax.scatter(list(my_order).index(algos_acronym[index]), [mn], marker='D', color='red', edgecolors='b', s=[200], alpha=0.7 ,zorder=4)


        ax.scatter(list(my_order).index(algos_acronym[index]), [mdn], marker='P', color='red', edgecolors='b', s=[200], alpha=0.7, zorder=4)
        i+=1

    ax.set_xlabel(ax.get_xlabel(), fontdict={"size" : 22})
    ax.set_ylabel(ax.get_ylabel(), fontdict={"size": 22})
    ax.legend(handles=patches, loc='upper left', fontsize=22)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_dot_{}.png".format(grid_type)))



def main():

    fig_4, axs_4 = plt.subplots(2,2, figsize=(20, 20))
    fig_5, axs_5 = plt.subplots(1,2, figsize=(20, 10))

    main_path = "/home/hag007/Desktop/aggregate_report/venn"
    df_measurements_counts=pd.read_csv(os.path.join(main_path, "count_matrix.tsv"), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
    # df_measurements_counts=df_measurements_counts.drop(labels=["SHERA"], axis=1)
    dot_grid(df_measurements_counts, "counts", "term count", ax=axs_4[0][0])


    df_measurements_ratio = pd.read_csv(os.path.join(main_path, "ratio_matrix.tsv"), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
    dot_grid(df_measurements_ratio, "ratio", "EHR", ax=axs_4[1][0])

    plot_scatter_with_size(df_measurements_counts, df_measurements_ratio, "size", ax=axs_5[0])

    main_path = "/home/hag007/Desktop/aggregate_gwas_report/venn"
    df_measurements_counts=pd.read_csv(os.path.join(main_path, "count_matrix.tsv"), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
    # df_measurements_counts=df_measurements_counts.drop(labels=["SHERA"], axis=1)
    dot_grid(df_measurements_counts, "counts", "term count", ax=axs_4[0][1])


    df_measurements_ratio = pd.read_csv(os.path.join(main_path, "ratio_matrix.tsv"), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
    dot_grid(df_measurements_ratio, "ratio", "EHR", ax=axs_4[1][1])

    plot_scatter_with_size(df_measurements_counts, df_measurements_ratio, "size", ax=axs_5[1])


    fig_4.text(0.01,0.98, "A:", weight='bold',fontsize=22)
    fig_4.text(0.5, 0.98, "B:", weight='bold',fontsize=22)
    fig_4.text(0.01, 0.49, "C:", weight='bold',fontsize=22)
    fig_4.text(0.5, 0.49, "D:", weight='bold',fontsize=22)
    fig_4.tight_layout()
    fig_4.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_4.png"))

    fig_5.text(0.01, 0.97, "A:", weight='bold', fontsize=22)
    fig_5.text(0.5, 0.97, "B:", weight='bold', fontsize=22)
    fig_5.tight_layout()
    fig_5.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_5.png"))


    # df_counts = pd.read_csv(os.path.join(main_path, "empirical_terms_counts.tsv"), sep='\t', index_col=0)
    # plot_scatter_with_size(df_measurements_counts, df_measurements_ratio, "size", ax=axs[1][0])

if __name__=="__main__":
    main()