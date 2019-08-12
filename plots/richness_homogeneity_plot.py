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
               "my_netbox_td": "netbox_td",
               "dcem": "dcem"
               }


def plot_richness_homogeneity(df_richness, df_homogeneity, grid_type, ax=None):
    df_richness=df_richness.loc[np.sort(df_richness.index), np.sort(df_richness.columns)]
    df_homogeneity=df_homogeneity.loc[np.sort(df_homogeneity.index), np.sort(df_homogeneity.columns)]
    algos = list(df_richness.index)
    sns.set_palette("husl", len(algos))

    if ax is None:
        fig, ax = plt.subplots(figsize=(13, 13))

    ax.set_facecolor('#fffde3')
    ax.grid(color='gray')

    i_x = np.arange(len(df_richness.columns))
    i_y = np.arange(len(df_richness.index))
    size = np.array([df_richness.iloc[b, a] for a in i_x for b in i_y])
    x = np.array([df_homogeneity.iloc[b, a] for a in i_x for b in i_y])
    y = np.array([df_richness.iloc[b, a] for a in i_x for b in i_y])
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
    ax.set_xlabel("homogeneity", fontdict={"size":22})
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)
    ax.set_ylabel("richness", fontdict={"size":22})



def dot_grid(df_measurements, grid_type, y_label, filter_zeros=False, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    ax.set_facecolor('#fffde3')
    ax.grid(color='gray')

    algos=df_measurements.index.values
    sns.set_palette("husl", algos.shape[0])

    df_new = pd.DataFrame()
    for i, cur_row in df_measurements.iterrows():
        for cur_entry in cur_row:
            df_new=df_new.append({"algo": algos_acronym[i], y_label : cur_entry}, ignore_index=True)

    ax = sns.stripplot(x='algo', y=y_label, data=df_new, jitter=True, size=15, ax=ax)
    ax.set_title("{}".format(y_label), fontdict={"size":22})
    colorlist = [sns.color_palette()[i] for i in
                 np.array(list(range(len(algos))))] #  / float(len(algos) - 1)]
    patches = [Line2D([0], [0], marker='o', color='gray', label=algos_acronym[a], markersize=12,
                      markerfacecolor=c, alpha=0.7) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]

    patches += [Line2D([0], [0], marker='D', color='gray', markersize=0)]
    patches += [Line2D([0], [0], marker='D', color='gray', label='mean', markersize=12, markerfacecolor='red', alpha=0.7)]
    patches += [Line2D([0], [0], marker='P', color='gray', label='median', markersize=12, markerfacecolor='red', alpha=0.7)]
    i=0
    for index, row in df_measurements.iterrows():

        if filter_zeros:
            mn = row.values[row.values != 0].mean()
            mdn = np.median(row.values[row.values != 0])
        else:
            mn = row.values.mean()
            mdn = np.median(row.values)


        ax.scatter([i], [mn], marker='D', color='red', edgecolors='b', s=[200], alpha=0.7 ,zorder=4)


        ax.scatter([i], [mdn], marker='P', color='red', edgecolors='b', s=[200], alpha=0.7, zorder=4)
        i+=1

    ax.set_xlabel("algos", fontdict={"size" : 22})
    ax.set_ylabel(y_label, fontdict={"size": 22})
    ax.legend(handles=patches, loc='upper left', fontsize=22)
    ax.set_xticklabels(tuple([algos_acronym[x[3:]] if x.startswith("GE_") else algos_acronym[x] for x in algos]), rotation=45)

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_dot_{}.png".format(grid_type)))



def main():

    cutoffs=[1.0,2.0,3.0,4.0,5.0]

    for cutoff in cutoffs:


        fig, axs = plt.subplots(1,2, figsize=(20, 10))

        suffix="GE_{}".format(cutoff)
        df_richness=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "solution_richness_matrix_{}.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
        df_homogeneity = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "homogeneity_avg_matrix_{}.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
        plot_richness_homogeneity(df_richness, df_homogeneity, "size", ax=axs[0])

        suffix="PASCAL_SUM_{}".format(cutoff)
        df_richness=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "solution_richness_matrix_{}.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
        df_homogeneity = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "homogeneity_avg_matrix_{}.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
        plot_richness_homogeneity(df_richness, df_homogeneity, "size", ax=axs[1])

        fig.text(0.01,0.97, "A:", weight='bold',fontsize=22)
        fig.text(0.5, 0.97, "B:", weight='bold',fontsize=22)
        fig.tight_layout()
        fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_16_{}.png".format(cutoff)))


if __name__=="__main__":
    main()