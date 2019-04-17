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
               "hotnet2": "hotnet2",
               "bionet": "bionet",
               "keypathwayminer_INES_GREEDY": "KPM"}


def dot_grid(df_measurements, grid_type, y_label):

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_facecolor('#fffde3')
    ax.grid(color='gray')



    colormap = cm.jet

    algos=df_measurements.index.values
    sns.set_palette("husl", algos.shape[0])

    df_new = pd.DataFrame()
    for i, cur_row in df_measurements.iterrows():
        for cur_entry in cur_row:
            df_new=df_new.append({"algo": algos_acronym[i], y_label : cur_entry}, ignore_index=True)

    ax = sns.stripplot(x='algo', y=y_label, data=df_new, jitter=True, size=10)

    colorlist = [sns.color_palette()[i] for i in
                 np.array(list(range(len(algos))))] #  / float(len(algos) - 1)]
    patches = [Line2D([0], [0], marker='o', color='gray', label=algos_acronym[a], markersize=12,
                      markerfacecolor=c, alpha=0.7) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]

    patches += [Line2D([0], [0], marker='D', color='gray', markersize=0)]
    patches += [Line2D([0], [0], marker='D', color='gray', label='mean', markersize=12, markerfacecolor='red', alpha=0.7)]
    patches += [Line2D([0], [0], marker='P', color='gray', label='median', markersize=12, markerfacecolor='red', alpha=0.7)]
    i=0
    for index, row in df_measurements.iterrows():
        # ax.scatter([i for a in np.arange(row.shape[0])], row.values, marker='o', c=colorlist[i], s=[400 for a in np.arange(row.shape[0])], alpha=0.7)
        ax.scatter([i], [row.values.mean()], marker='D', color='red', edgecolors='b', s=[200], alpha=0.7 ,zorder=4)
    #     # ax.errorbar([i], [row.values.mean()], yerr=np.std(row.values))
        ax.scatter([i], [np.median(row.values)], marker='P', color='red', edgecolors='b', s=[200], alpha=0.7, zorder=4)
        i+=1
    #
    #
    ax.set_xlabel("algos", fontdict={"size" : 12})
    ax.set_ylabel(y_label, fontdict={"size": 12})
    ax.legend(handles=patches, loc='upper left', fontsize=12)
    plt.xticks(np.arange(len(algos)),
               tuple([algos_acronym[x[3:]] if x.startswith("GE_") else algos_acronym[x] for x in algos]),
               rotation='vertical', size=12)

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_dot_{}.png".format(grid_type)))
    plt.clf()



def main():
    main_path = "/home/hag007/Desktop/aggregate_report/visual"
    df_measurements=pd.read_csv(os.path.join(main_path, "true_positive_counts.tsv"), sep='\t', index_col=0)
    dot_grid(df_measurements, "counts", "n_terms")

    main_path = "/home/hag007/Desktop/aggregate_report/visual"
    df_measurements = pd.read_csv(os.path.join(main_path, "true_positive_ratio.tsv"), sep='\t', index_col=0)
    dot_grid(df_measurements, "ratio", "EHR")

    main_path = "/home/hag007/Desktop/aggregate_report/visual"
    df_measurements = pd.read_csv(os.path.join(main_path, "empirical_terms_variability.tsv"), sep='\t', index_col=0)
    i_zeros = df_measurements==0
    df_measurements=df_measurements+3
    df_measurements[i_zeros]=0
    dot_grid(df_measurements, "variability", "heterogeneity")


if __name__=="__main__":
    main()