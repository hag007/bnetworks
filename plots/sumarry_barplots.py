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


from grid_plot_criteria import plot_scatter_with_size

ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.01
SIM_TH= 0.4


def barplot(df_summary_statistic, algos, ax):
    for i, dataset in enumerate(algos):
        barplot(df_summary_statistic, ax)

def barplots(df_summary_statistic, axs):

    algos = np.sort(np.unique(df_summary_statistic['Algorithm'].values))
    datasets = np.sort(df_summary_statistic['dataset'].values)


    for i, dataset in enumerate(datasets):
        barplot(df_summary_statistic, algos, axs[i])




def main():

    fig, axs = plt.subplots(2,2, figsize=(20, 20))

    main_path = constants.OUTPUT_GLOBAL_DIR
    prefix="GE"
    algos=[]
    datasets = []
    df_statistic=pd.read_csv(os.path.join(main_path, "summary_statistic_{}.tsv".format(prefix)), sep='\t', index_col=0)
    # df_measurements_counts=df_measurements_counts.drop(labels=["SHERA"], axis=1)
    barplots(df_statistic, ax=axs)

    prefix="PASCAL_SUM"
    algos=[]
    datasets = []
    df_statistic=pd.read_csv(os.path.join(main_path, "summary_statistic_{}.tsv".format(prefix)), sep='\t', index_col=0)
    # df_measurements_counts=df_measurements_counts.drop(labels=["SHERA"], axis=1)
    barplots(df_statistic, ax=axs)

    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_18_{}.png".format("prefix")))

if __name__=="__main__":
    main()