import pandas as pd

from fastsemsim.SemSim import *
from infra import *

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


algos_acronym={"jactivemodules_greedy":"jAM_greedy",
               "jactivemodules_sa": "jAM_SA",
               "netbox": "netbox",
               "my_netbox_td": "my_netbox_td",
               "bionet": "bionet",
               "hotnet2": "hotnet2",
               "keypathwayminer_INES_GREEDY": "KPM"
}


def plot_itegrated_recovery(ss_ratios, base_folder, average_file_format, auc_file_format, p_file_format, r_file_format, f1_file_format, empty_file_format, suffix, axs):
      #

    ps = pd.DataFrame()
    rs = pd.DataFrame()
    f1s = pd.DataFrame()
    prs = pd.DataFrame()
    empties = pd.DataFrame()

    for cur_ss_ratio in ss_ratios:
        n_iterations = len(
            pd.read_csv(os.path.join(base_folder, average_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                        index_col=0)["precisions"][0][1:-1].split(", "))
        print "n_teration for ss_ratio={}: {}".format(cur_ss_ratio, n_iterations)
        df_cur_pr = pd.read_csv(os.path.join(base_folder, auc_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                                index_col=0).mean(axis=1)
        df_cur_p = pd.read_csv(os.path.join(base_folder, p_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                               index_col=0).mean(axis=1)
        df_cur_r = pd.read_csv(os.path.join(base_folder, r_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                               index_col=0).mean(axis=1)
        df_cur_f1 = pd.read_csv(os.path.join(base_folder, f1_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                                index_col=0).mean(axis=1)
        df_cur_empty = pd.read_csv(os.path.join(base_folder, empty_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                                   index_col=0).mean(axis=1) / float(n_iterations)

        prs = pd.concat([prs, df_cur_pr], axis=1)
        ps = pd.concat([ps, df_cur_p], axis=1)
        rs = pd.concat([rs, df_cur_r], axis=1)
        f1s = pd.concat([f1s, df_cur_f1], axis=1)
        empties = pd.concat([empties, df_cur_empty], axis=1)

    prs = prs.loc[np.sort(prs.index.values)]
    ps = ps.loc[np.sort(prs.index.values)]
    rs = rs.loc[np.sort(prs.index.values)]
    f1s = f1s.loc[np.sort(prs.index.values)]
    empties = empties.loc[np.sort(prs.index.values)]

    colorlist = [sns.color_palette()[i] for i in
                 np.array(list(range(prs.shape[0])))]  # / float(len(algos) - 1)]
    patches = [Line2D([0], [0], marker='o', color='gray', label=algos_acronym[a], markersize=12,
                      markerfacecolor=c, alpha=0.7) for i, a, c in
               zip(list(range(prs.shape[0])), prs.index.values, colorlist)]

    patches_0 = patches + \
                [Line2D([0], [0], linestyle='-', color='black', label='PR-AUC', markersize=12, markerfacecolor='gray',
               alpha=0.7)] + \
                [Line2D([0], [0], linestyle='dotted', color='black', label='fraction of\nnon-empty solutions', markersize=12,
               markerfacecolor='gray', alpha=0.7)]

    patches_1 = patches + [
        Line2D([0], [0], linestyle='-', linewidth=3.0, color='black', label='f1', markersize=12, markerfacecolor='gray',
               alpha=0.7)]+ \
                [Line2D([0], [0], linestyle='dashed', color='black', label='precision', markersize=12, markerfacecolor='gray',
               alpha=0.7)] + \
                [Line2D([0], [0], linestyle='dotted', color='black', label='recall', markersize=12,
                                  markerfacecolor='gray', alpha=0.7)]

    ss_ratios = 1 - np.array(ss_ratios)
    for cur_i in np.arange(prs.shape[0]):
        if prs.index[cur_i] == "my_netbox_td": continue
        axs[0].plot(ss_ratios, prs.iloc[cur_i], c=colorlist[cur_i])
        axs[0].plot(ss_ratios, empties.iloc[cur_i], c=colorlist[cur_i], linestyle='dotted')
        # plt.plot(ss_ratios, np.multiply(empties.iloc[cur_i] / 100.0, prs.iloc[cur_i]) , c=colorlist[cur_i], linestyle='dashed')
        axs[1].plot(ss_ratios, ps.iloc[cur_i], c=colorlist[cur_i], linestyle='dashed')
        axs[1].plot(ss_ratios, rs.iloc[cur_i], c=colorlist[cur_i], linestyle='dotted')
        axs[1].plot(ss_ratios, f1s.iloc[cur_i], c=colorlist[cur_i], linestyle='-', linewidth=3.0)

    axs[0].set_xlabel("subsample fraction", fontsize=22)
    axs[0].set_ylabel("PR-AUC/\nfration of non-empty iterations", fontsize=22)  #

    axs[1].set_xlabel("subsample fraction", fontsize=22)
    axs[1].set_ylabel("precision/\nrecall/\nf1", fontsize=22)

    axs[0].legend(handles=patches_0, loc='lower right', prop={"size": 20})
    axs[1].legend(handles=patches_1, loc='lower right', prop={"size": 20})



if __name__=="__main__":

    base_folder="/media/hag007/Data/bnet/output/emp_fdr/MAX"
    auc_file_format = "pr_auc_recovery_summary_{}_{}.tsv"
    p_file_format = "recovery_results_{}_{}_matrix_p.tsv"
    r_file_format = "recovery_results_{}_{}_matrix_r.tsv"
    f1_file_format = "recovery_results_{}_{}_matrix_f1.tsv"
    empty_file_format = "recovery_results_{}_{}_matrix_empty.tsv"
    average_file_format="recovery_results_{}_{}.tsv"
    ss_ratios = [0.1, 0.2, 0.3, 0.4]

    fig,axs=plt.subplots(2,2,figsize=(20,16))
    suffix = "GE_100"
    plot_itegrated_recovery(ss_ratios, base_folder, average_file_format, auc_file_format, p_file_format, r_file_format, f1_file_format, empty_file_format, suffix, axs=axs[:,0])
    suffix = "PASCAL_SUM_100"
    plot_itegrated_recovery(ss_ratios, base_folder, average_file_format, auc_file_format, p_file_format, r_file_format,
                            f1_file_format, empty_file_format, suffix, axs=axs[:,1])


    plt.tight_layout()
    plt.figtext(0.01, 0.97, "A:", weight='bold', fontsize=22)
    plt.figtext(0.01, 0.5, "B:", weight='bold', fontsize=22)
    plt.figtext(0.5, 0.97, "C:", weight='bold', fontsize=22)
    plt.figtext(0.5, 0.5, "D:", weight='bold', fontsize=22)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_14.png"))




