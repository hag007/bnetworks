import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.SemSim.SetSemSim import SetSemSim
import matplotlib
matplotlib.use("Agg")

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star

import argparse

import math
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

import scipy.spatial.distance as ssd

METHOD='weighted'
def main():


    df_summary=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "cancer_type_go_distance.tsv"), sep='\t', index_col=0)
    df_summary-=1

    df_summary+=np.abs(np.nanmin(df_summary.values))+1
    df_summary = df_summary.where(~np.isnan(df_summary.values), other=1)
    for cur in df_summary.columns:
        df_summary.loc[cur,cur]=0

    distArray = ssd.squareform(df_summary.values )

    linked = linkage(distArray, method=METHOD, metric='euclidean')
    # algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "hotnet2", "keypathwayminer_INES_GREEDY"]
    datasets = ["TNFa_2", "HC12", "SHERA", "ROR_1", "SHEZH_1", "ERS_1", "IEM"]
    color_strategy=datasets
    my_color=["red", "green", "blue", "black", "brown", "orange", "purple"]
    # colors = [[i_y for i_y, y in enumerate(algos) if y in x][0] for x in list(df_summary.index)]
    plt.figure(figsize=(10, 7))
    ax = plt.gca()
    dendrogram(linked,
               orientation='top',
               labels=df_summary.index.values,
               distance_sort='descending',
               show_leaf_counts=True,
               leaf_rotation=90)

    # Apply the right color to each label

    xlbls = ax.get_xmajorticklabels()
    num = -1
    for lbl in xlbls:
        val = my_color[[i_y for i_y, y in enumerate(color_strategy) if y in lbl._text][0]]
        lbl.set_color(val)

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hierarchical_go_clustering.png"))








    algos=["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "hotnet2", "keypathwayminer_INES_GREEDY"]
    datasets=["TNFa_2", "HC12", "SHERA", "ROR_1", "SHEZH_1", "ERS_1", "IEM"]
    df_summary3d = np.empty((len(algos), len(algos), len(datasets)),dtype=np.object)
    for cur_col in df_summary.columns:
        for cur_index in df_summary.index:
            cur_ds_col = [a for a in datasets if a in cur_col][0]
            cur_algo_col = [a for a in algos if a in cur_col][0]
            cur_ds_index = [a for a in datasets if a in cur_index][0]
            cur_algo_index = [a for a in algos if a in cur_index][0]
            if cur_ds_col != cur_ds_index: continue

            df_summary3d[algos.index(cur_algo_index), algos.index(cur_algo_col)
            , datasets.index(cur_ds_index)]=df_summary.loc[cur_index, cur_col]

    df_summary3d[df_summary3d == None] = 5
    for cur in range(len(algos)):
        for cur_z in range(len(datasets)):
            df_summary3d[cur, cur,cur_z]=0

    df_summary_avg=np.nanmean(df_summary3d,axis=2)
    df_summary_max = np.nanmax(df_summary3d, axis=2)
    df_summary_min = np.nanmin(df_summary3d, axis=2)



    np_similarity = df_summary_avg # df_summary.values
    labels=algos # [x.split(".")[0] for x in df_summary.index.values]

    for i, cur_ds in enumerate(datasets):
        # plot_hierarchical_algos(df_summary3d[:,:,i], labels, cur_ds)
        pd.DataFrame(data=df_summary3d[:,:,i], index=algos, columns=algos).to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "similarity_matrix_{}.tsv".format(datasets[i])),sep='\t')

    plot_hierarchical_algos(df_summary_avg, labels, "avg")
    plot_hierarchical_algos(df_summary_min, labels, "min")
    plot_hierarchical_algos(df_summary_max, labels, "max")


def plot_hierarchical_algos(df_summary_avg, labels, cur_ds):
    distArray = ssd.squareform(df_summary_avg)
    linked = linkage(distArray, method=METHOD, metric='euclidean')
    plt.figure(figsize=(10, 7))
    dendrogram(linked,
               orientation='top',
               labels=labels,
               leaf_rotation=90,
               distance_sort='descending',
               show_leaf_counts=True)
    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hierarchical_go_clustering_{}.png".format(cur_ds)))


if __name__ == "__main__":
    main()
