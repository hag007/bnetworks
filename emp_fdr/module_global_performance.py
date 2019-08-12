import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as ml_colors
from matplotlib.lines import Line2D

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import argparse

from utils.go_similarity import calc_similarity_matrix, calc_intra_similarity
from variability_between_go_sets import cc_to_v_ratio_analysis, community_performance_analysis

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50") #  TNFa_2,HC12,ROR_1,ERS_1,IEM,SHERA,SHEZH_1 "Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM")  # GE  PASCAL_SUM
    parser.add_argument('--algos', dest='algos', default="dcem,jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td,hotnet2,keypathwayminer_INES_GREEDY") # my_netbox_td ,hotnet2,keypathwayminer_INES_GREEDY
    # parser.add_argument('--module_indices', dest='module_indices',
    #                     default="0,1,2")  # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    parser.add_argument('--pf', dest='pf', default=3)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX"))
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--cutoffs', dest='cutoffs', default="0.0, 1.0,2.0,3.0,4.0") # ,5.0,6.0,7.0,8.0,9.0,10.0,11.0

    args = parser.parse_args()

    prefix = args.prefix
    base_folder = args.base_folder
    sim_method= args.sim_method
    file_format=args.file_format
    # semsim = SetSemSim(ontology, ac, TSS=sim_method, MSS="BMA")
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    # module_indices=args.module_indices.split(",")
    algos = args.algos.split(",")
    cutoffs = np.array(args.cutoffs.split(','),dtype=float)
    pf=int(args.pf)
    print "test"
    h_scores = pd.DataFrame()
    df_full_data = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_full_data_{}.tsv".format(prefix)), index_col=0, sep='\t')

    for cutoff in cutoffs:
        suffix=str(cutoff)
        m_homogeneities = []
        m_heterogeneities = []
        m_cs = []
        m_txts = []
        m_scores={}
        for cur_alg in algos:
            m_scores[cur_alg]=[]
        for cur_ds in datasets:
            for cur_alg in algos:

                all_m_scores=[]
                ####

                try:

                    emp_results = pd.read_csv(
                        os.path.join(base_folder,
                                     file_format.format(cur_ds, cur_alg)), sep='\t', index_col=0)

                except:
                    print "could not find {}".format(os.path.join(base_folder,
                                                                  file_format.format(cur_ds, cur_alg)), cur_ds, cur_alg)
                    continue

                emp_results = emp_results.sort_values(by='emp_rank')
                emp_results_fdr = emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(
                    lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :]

                #####


                module_indices=sorted([int(a.split("_")[-1]) for a in df_full_data.loc[(df_full_data['algo']== cur_alg).values & (df_full_data['dataset']== cur_ds).values & (df_full_data['EHR']>0.2).values].index.values])
                for cur_module_index_0 in module_indices:
                    set_0=df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index_0), "tp"]
                    if not type(set_0) is str:
                        set_0=[]
                    else:
                        set_0=[a.split(": ")[0] for a in df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index_0), "tp"].split("\n")]

                    score_p=None
                    for cur_module_index_1 in module_indices:

                        if cur_module_index_0 > cur_module_index_1:
                            continue

                        set_1 = df_full_data.loc["{}_{}_module_{}".format(cur_ds, cur_alg, cur_module_index_1), "tp"]
                        if not type(set_1) is str:
                            set_1 = []
                        else:
                            set_1 = [a.split(": ")[0] for a in df_full_data.loc["{}_{}_module_{}".format(cur_ds, cur_alg, cur_module_index_1), "tp"].split("\n")]

                        print "current dataset: {} {} {} {}".format(cur_ds, cur_alg, cur_module_index_0, cur_module_index_1)
                        cache_file=os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_emp_{}_{}_{}.npy".format(cur_ds, cur_alg, sim_method))


                        if set_1==set_0:
                            score=np.mean(community_performance_analysis(cache_file, intervals=np.linspace(cutoff,cutoff,1), modules=[set_0, set_1], dataset=cur_ds, algo=cur_alg, base_folder=base_folder, file_format=file_format))
                        else:
                            score=0
                            # if score_p==None:
                            #     score_p=np.mean(community_performance_analysis(cache_file, modules=[set_0, set_1]))
                            # score=score_p

                        print "score: {}".format(score)

                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0, cur_module_index_1), "homogeneity_score"]=score
                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0,
                                                                  cur_module_index_1), "dataset"] = cur_ds
                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0,
                                                                  cur_module_index_1), "algo"] = cur_alg
                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0,
                                                                  cur_module_index_1), "module_0"] = cur_module_index_0
                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0,
                                                                  cur_module_index_1), "module_1"] = cur_module_index_1


                    if len(module_indices) > 1:

                        m_i0=cur_module_index_0
                        m_i1=module_indices[module_indices.index(cur_module_index_0)-1]
                        if m_i0 > m_i1:
                            tmp=m_i0
                            m_i0=m_i1
                            m_i1=tmp

                        m_homogeneities.append(h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0, cur_module_index_0
                                                                  ), "homogeneity_score"])

                        m_heterogeneities.append(h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, m_i0,m_i1), "homogeneity_score"])
                        m_cs.append(algos.index(cur_alg)/float(len(algos)))
                        m_txts.append(cur_ds)
                        all_m_scores.append((1 if m_homogeneities[-1] > 1 and m_heterogeneities[-1] < 1 else 0))

                if len(all_m_scores) != 0:
                    m_scores[cur_alg].append(np.mean(all_m_scores))

        for k,v in m_scores.iteritems():
            print k, "n={}".format(len(v)), np.mean(v)
        fig, ax = plt.subplots()
        ax.scatter(m_heterogeneities, m_homogeneities, c=m_cs, cmap='jet')
        for i, data in enumerate(zip(m_heterogeneities, m_homogeneities)):
            x, y = data
            ax.annotate(m_txts[i], (x, y))
        ax.set_xlabel("m_heterogeneity")
        ax.set_ylabel("m_homogeneity")

        cmap = plt.cm.jet
        colorlist = [ml_colors.rgb2hex(cmap(a / float(np.size(algos) - 1))) for a in np.arange(np.size(algos))]
        patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                          markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
        ax.legend(handles=list(reversed(patches)), loc='lower right', framealpha=0.5)

        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"m_h_plot_{}_{}.png".format(prefix,suffix)))


        h_scores.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "h_scores_{}_{}.tsv".format(prefix,suffix)), sep='\t')

        h_scores=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "h_scores_{}_{}.tsv".format(prefix,suffix)), sep='\t',index_col=0)

        std_homogeneities=[]
        averaged_homogeneities=[]
        cs=[]
        txts=[]
        fig, ax = plt.subplots()
        df_homogeneity_avg=pd.DataFrame()
        df_homogeneity_std= pd.DataFrame()
        for cur_ds in datasets:
            for cur_alg in algos:
                ds_filtered_scores=h_scores.loc[(h_scores['algo']== cur_alg).values & (h_scores['dataset']== cur_ds).values]
                std_homogeneities.append(ds_filtered_scores[ds_filtered_scores["module_0"] == ds_filtered_scores["module_1"]]['homogeneity_score'].values.std())
                averaged_homogeneities.append(ds_filtered_scores[ds_filtered_scores["module_0"] == ds_filtered_scores["module_1"]]['homogeneity_score'].values.mean())
                df_homogeneity_avg.loc[cur_alg, cur_ds]=averaged_homogeneities[-1]
                df_homogeneity_std.loc[cur_alg, cur_ds] = std_homogeneities[-1]
                cs.append(algos.index(cur_alg)/float(len(algos)))
                txts.append(cur_ds)
                # labels.append()

        df_homogeneity_avg.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "homogeneity_avg_matrix_{}_{}.tsv".format(prefix,suffix)), sep='\t')
        df_homogeneity_std.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "homogeneity_std_matrix_{}_{}.tsv".format(prefix,suffix)), sep='\t')

        averaged_homogeneities=np.array(averaged_homogeneities)
        std_homogeneities= np.array(std_homogeneities)
        cs = np.array(cs)
        txts=np.array(txts)
        # cs=cs[~np.isnan(averaged_heterogeneities)]
        # averaged_homogeneities=averaged_homogeneities[~np.isnan(averaged_heterogeneities)]
        # txts=txts[~np.isnan(averaged_heterogeneities)]
        # averaged_heterogeneities=averaged_heterogeneities[~np.isnan(averaged_heterogeneities)]
        DEFAULT_VAL=0
        ax.scatter([a if not np.isnan(a) else DEFAULT_VAL for a in averaged_homogeneities], [a if not np.isnan(a) else DEFAULT_VAL for a in averaged_homogeneities], c=cs, cmap='jet')
        for i,data in enumerate(zip([a if not np.isnan(a) else DEFAULT_VAL for a in averaged_homogeneities], [a if not np.isnan(a) else DEFAULT_VAL for a in averaged_homogeneities])):
            x, y = data
            ax.annotate(txts[i], (x,y))
        ax.set_xlabel("avg_homogeneity (var={})".format(round(np.var([filter(lambda a : not np.isnan(a), std_homogeneities)]),3)))
        ax.set_ylabel("std_homogeneity (var={})".format(round(np.var([filter(lambda a : not np.isnan(a), averaged_homogeneities)]),3)))

        cmap = plt.cm.jet
        colorlist = [ml_colors.rgb2hex(cmap(a / float(np.size(algos) - 1))) for a in np.arange(np.size(algos))]
        patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                          markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
        ax.legend(handles=list(reversed(patches)), loc='lower right', framealpha=0.5)


        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"h_plot_{}_{}.png".format(prefix,suffix)))



