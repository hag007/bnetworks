#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import sys
sys.path.insert(0, '../')

import pandas as pd
from pandas.errors import EmptyDataError

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import utils.go_hierarcies

import matplotlib.cm as cm
import matplotlib.colors as ml_colors

import multiprocessing

from utils.daemon_multiprocessing import func_star
ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = False

if TERMS_SIMILARITY_TO_NUM_OF_TERMS:
    ontology_type = 'GeneOntology'
    ignore_parameters = {'ignore': {}}
    source_type = 'obo'
    source = os.path.join(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

    print "\n######################"
    print "# Loading ontology... #"
    print "######################\n"

    ontology = ontologies.load(source=source, source_type=source_type, ontology_type=ontology_type,
                               parameters=ignore_parameters)

    print "\n######################"
    print "# Loading Annotation Corpus... #"
    print "######################\n"
    ac = AnnotationCorpus.AnnotationCorpus(ontology)
    ac.parse(os.path.join(constants.GO_DIR, "goa_human.gaf"), "gaf-2.0")
    ac.isConsistent()

    print "\n#################################"
    print "# Annotation corpus successfully loaded."
    print "#################################\n"

    semsim = GSESAMESemSim(ontology, ac)  # maxSemSim(ontology, ac) #
#################
if ENABLE_GO_GRAPH:
    dict_result, go2geneids, geneids2go, entrez2ensembl = utils.go_hierarcies.build_hierarcy(
        roots=['GO:0008150', 'GO:0005575', 'GO:0003674'])


#################

def calc_similarity(mat_adj, i_x, i_y, x, y):
    key="{}_{}".format(i_x,i_y)
    key_inv="{}_{}".format(i_y,i_x)
    if mat_adj[key] != -2: return
    mat_adj[key] = semsim.SemSim(x, y)  # , ResnikSemSim(ontology,ac))
    print mat_adj[key]
    if np.isnan(mat_adj[key]):
        mat_adj[key] = -1
    mat_adj[key_inv] = mat_adj[key]

def main(dataset_name, algos = ["keypathwayminer_INES_GREEDY", "netbox", "hotnet2", "jactivemodules_greedy", "bionet",
             "jactivemodules_sa"]):
    global df_go_metadata, all_hg_score_modules, df_hg_output
    colormap = cm.rainbow
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    GO_RANK_CUTOFF = 150
    ##########################

    all_homogeneity = []
    all_separability = []
    agg_homogeneity = []
    agg_separability = []
    algo_go_sim_score = []
    colors = []
    df_all_hg_pval = pd.DataFrame()
    df_go_metadata = pd.DataFrame()
    all_hg_score_labels = []
    all_hg_score_modules = []
    labels_by_sample = []
    total_num_genes = []
    algos_signals = []
    modules_signals = []
    df_all_hg_qval = pd.DataFrame()
    df_module2best_rows = []
    df_module2avg_rows = []
    df_algo2best_rows = []
    df_algo2avg_rows = []
    #
    # n_terms_summary = pd.read_csv(
    #     os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "n_terms_summary.tsv"),
    #     sep='\t', index_col=0)

    
    for i_algo, algo in enumerate(algos):
        df_all_hg_qval = pd.DataFrame()
        print "current algo: {}".format(algo)
        try:
            emp_results = pd.read_csv(
                os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr",
                             "emp_diff_{}_{}_passed_oob.tsv".format(dataset_name[dataset_name.index("_") + 1:], algo)), sep='\t', index_col=0)
            # emp_results = pd.read_csv(
            #     os.path.join("/home/hag007/Desktop/fdr_terms/fdr_005_i_1000/terms/emp_diff_{}_{}_passed_oob.tsv"
            #                  .format(dataset_name[dataset_name.index("_") + 1:], algo)),
            #     sep='\t', index_col=0)
        except:
            total_num_genes.append(0)
            algos_signals.append(0)
            algo_go_sim_score.append(1)
            continue

        emp_results=emp_results.sort_values(by='emp_rank')
        emp_results_fdr=emp_results.loc[emp_results["passed_oob_permutation_test"].values,:]["GO name"]

        homogeneity = []
        separability = []
        algos_signals.append(len(emp_results_fdr.index))
        all_go_terms = emp_results_fdr.index.values

        try:
            total_num_genes.append(pd.read_csv(
                os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, algo, "all_modules_general.tsv"),
                sep="\t")["total_num_genes"][0])
        except:
            total_num_genes.append(0)

        algo2best_go_ratio = 0
        algo2avg_go_ratio = 0
        module2best_go_ratio=[]
        module2avg_go_ratio = []

        # df_modules_summary = pd.read_csv(
        #     os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, algo, "modules_summary.tsv"),
        #     sep='\t')



        # if RATIO_TO_GO_TERM:
        #     df_algo2best_rows.append(
        #         {'name': '{}_total_avg'.format(algo), 'score': algo2best_go_ratio / max(i, 1)})
        #     df_algo2avg_rows.append({'name': '{}_total_avg'.format(algo), 'score': algo2avg_go_ratio / max(i, 1)})


        # adj = np.ones((len(all_go_terms), len(all_go_terms))) * (-2)
        manager=multiprocessing.Manager()
        adj=manager.dict()
        for x in range(len(all_go_terms)): 
            for y in range(len(all_go_terms)):
                adj["{}_{}".format(x,y)]=-2        
         
        if TERMS_SIMILARITY_TO_NUM_OF_TERMS:
            params=[]
            for i_x, x in enumerate(all_go_terms):
                # print "calc distance between terms {}/ {}".format(i_x, len(all_go_terms))
                for i_y, y in enumerate(all_go_terms):
                    params.append([calc_similarity, [adj, i_x, i_y, x, y]])
     
            p = multiprocessing.Pool(10)
            p.map(func_star,params)
            
            adj_sum=sum([adj["{}_{}".format(x,y)] for x in range(len(all_go_terms)) for y in range(len(all_go_terms)) if adj["{}_{}".format(x,y)]!=-1])
            adj_count=float(len([adj["{}_{}".format(x,y)] for x in range(len(all_go_terms)) for y in range(len(all_go_terms)) if adj["{}_{}".format(x,y)]!=-1]))
            print "adj_sum: ",adj_sum
            print "adj_count: ",adj_count
            algo_go_sim_score.append(1-adj_sum /adj_count)
    
            p.close()
            p.join()
            
        
         
    if RATIO_TO_GO_TERM:
        pd.DataFrame(df_module2best_rows).set_index('name').to_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "GO_terms_per_module_ratio_best.tsv"),
            sep='\t')
        pd.DataFrame(df_module2avg_rows).set_index('name').to_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "GO_terms_per_module_ratio_avg.tsv"),
            sep='\t')
        pd.DataFrame(df_algo2best_rows).set_index('name').to_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "GO_terms_ratio_best.tsv"), sep='\t')
        pd.DataFrame(df_algo2avg_rows).set_index('name').to_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "GO_terms_ratio_avg.tsv"), sep='\t')

    if TERMS_SIMILARITY_TO_NUM_OF_TERMS:

        colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                     np.array(list(range(len(algos)))) / float(len(algos) - 1)]

        fig, ax = plt.subplots(figsize=(10, 10))
        df_summary = pd.DataFrame()
        for h, s, c, a, gene_size in zip(algos_signals, algo_go_sim_score, colorlist, algos,
                                         total_num_genes):  # [0 for x in range(len(algo_go_sim_score))]
            print (h, s)
            ax.scatter(h, s, s=(50 + 2000 * (float(gene_size) / (1+np.max(total_num_genes)))),
                       c=c, cmap='jet', label=a)
            df_series=pd.Series({"algo": a, "dataset": dataset_name, "sig_terms": h,
                       "sig_terms_rank": pd.Series(np.array(algos_signals)).rank(ascending=0).values[
                           np.where(np.array(algos_signals) == h)[0][0]], "variability": s,
                       "variability_rank": pd.Series(np.array(algo_go_sim_score)).rank(ascending=0).values[
                           np.where((np.array(algo_go_sim_score)) == s)[0][0]], "n_genes": gene_size})
            df_series.name = "{}_{}".format(dataset_name, algo)
            df_summary=df_summary.append(df_series)
            colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                         np.array(list(range(len(algos)))) / float(len(algos) - 1)]
            patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                              markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
            ax.set_xlabel("# GO terms (-log10(qval)) above threshold")
            ax.set_ylabel("GO terms variability")
            ax.legend(handles=patches)
            ax.grid(True)
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                                 "hs_plot_terms_signal_algo_{}.png".format(constants.DATASET_NAME)))
        return df_summary


if __name__ == "__main__":

    go_ratio_ds_summary = pd.DataFrame()
    # datasets = [name for name in os.listdir(constants.DATASETS_DIR) if
    #             os.path.isdir(os.path.join(constants.DATASETS_DIR, name)) and name.startswith("GWAS_")]
    datasets=["GE_IEM", "GE_TNFa_2",  "GE_HC12", "GE_SHERA", "GE_ROR_1", "GE_SHEZH_1", "GE_ERS_1"]
    ds_summary=pd.DataFrame()
    for cur_ds in datasets:
        print "current dataset: {}".format(cur_ds)
        # ["jactivemodules_sa", "jactivemodules_greedy", "hotnet2", "bionet"]
        df_results=main(dataset_name=cur_ds, algos =["jactivemodules_sa", "jactivemodules_greedy", "hotnet2", "bionet"])
        ds_summary=pd.concat([ds_summary, df_results])
    ds_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "ds_go_rank_summary.tsv"), sep='\t')


