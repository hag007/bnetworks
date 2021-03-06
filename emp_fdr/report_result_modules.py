import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.SemSim.SetSemSim import SetSemSim

import matplotlib
matplotlib.use("Agg")

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star

import argparse

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

semsim = SetSemSim(ontology, ac, TSS="Resnik", MSS="BMA")

def calc_similarity(mat_adj, i_x, i_y, x, y):
    key="{}_{}".format(i_x,i_y)
    key_inv="{}_{}".format(i_y,i_x)
    if mat_adj[key] != -200: return
    if(x==y):
       mat_adj[key]=-100
       return 

    mat_adj[key] = semsim.SemSim(x, y)

    # if mat_adj[key]<0:
    # print x,y,mat_adj[key]

    if np.isnan(mat_adj[key]):
        mat_adj[key] = -100
    mat_adj[key_inv] = mat_adj[key]

def main(datasets, algos, prefix, pf=10):


    if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores")):
        try:
            os.makedirs(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores"))
        except Exception, e:
            print "error while creating ds_2_alg_scores folder: {}".format(e)

    for cur_ds in datasets:
        print "cur ds: {}".format(cur_ds) 
        constants.update_dirs(DATASET_NAME_u=cur_ds)
        algo_go_sim_score = []
        total_num_genes = []
        algos_signals = []

        for i_algo, cur_algo in enumerate(algos):
            print "current cur_algo: {}".format(cur_algo)
            try:
                emp_results = pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", 
                                 "emp_diff_modules_{}_{}_passed_oob.tsv".format(cur_ds[len(prefix) + 1:], cur_algo)), sep='\t', index_col=0)

            except:
                print "file {} wasn't found".format("emp_diff_modules_{}_{}_passed_oob.tsv".format(cur_ds[len(prefix) + 1:], cur_algo))
                total_num_genes.append(0)
                algos_signals.append(0)
                algo_go_sim_score.append(1)
                file(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores", "{}_{}_{}".format(cur_ds,cur_algo, "n_sig.txt")), 'w+').write(str(0))
                file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores", "{}_{}_{}".format(cur_ds, cur_algo, "var.txt")), 'w+').write(str(0))
                continue

            emp_results=emp_results.sort_values(by='emp_rank')
            emp_results_fdr=emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values,:]["GO name"]            
            

            algos_signals.append(len(emp_results_fdr.index))
            all_go_terms = emp_results_fdr.index.values

            try:
                total_num_genes.append(pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, "all_modules_general.tsv"),
                    sep="\t")["total_num_genes"][0])
            except:
                total_num_genes.append(0)
                print "cannot read file {}".format(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, "all_modules_general.tsv"))

            adj_sum, adj_count = calc_intra_similarity(all_go_terms, pf)
            print "about to save files..."
            file(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores", "{}_{}_{}".format(cur_ds,cur_algo, "n_sig.txt")), 'w+').write(str(len(all_go_terms)))
            file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores", "{}_{}_{}".format(cur_ds, cur_algo, "var.txt")), 'w+').write(str( - adj_sum / adj_count) if adj_count>0 else str(0))

            # manager=multiprocessing.Manager()
            # adj=manager.dict()
            # for x in range(len(all_go_terms)):
            #     for y in range(len(all_go_terms)):
            #         adj["{}_{}".format(x,y)]=-200

def calc_intra_similarity(all_go_terms, pf):
    manager = multiprocessing.Manager()
    adj = manager.dict()
    for x in range(len(all_go_terms)):
        for y in range(len(all_go_terms)):
            adj["{}_{}".format(x, y)] = -200
    params = []
    for i_x, x in enumerate(all_go_terms):
        for i_y, y in enumerate(all_go_terms):
            # calc_similarity(adj, i_x, i_y, x, y)  
            params.append([calc_similarity, [adj, i_x, i_y, x, y]])
    print "len(params): {}".format(len(params))

    
    p = multiprocessing.Pool(pf)
    p.map(func_star, params)
    p.close()
    p.join()

    adj_sum = sum([adj["{}_{}".format(x, y)] for x in range(len(all_go_terms)) for y in range(len(all_go_terms)) if
                   adj["{}_{}".format(x, y)] != -100])
    adj_count = float(len(
        [adj["{}_{}".format(x, y)] for x in range(len(all_go_terms)) for y in range(len(all_go_terms)) if
         adj["{}_{}".format(x, y)] != -100]))
    print "adj_sum: ", adj_sum
    print "adj_count: ", adj_count
    print "adj_avg: ", adj_sum/ max(adj_count, 1)

    return adj_sum, adj_count



            # adj_sum=sum([adj["{}_{}".format(x,y)] for x in range(len(all_go_terms)) for y in range(len(all_go_terms)) if adj["{}_{}".format(x,y)]!=-100])
            # adj_count=float(len([adj["{}_{}".format(x,y)] for x in range(len(all_go_terms)) for y in range(len(all_go_terms)) if adj["{}_{}".format(x,y)]!=-100]))
            # print "adj_sum: ",adj_sum
            # print "adj_count: ",adj_count
            #
            # p.close()
            # p.join()
            # print "n sig_terms: {}, variability: {}".format(len(all_go_terms), str(1 - adj_sum / adj_count))
            # file(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores", "{}_{}_{}".format(cur_ds,cur_algo, "n_sig.txt")), 'w+').write(str(len(all_go_terms)))
            # file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores", "{}_{}_{}".format(cur_ds, cur_algo, "var.txt")), 'w+').write(str(1 - adj_sum / adj_count) if adj_count>0 else str(0) )



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--pf', dest='pf', default=5)
    args = parser.parse_args()

    prefix = args.prefix
    datasets=["{}_{}".format(prefix,x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    pf=int(args.pf)
    print "test" 
    ds_summary=pd.DataFrame()
    main(datasets=datasets, algos=algos, prefix=prefix, pf=pf)
