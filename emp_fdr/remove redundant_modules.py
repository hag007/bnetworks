import math
import random
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

import utils.go_hierarcies as go_hierarcies

import simplejson as json

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


ROOT='GO:0008150'
dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
    roots=[ROOT])
vertices = dict_result.values()[0]['vertices']

sim_metric="Resnik"
semsim = SetSemSim(ontology, ac, TSS=sim_metric, MSS="BMA")
go_hierarcies
CUTOFF=7

def calc_similarity(mat_adj, i_x, i_y, x, y):
    key="{}_{}".format(i_x,i_y)
    key_inv="{}_{}".format(i_y,i_x)
    if mat_adj[key] != -200: return
    if(x==y):
       mat_adj[key]=-100
       return 

    mat_adj[key] = semsim.SemSim(x, y)

    # if mat_adj[key]<0:
    print i_x,i_y,mat_adj[key]

    if np.isnan(mat_adj[key]):
        mat_adj[key] = -100
    mat_adj[key_inv] = mat_adj[key]

def main(datasets, algos, pf=10, base_folder='/home/hag007/Desktop/aggregate_gwas_report/oob', file_format="emp_diff_modules_{}_{}_passed_oob.tsv"):


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
                    os.path.join(base_folder,
                                 file_format.format(cur_ds, cur_algo)), sep='\t', index_col=0)

            except:
                total_num_genes.append(0)
                algos_signals.append(0)
                algo_go_sim_score.append(1)
                continue

            emp_results=emp_results.sort_values(by='emp_rank')
            emp_results_fdr=emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values,:]
            

            algos_signals.append(len(emp_results_fdr.index))
            all_go_terms = emp_results_fdr.index.values

            try:
                total_num_genes.append(pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix,cur_ds), cur_algo, "all_modules_general.tsv"),
                    sep="\t")["total_num_genes"][0])
            except:
                total_num_genes.append(0)

            cache_file=os.path.join(base_folder, "similarity_cache_{}_{}_{}.npy".format(sim_metric, cur_ds,cur_algo))
            go_terms_reduced, adj = calc_intra_similarity(all_go_terms, pf, emp_results_fdr, cache_file)

            emp_results_fdr[emp_results_fdr.index.isin(go_terms_reduced)].to_csv(os.path.join("/home/hag007/Desktop/aggregate_report/oob", file_format.format(cur_ds,cur_algo)[:-4]+"_reduced.tsv"), sep='\t')
            open(cache_file, 'w+').write(json.dumps(dict(adj)))

            print "total redundancy removal: {}/{}".format(emp_results_fdr[emp_results_fdr.index.isin(go_terms_reduced)].shape[0], emp_results_fdr.shape[0])

def calc_intra_similarity(all_go_terms, pf, enrichment_scores, cache_file):
    all_go_terms = list(all_go_terms)
    if os.path.exists(cache_file):
        adj=json.load(open(cache_file,'r'))
    else:
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

    for cur_go_1, cur_go_id_1 in enumerate(list(all_go_terms)):
        for cur_go_2, cur_go_id_2 in enumerate(list(all_go_terms)):
            if cur_go_id_1 not in all_go_terms or cur_go_id_2 not in all_go_terms or cur_go_id_2 == cur_go_id_1:
                continue
            if adj["{}_{}".format(cur_go_1, cur_go_2)] > CUTOFF:
                if enrichment_scores.loc[cur_go_id_1,'hg_pval_max'] > enrichment_scores.loc[cur_go_id_2,'hg_pval_max'] + 1:
                    all_go_terms.remove(cur_go_id_2)
                elif enrichment_scores.loc[cur_go_id_1,'hg_pval_max'] + 1 < enrichment_scores.loc[cur_go_id_2,'hg_pval_max'] :
                    all_go_terms.remove(cur_go_id_1)
                    break
                elif is_parent(cur_go_id_1,cur_go_id_2):
                    all_go_terms.remove(cur_go_id_2)
                elif is_parent(cur_go_id_2,cur_go_id_1):
                    all_go_terms.remove(cur_go_id_1)
                else:
                    all_go_terms.remove([cur_go_id_1, cur_go_id_2][int(math.floor(random.random()*2))])


    print all_go_terms
    return all_go_terms, adj


def is_parent(parent,child):
    vertices = dict_result.values()[0]['vertices']
    c_parents=vertices[child]["obj"].parents
    if child==ROOT:
        return False

    elif len(c_parents)==0:
        return True

    else:
        return any([is_parent(parent, a.id) for a in list(c_parents)])








    return adj, all_go_terms



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
    parser.add_argument('--datasets', dest='datasets', default="Breast_Cancer2.G50") # ",Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM")
    parser.add_argument('--algos', dest='algos', default="my_netbox_td") # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    parser.add_argument('--pf', dest='pf', default=3)
    args = parser.parse_args()

    prefix = args.prefix
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    pf=int(args.pf)
    print "test" 
    ds_summary=pd.DataFrame()
    for cur_ds in datasets:
        print "current dataset: {}".format(cur_ds)
        main(datasets=datasets, algos=algos, pf=pf)
