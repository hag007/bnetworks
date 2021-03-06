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

semsim = SetSemSim(ontology, ac, TSS="Resnik", MSS="BMA") # SetSemSim(ontology, ac, TSS="GSESAME", MSS="BMA")# GSESAMESemSim(ontology, ac)  # maxSemSim(ontology, ac) #

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


PVAL="emp_pval" # "pval"# "P-value"
QVAL="qval"# "FDR q-value"
METHOD='average'
def main(base_folder="/home/hag007/Desktop/aggregate_report/oob", pf=5):

    go_terms_tables={}
    go_terms = {}
    norm=0
    for cur in [x for x in os.listdir(base_folder) if not os.path.isdir(x) and x and x.endswith("oob.tsv")]:

        go_terms_tables[cur[9:-15]]=pd.read_csv(os.path.join(base_folder, cur),sep='\t', index_col=0)
        norm=max(-np.log10(go_terms_tables[cur[9:-15]][PVAL].min()),norm)

    for k,v in go_terms_tables.iteritems():
        go_terms[k]= v[v["passed_oob_permutation_test"]==True ] # v[v[QVAL]<=0.05][v["B"]<=500][v["B"]>=10]

    df_summary=pd.DataFrame()
    for cur_x in go_terms.keys():
        for cur_y in go_terms.keys():
            if go_terms[cur_y].index.size == 0 or go_terms[cur_x].index.size == 0: continue
            if cur_y==cur_x:
                df_summary.loc[cur_x, cur_y] = 0

            if go_terms.keys().index(cur_y) > go_terms.keys().index(cur_x):

                manager=multiprocessing.Manager()
                adj=manager.dict()
                for x in range(len(go_terms[cur_x].index)):
                    for y in range(len(go_terms[cur_y].index)):
                        adj["{}_{}".format(x,y)]=-200

                params=[]
                for i_x, x in enumerate(go_terms[cur_x].index):
                    for i_y, y in enumerate(go_terms[cur_y].index):
                        params.append([calc_similarity, [adj, i_x, i_y, x, y]]) # , -np.log10(go_terms[cur_x].loc[x,PVAL]), -np.log10(go_terms[cur_y].loc[y,PVAL]), norm]])

                p = multiprocessing.Pool(pf)
                p.map(func_star,params)

                adj_sum = sum(
                    [adj["{}_{}".format(x, y)] if adj["{}_{}".format(x, y)] != -100 else 10 for x in range(len(go_terms[cur_x])) for y in range(len(go_terms[cur_y]))])
                adj_count = float(len(
                    [adj["{}_{}".format(x, y)] if adj["{}_{}".format(x, y)] != -100 else 10 for x in range(len(go_terms[cur_x])) for y in range(len(go_terms[cur_y]))]))

                # adj_sum_x_max=[[adj["{}_{}".format(x,y)]  if adj["{}_{}".format(x,y)]!=-100 else 10 for y in range(len(go_terms[cur_y]))] for x in range(len(go_terms[cur_x]))]
                # adj_sum_x_max=[max(x) for x in adj_sum_x_max]
                # adj_sum_y_max = [[adj["{}_{}".format(x, y)] if adj["{}_{}".format(x, y)] != -100 else 10 for x in range(len(go_terms[cur_x]))] for y in range(len(go_terms[cur_y]))]
                # adj_sum_y_max = [max(x) for x in adj_sum_y_max]

                # adj_sum_max=adj_sum_x_max+adj_sum_y_max
                # adj_sum=sum(adj_sum_max) # - len(set(go_terms[cur_x]).intersection(go_terms[cur_y]))
                # adj_count= len(list(np.append(go_terms[cur_x], go_terms[cur_y])))

                print "adj_sum: ",adj_sum
                print "adj_count: ",adj_count

                if adj_sum >0:
                    df_summary.loc[cur_x,cur_y]=- adj_sum/adj_count
                    df_summary.loc[cur_y, cur_x]=df_summary.loc[cur_x, cur_y]
                print cur_x, cur_y, go_terms.keys().index(cur_x), go_terms.keys().index(cur_y)
                p.close()
                p.join()

    print df_summary
    df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "cancer_type_go_distance_original_1.tsv"), sep='\t')

    distArray = ssd.squareform(df_summary[df_summary.index.values].values+1)

    linked = linkage(distArray, method=METHOD, metric='euclidean')
    colors=[[i_y for i_y, y in enumerate(algos) if y in x][0] for x in list(df_summary.index)]
    plt.figure(figsize=(10, 7))
    dendrogram(linked,
               orientation='top',
               labels=[x.split(".")[0] for x in df_summary.index.values],
               distance_sort='descending',
               colors=colors,
               show_leaf_counts=True)
    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hierarchical_go_clustering.png"))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,netbox,bionet,hotnet2,keypathwayminer_INES_GREEDY")
    parser.add_argument('--pf', dest='pf', default=40)
    args = parser.parse_args()

    prefix = args.prefix
    datasets=["{}_{}".format(prefix,x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    pf=int(args.pf)
    print "test" 
    ds_summary=pd.DataFrame()
    for cur_ds in datasets:
        print "current dataset: {}".format(cur_ds)
        main(pf=pf, base_folder=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr/MAX"))

