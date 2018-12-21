import sys
sys.path.insert(0, '../')
import constants
import go_hierarcies
import pandas as pd
import numpy as np
import os


dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
            roots=['GO:0008150'])
vertices=dict_result.values()[0]['vertices']


def get_all_genes_for_term(vertices, cur_root, term, in_subtree):

    in_subtree= in_subtree or term==cur_root
    all_genes = set()
    if in_subtree:
        all_genes.update(go2geneids[cur_root])
    for cur_child in vertices[cur_root]["obj"].children:
        all_genes.update(get_all_genes_for_term(vertices, cur_child.id, term, in_subtree))
    return all_genes

csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR,  "emp_diff.tsv")
df=pd.read_csv(csv_file_name,sep='\t',index_col=0)

n_genes=[]
depth=[]
for i, cur_go_id in enumerate(df.index.values):
    print "current i: {}/{}".format(i, len(df.index.values))
    n_genes.append(get_all_genes_for_term(vertices, 'GO:0008150', cur_go_id, 'GO:0008150'==cur_go_id))
    depth.append(dict_result.values()[0]['vertices'][cur_go_id]['D'])

df["n_genes"]=pd.Series(n_genes, index=df.index)
df["depth"]=pd.Series(depth, index=df.index)

df.to_csv(csv_file_name[:-4]+"_md.tsv",sep='\t')
