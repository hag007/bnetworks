import sys
sys.path.insert(0, "../")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.community.quality import modularity

from networkx.algorithms.community.centrality import girvan_newman

from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

import os
import constants

def read_modules(file_path):
    passed_fdr=[]
    modules=[]

    with open(file_path, 'r') as f:
        line = f.readline()
        while line!="":
            line = f.readline()
            if line.startswith("cc"):
                modules.append(f.readline().strip()[1:-1].split(', '))
            if line.startswith("fdr"):
                passed_fdr=f.readline()[1:-1].split(' ')

        # i=0
        # for cur_module, cur_fdr_status in zip(modules, passed_fdr):
        #     if cur_fdr_status=='True':
        #         print "module number {} (n={}):\n{}".format(i, len(cur_module), cur_module)
        #     i+=1

        return modules


if __name__=="__main__":
    read_modules(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules.txt"))



