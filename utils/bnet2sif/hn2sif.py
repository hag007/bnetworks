import sys
import subprocess
import constants
import os
sys.path.insert(0, '../')

from utils import ensembl2entrez

entrez_dict = ensembl2entrez.get_entrez2ensembl_dictionary()



with file(os.path.join(constants.DICTIONARIES_DIR, "hn_ppi.txt")) as f1:
    global entrez_dict

    lns = f1.readlines()
    headers = lns[0]
    lns = lns[1:]

    lns_formatted = [[x.split()[0], x.split()[1]] for i, x in enumerate(lns)]
    lns_formatted = [x for i, x in enumerate(lns_formatted) if x[0] in entrez_dict and x[1] in entrez_dict]

    gene_list = list(set([entrez_dict[y] for x in lns_formatted for y in x]))
    file(os.path.join(constants.LIST_DIR, "hn_bg.txt"), "w+").write("\n".join(gene_list))

    lns_formatted = ["\t".join([entrez_dict[x[0]], "ppi", entrez_dict[x[1]]]) for i, x in enumerate(lns_formatted)]

    print len(lns_formatted)
    with file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hn_out.sif"), 'w+') as f2:
        f2.write("\n".join(lns_formatted))



