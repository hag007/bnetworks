import constants
import os
from utils import uniprot_dictionary

uniport_dict = uniprot_dictionary.get_ensg_dict()


with file(os.path.join(constants.DICTIONARIES_DIR, "dip_ppi.txt")) as f1:
    lns = f1.readlines()

    lns_formatted = [[x.split("\t")[0], x.split("\t")[1]] for x in lns]
    lns_formatted = [x for i, x in enumerate(lns_formatted) if i==0 or ("uniprotkb" in x[0] and "uniprotkb" in x[1])]

    lns_formatted = [[x[0].split(":")[-1], x[1].split(":")[-1]] if i>0 else x for i, x in enumerate(lns_formatted)]
    lns_formatted = [[x[0].split(":")[-1], x[1].split(":")[-1]] for i, x in enumerate(lns_formatted) if i == 0 or (uniport_dict.has_key(x[0]) and uniport_dict.has_key(x[1]))]
    lns_formatted = [[list(uniport_dict[x[0]]['ENSG'])[0], list(uniport_dict[x[1]]['ENSG'])[0]] if i>0 else x for i, x in enumerate(lns_formatted)]

    gene_list = list(set([y for x in lns_formatted for y in x]))
    file(os.path.join(constants.LIST_DIR, "dip_bg.txt"), "w+").write("\n".join(gene_list))

    lns_formatted = ["\t".join([x[0], "ppi", x[1]]) for x in lns_formatted]

    print len(lns_formatted)
    with file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "dip_out.sif"), 'w+') as f2:
        f2.write("\n".join(lns_formatted))





