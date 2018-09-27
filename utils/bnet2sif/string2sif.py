import constants
import os
from utils import ensp_dictionary

ensp_dict = ensp_dictionary.get_ensp_dict()



with file(os.path.join(constants.DICTIONARIES_DIR, "string_ppi_small.txt")) as f1:
    global ensp_dict

    lns = f1.readlines()
    headers = lns[0]
    lns = lns[1:]

    lns_formatted = [[x.split()[0][5:], x.split()[1][5:]] for i, x in enumerate(lns)]
    lns_formatted = [x for i, x in enumerate(lns_formatted) if x[0] in ensp_dict and x[1] in ensp_dict]
    lns_formatted = ["\t".join([list(ensp_dict[x[0]]["ENSG"])[0], "ppi", list(ensp_dict[x[1]])[0]]) for i, x in enumerate(lns_formatted)]

    print len(lns_formatted)
    with file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "string_out_small.sif"), 'w+') as f2:
        f2.write("\n".join(lns_formatted))



