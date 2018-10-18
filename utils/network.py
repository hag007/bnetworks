import pandas as pd
import numpy as np
import os
import constants


def get_bg_genes(network_name="dip", h_src="ID_interactor_A", h_dst="ID_interactor_B"):
    network_df = pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_name+".sif"), sep="\t")
    src = np.array(network_df[h_src])
    dst = np.array(network_df[h_dst])
    vertices = list(set(np.append(src, dst)))
    return vertices