import os
from multiprocessing import Process

import constants
import utils.aggregate_reports as aggregate_reports
from runners import jactivemodules_greedy_runner
from runners import bionet_runner
from runners import hotnet2_runner
from runners import keypathwayminer_ines_runner
from runners import reactomefi_runner
from runners import netbox_runner
from runners import jactivemodules_sa_runner


def create_ds_folders(dataset_name):
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "data")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "cache")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "output")))




def run_dataset(dataset_name, expected_genes=None, disease_name=None, fdr=0.05):
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    prcs = [
            Process(target=bionet_runner.main, args=[dataset_name, disease_name, expected_genes, fdr]),
            Process(target=hotnet2_runner.main, args=[dataset_name, disease_name, expected_genes]),
            Process(target=netbox_runner.main, args=[dataset_name, disease_name, expected_genes]),
            Process(target=keypathwayminer_ines_runner.main, args=[dataset_name, disease_name, expected_genes]),
            Process(target=reactomefi_runner.main, args=[dataset_name, disease_name, expected_genes])
    ]
    for cur in prcs:
        cur.start()
    jac_g = Process(target=jactivemodules_greedy_runner.main, args=[dataset_name, disease_name, expected_genes])
    jac_g.start()
    jac_g.join()
    jac_s = Process(target=jactivemodules_sa_runner.main, args=[dataset_name, disease_name, expected_genes])
    jac_s.start()
    jac_s.join()
    for cur in prcs:
        cur.join()


if __name__ == "__main__":
    run_dataset("MCF7_2")
    aggregate_reports.aggregate_datasets("MCF7_2")

    # run_dataset("MCF7_2")
    # aggregate_reports.aggregate_datasets("MCF7_2")

