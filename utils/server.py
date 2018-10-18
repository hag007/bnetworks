import sys
sys.path.insert(0, '../')

import argparse
import constants
import os
import time
import shutil


def get_parameters():
    parser = argparse.ArgumentParser(description='args from wev client')
    parser.add_argument('--gene_expression', dest='ge')
    parser.add_argument('--score', dest='score')
    parser.add_argument('--network', dest='nw')
    parser.add_argument('--mode', dest='mode', default=0)
    args = parser.parse_args()
    if args.mode==0:
        return None

    NETWORK_NAME = os.path.splitext(os.path.basename(args.nw))[0]
    dataset_name = "user" + str(time.time())
    constants.update_dirs(DATASET_NAME=dataset_name)
    os.makedirs(os.path.join(constants.DATA_DIR, "data"))
    os.makedirs(os.path.join(constants.OUTPUT_DIR, "output"))
    os.makedirs(os.path.join(constants.CACHE_DIR, "cache"))

    if args.ge != "" and args.ge != os.path.join(constants.DATA_DIR, os.path.basename(args.ge)):
        shutil.copy(args.ge, os.path.join(constants.DATA_DIR, os.path.basename(args.ge)))

    if args.score != "" and  args.score != os.path.join(constants.DATA_DIR, os.path.basename(args.score)):
        shutil.copy(args.score, os.path.join(constants.DATA_DIR, os.path.basename(args.score)))

    if args.nw != "" and args.nw != os.path.join(constants.NETWORKS_DIR, NETWORK_NAME + ".sif"):
        shutil.copy(args.nw, os.path.join(constants.NETWORKS_DIR, NETWORK_NAME + ".sif"))

    return args, NETWORK_NAME, dataset_name