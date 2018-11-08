import os
import constants
import shutil
import json
import numpy as np
import pandas as pd
import ebm
import network
from utils.scripts import format_script
from pandas.io.common import EmptyDataError

hg_report_header = ["GO id", "value", "pval"]
emb_headers= ["algo", "empirical_brown_score", "enriched_terms"]

def emb_algos(root_path, hg_report_name, general_report_name):

    reports = [os.path.join(root_path, name,hg_report_name)+".tsv" for name in os.listdir(root_path) if os.path.isdir(os.path.join(root_path, name))]
    df_rv = pd.DataFrame() # pd.read_csv(reports[0], sep="\t")[["value"]].rename(columns={'value': reports[0].split("/")[-2]})
    df_pval = pd.DataFrame() # pd.read_csv(reports[0], sep="\t")[["pval"]].rename(columns={'pval': reports[0].split("/")[-2]})
    for hg_report in reports:
        df_rv = pd.concat([df_rv, pd.read_csv(hg_report, sep="\t")["value"]], axis=1, join='outer').rename(columns={'value': hg_report.split("/")[-2]})
        df_rv = df_rv.fillna(0)
        df_pval = pd.concat([df_pval, pd.read_csv(hg_report, sep="\t")["pval"]], axis=1, join='outer').rename(columns={'pval': hg_report.split("/")[-2]})
        df_pval = df_pval.fillna(1)

    emb_report=[]
    for cur_col in df_pval.columns:
        print "current algo: {}".format(cur_col)
        general_report = pd.read_csv(os.path.join(root_path, cur_col, general_report_name) + ".tsv",sep="\t")
        emb_report.append({
            "algo" : str(cur_col),
            "empirical_brown_score" : str(ebm.KostsMethod(np.array(df_rv, dtype=np.float), np.array(df_pval[cur_col], dtype=np.float))),
            "total_num_genes" : general_report["total_num_genes"][0],
            "enriched_terms" : str((np.sum(df_pval[cur_col]!=1)))})


    pd.DataFrame(emb_report).set_index("algo").to_csv(os.path.join(root_path, hg_report_name+"_aggregated.tsv"), sep="\t")

def emb_modules(root_path, hg_report_name, modules_hg_report_name, algo_to_modules, modules_summary_file_name):

    reports = [os.path.join(root_path, algo,"module_{}_{}.tsv".format(module, modules_hg_report_name)) for algo, modules in algo_to_modules.iteritems() for module in modules]
    df_rv = pd.DataFrame()
    df_pval = pd.DataFrame()
    for hg_report in reports:
        new_col_name = hg_report.split("/")[-2] +"_"+ os.path.basename(hg_report).split("_")[1]
        try:
            df_rv = pd.concat([df_rv, pd.read_csv(hg_report, sep="\t")["value"]], axis=1, join='outer').rename(
                columns={'value': new_col_name})
            df_rv = df_rv.fillna(0)
            df_pval = pd.concat([df_pval, pd.read_csv(hg_report, sep="\t")["pval"]], axis=1, join='outer').rename(
                columns={'pval': new_col_name})
            df_pval = df_pval.fillna(1)
        except EmptyDataError:
            df_rv[new_col_name] = 0
            df_pval[new_col_name] = 1



    emb_report=[]

    for algo, modules in algo_to_modules.iteritems():
        modules_summary = pd.read_csv(os.path.join(root_path, algo, modules_summary_file_name) + ".tsv", sep="\t").set_index(
            "module")
        for module in modules:
            print "emb algo: {}, module: {}".format(algo, module)
            cur_algo_col = algo+"_"+str(module)
            emb_report.append({
                "algo_module": str(algo)+"_"+str(module),
                "algo" : str(algo),
                "module": str(module),
                "empirical_brown_score" : str(ebm.KostsMethod(np.array(df_rv, dtype=np.float), np.array(df_pval[cur_algo_col], dtype=np.float))),
                "total_num_genes" : modules_summary.loc[[module],:]["#_genes"].iloc[0],
                "enriched_terms" : str((np.sum(df_pval[cur_algo_col]!=1)))})


    pd.DataFrame(emb_report).set_index("algo_module").to_csv(os.path.join(root_path, hg_report_name+"_aggregated.tsv"), sep="\t")



if __name__ == "__main__":
    print "TNFa_2"
    constants.update_dirs(DATASET_NAME_u="TNFa_2")
    # calculate_emb("/home/hag007/bnet/output/TNFa_2", "k_{}_modules_hg_samples".format(network.MODULE_TH), "k_{}_modules_general".format(network.MODULE_TH))
    # emb_algos("/home/hag007/bnet/output/TNFa_2", "all_modules_hg_samples", "all_modules_general")
    # print "MCF7"
    # aggregate_variable("/home/hag007/bnet/output/MCF7_2", "k_10_modules_summary_hg_samples")
    # aggregate_variable("/home/hag007/bnet/output/MCF7_2", "all_modules_summary_hg_samples")

    emb_modules("/home/hag007/bnet/output/TNFa_2", "all_separated_modules_hg_samples", "separated_modules_hg_samples",{"jactivemodules_greedy" : [0,1], "netbox" : [0,1]},"modules_summary")

    emb_modules("/home/hag007/bnet/output/TNFa_2", "k_{}_separated_modules_hg_samples".format(network.MODULE_TH),"separated_modules_hg_samples", {"jactivemodules_greedy": [0,1], "netbox" : [0,1]},"modules_summary")


