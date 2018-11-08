import pandas as pd
import numpy as np
import os
import constants
import time
import shutil
import sys
import json
import pandas as pd

from df_helpers import to_full_list
from df_helpers import to_full_np
from utils.scripts import format_script
import zipfile

from utils.go import check_group_enrichment
import go

SH_MODULE_NAME = "module"
SH_NUM_GENES = "#_genes"
SH_ENRICHED = "enriched_groups"
SH_DETAILS = "more_details"
SH_TABLE_HEADERS = [SH_MODULE_NAME, SH_NUM_GENES, SH_ENRICHED, SH_DETAILS]

MODULE_TH = 10

def zipdir(path_to_zip, zip_file_path):
    ziph = zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(path_to_zip):
        for file in files:
            ziph.write(os.path.join(root, file))

def get_network_genes(network_name="dip", h_src="ID_interactor_A", h_dst="ID_interactor_B"):
    network_df = pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_name+".sif"), sep="\t")
    src = np.array(network_df[h_src])
    dst = np.array(network_df[h_dst])
    vertices = list(set(np.append(src, dst)))
    return vertices

def remove_subgraph_self_loops(nodes_to_remove, network_file_name=os.path.join(constants.NETWORKS_DIR,"dip.sif"), h_src="ID_interactor_A", h_dst="ID_interactor_B"):
    if len(nodes_to_remove) == 0:
        return network_file_name
    network_df = pd.read_csv(network_file_name, sep="\t")
    filtered_network = network_df[network_df[h_src]!=network_df[h_dst.isin(nodes_to_remove)]]
    new_file_name = os.path.splitext(network_file_name) + "_no_loops" +".sif"
    filtered_network.to_csv(new_file_name, sep="\t", index=False)
    return filtered_network

def remove_subgraph_by_nodes(nodes_to_remove, network_file_name=os.path.join(constants.NETWORKS_DIR,"dip.sif"), h_src="ID_interactor_A", h_dst="ID_interactor_B", ts=str(time.time())):
    if len(nodes_to_remove) == 0:
        return network_file_name
    network_df = pd.read_csv(network_file_name, sep="\t")
    filtered_network = network_df[~(network_df[h_src].isin(nodes_to_remove) | network_df[h_dst].isin(nodes_to_remove))]
    new_file_name = os.path.splitext(network_file_name)[0] + ts +".sif"
    filtered_network.to_csv(new_file_name, sep="\t", index=False)
    return new_file_name



def summary_intergrative_reports(all_hg_reports, modules_summary, total_hg_report, algo_name, module_genes, disease_name, expected_genes, report_file_name):



    general_algo_report(algo_name, all_hg_reports, module_genes, modules_summary, report_file_name, total_hg_report)

    if constants.DISEASE_MODE:
        disease_algo_report(algo_name, disease_name, expected_genes, module_genes, modules_summary, report_file_name)

    if constants.EMB_MODE:
        emb_score_report(algo_name, report_file_name, "hg_samples", total_hg_report)


def emb_score_report(algo_name, report_file_name, hg_sample_file_name, hg_report):
    samples = [{go.HG_GO_ID : cur_term[go.HG_GO_ID], go.HG_VALUE : cur_term[go.HG_VALUE], go.HG_PVAL : cur_term[go.HG_PVAL]} for cur_term in hg_report]
    df_emb = pd.DataFrame(samples)
    df_emb.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, algo_name,
                     "{}_{}.tsv".format(report_file_name, hg_sample_file_name)),sep="\t", index=False)
    return df_emb


def disease_algo_report(algo_name, disease_name, expected_genes, module_genes, modules_summary, report_file_name):

    disease_data = {
        "disease_name": disease_name,
        "num_of_modules": len(modules_summary),
        "TP+FN_(_true_)": len(expected_genes),
        "TP+TN_(_retrieved_)": len(module_genes),
        "TP/(TP+TN)_(_precision_)": 0,
        "TP/(TP+FN)_(_recall_)": 0,
        "F1": 0,
        "TP": 0,
        "module_size_avg" : 0,
        "module_size_std" :0
    }
    if len(modules_summary) > 0:
        modules_summary = pd.DataFrame(modules_summary)
        disease_genes_extracted = float(len(set(module_genes).intersection(expected_genes)))
        disease_data["TP"] = disease_genes_extracted
        disease_data["TP/(TP+TN)_(_precision_)"] = disease_genes_extracted / len(module_genes)
        disease_data["TP/(TP+FN)_(_recall_)"] = disease_genes_extracted / len(expected_genes)
        disease_data["F1"] = 2 * ((disease_data["TP/(TP+TN)_(_precision_)"] * disease_data["TP/(TP+FN)_(_recall_)"]) /
                                  (disease_data["TP/(TP+TN)_(_precision_)"] + disease_data["TP/(TP+FN)_(_recall_)"]))

        disease_data["module_size_avg"] = modules_summary[SH_NUM_GENES].mean()
        disease_data["module_size_std"] = modules_summary[SH_NUM_GENES].std()


    pd.DataFrame([disease_data]).to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, algo_name, "{}_disease.tsv".format(report_file_name)),sep="\t", index=False)


def general_algo_report(algo_name, all_hg_reports, module_genes, modules_summary, report_file_name, total_hg_report):
    data = {}
    if len(modules_summary) > 0 :
        df_summary = pd.DataFrame(modules_summary)
        data = {"num_of_modules": df_summary.index.size,
                "module_size_avg": df_summary[SH_NUM_GENES].mean(),
                "module_size_std": df_summary[SH_NUM_GENES].std(),
                "total_num_genes": len(module_genes)
                }


    if len(all_hg_reports) > 0:
        df_all_hg = [pd.DataFrame(x) for x in all_hg_reports]
        enrichment_dist = [x.index.size for x in df_all_hg]
        pval_dist = [np.array(x[go.HG_PVAL].astype(np.float).apply(lambda x: -np.log10(x))) if x.index.size > 0 else np.array([]) for x in df_all_hg]
        modules_enrichment_data = {"module_enriched_terms_avg": np.average(enrichment_dist),
                                   "module_enriched_terms_std": np.std(enrichment_dist),
                                   "module_enriched_terms_signal_avg_avg": np.average([np.average(x) if len(x) > 1 else 0
                                                                                        for x in pval_dist]),
                                   "module_enriched_terms_signal_avg_std": np.average([np.std(x) if len(x) > 1 else 0
                                        for x in pval_dist]),
                                   "module_enriched_terms_signal_std_avg": np.std([np.average(x) if len(x) > 1 else 0
                                        for x in pval_dist]),
                                   "module_enriched_terms_signal_std_std": np.std([np.std(x) if len(x) > 1 else 0
                                        for x in pval_dist])}


        data.update(modules_enrichment_data)
        data["module_enriched_terms_signal_score"] = \
            data['module_enriched_terms_signal_avg_avg'] / ((data[
                                                                 'module_enriched_terms_signal_avg_std'] +
                                                             data[
                                                                 'module_enriched_terms_signal_std_avg'] +
                                                             data[
                                                                 'module_enriched_terms_signal_std_std']) * \
                                                            (data[
                                                                 "total_num_genes"] /
                                                             data[
                                                                 "module_size_avg"] *
                                                             data[
                                                                 "num_of_modules"]))

    if len(total_hg_report) > 0:
        df_total_hg = pd.DataFrame(total_hg_report)
        all_enrichment_data = {
            "total_enriched_terms_avg" : df_total_hg[go.HG_PVAL].astype(np.float).apply(lambda x: -np.log10(x)).mean(),
            "total_enriched_terms_std" : df_total_hg[go.HG_PVAL].astype(np.float).apply(lambda x: -np.log10(x)).std(),
            "total_num_enriched_terms": len(total_hg_report)
        }
        data.update(all_enrichment_data)

    df = pd.DataFrame()
    if len(data) >0:
        df = pd.DataFrame([data])

    df.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, algo_name,
                     "{}_general.tsv".format(report_file_name)), sep="\t", index=False)


def draw_network(active_genes, network_name="dip", h_src="ID_interactor_A", h_dst="ID_interactor_B"):
    # type: (object, object, object, object) -> object

    reduced_nodes = [[x[h_src], x[h_dst]] for i, x in pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_name+".sif"), sep="\t").iterrows() if x[h_src] in active_genes or x[h_dst] in active_genes]
    reduced_nodes = [y for x in reduced_nodes for y in x]
    reduced_nodes = list(set(reduced_nodes))
    cy_edges_all = [{"data": { "id": x[h_src]+x[h_dst], "source": x[h_src], "target": x[h_dst], "label" : "-"}, "selected" : True} if x[h_src] in active_genes and x[h_dst] in active_genes
                    else {"data": { "id": x[h_src]+x[h_dst], "source": x[h_src], "target": x[h_dst]}}
                    for i, x in pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_name+".sif"), sep="\t").iterrows() if x[h_src] in active_genes or x[h_dst] in active_genes]

    cy_nodes_all = [{"data": {"id" : x, "label": x}, "selected" : True} if x in active_genes else {"data": {"id" : x}} for x in reduced_nodes]



    cy_nodes_active = [{"data": {"id" : x, "label": x}, "selected" : True} for x in active_genes]
    cy_edges_active = [{"data": { "id": x[h_src]+x[h_dst], "source": x[h_src], "target": x[h_dst], "label" : "-"}, "selected" : True}
                       for i, x in pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_name+".sif"), sep="\t").iterrows() if x[h_src] in active_genes and x[h_dst] in active_genes]

    return (cy_edges_active + cy_nodes_active), (cy_edges_all + cy_nodes_all)



def generate_report_from_template(output_file_name, cy_active, cy_all, algo_name="", hg_report=[],
                                  disease_genes_statistics=[], modules_summary=[]):

    hg_report = to_full_list(hg_report, "#")
    disease_genes_statistics = to_full_list(disease_genes_statistics, "#")
    modules_summary = to_full_list(modules_summary, "#")

    format_script(os.path.join(constants.TEMPLATES_DIR, "graph.html"), DATA_ALL=json.dumps(cy_all),
                  DATA_ACTIVE=json.dumps(cy_active), HG_REPORT=json.dumps(hg_report),
                  MODULES_SUMMARY=json.dumps(modules_summary), NUM_OF_GENES=len([x for x in cy_active if not x["data"].has_key("source")]),
                  DISEASE_GENES=json.dumps(disease_genes_statistics))
    output_dir = os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, algo_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    shutil.move(os.path.join(constants.TEMPLATES_DIR, "graph.html"),
                os.path.join(output_dir, "graph_{}.html".format(output_file_name)))
    return "graph_{}.html".format(output_file_name)


def build_all_reports(algo_name, modules, all_bg_genes, network_file_name, disease_name=None, expected_genes=None):
    output_base_dir = os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, algo_name)
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)

    all_hg_reports = []
    modules_summary = []

    for i, module in enumerate(modules):
        hg_row, summary_row = module_report(algo_name, i, module, all_bg_genes[i])
        all_hg_reports.append(hg_row)
        modules_summary.append(summary_row)

    module_genes = list(set(reduce((lambda x, y: x + y), modules)))
    modules_larger_than_k, module_larger_than_k_genes, k_hg_reports, k_modules_summary = \
        get_k_threshold_modules(modules, all_hg_reports, modules_summary)

    pd.DataFrame(modules_summary).set_index("module").to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, algo_name, "modules_summary.tsv"), sep="\t")
    generate_algo_report(algo_name, all_bg_genes[0], all_hg_reports, disease_name, expected_genes, module_genes,
                         modules_summary, "all_modules")

    generate_algo_report(algo_name, all_bg_genes[0], k_hg_reports, disease_name, expected_genes, module_larger_than_k_genes,
                         k_modules_summary, "k_{}_modules".format(MODULE_TH))

    # zipdir(output_main_dir, os.path.join(output_base_dir, "report.zip"))
    return output_base_dir






def get_k_threshold_modules(modules, all_hg_reports, modules_summary):
    modules_larger_than_k = [cur for cur in modules if len(cur) >= MODULE_TH]
    module_larger_than_k_genes = list(set(reduce((lambda x, y: x + y), modules_larger_than_k, [])))
    k_modules_summary = [modules_summary[i] for i, cur in enumerate(modules) if len(cur) >= MODULE_TH]
    k_hg_reports = []
    if constants.HG_MODE:
        k_hg_reports = [all_hg_reports[i] for i, cur in enumerate(modules) if len(cur) >= MODULE_TH]
    return modules_larger_than_k, module_larger_than_k_genes, k_hg_reports, k_modules_summary


def generate_algo_report(algo_name, bg_genes, all_hg_reports, disease_name, expected_genes, module_genes,
                         modules_summary, report_name):
    hg_report = []
    if constants.HG_MODE or constants.EMB_MODE:
        hg_report = check_group_enrichment(module_genes, bg_genes)
    active_cy, all_cy = draw_network(active_genes=module_genes)
    generate_report_from_template(report_name, active_cy, all_cy, algo_name, hg_report, [], modules_summary)
    summary_intergrative_reports(all_hg_reports, modules_summary, hg_report, algo_name, module_genes, disease_name,
                                 expected_genes, report_name)


def module_report(algo_name, i, module, bg_genes):
    file(os.path.join(constants.OUTPUT_DIR, "{}_module_genes_{}.txt".format(algo_name, i)), "w+").write(
        "\n".join(module))
    file(os.path.join(constants.OUTPUT_DIR, "{}_bg_genes_{}.txt".format(algo_name, i)), "w+").write(
        "\n".join(bg_genes))
    modules_summary_row = {SH_MODULE_NAME: i, SH_NUM_GENES: len(module)}
    hg_report = []
    if constants.HG_MODE:
        hg_report = check_group_enrichment(list(module), list(bg_genes))
        modules_summary_row[SH_ENRICHED] = len(hg_report) - 1
        if constants.EMB_MODE:
            hg_report = emb_score_report(algo_name, "module_" + str(i), "separated_modules_hg_samples", hg_report)
    active_cy, all_cy = draw_network(active_genes=module)
    report_output_file_name = generate_report_from_template(algo_name + str(i), active_cy, all_cy, algo_name, hg_report)
    modules_summary_row[SH_DETAILS] = report_output_file_name
    return hg_report, modules_summary_row


