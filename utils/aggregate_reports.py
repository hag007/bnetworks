import os
import constants
import shutil
import json
import numpy as np
import pandas as pd
import network
from utils.df_helpers import to_full_list
from utils.scripts import format_script
from utils.integrate_pvals_ebm import emb_algos
from utils.integrate_pvals_ebm import emb_modules
from scipy.stats import mannwhitneyu
# score_headers = ["overlap_score_(_OS_)", "enriched_terms_score_(_ETS_)", "pval_signal_score_(_PSS_)",
#                  "PSS_OS_combined_(_PSOS_)", "PSOS_*_ETS_(_PSOSETS_)", ]


def aggregate_all(root_path, report_name):

    aggregated_report = aggregate_report(report_name, root_path)
    aggregated_score = pd.DataFrame()
    if constants.HG_MODE:
        aggregated_score = aggregate_scores(aggregated_report, report_name, root_path)

    return aggregated_report, aggregated_score


def aggregate_scores(aggregated_report, report_name, root_path):
    aggregated_score = []
    C1 = 2 * aggregated_report[["module_enriched_terms_signal_std_avg", "module_enriched_terms_signal_avg_std"]].mean().sum()
    C2 = aggregated_report["module_enriched_terms_std"].mean()
    for i, report in aggregated_report.iterrows():
        aggregated_score_row = calc_row_score(C1, C2, i, report)
        aggregated_score.append(aggregated_score_row)

    df_aggregated_score = pd.DataFrame(aggregated_score).set_index("algo")
    df_aggregated_score.to_csv(os.path.join(root_path, report_name + "_score_aggregated.tsv"), sep="\t")

    return df_aggregated_score


def aggregate_report(report_name, root_path):
    aggregated_report = None
    reports = [os.path.join(root_path, name, report_name) + ".tsv" for name in os.listdir(root_path) if
               os.path.isdir(os.path.join(root_path, name))]

    for report in reports:
        report_data = pd.read_csv(report, sep="\t")
        report_data["algo"] = report.strip().split("/")[-2]
        if aggregated_report is None:
            aggregated_report = report_data
        else:
            aggregated_report = pd.concat([aggregated_report, report_data])

    aggregated_report =aggregated_report.set_index('algo')

    aggregated_report.to_csv(os.path.join(root_path, report_name + "_aggregated.tsv"), sep="\t")
    return aggregated_report


def calc_row_score(C1, C2, algo, report_data):

    return {"algo" : algo,
                                 "overlap_score_(_OS_)" : report_data["total_num_genes"] / (report_data["module_size_avg"] * report_data["num_of_modules"]),
                                 "enriched_terms_score_(_ETS_)" : report_data['module_enriched_terms_avg'] / (report_data['module_enriched_terms_std'] + C2),
                                 "pval_signal_score_(_PSS_)" : report_data['module_enriched_terms_signal_avg_avg'] /
                                                               ((report_data["module_size_avg"] / report_data["total_num_genes"]) * C1 +
                                                                report_data['module_enriched_terms_signal_std_avg'] + report_data[
                                                                    'module_enriched_terms_signal_avg_std']),
                                 "PSS_OS_combined_(_PSOS_)" : report_data['module_enriched_terms_signal_avg_avg'] /
                                                              ((report_data["module_size_avg"] * report_data["num_of_modules"] / report_data[
                                                                  "total_num_genes"]) * C1 +
                                                               report_data['module_enriched_terms_signal_std_avg'] + report_data[
                                                                   'module_enriched_terms_signal_avg_std']),
                                 "PSOS_*_ETS_(_PSOSETS_)" : report_data['module_enriched_terms_signal_avg_avg'] /
                                                            ((report_data["module_size_avg"] * report_data["num_of_modules"] / report_data[
                                                                "total_num_genes"]) * C1 +
                                                             report_data['module_enriched_terms_signal_std_avg'] + report_data[
                                                                 'module_enriched_terms_signal_avg_std']) *
                                                            report_data['module_enriched_terms_avg'] / (
                                                                    report_data['module_enriched_terms_std'] + C2)}

def aggregate_datasets(DATASET_NAME = constants.DATASET_NAME):

    if constants.EMB_MODE:
        root_path = os.path.join(constants.OUTPUT_GLOBAL_DIR, DATASET_NAME)
        emb_algos(root_path, "all_modules_hg_samples", "all_modules_general")
        emb_algos(root_path, "k_{}_modules_hg_samples".format(network.MODULE_TH), "k_{}_modules_general".format(network.MODULE_TH))
        if constants.HG_MODE:


            all_algo_modules = {}
            k_algo_modules = {}
            for name in os.listdir(root_path):
                if os.path.isdir(os.path.join(root_path, name)):
                    modules_summary = pd.read_csv(os.path.join(root_path, name,"modules_summary.tsv"), sep="\t").set_index(
                        "module")
                    all_algo_modules[name] = np.array(
                        modules_summary.index)
                    k_algo_modules[name] = np.array(
                        modules_summary[modules_summary["#_genes"] >= 10].index)


            emb_modules(os.path.join(constants.OUTPUT_GLOBAL_DIR, DATASET_NAME), "all_separated_modules_hg_samples",
                        "separated_modules_hg_samples", all_algo_modules,
                        "modules_summary")
            emb_modules(os.path.join(constants.OUTPUT_GLOBAL_DIR, DATASET_NAME), "k_{}_separated_modules_hg_samples".format(network.MODULE_TH),
                        "separated_modules_hg_samples", k_algo_modules,
                        "modules_summary")

    constants.update_dirs(DATASET_NAME_u=DATASET_NAME)

    reports_metadata = {"all": [aggregate_all(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME),
                                              "all_modules_general"),
                                pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME,
                                                         "all_modules_hg_samples_aggregated.tsv"), sep="\t", index_col=0, dtype=np.str),
                                pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME,
                                                         "all_separated_modules_hg_samples_aggregated.tsv"), sep="\t",
                                            index_col=0, dtype=np.str)
                                ],
                        "k_{}".format(network.MODULE_TH): [aggregate_all(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME),
                                                                         "k_{}_modules_general".format(network.MODULE_TH)),
                                                           pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME,
                                                                                    "k_{}_modules_hg_samples_aggregated.tsv".format(network.MODULE_TH)), sep="\t", index_col=0, dtype=np.str),
                                                           pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                                                                                    constants.DATASET_NAME,
                                                                                    "k_{}_separated_modules_hg_samples_aggregated.tsv".format(
                                                                                        network.MODULE_TH)), sep="\t",
                                                                       index_col=0, dtype=np.str)
                                                           ]}


    for k, v in reports_metadata.iteritems():
        mannwhitneyu_rank = calc_mann_whitney_U(root_path, v)



        concated = to_full_list(pd.concat([v[0][1], v[1]], axis=1), "algo")
        format_script(os.path.join(constants.TEMPLATES_DIR, "report.html"), REPORT=json.dumps(to_full_list(v[0][0], "algo")),
                      SCORE=json.dumps(concated), DISEASE_GENES=[], DISEASE_GENES_SUMMARY=[], MODULES_SCORE=json.dumps(to_full_list(v[2], "algo_module")), EMB_WU=json.dumps(to_full_list(mannwhitneyu_rank, "algo")), MODULE_FILTER=(k+"_modules"))
        output_dir = os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        shutil.move(os.path.join(constants.TEMPLATES_DIR, "report.html"),
                    os.path.join(output_dir, "report_{}.html".format(k)))


def calc_mann_whitney_U(root_path, v):
    mannwhitneyu_rank = []
    for name in os.listdir(root_path):
        if os.path.isdir(os.path.join(root_path, name)):
            v_temp = v[2].copy()
            v_temp["empirical_brown_score"] = v_temp["empirical_brown_score"].astype(np.float).apply(
                lambda x: -np.log10(x))
            x = np.array(v_temp[v_temp["algo"] == name]["empirical_brown_score"])
            y = np.array(v_temp[v_temp["algo"] != name]["empirical_brown_score"])
            mannwhitneyu_rank.append({"algo": name, "MWU_pval": mannwhitneyu(x, y).pvalue})
            print "{}".format(mannwhitneyu_rank[-1])
    return mannwhitneyu_rank


def aggregate_disease(DATASET_NAME):
    counter = 1
    report_all = pd.DataFrame()
    report_k = pd.DataFrame()
    while True:
        constants.update_dirs(DATASET_NAME_u=DATASET_NAME + "_{}".format(counter))
        if not os.path.exists(constants.DATASET_DIR): break

        reports_metadata = {"all": aggregate_all(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME),
                                                 "all_modules_disease"),
                            "k": aggregate_all(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME),
                                               "k_{}_modules_disease".format(network.MODULE_TH))}
        report_all = pd.concat([report_all, reports_metadata["all"][0]])
        report_k = pd.concat([report_k, reports_metadata["k"][0]])
        counter += 1
    report_all = report_all.set_index("algo")
    report_k = report_k.set_index("algo")
    df_summary_all = calc_p_r_f_scores(report_all)
    df_summary_k = calc_p_r_f_scores(report_k)
    for k, v in {"all": [report_all, df_summary_all], "k_{}".format(network.MODULE_TH): [report_k, df_summary_k]}.iteritems():
        format_script(os.path.join(constants.TEMPLATES_DIR, "report.html"), REPORT=[],
                      SCORE=[], DISEASE_GENES=json.dumps(to_full_list(v[0], "algo")),
                      DISEASE_GENES_SUMMARY=json.dumps(to_full_list(v[1], "algo")), MODULE_FILTER=(k+"_modules"))
        output_dir = os.path.join(constants.OUTPUT_GLOBAL_DIR, DATASET_NAME)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        shutil.move(os.path.join(constants.TEMPLATES_DIR, "report.html"),
                    os.path.join(output_dir, "report_{}.html".format(k)))


def calc_p_r_f_scores(df):
    algos = df["algo"].unique()
    disease_summary = []
    for cur_algo in algos:
        df_2 = df[df["algo"] == cur_algo]
        disease_summary.append({
            "algo": str(cur_algo),
            "TP_mean": df_2["TP"].astype(np.int).mean(),
            "TP_std": df_2["TP"].astype(np.int).std(ddof=0),
            "TP+FN_mean": df_2["TP+FN_(_true_)"].astype(np.int).mean(),
            "TP+FN_std": df_2["TP+FN_(_true_)"].astype(np.int).std(ddof=0),
            "TP/(TP+FN)_mean": df_2["TP/(TP+FN)_(_recall_)"].astype(np.float).mean(),
            "TP/(TP+FN)_std": df_2["TP/(TP+FN)_(_recall_)"].astype(np.float).std(ddof=0),
            "TP+TN_mean": df_2["TP+TN_(_retrieved_)"].astype(np.int).mean(),
            "TP+TN_std": df_2["TP+TN_(_retrieved_)"].astype(np.int).std(ddof=0),
            "TP/(TP+TN)_mean": df_2["TP/(TP+TN)_(_precision_)"].astype(np.float).mean(),
            "TP/(TP+TN)_std": df_2["TP/(TP+TN)_(_precision_)"].astype(np.float).std(ddof=0),
            "F1_mean": df_2["F1"].astype(np.float).mean(),
            "F1_std": df_2["F1"].astype(np.float).std(ddof=0)})

    return pd.DataFrame(disease_summary)


if __name__ == "__main__":
    aggregate_datasets(DATASET_NAME="TNFa_2")
    # aggregate_datasets(DATASET_NAME="MCF7_2")



