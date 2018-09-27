import re
import os
import math
from numpy import *
from decimal import Decimal, ROUND_HALF_EVEN
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction
import numpy.random
from sklearn.datasets import fetch_mldata
import sklearn.preprocessing
import numpy as np
# import matplotlib as mpl
#mpl.use('Agg')
import scipy.special
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
# import matplotlib.pyplot as plt
# from matplotlib import style
# style.use("ggplot")
import time
# from matplotlib.ticker import FormatStrFormatter
import math
import logging
from ge_pca_by_samples import genes_pca_by_samples
from patient_sets_distribution_differences import patient_sets_distribution_differences
import json
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
#from clusters_similarity import *
from significant_enrichment import *
from differentially_expressed_genes import *
from gene_sets_correlation import *
from genes_correlation import find_genes_correlations
from gene_sets_overlapping import *
#from svm import *
import fre
from fre import *
from feature_selection import *
from utils.param_builder import *
from genes_clustering import *
from patients_clustering import *
from cox_penotype_analysis import *
from cox_gene_analysis import *
from genes_statistic import *
from ge_prediction_by_mutations import *
from gene_rank_correlation import gene_correlation_scores
from genes_pca import genes_pca
from mutations_pca import mutation_pca
from mutations_pca_by_samples import mutation_pca_by_samples
from patients_clustring_filtered_by_mutations import cluster_patients_filtered_by_mutation
from mutations_pca_by_samples import plot_pca_by_samples

import sys
sys.setrecursionlimit(50000)

######## mHGT ###############

# dataset = "UVM"
# data_normalizaton = "fpkm"
# tested_gene_file_name = "mito.txt"
# total_gene_file_name = "protein_coding.txt"
# gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# constants.update_dirs(CANCER_TYPE_u=dataset)
# groups_name = "mutation_binary"
# groups = json.load(file("groups/{}.json".format(groups_name)))
# N=len(load_gene_list(total_gene_file_name))
# B = len(load_gene_list(tested_gene_file_name))
# pval_preprocessing_file_name = "pvals_{}_{}_{}.txt".format(total_gene_file_name.split(".")[0], data_normalizaton, groups_name)
# hgt_preprocessing_file_name = "HGTs_out_{}_{}.npy".format(N,B)
# find_expression_significance(tested_gene_file_name=tested_gene_file_name, total_gene_file_name=total_gene_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, groups=groups, hgt_preprocessing_file_name = hgt_preprocessing_file_name, pval_preprocessing_file_name = pval_preprocessing_file_name, N=N, B = B)

######## DEG ###############

# dataset = "UVM"
# data_normalizaton = "fpkm"
# tested_gene_file_name = "protein_coding.txt"
# total_gene_file_name = "protein_coding.txt"
# gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# constants.update_dirs(CANCER_TYPE_u=dataset)
# groups_name = "additional_metastases_after_traetment"
# groups = json.load(file("groups/{}.json".format(groups_name)))
#
# deg(tested_gene_file_name=tested_gene_file_name, total_gene_file_name=total_gene_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, groups=groups, groups_name=groups_name)


######## RFE ###############

## randomized by size
#
#  dataset = "SKCM"
# data_normalizaton = "counts"
# gene_expression_file_name, phenotype_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#
# gene_list_file_name = "protein_coding.txt"
# rounds=100
# rank_method = DISTANCE
# recursion_number_of_steps=1
# permutation=RANDOMIZED
# gene_sets_sizes = [13, 6, 7, 13]
# results = []
# for i in range(5):
#     for cur_size in gene_sets_sizes:
#         results.append(RFE(gene_list_file_name, gene_expression_file_name, phenotype_file_name, rank_method = rank_method, rounds=rounds, recursion_step_size=cur_size, recursion_number_of_steps=recursion_number_of_steps, pval_preprocessing_file_name = pval_preprocessing_file_name, permutation=permutation))
#     fre.print_to_excel(results=results, gene_sets_sizes=gene_sets_sizes,rank_method = rank_method, permutation=permutation)



##################### FEATURE SELECTION ################

# dataset = "UVM"
# data_normalizaton = "fpkm"
# gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# constants.update_dirs(CANCER_TYPE_u=dataset)
# # gene_list_file_names = ["warburg.txt", "warburg_high.txt", "warburg_low.txt", "mito_warburg.txt", "mito.txt", "test.txt", "test_2.txt", "review_high.txt", "review_low.txt", "review_low.txt", "review.txt"]
# # gene_list_file_names = [ "warburg_high.txt", "warburg_low.txt", "ldha_singular.txt", "glycolysis_go.txt"]
# gene_list_file_names = ["fre_pr_uvm_additional_metastatses_1124.txt"]
# gene_filter_file_name = "protein_coding.txt"
# rounds=30
# recursion_step_size = 10
# feature_selection_method = "rfe"
# rank_method = DISTANCE
# score_method = "average_precision" # "roc_auc"
# target_genes_subset = "uvm_mito_part.txt"
# groups=json.load(file("groups/additional_metastases_after_traetment.json"))
# feature_selection(gene_list_file_name=gene_list_file_names, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, groups = groups, gene_filter_file_name=gene_filter_file_name, rounds=rounds, target_genes_subset = target_genes_subset, recursion_step_size=recursion_step_size, feature_selection_method=feature_selection_method, score_method=score_method)

##################### SETS CORRELATION ################

# tested_gene_list_file_name = ["mito.txt", "glycolysis_go.txt", "warburg_high.txt", "warburg_low.txt", "ldha_singular.txt"]
#
# dataset = "UVM"
# constants.update_dirs(CANCER_TYPE_u=dataset)
# data_normalizaton = "fpkm_normalized_by_genes_standardization"
# gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name,  pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# find_sets_correlations(tested_gene_list_file_name=tested_gene_list_file_name,
#                        total_gene_list_file_name="protein_coding.txt", gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name)

##################### GENES CORRELATION ################

# dataset = "UVM"
# constants.update_dirs(CANCER_TYPE_u=dataset)
# data_normalizaton = "fpkm"
# gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# filter_expression = None
# # filter_expression = json.load(file("filters/mutation_binary_1.json"))
# intersection_gene_file_names= ["uvm_mito_part"] # ["uvm_mito_inner_membrane.txt", "uvm_mito_membrane.txt", "uvm_mito_protein_complex.txt", "uvm_mito_matrix.txt", "uvm_mitochondrion.txt", "uvm_mito_part.txt"]
# total_gene_list_file_name = "protein_coding_included.txt"
# output_f = open(os.path.join(constants.OUTPUT_DIR,"mir_HG_{}.txt".format(time.time())), "w+")
# counter = 0
# # included_mirs = ["hsa-mir-2116", "hsa-mir-3156-1", "hsa-mir-513a-1", "hsa-mir-513a-2", "hsa-mir-4763", "hsa-mir-513c", "hsa-mir-5003", "hsa-mir-548s", "hsa-mir-6892", "hsa-mir-1468", "hsa-mir-548b", "hsa-mir-513b", "hsa-mir-940", "hsa-mir-4797", "hsa-mir-4725", "hsa-mir-125b-2", "hsa-mir-224", "hsa-mir-452", "hsa-mir-887", "hsa-mir-5588", "hsa-mir-6783", "hsa-mir-510", "hsa-mir-143", "hsa-mir-193b", "hsa-mir-4781", "hsa-mir-766", "hsa-mir-4709", "hsa-mir-1228", "hsa-mir-7974", "hsa-mir-1293", "hsa-mir-6764", "hsa-mir-1301", "hsa-mir-935", "hsa-mir-6509", "hsa-mir-514b", "hsa-mir-4758", "hsa-mir-508", "hsa-mir-937", "hsa-mir-103a-2", "hsa-mir-6728", "hsa-mir-155", "hsa-mir-148b", "hsa-mir-877", "hsa-mir-6811", "hsa-mir-3156-2", "hsa-mir-873", "hsa-mir-324", "hsa-mir-200a", "hsa-mir-1910", "hsa-mir-200b", "hsa-mir-181b-2", "hsa-mir-4661", "hsa-mir-181a-2", "hsa-let-7c", "hsa-mir-3682", "hsa-mir-181b-1", "hsa-mir-145", "hsa-mir-514a-3", "hsa-mir-192", "hsa-mir-652", "hsa-mir-509-1", "hsa-mir-29c", "hsa-mir-509-2", "hsa-mir-509-3", "hsa-mir-5579", "hsa-mir-4740", "hsa-mir-3117", "hsa-mir-211", "hsa-mir-501", "hsa-mir-3940", "hsa-mir-514a-1", "hsa-mir-10b", "hsa-mir-514a-2", "hsa-mir-6507", "hsa-mir-30a", "hsa-mir-3186", "hsa-mir-24-1", "hsa-mir-24-2", "hsa-mir-1249", "hsa-mir-939", "hsa-mir-365b", "hsa-mir-197", "hsa-mir-6837", "hsa-mir-506", "hsa-mir-335", "hsa-mir-365a", "hsa-mir-6516", "hsa-mir-943", "hsa-mir-27b", "hsa-mir-6715a", "hsa-mir-3619", "hsa-mir-200c", "hsa-mir-548v", "hsa-mir-194-2", "hsa-mir-1307", "hsa-mir-6802", "hsa-mir-503", "hsa-mir-5683", "hsa-mir-708", "hsa-mir-4645", "hsa-mir-1231", "hsa-mir-1296", "hsa-mir-378g", "hsa-mir-561", "hsa-mir-221", "hsa-mir-3190", "hsa-mir-4665", "hsa-mir-218-2", "hsa-mir-330", "hsa-mir-222", "hsa-mir-507", "hsa-mir-3679", "hsa-mir-484", "hsa-mir-412", "hsa-mir-342", "hsa-mir-876", "hsa-mir-4638", "hsa-mir-188", "hsa-mir-4470", "hsa-mir-30c-2", "hsa-mir-4733", "hsa-mir-363", "hsa-mir-130a", "hsa-mir-592", "hsa-mir-132", "hsa-let-7g", "hsa-mir-99a", "hsa-mir-6882", "hsa-mir-214", "hsa-mir-4422", "hsa-mir-504", "hsa-mir-26a-2", "hsa-mir-219a-1", "hsa-mir-670", "hsa-mir-26a-1", "hsa-mir-4723", "hsa-mir-4726", "hsa-mir-942", "hsa-mir-15a", "hsa-mir-3189", "hsa-mir-497", "hsa-mir-1262", "hsa-mir-3187", "hsa-mir-218-1", "hsa-mir-4647", "hsa-mir-466", "hsa-mir-6785", "hsa-mir-181c", "hsa-mir-6503", "hsa-mir-6850", "hsa-mir-6858", "hsa-mir-146b", "hsa-mir-3153"]
# included_mirs = ["hsa-mir-548b"]
# mir_total_list = load_dictionary("mir_to_mrna.txt")
#
#
# #####
#
# # total_gene_list_file_name = load_gene_list(total_gene_list_file_name)
# # intersection_gene_sets = []
# # if intersection_gene_file_names is not None:
# #     intersection_gene_file_names = [
# #         np.array([y.split(".")[0] for y in load_gene_list(x)]) if type(x) == str else [y.split(".")[0] for y in x] for x
# #         in intersection_gene_file_names]
#
# #####
#
# for cur_mir in mir_total_list:
#     if not cur_mir[0] in  included_mirs: continue
#     print "current mir: {}".format(cur_mir[0])
#     # with file(os.path.join(constants.LIST_DIR, "dummy_mir_1.txt"),"w+") as f:
#     #     f.write(cur_mir[0])
#     # with file(os.path.join(constants.LIST_DIR, "dummy_genes_1.txt"), "w+") as f:
#     #     f.write("\r\n".join(cur_mir[1:]))
#
#     tested_gene_list_file_names=[cur_mir[1:], [cur_mir[0]]]
#     var_th_index=None
#     gene_expression_file_names = [gene_expression_file_name, mirna_file_name]
#
#     # output_f.write(load_gene_list("dummy_mir_1.txt")[0]+ "\t")
#     data_output, hg_output = find_genes_correlations(tested_gene_list_file_names=tested_gene_list_file_names, total_gene_list_file_name=total_gene_list_file_name,
# gene_expression_file_names=gene_expression_file_names, intersection_gene_file_names=intersection_gene_file_names, phenotype_file_name=phenotype_file_name, filter_expression=filter_expression, var_th_index=var_th_index)
#     output_f.write("\t".join([str(x) for x in hg_output])+"\r\n")
#     counter +=1
#     # if counter ==2:
#     #     break
# output_f.close()
#


##################### HG ENRICHMENT #########################

# dataset = "UVM"
# tested_gene_file_name = "rfe_uvm_10_3014.txt" # "warburg_low.txt"
# total_gene_file_name = "protein_coding.txt" # "rfe_uvm_10_3014.txt"
# constants.update_dirs(CANCER_TYPE_u=dataset)
# check_group_enrichment(tested_gene_file_name,total_gene_file_name)


##################### GO ENRICHMENT #########################

# dataset = "UVM"
# tested_gene_file_name = "glycolysis_go.txt" # "rfe_uvm_10_3014.txt"
# total_gene_file_name = "rfe_uvm_10_3014.txt" # "protein_coding.txt" #
# constants.update_dirs(CANCER_TYPE_u=dataset)
# check_group_enrichment(tested_gene_file_name,total_gene_file_name)

##################### GENE CORRELATION SCORE ################

# gene_correlation_scores(tested_gene_list_file_name="mito.txt",
#                        total_gene_list_file_name="protein_coding.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", top_n=20)

##################### OVERLAPING ################

# find_sets_overlaps(tested_gene_list_file_name=["tca_pathcards.txt", "oxidative_phosphorylation_pathcards.txt", "glycolysis_pathcards.txt", "mito.txt", "oxidative.txt", "review.txt", "review_high.txt", "review_low.txt", "warburg_high.txt", "warburg_low.txt", "warburg.txt", "oxidative_HIF.txt", "pyruvate.txt", "ldha_singular.txt"],
#                        total_gene_list_file_name="protein_coding.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", phenotype_file_name="TCGA-SKCM.GDC_phenotype.tsv")


##################### GENES CLUSTERING AND ENRICHMENT ################

# tested_gene_list_file_name="proliferation_invasion.txt"
# total_gene_list_file_name="proliferation_invasion.txt"
# var_th_index=None
# start_k=2
# end_k=4
# dataset="SKCM"
# data_normalizaton="fpkm"
# constants.update_dirs(CANCER_TYPE_u=dataset)
# calc_go=False
# enrichment_list_file_names = ["invasion.txt", "proliferation.txt", "proliferation_invasion.txt"]
# gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# find_clusters_and_gene_enrichment(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, var_th_index=var_th_index, start_k=start_k, end_k=end_k, calc_go=calc_go, enrichment_list_file_names = enrichment_list_file_names)

##################### CLUSTERS AND SURVIVAL ################


# # for x in range(10):
# #     random_set_file_name = generate_random_set(random_size=13, meta_gene_set="protein_coding_long.txt")
#
# for dataset in ["ESCA",  "ACC", "CHOL", "BLCA", "BRCA", "CESC", "COAD", "UCEC", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "DLBC", "LIHC", "LGG", "LUAD", "LUSC", "SKCM", "MESO", "UVM", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "STAD", "TGCT", "THYM", "THCA", "UCS"]:
# for dataset in ["SKCM"]:
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=3
#     end_k=3
#     #[{"gender.demographic": {"type": "string", "value": ["male"]}},
#                                        # {"gender.demographic": {"type": "string", "value": ["female"]}}]
#     filter_expression = None
#     meta_groups = []
#     meta_groups=[
#         [
#             {"sample_type.samples": {"type": "string", "value": ["Primary Tumor", "Metastatic"]
#                                      },
#              "tumor_stage.diagnoses": {"type": "string", "value": ["stage 0", "stage i", "stage ia", "stage ib", "stage ic"]
#                                      }
#              },
#             {"sample_type.samples": {"type": "string", "value": ["Primary Tumor", "Metastatic"]
#                                      },
#              "tumor_stage.diagnoses": {"type": "string",
#                                        "value": ["stage iv", "stage iii", "stage iiia", "stage iiib", "stage iiic"]
#                                        }
#              }
#             ]
#     ]
#
#     # , "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]}
#
#     # for cur_tt in ["Primary Tumor"]:
#     filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor", "Metastatic"]},
#                           "tumor_stage.diagnoses": {"type": "string",
#                                                     "value": ["stage 0", "stage i", "stage ia", "stage ib", "stage ic", "stage iv", "stage iii", "stage iiia", "stage iiib",
#                                                               "stage iiic"]
#                                                     }
#                           }]
#     print "process {}".format(dataset)
#     find_clusters_and_survival(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups)

#################### SVM PREDICTION ##########################

# dataset = "UVM"
# data_normalizaton = "fpkm"
# gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# constants.update_dirs(CANCER_TYPE_u=dataset)
# # gene_list_file_names = ["warburg.txt", "warburg_high.txt", "warburg_low.txt", "mito_warburg.txt", "mito.txt", "test.txt", "test_2.txt", "review_high.txt", "review_low.txt", "review_low.txt", "review.txt"]
# # gene_list_file_names = ["uvm_mito_matrix.txt", "uvm_mito_membrane.txt", "uvm_mito_part.txt", "uvm_mito_protein_complex.txt", "uvm_mitochondrion.txt"]
# gene_list_file_names = ["mito.txt"]
# gene_filter_file_name = "protein_coding.txt"
# rounds=100
# rank_method = DISTANCE
# results = []
# groups=json.load(file("groups/additional_metastases_after_traetment.json"))
# prediction_by_gene_expression(gene_list_file_names=gene_list_file_names, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, gene_filter_file_name=gene_filter_file_name, rounds=rounds, rank_method=rank_method, labels_permutation=constants.LABELS_NORMAL, compare_to_random = True,
#                               groups=groups)



#################### COX PHENOTYPE ANALYSIS ##########################

# dataset = "SKCM"
# data_normalizaton = "counts"
# gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# # gene_list_file_names = ["warburg.txt", "warburg_high.txt", "warburg_low.txt", "mito_warburg.txt", "mito.txt", "test.txt", "test_2.txt", "review_high.txt", "review_low.txt", "review_low.txt", "review.txt"]
# test_independently = False
# filtered_out = ['bcr', 'days_to_initial_pathologic_diagnosis', 'disease_code', 'informed_consent_verified', 'project_code', 'withdrawn', 'year_of_dcc_upload', "therapy_type_notes", "pathology_report_file_name", "tx_on_clinical_trial", "units", "year_of_form_completion", "vial_number", "ethnicity.demographic", "is_ffpe.samples", "project.tissue_source_site", "name.tissue_source_site","bcr_id.tissue_source_site", "project_id.project", "submitter_id", "therapy_type", "therapy_type_notes", 'classification_of_tumor.diagnoses', 'last_known_disease_status.diagnoses', 'prior_malignancy.diagnoses', 'progression_or_recurrence.diagnoses', 'tumor_grade.diagnoses', 'disease_type', 'primary_site', 'disease_type.project', 'name.project', 'primary_site.project', "vital_status.diagnoses", "day_of_dcc_upload", "month_of_dcc_upload", "tissue_retrospective_collection_indicator", "tissue_prospective_collection_indicator", "code.tissue_source_site", "tissue_source_site", "site_of_resection_or_biopsy.diagnoses", "disease_type", "primary_site", "disease_type.project", "name.project", "primary_site.project", "project_id.project", "submitter_id", "bcr_id.tissue_source_site", "code.tissue_source_site", "name.tissue_source_site", "project.tissue_source_site", "is_ffpe.samples", "oct_embedded.samples", "sample_type_id.samples", "state.samples", "system_version",  "age_at_diagnosis.diagnoses", "day_of_form_completion", "patient_id", "batch_number"]
# filtered_in = ["sample_type.samples", 'gender.demographic', "pathologic_M","pathologic_T","pathologic_N", "tumor_stage.diagnoses", "melanoma_clark_level_value", "days_to_birth.diagnoses" ,"weight.exposures","height.exposures", "person_neoplasm_cancer_status"]
# filter_type = FILTER_IN
# filter_na_by_rows = True
# f = file(os.path.join(constants.OUTPUT_DIR,"output_expression"), 'w+')
# gene_expression_file_name, phenotype_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# cox_phenotype(test_independently=test_independently, filtered_out=filtered_out, filtered_in=filtered_in, filter_type=filter_type, filter_na_by_rows=filter_na_by_rows, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name)
# f.close()

#################### COX GENE ANALYSIS ##########################
# dataset = "UVM"
# data_normalizaton = "fpkm"
# gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# constants.update_dirs(CANCER_TYPE_u=dataset)
# test_independently = False
# filtered_out = []
# filtered_in = []
# tested_genes_file_name = "warburg_high.txt"
# filter_type = FILTER_IN
# var_th_index = None
# filter_na_by_rows = True
# filter_expression = None
# f = file(os.path.join(constants.BASE_PROFILE, "output","cox_genes_{}_{}.txt".format(dataset, time.time())), 'w+')
# # filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor"]}, "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR", "TUMOR FREE", ""]}}]
# results = cox_gene(test_independently, filtered_out, filtered_in, filter_type, filter_na_by_rows,tested_genes_file_name = tested_genes_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, filter_expression=filter_expression, var_th_index=var_th_index)
# f.write(results)
# f.close()


##################### CLUSTERS AND SURVIVAL ################

# for dataset in ["PAAD"]:
#     meta_groups=[]
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=3
#     end_k=3
#     filter_expression = None
#     d_codes = []
#     # for cur in meta_groups[0]:
#     #     d_codes.append(cur["disease_code"]["value"][0])
#
#     # for cur_tt in ["Primary Tumor"]:
#     filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor","Metastatic"]},
#                           # "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]} ,
#                           # "disease_code": {"type": "string",
#                           #              "value": d_codes
#                           #              }
#                           }]
#     print "process {}".format(dataset)
#     find_clusters_and_survival(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups)

############# PCA for gene expression ##################

# for dataset in ["PANCAN"]:
#     meta_groups= [] # json.load("package_utils/tumor_type_groups.json")
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=2
#     end_k=2
#     #[{"gender.demographic": {"type": "string", "value": ["male"]}},
#                                        # {"gender.demographic": {"type": "string", "value": ["female"]}}]
#     filter_expression = None
#
#     d_codes = []
#     for cur in meta_groups[0]:
#         d_codes.append(cur["disease_code"]["value"][0])
#
#     # , "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]}
#
#     # for cur_tt in ["Primary Tumor"]:
#     filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor","Metastatic"]},
#                           # "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]} ,
#                           # "disease_code": {"type": "string",
#                           #              "value": d_codes
#                           #              }
#                           }]
#     print "process {}".format(dataset)
#     genes_pca(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups)


############ PREDICTION_BY_MUTATION ###################
# for cur_suffix in ["gt","gly","lac","tca"]:
#     for cur_dir in ["high","low"]:

for cur_unsupervised in [True]:
    # for cur_tested_mir in ["hsa-mir-513a-1", "hsa-mir-513a-2", "hsa-mir-4781", "hsa-mir-4797", "hsa-mir-548b", "hsa-mir-2116", "hsa-mir-510", "hsa-mir-5588", "hsa-mir-548s", "hsa-mir-1468", "hsa-mir-561", "hsa-mir-513b", "hsa-mir-125b-2", "hsa-mir-1262", "hsa-mir-6509", "hsa-mir-6507", "hsa-mir-181b-2", "hsa-mir-181b-1", "hsa-mir-873", "hsa-mir-221", "hsa-mir-181a-2", "hsa-mir-508", "hsa-mir-192", "hsa-mir-1910", "hsa-mir-548y", "hsa-mir-506", "hsa-mir-197", "hsa-mir-5683", "hsa-mir-363", "hsa-mir-514a-3", "hsa-mir-876", "hsa-mir-514a-1", "hsa-mir-378g", "hsa-mir-514a-2", "hsa-mir-514b", "hsa-mir-548v", "hsa-mir-509-1", "hsa-mir-30a", "hsa-mir-509-2", "hsa-mir-509-3", "hsa-mir-935", "hsa-mir-507", "hsa-mir-942", "hsa-mir-211", "hsa-let-7c", "hsa-mir-30c-2", "hsa-mir-99a", "hsa-mir-29c", "hsa-mir-194-2", "hsa-mir-651", "hsa-mir-335", "hsa-mir-219a-1", "hsa-mir-222", "hsa-mir-1249", "hsa-mir-203a", "hsa-mir-4733"]:
    for cur_tested_file in ["mir-30c-2.txt"]: # ["uvm_mito_inner_membrane.txt", "uvm_mito_matrix.txt", "uvm_mito_membrane.txt", "uvm_mito_part.txt", "uvm_mito_protein_complex.txt", "uvm_mitochondrion.txt"]: #  ["hsa-mir-181b-1", "hsa-mir-181b-2", "hsa-mir-1468" ]:
    # for cur_tested_mir in ["hsa-mir-194-2"]:
    #     cur_tested_file = cur_tested_mir[cur_tested_mir.index("-")+1:]+".txt"
    #     with file(os.path.join(constants.LIST_DIR, cur_tested_file), "w") as f:
    #         f.write(cur_tested_mir)
        for cur_json in ["additional_metastases_after_traetment"]:
            for dataset in ["UVM"]:
                # if dataset == "PANCAN": continue
                meta_groups = None
                meta_groups=[json.load(file("groups/{}.json".format(cur_json)))]

                constants.update_dirs(CANCER_TYPE_u=dataset)
                data_normalizaton = "fpkm"
                gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
                random_set_file_name = generate_random_set(random_size=1, meta_gene_set="mir_total.txt")
                tested_gene_list_file_name=  cur_tested_file # ""mir_warburg_{}_{}.txt".format(cur_dir, cur_suffix) # random_set_file_name #
                total_gene_list_file_name= None # "protein_coding_long.txt"
                var_th_index= None# 70
                is_unsupervised=cur_unsupervised
                integ=False
                start_k=2
                end_k=4
                min_ratio = 0.01
                excluded_mutation_gene_list = None # ["TTN", "BRAF", "GNAQ", "GNA11"]
                included_mutation_gene_list = "uvm_mutations.txt"

                d_codes = []
                # if meta_groups is not None:
                #     for cur in meta_groups[0]:
                #         d_codes.append(cur["tumor_stage.diagnoses"]["value"][0])
                # , "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]}

                # for cur_tt in ["Primary Tumor"]:
                filter_expression = None
                filter_expression =  json.load(file("filters/{}.json".format(cur_json)))
                print "process {}".format(dataset)
                phenotype_labels_heatmap = None#["breslow_depth_value"]
                clustering_algorithm = "euclidean"
                # try:
                # find_clusters_and_gene_enrichment(tested_gene_list_file_name=tested_gene_list_file_name,
                #                                   total_gene_list_file_name=tested_gene_list_file_name,
                #                                   gene_expression_file_name=mirna_file_name,
                #                                   phenotype_file_name=phenotype_file_name,
                #                                   var_th_index=var_th_index, start_k=start_k, end_k=end_k,
                #                                   calc_go=False,
                #                                   enrichment_list_file_names=None, meta_groups=meta_groups, filter_expression=filter_expression,
                #                                   cluster_algorithm = "hierarchical")
                find_clusters_and_survival(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=mirna_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, clustering_algorithm = clustering_algorithm)
                # for k in range(10):
                # predict_ge_by_mutation(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, phenotype_labels_heatmap=phenotype_labels_heatmap, integ=integ, min_ratio=min_ratio, included_mutation_gene_list=included_mutation_gene_list, excluded_mutation_gene_list=excluded_mutation_gene_list)
                # except:
                #     pass
                # patient_sets_distribution_differences(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, clustering_algorithm = clustering_algorithm)
                #cluster_patients_filtered_by_mutation(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, integ=True ,min_ratio=0.4)

############ mutation PCA ##################

# for dataset in ["PANCAN"]:
#     meta_groups = json.load("package_utils/tumor_type_groups.json")
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=2
#     end_k=2
#     filter_expression = None
#     d_codes = []
#     # for cur in meta_groups[0]:
#     #     d_codes.append(cur["disease_code"]["value"][0])
#
#     filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor","Metastatic"]},
#                           # "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]} ,
#                           # "disease_code": {"type": "string",
#                           #              "value": d_codes
#                           #              }
#                           }]
#     print "process {}".format(dataset)
#     mutation_pca(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups)

############ mutation PCA ##################

# lists =  ["protein_coding.txt"]# ["mito.txt", "uvm_mito_inner_membrane.txt", "uvm_mito_membrane.txt", "uvm_mito_protein_complex.txt", "uvm_mito_matrix.txt", "uvm_mitochondrion.txt", ]
# json_file = "tss_tcga_skcm"
#
# for dataset in ["SKCM"]:
#     for cur_list in lists:
#         # meta_groups = json.load("package_utils/tumor_type_groups.json")
#         constants.update_dirs(CANCER_TYPE_u=dataset)
#         data_normalizaton = "fpkm"
#         gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mutation_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#         tested_gene_list_file_name=cur_list
#         total_gene_list_file_name="protein_coding_long.txt"
#         var_th_index=2000
#         is_unsupervised=True
#         start_k=2
#         end_k=2
#         n_components = 2
#         filter_expression = None
#         d_codes = []
#         # for cur in meta_groups[0]:
#         #     d_codes.append(cur["disease_code"]["value"][0])
#
#         filter_expression = json.load(file("filters/{}.json".format(json_file)))
#         groups = [json.load(file("groups/{}.json".format(json_file)))]
#         print "process {}".format(dataset)
#         genes_pca_by_samples(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name,  var_th_index=var_th_index, filter_expression= filter_expression, meta_groups = groups, n_components=n_components)

############ PREDICTION_BY_MUTATION ###################

# for dataset in ["UVM"]:
#     if dataset == "PANCAN": continue
#     meta_groups = None
#     meta_groups= [json.load(file("groups/mutation_cluster.json"))]
#
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     random_set_file_name = generate_random_set(random_size=25, meta_gene_set="protein_coding_long.txt")
#     tested_gene_list_file_name= "warburg_high.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     filter_expression = None
#     print "process {}".format(dataset)
#     genes_pca_by_samples(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, filter_expression= filter_expression, meta_groups = meta_groups)
#     # predict_ge_by_mutation(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, phenotype_labels_heatmap=phenotype_labels_heatmap, integ=integ, min_ratio=0.9, omitted_genes=omitted_genes)
#     # cluster_patients_filtered_by_mutation(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, integ=True ,min_ratio=0.4)
