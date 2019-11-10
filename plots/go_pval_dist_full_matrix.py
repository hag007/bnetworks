from matplotlib import style



style.use("ggplot")
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt

from statsmodels.sandbox.stats.multicomp import fdrcorrection0

from matplotlib_venn import venn2

fontsize=30

font = {'size'   : fontsize}
mpl.rc('xtick', labelsize=fontsize)    # fontsize of the tick labels
mpl.rc('ytick', labelsize=fontsize)


def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def false_positive_example(algo="jactivemodules_greedy",dataset="TNFa_2", df_summary=pd.DataFrame()):

    output=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX","emp_diff_modules_{}_{}.tsv".format(dataset, algo)), sep='\t', error_bad_lines=False)
    output_md = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_modules_{}_{}_md.tsv".format(dataset, algo)), sep='\t', error_bad_lines=False)
    output_md_ids=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500])]['GO id']
    output=output.loc[output["GO id"].isin(output_md_ids.values)]
    output_md=output_md[output_md["GO id"].isin(output_md_ids)]

    n_of_bg=1
    for cur_bg in np.arange(n_of_bg):

        output["cur_bg_hg_pval"]=output.apply(lambda a: float(a["dist_n_samples"][1:-1].split(", ")[cur_bg]) if not pd.isnull(a["GO name"]) and a["dist_n_samples"].startswith("[") else np.nan ,axis=1)
        bg_pvals=output[output["cur_bg_hg_pval"] != -1].dropna()


    bg_sig_ids=bg_pvals[fdrcorrection0([10**-a for a in  bg_pvals["cur_bg_hg_pval"]], alpha=0.05, method='indep', is_sorted=False)[0]]["GO id"]

    real_pvals=output_md[output_md["hg_pval_max"].apply(lambda a: isFloat(a))].dropna()
    real_pvals["hg_pval_max"]=real_pvals["hg_pval_max"].values.astype(np.float)

    values=real_pvals["hg_pval_max"][real_pvals["hg_pval_max"]!=-1]
    if is_corrected:
        values=-np.log10(fdrcorrection0(10**-values, alpha=0.05, method='indep', is_sorted=False)[1])


    real_sig_ids = real_pvals[fdrcorrection0([10**-a for a in real_pvals["hg_pval_max"]], alpha=0.05, method='indep', is_sorted=False)[0]]["GO id"]
    # out=venn2([set(real_sig_ids), set(bg_sig_ids)], set_labels=('Original\ndataset', 'Permuted\ndataset'), ax=axs[1])
    n_intersection = len(set(real_sig_ids).intersection(set(bg_sig_ids)))
    n_real = len(real_sig_ids)
    df_summary.loc[dataset,algo]="{}\{}={}".format(n_intersection,n_real,float(n_intersection)/round(n_real) )



if __name__ == "__main__":

    for is_corrected in [False, True]:
        suffix = ("qval" if is_corrected else "pval")
        datasets=["TNFa_2", "Schizophrenia.G50"]  #, "Schizophrenia.G50"
        algos=["jactivemodules_greedy","netbox"] #
        i=0
        fig,axs=plt.subplots(len(datasets)*len(algos), 2, figsize=(26, 8*len(datasets)*len(algos)))
        for dataset in datasets:
            for algo in algos:
                df_summary=pd.DataFrame()
                false_positive_example(algo=algo, dataset=dataset,df_summary=df_summary) # "Breast_Cancer.G50"
                i += 1

        df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_12_matrix.tsv"), sep='\t', index_label="algo")
