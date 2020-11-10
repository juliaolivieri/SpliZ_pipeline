import argparse
import math
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import time

def get_args():
  parser = argparse.ArgumentParser(description="estimate FDRs for anova significance")
  parser.add_argument("--dataname",help="name of dataset to use")
  parser.add_argument("--suffix",help="suffix to read in")
  parser.add_argument("--all_datanames",nargs="+",help="datanames to compare to")
  args = parser.parse_args()
  return args

def get_out_df_factor(df_sig,df_coeff):
  out_dict = {"gene" : [], "variable" : [], "sig" : [], "coeff" : []}
  for gene in df_sig.index:
    sig_val = df_sig.loc[gene,:]["compartment"]
    if not math.isnan(sig_val):
      coeff = df_coeff.loc[gene,:][~df_coeff.loc[gene,:].isna()]

      for j, coeff_val in coeff.iteritems():
        if j.startswith("compartment"):
          out_dict["gene"].append(gene)
          out_dict["variable"].append(j[11:])
          out_dict["sig"].append(sig_val)
          out_dict["coeff"].append(coeff_val)
  out_df = pd.DataFrame.from_dict(out_dict)
  out_df["gene_variable"] = out_df["gene"] + out_df["variable"]

  return out_df

def calc_FDR_anova(eps, eff_size,pre_dfs,datanames):
  out_dfs = {}
  for dataname, out_df in pre_dfs.items():
    out_df["diff"] = (out_df["sig"] < eps) & (abs(out_df["coeff"]) > eff_size)
    out_dfs[dataname] = out_df
  merged = out_dfs[datanames[0]].merge(out_dfs[datanames[1]][["gene_variable","sig","coeff","diff"]],on="gene_variable",suffixes=["_" + x for x in datanames])
  merged_sig = merged[merged["diff_" + datanames[0]] & merged["diff_" + datanames[1]]]
  merged_sig["concord"] = merged_sig["coeff_" + datanames[0]]*merged_sig["coeff_" + datanames[1]] > 0
  try:
    return min(1,2*(1 - merged_sig["concord"].value_counts()[True]/merged_sig.shape[0]))
  except:
    return 0

def main():
  args = get_args()
  t0 = time.time()
  in_path = "scripts/output/anova_factor/"
  resid_suff = ""
  model_types = ["","_unweight"]

  for mt in model_types:
    out_cols = []
    sig = pd.read_csv("scripts/output/significant_genes/{}-scZ_allp{}.tsv".format(args.dataname,args.suffix),sep="\t")
  
    outpath = "scripts/output/final_FDRs_anova_factor/"
  #  all_datanames = ["TS_10x_redo","TSP2_10x_rerun_3prime","TSP1_SS2"]
  
  #  dn = "TS_10x_redo"
  
    df_sig = pd.read_csv("{}{}_scZ_anova_out{}{}{}.tsv".format(in_path,args.dataname,mt,args.suffix,resid_suff),sep="\t",index_col=0)
    df_coeff = pd.read_csv("{}{}_scZ_anova_coeff{}{}{}.tsv".format(in_path,args.dataname,mt,args.suffix,resid_suff),sep="\t", index_col = 0)
    df = get_out_df_factor(df_sig, df_coeff)
    for i in range(len(args.all_datanames)):
      for j in range(i + 1, len(args.all_datanames)):
        pre_dfs = {}
        datanames = [args.all_datanames[i],args.all_datanames[j]]
        print("datanames",datanames)
        print("time",time.time() - t0)
  
        for dataname in datanames:
  
          df_sig = pd.read_csv("{}{}_scZ_anova_out{}{}{}.tsv".format(in_path,dataname,mt,args.suffix,resid_suff),sep="\t",index_col=0)
          df_coeff = pd.read_csv("{}{}_scZ_anova_coeff{}{}{}.tsv".format(in_path,dataname,mt,args.suffix,resid_suff),sep="\t", index_col = 0)
          out_df = get_out_df_factor(df_sig, df_coeff)
          pre_dfs[dataname] = out_df
      
      
        df["FDR_{}_{}".format(*datanames)] = df.apply(lambda row : calc_FDR_anova(row["sig"], row["coeff"],pre_dfs,datanames),axis=1)
        out_cols.append("FDR_{}_{}".format(*datanames))
  
  #  df.to_csv("{}{}_FDR_test.tsv".format(outpath,args.dataname),sep="\t",index=False)
    df = df.sort_values(by=out_cols)
    df["pval_adj"] = multipletests(df["sig"],method="fdr_bh")[1]
    df.to_csv("{}{}_FDR{}{}{}.tsv".format(outpath,args.dataname,mt,args.suffix,resid_suff),sep="\t",index=False)
    print("wrote file",time.time() - t0)
main()
