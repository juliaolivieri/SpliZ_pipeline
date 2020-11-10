import argparse
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import time

def get_args():
  parser = argparse.ArgumentParser(description="estimate FDRs for outlying mzs")
  parser.add_argument("--dataname",help="name of dataset to use")
  parser.add_argument("--suffix",help="suffix to read in")
  parser.add_argument("--all_datanames",nargs="+",help="datanames to compare to")
  args = parser.parse_args()
  return args

def calc_FDR_mz(q,mz, pre_dfs,p, datanames):
  dfs = []
  for dataname, df in pre_dfs.items():
    df["sig"] = ((df[p] <= q) | (df[p] >= 1 - q)) & ((df["mz"] <= -mz) | (df["mz"] >= mz))
    df["plus"] = df["mz"] > 0
#     df["sig_plus"] = df["sig"] & df["plus"]
#     df["sig_minus"] = df["sig"] & ~(df["plus"])
    df = df[df["n_cells_ont"] > 20].drop_duplicates("ontology_gene")
    dfs.append(df)
  merge_df = dfs[0].merge(dfs[1][["mz","ontology_gene",p,"sig","plus","n_cells_ont"]],on="ontology_gene",suffixes=["_" + x for x in datanames[:2]])
#   full_merge = merge_df
#   for dataname in datanames:
#     for s in ["plus","minus"]:
#       merge_df["num_" + s + "_" + dataname] = merge_df["geneR1A_uniq"].map(merge_df.groupby("geneR1A_uniq")["sig_plus_" + dataname].sum())
  num_sig_both_all = merge_df[merge_df["sig_" + datanames[0]] & merge_df["sig_" + datanames[1]] ].shape[0]
  num_sig_concord_all = merge_df[merge_df["sig_" + datanames[0]] & merge_df["sig_" + datanames[1]] & (merge_df["plus_" + datanames[1]] == merge_df["plus_" + datanames[0]])].shape[0]
  try:
    return min(1,2*(1 - num_sig_concord_all/num_sig_both_all))
  except:
    
    # happens when this is more extreme than everything shared between the datasets
    return 0

def main():
  t0 = time.time()
  args = get_args()
  in_path = "scripts/output/significant_genes/"
  outpath = "scripts/output/final_FDRs_mz/"

  #[["TS_10x_redo","TSP2_10x_rerun_3prime","TSP1_SS2"],["TS_10x_redo_bestrefseq","TSP2_10x_rerun_3prime_bestrefseq","TSP1_SS2_bestrefseq"],["lemur_ss2","lemur_Antoine_4","lemur_Stumpy_4"]]
#  all_datanames = ["TS_10x_redo","TSP2_10x_rerun_3prime","TSP1_SS2"]
  p = "allp"

  out_cols = ["geneR1A_uniq","tissue","compartment","free_annotation","ontology",p,"mz","pval_adj"]

  df = pd.read_csv("{}{}-scZ_allp{}.tsv".format(in_path, args.dataname,args.suffix),usecols=["geneR1A_uniq","tissue","compartment","free_annotation",p,"mz","ontology","ontology_gene","cell_gene","n_cells_ont"],sep="\t")
  print("read in","{}{}-scZ_allp{}.tsv".format(in_path, args.dataname,args.suffix))


  for i in range(len(args.all_datanames)):
    print("i",i)
    for j in range(i + 1, len(args.all_datanames)):
      print("j",j)
      pre_dfs = {}
      datanames = [args.all_datanames[i],args.all_datanames[j]]
      print("datanames",datanames)
      print("time",time.time() - t0)

      for dataname in datanames:
        print("dataname",dataname)
        temp = pd.read_csv("{}{}-scZ_allp{}.tsv".format(in_path, dataname,args.suffix),usecols=["geneR1A_uniq","tissue","compartment","free_annotation",p,"mz","ontology","ontology_gene","cell_gene","n_cells_ont"],sep="\t")
        pre_dfs[dataname] = temp
    
    
      df["rev_q"] = 1 - df[p]
      df["q"] = df[["rev_q",p]].min(axis=1)
      df["eff_size"] = abs(df["mz"])  
    
      df["FDR_{}_{}".format(*datanames)] = df.apply(lambda row : calc_FDR_mz(row["q"], row["eff_size"],pre_dfs,p,datanames),axis=1)

      out_cols.append("FDR_{}_{}".format(*datanames))

#  df.to_csv("{}{}_FDR_test.tsv".format(outpath,dn),sep="\t",index=False)
  df["rev"] = 1 - df["allp"]
  df["pval"] = df[["rev","allp"]].min(axis=1)
  df["pval_adj"] = multipletests(df["pval"],method="fdr_bh")[1]

  df = df[out_cols].sort_values(by=out_cols[:-len(args.all_datanames)])
  print("saved at {}{}_FDR{}.tsv".format(outpath,args.dataname,args.suffix))

  df.to_csv("{}{}_FDR{}.tsv".format(outpath,args.dataname,args.suffix),sep="\t",index=False)
  print("wrote file",time.time() - t0)
main()
