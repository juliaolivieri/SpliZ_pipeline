import argparse
import numpy as np
import pandas as pd

def get_args():
  parser = argparse.ArgumentParser(description="Create final summary file")
  parser.add_argument("--dataname",help="name of dataset to use")
  parser.add_argument("--suffix",help="suffix to save with")
  parser.add_argument("--group_col",help="column to group the data by (e.g. ontology, compartment, tissue)",default="ontology")
  parser.add_argument("--sub_col",help="subset data by this column before checking for differences (e.g. tissue, compartment)",default="dummy")
  parser.add_argument("--num_perms",type=int,help="number of permutations to run for")

  args = parser.parse_args()
  return args


def main():
  args = get_args()
  # load in data
  in_dir = "scripts/output/"

  pval_df = pd.read_csv("{}variance_adjusted_permutations/{}_pvals_{}-{}_{}{}.tsv".format(in_dir, args.dataname,args.group_col,args.sub_col,args.num_perms, args.suffix), sep = "\t")
  
  splizsite_dfs = []
  for prefix in ["","second_evec_","third_evec_"]:
    splizsite_dfs.append(pd.read_csv("{}SpliZsites/{}{}_{}-{}_{}{}.tsv".format(in_dir,prefix,args.dataname,args.group_col,args.sub_col,args.num_perms, args.suffix),sep="\t"))
  splizsite_df = pd.concat(splizsite_dfs,axis=0).drop_duplicates()
  
  df = pd.read_csv("{}rijk_zscore/{}_sym_SVD_normdonor{}_subcol.tsv".format(in_dir,args.dataname,args.suffix),sep="\t")
  if (args.sub_col == "tiss_comp") & (args.sub_col not in df.columns):
    df["tiss_comp"] = df["tissue"] + df["compartment"]
  elif args.sub_col == "dummy":
    df["dummy"] = "dummy"

  # combine outputs
  out_dict = {"gene" : [],"sub_col" : [], "group_col" : [],  "SpliZsites" : []}
  z_cols = ["scZ","svd_z0","svd_z1","svd_z2"]
  
  for z_col in z_cols:
    out_dict["{}_median".format(z_col)] = []
    out_dict["{}_pval".format(z_col)] = []
  
  for gene, gene_df in df.groupby("geneR1A_uniq"):
    for tiss, tiss_df in gene_df.groupby(args.sub_col):
      for ont, ont_df in tiss_df.groupby(args.group_col):
        out_dict["gene"].append(gene)
        out_dict["sub_col"].append(tiss)
        out_dict["group_col"].append(ont)
        out_dict["SpliZsites"].append(",".join([str(x) for x in splizsite_df[splizsite_df["gene"] == gene]["end"]]))
        
        
        for z_col in z_cols:
  
          out_dict["{}_median".format(z_col)].append(ont_df[z_col].median())
          try:
            pval = pval_df[(pval_df["geneR1A_uniq"] == gene) & (pval_df["sub_col"] == tiss)]["perm_pval_adj_{}".format(z_col)].iloc[0]
          except:
            pval = np.nan
          out_dict["{}_pval".format(z_col)].append(pval)
  out_df = pd.DataFrame.from_dict(out_dict)
  out_df = out_df.sort_values(["gene","sub_col","scZ_median"])
  out_df.to_csv("{}/final_summary/summary_{}_{}-{}_{}{}.tsv".format(in_dir,args.dataname,args.group_col,args.sub_col,args.num_perms, args.suffix),sep="\t",index=False)

main()
