import argparse
import pandas as pd

def get_args():
  parser = argparse.ArgumentParser(description="calculate splicing scores per gene/cell")
  parser.add_argument("--dataname",help="name of dataset to use")
  parser.add_argument("--pinning_S",type=float,help="pinning level for S_ijks")
  parser.add_argument("--pinning_z",type=float,help="pinning level for zs")

#  parser.add_argument("--bound_lower",action="store_true",help="include lower bound on number of junctional reads a cell + gene needs to have in order to get a z score")

  parser.add_argument("--lower_bound",type=int,help="only include cell/gene pairs the have at least this many junctional reads for the gene")


  args = parser.parse_args()
  return args


def main():
  outpath = "scripts/output/rijk_zscore/"
  suff = ""
  args = get_args()
  sub_cols = ["cell","geneR1A_uniq","tissue","compartment","free_annotation","ontology","scZ","n.g_A","n.g_B","svd_z0","svd_z1","svd_z2","cell_gene"]
  df = pd.read_parquet("{}{}_sym_SVD_normdonor_S_{}_z_{}_b_{}{}.pq".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff),columns=sub_cols)
  df.drop_duplicates("cell_gene")[sub_cols].to_csv("{}{}_sym_SVD_normdonor_S_{}_z_{}_b_{}{}_subcol.tsv".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff),index=False,sep="\t")
main()

