import argparse
import datetime 
import numpy as np
import os
import pandas as pd
from scipy import linalg
import time
from tqdm import tqdm
import warnings

def get_args():
  parser = argparse.ArgumentParser(description="calculate splicing scores per gene/cell")
  parser.add_argument("--dataname",help="name of dataset to use")
  parser.add_argument("--pinning_S",type=float,help="pinning level for S_ijks")
  parser.add_argument("--pinning_z",type=float,help="pinning level for zs")
  parser.add_argument("--light",action="store_true",help="if included, don't calculate extra columns (saves time)")
  parser.add_argument("--verbose",action="store_true",help="print times")
  parser.add_argument("--unfilt",action="store_true",help="don't filter by SICILIAN")
  parser.add_argument("--v2",action="store_true",help="filter by SICILIAN v2")
  parser.add_argument("--temp",action="store_true",help="overwrite temp file")
  parser.add_argument("--svd_type",choices=["normgene","normdonor"],help="method of calculating matrix before SVD")
  parser.add_argument("--test",help="just run for MYL6 for testing purposes")





#  parser.add_argument("--bound_lower",action="store_true",help="include lower bound on number of junctional reads a cell + gene needs to have in order to get a z score")

  parser.add_argument("--lower_bound",type=int,help="only include cell/gene pairs the have at least this many junctional reads for the gene")


  args = parser.parse_args()
  return args

def main():
  args = get_args()

  outpath = "scripts/output/rijk_zscore/"

  if not os.path.exists("{}SVD_{}/".format(outpath,args.svd_type)):
    os.makedirs('{}SVD_{}/'.format(outpath,args.svd_type))

  suff = ""
  if args.light:
    suff += "_light"
  if args.unfilt:
    suff += "_unfilt"


  usecols = ["juncPosR1A","juncPosR1B","cell_gene","numReads","nSijkA","nSijkB","refName_newR1","cell","geneR1A_uniq"] 
  if args.temp:
    df = pd.read_parquet("{}{}_sym_temp_S_{}_z_{}_b_{}{}.pq".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff))
  else:
    df = pd.read_parquet("{}{}_sym_S_{}_z_{}_b_{}{}.pq".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff))


  ##### PERFORM SVD ZSCORE CALCULATION #####

  if args.svd_type == "normgene":
    zcontrib_col = "zcontrib"
  elif args.svd_type == "normdonor":
    
    for let in ["A","B"]:
      
      # find number of reads per donor (or acceptor) per cell
      df["cell_gene_pos" + let] = df["cell_gene"] + df["juncPosR1" + let].astype(str)
      df["n.g_pos" + let] = df.groupby("cell_gene_pos" + let)["numReads"].transform("sum")
      
      # normalize on a donor/acceptor rather than a gene basis
      # TRY OUT NOT SQRT-ING denominator as normalization
      df["zcontrib_posnorm" + let] = df["numReads"] * df["nSijk" + let] / df["n.g_pos" + let]
    zcontrib_col = "zcontrib_posnorm"

  letters = ["A","B"]
  for let in letters:

    # replace NANs with zeros
    df["zcontrib{}_rep".format(let)] = df[zcontrib_col + let].fillna(0)

    # create label for each junction + donor/acceptor
    df["str_juncPosR1" + let] = df["juncPosR1" + let].astype(int).astype(str) + "_" + let
    df["cell_gene_pos" + let] = df["cell"] + df["geneR1A_uniq"] + df["juncPosR1" + let].astype(str)

    # get sum of zcontribs for the given cell and splice site
    df["summed_zcontrib" + let] = df.groupby("cell_gene_pos" + let)["zcontrib{}_rep".format(let)].transform('sum')

  k = 3 # number of components to include
  loads = {"f{}".format(i) : {} for i in range(k)}
  zs = {"svd_z{}".format(i) : {} for i in range(k)}
  
  for gene, gene_df in tqdm(df.groupby("geneR1A_uniq")):
    
    # get zcontrib matrix
    gene_mats = []
    for let in letters:
      gene_mat = gene_df.drop_duplicates("cell_gene_pos" + let).pivot_table(index="cell_gene",columns="str_juncPosR1{}".format(let),values="summed_zcontrib" + let,fill_value=0)

      gene_mats.append(gene_mat)
    gene_mat = gene_mats[0].merge(gene_mats[1],on="cell_gene")

    # remove columns that are only zeros (could save space)
#    gene_mat = gene_mat.loc[:, (gene_mat != 0).any(axis=0)]

    # mean-normalize the rows
    gene_mat = gene_mat.subtract(gene_mat.mean(axis=1),axis=0)
    
#    print("gene mat shape",gene_mat.shape)
    # calculate svd
    u, s, vh = linalg.svd(gene_mat,check_finite=False,full_matrices=False)
    
    if len(s) >= k:
      # calculate new z scores based on svd
      new_zs = gene_mat.dot(np.transpose(vh[:k,:]))
  
      # calculate load on each component
      load = np.square(s)/sum(np.square(s))
  
      # save new zs and fs in dictionaries to save later
      for i in range(k):
        loads["f{}".format(i)][gene] = load[i]
        zs["svd_z{}".format(i)].update(pd.Series(new_zs[i].values,index=new_zs.index).to_dict())
  
      # save loadings
      v_out = pd.DataFrame(vh,columns=gene_mat.columns)
      v_out.to_csv("{}SVD_{}/{}_{}_S_{}_z_{}_b_{}{}.tsv".format(outpath, args.svd_type,gene,args.dataname, args.pinning_S, args.pinning_z, args.lower_bound, suff), index=False, sep = "\t")
      
  for i in range(k):
    df["f{}".format(i)] = df["geneR1A_uniq"].map(loads["f{}".format(i)])
    df["svd_z{}".format(i)] = df["cell_gene"].map(zs["svd_z{}".format(i)])
  df["svd_z_sumsq"] = (df[["svd_z{}".format(i) for i in range(k)]]**2).sum(axis=1)

  sub_cols = ["cell","geneR1A_uniq","scZ","svd_z_sumsq","n.g_A","n.g_B"] + ["f{}".format(i) for i in range(k)] + ["svd_z{}".format(i) for i in range(k)] #+ velocity_cols
  if "ontology" in df.columns:
    sub_cols = sub_cols + ["tissue","compartment","free_annotation","ontology"]
  if args.temp:
    df.drop_duplicates("cell_gene")[sub_cols].to_csv("{}{}_sym_S_{}_z_{}_b_{}{}_subcol.tsv".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff),index=False,sep="\t")

  sub_cols = ["cell","geneR1A_uniq","tissue","compartment","free_annotation","ontology","scZ","n.g_A","n.g_B","svd_z0","svd_z1","svd_z2"] 
  print("{}{}_sym_SVD_S_{}_z_{}_b_{}{}_subcol.tsv".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff))
  df.drop_duplicates("cell_gene")[sub_cols].to_csv("{}{}_sym_SVD_normdonor_S_{}_z_{}_b_{}{}_subcol.tsv".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff),index=False,sep="\t")

  df.to_parquet("{}{}_sym_S_{}_z_{}_b_{}{}.pq".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff))

  if args.temp:
    df.to_parquet("{}{}_sym_S_{}_z_{}_b_{}{}.pq".format(outpath,args.dataname,args.pinning_S, args.pinning_z, args.lower_bound, suff))
  else:
    df.to_parquet("{}{}_sym_SVD_{}_S_{}_z_{}_b_{}{}.pq".format(outpath,args.dataname,args.svd_type,args.pinning_S, args.pinning_z, args.lower_bound, suff))

main()
