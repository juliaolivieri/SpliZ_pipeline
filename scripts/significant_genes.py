import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow
from scipy.stats import norm
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

def get_args():
  parser = argparse.ArgumentParser(description="calculate splicing scores per gene/cell")
  parser.add_argument("--dataname",help="name of dataset to use")
  parser.add_argument("--z_col", choices=["z","scaled_z","scZ"], help="column to use as the z score")
  parser.add_argument("--pinning_S",type=float,help="level of pinning for S_ijks")
  parser.add_argument("--pinning_z",type=float,help="level of pinning for zs")
  parser.add_argument("--bound_lower",action="store_true",help="use values based on lower bound cutoff")
  parser.add_argument("--unfilt",action="store_true",help="don't filter by SICILIAN")
  parser.add_argument("--lower_bound",type=int,help="only include cell/gene pairs the have at least this many junctional reads for the gene")
  parser.add_argument("--bootstrap",action="store_true",help="get a bootstrap q value")


  args = parser.parse_args()
  return args

def main():
  args = get_args()
  print("z_col",args.z_col)
  outpath = "scripts/output/significant_genes/"
  suff = ""
  if args.unfilt:
    suff += "_unfilt"
  df = pd.read_parquet("scripts/output/rijk_zscore/{}_sym_S_{}_z_{}_b_{}{}.pq".format(args.dataname,args.pinning_S, args.pinning_z,args.lower_bound, suff),columns=["free_annotation","tissue","compartment","geneR1A_uniq","z","cell","numReads", "cell_gene","scaled_z","scZ","n.g"])

#  df["n.g"] = df["cell_gene"].map(df.groupby("cell_gene")["numReads"].sum())
#  df["scaled_z"] = df["z"] / np.sqrt(df["n.g"])

  # remove multiple z scores from the same cell + gene
#  df["cell_gene"] = df["cell"] + df["geneR1A_uniq"]
  df = df.drop_duplicates("cell_gene")

  # remove those that lack annotation or gene name
  df = df[(df["free_annotation"] != "") & (~df["geneR1A_uniq"].isin(["","unknown"]))]

  # get median by ontology
  df["ontology"] = df["tissue"] + df["compartment"] + df["free_annotation"]
  df["ontology_gene"] = df["ontology"] + df["geneR1A_uniq"]
  df = df[~df["ontology_gene"].isna()]
  df["mz"] = df["ontology_gene"].map(df.groupby("ontology_gene")[args.z_col].median())
  df["meanz"] = df["ontology_gene"].map(df.groupby("ontology_gene")[args.z_col].mean())

  # Find number of cells in this gene + ontology
  df["n_cells_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell"].nunique())

  if args.bootstrap:
    # run bootstrap
    num_bootstrap = 100
    boot_dict = {}
    for gene, gene_df in tqdm(df.groupby("geneR1A_uniq")):
    #   gene_df = df[df["geneR1A_uniq"] == "MYL6"]
      
    
      for ontology in gene_df["ontology"].unique():
        mzs = []
        num_cells = gene_df[gene_df["ontology"] == ontology].shape[0]
        for i in range(num_bootstrap):
          mzs.append(np.median(np.random.choice(gene_df[args.z_col],num_cells)))
        boot_dict[ontology + gene] = len([x for x in mzs if x < gene_df[gene_df["ontology"] == ontology]["mz"].iloc[0]])/num_bootstrap
    df["boot_q"] = df["ontology_gene"].map(boot_dict)  

  # get standard deviation of z "globally"
  df["allsig"] = df[args.z_col].std()

  df["sig_bygene"] = df["geneR1A_uniq"].map(df.groupby("geneR1A_uniq")[args.z_col].std())

  df["n.g_med"]  = df["ontology"].map(df.groupby("ontology")["n.g"].median())

  # calculate the median "globally"
  df["mz_all"] = df[args.z_col].median()
  df["mz_bygene"] = df["geneR1A_uniq"].map(df.groupby("geneR1A_uniq")[args.z_col].median())

  df["meanz_all"] = df[args.z_col].mean()

  # calculate standard deviation of the median (based only on number of cells with the gene in the ontology)
#  df["allsd_med"] = np.sqrt(np.pi*df["allsig"]**2/(2*df["n_cells_ont"]))
#  df["sd_bygene_med"] = np.sqrt(np.pi*df["sig_bygene"]**2/(2*df["n_cells_ont"]))

  # assume all variances are 1
  df["allsd_med"] = np.sqrt(np.pi/(2*df["n_cells_ont"]))
  df["sd_bygene_med"] = np.sqrt(np.pi/(2*df["n_cells_ont"]))

  # use variance == 1/n.g**2
  df["sd_n.g"] = np.sqrt((np.pi*(1/df["n.g_med"]**2))/(2*df["n_cells_ont"]))


  df["allsd_mean"] = np.sqrt(df["allsig"]**2/df["n_cells_ont"])

  # calculate the cdf of the median
#  df["allp"] = norm.cdf(df["mz"],loc=df["mz_all"], scale=df["allsd_med"])
  df["allp"] = norm.cdf(df["mz"],loc=df["mz_all"], scale=df["allsd_med"])
  df["allp_n.g"] = norm.cdf(df["mz"],loc=df["mz_all"], scale=df["sd_n.g"])


  df["p"] = norm.cdf(df["mz"],loc=df["mz_bygene"], scale=df["sd_bygene_med"])

  df["allp_mean"] = norm.cdf(df["meanz"],loc=df["meanz_all"], scale=df["allsd_mean"])

  # get one value per gene/ontology pair
  m_df = df.drop_duplicates("ontology_gene")

  ep = 0.025

  m_df["diff"] = False
  m_df["diff_mean"] = False
  m_df["diff_p"] = False
  m_df.loc[(m_df["p"] < ep) | (m_df["p"] > 1 - ep), "diff_p"] = True
  m_df = m_df[m_df["n_cells_ont"] > 20]
  m_df.loc[(m_df["allp"] < ep) | (m_df["allp"] > 1 - ep), "diff"] = True

  m_df.loc[(m_df["allp_mean"] < ep) | (m_df["allp_mean"] > 1 - ep), "diff_mean"] = True

  # add condensed freeannotation names
#  ann_df = pd.read_csv("notebooks/output/condense_freeanns/TS_condense_freeanns.tsv",sep="\t")
#  ann_dict = pd.Series(ann_df.condensed_freeann.values,index=ann_df.free_annotation).to_dict()
#  m_df["condensed_freeann"] = m_df["free_annotation"].map(ann_dict)
#  m_df["is_stem"] = m_df[m_df["free_annotation"].str.contains("stem")]

  m_df.to_csv("{}{}-{}_allp_S_{}_z_{}_b_{}{}.tsv".format(outpath,args.dataname,args.z_col,args.pinning_S, args.pinning_z,args.lower_bound, suff),sep="\t",index=False)
main()
