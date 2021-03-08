import argparse
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from tqdm import tqdm

def get_args():
  parser = argparse.ArgumentParser(description="calculate fdr based on permutations")
  parser.add_argument("--dataname",help="name of dataset to use")
  parser.add_argument("--suffix",help="suffix to save with")

  parser.add_argument("--num_perms",type=int,help="number of permutations to run for")
  parser.add_argument("--z_col",choices=["scZ","svd_z0","svd_z1","svd_z2"],help="z column to calculate p values for")

  parser.add_argument("--temp",action="store_true",help="use temp rijk file")

  args = parser.parse_args()
  return args

def rand_ont_bygene(gene, gene_ont_dict):

  # randomly choose ontology that the gene is expressed in from list of all (with duplicates)
  return random.choice(gene_ont_dict[gene])

def rand_rijk(df,tiss_dict,comp_dict,freeann_dict):

  shuffle_onts = {}

  # shuffle independently for each ncells quantile
  for quant, quant_df in df.groupby("n_cells_gene_ont_quant"):
      
      # shuffle independently for each gene
      for gene, gene_df in quant_df.groupby("geneR1A_uniq"):

        # get all cells+genes and distribution of ontologies
        cell_genes = list(gene_df["cell_gene"].unique())
        ontologies = list(gene_df.drop_duplicates("cell_gene")["ontology"])

        # if there's only one ontology for this gene in this range don't assign (will be dropped)
        if len(set(ontologies)) > 1:

          # randomize ontologies (occurs in place)
          random.shuffle(ontologies)

          # map back to cell + gene
          so = {cell_gene : ont for cell_gene, ont in zip(cell_genes,ontologies)}

          # update dictionary in place
          shuffle_onts.update(so)

  # map ontologies in dataframe
  df["ontology"] = df["cell_gene"].map(shuffle_onts)
  df = df[~df["ontology"].isna()]
  df["ontology_gene"] = df["ontology"] + df["geneR1A_uniq"]

  # return shuffled dataframe
  return df


def sig_genes(df,z_cols):
  df = df.drop_duplicates("cell_gene")
  df = df[~df["ontology_gene"].isna()]

  for z_col in z_cols:
    # calculate the median SpliZ for each ontology
    df["median_" + z_col] = df["ontology_gene"].map(df.groupby("ontology_gene")[z_col].median())

  # Find number of cells in this gene + ontology
  df["n_cells_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell"].nunique())

  # only keep ontologies with at least 20 cells
  df = df[df["n_cells_ont"] > 20]
  return df.drop_duplicates("ontology_gene")


def main():
  mean_normalize = False

  outpath = "scripts/output/perm_pvals/"

  # load data
  args = get_args()
  num_quants = 5
  z_cols = ["scZ","svd_z0","svd_z1","svd_z2"]
  usecols = ["ontology","geneR1A_uniq","cell","cell_gene","n_cells_gene_ont_quant","ontology_gene"] + z_cols
  columns = ["geneR1A_uniq","tissue","compartment","free_annotation","ontology","cell","n.g"] + z_cols
  inpath = "scripts/output/rijk_zscore/"

  # read in z score file
  if args.temp:
    df = pd.read_parquet("{}{}_sym_temp_S_0.1_z_0.0_b_5.pq".format(inpath,args.dataname),columns=columns)
  else:
    df = pd.read_parquet("{}{}_sym_SVD_normdonor_S_0.1_z_0.0_b_5.pq".format(inpath,args.dataname),columns=columns)
#    df = pd.read_parquet("{}{}_sym_S_0.1_z_0.0_b_5.pq".format(inpath,args.dataname),columns=columns)


  df["cell_gene"] = df["cell"] + df["geneR1A_uniq"]
#  z_cols = ["scZ","svd_z0","svd_z1","svd_z2"]

  if mean_normalize:
    for z_col in z_cols:
      df["mean_" + z_col + "_bygene"] = df["geneR1A_uniq"].map(df.drop_duplicates("cell_gene").groupby("geneR1A_uniq")[z_col].mean())
      df["n" + z_col] = df[z_col] - df["mean_" + z_col + "_bygene"]
    usecols = usecols + ["n" + z_col for z_col in z_cols]
    columns = columns + ["n" + z_col for z_col in z_cols]
  
    z_cols = z_cols + ["n" + z_col for z_col in z_cols]
  
  # find ontologies, tissues, compartments, cell types
  df["ontology"] = df["tissue"] + df["compartment"] + df["free_annotation"]
  df["ontology_gene"] = df["ontology"] + df["geneR1A_uniq"]
#  df["numReads_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["numReads"].sum())
  # remove those that lack annotation or gene name
  df = df[(df["free_annotation"] != "") & (~df["geneR1A_uniq"].isin(["","unknown"]))]
  tiss_dict = pd.Series(df.tissue.values,index=df.ontology).to_dict()
  comp_dict = pd.Series(df.compartment.values,index=df.ontology).to_dict()
  freeann_dict = pd.Series(df.free_annotation.values,index=df.ontology).to_dict()
#  df["numReads_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["numReads"].sum())
  
  # mark which quantile each gene is in (by cell)
  df["num_cells_gene_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell"].nunique())
  temp = df.drop_duplicates("ontology_gene")
  temp["n_cells_gene_ont_quant"] = pd.qcut(temp['num_cells_gene_ont'], num_quants, labels=False,duplicates="drop") 
  quant_dict = pd.Series(temp.n_cells_gene_ont_quant.values,index=temp.ontology_gene).to_dict()
  df["n_cells_gene_ont_quant"] = df["ontology_gene"].map(quant_dict)

  # save this so we don't have to read it in at each permutation
  orig_df = df.copy()

  # get median z scores for original data
  fdr = sig_genes(df,z_cols)
#  fdr["abs_mz"] = abs(fdr["median_" + args.z_col])

  # perform permutations
  dfs = []
  for num in tqdm(range(args.num_perms)):
  
    print("round",num)

    # randomize ontologies
    df = rand_rijk(orig_df.copy()[usecols],tiss_dict,comp_dict,freeann_dict)
    print("rank_rijk")

    # find median z scores
    df = sig_genes(df,z_cols)
    print("sig")

    # add relevant info to list
    dfs.append(df[["n_cells_gene_ont_quant","geneR1A_uniq","ontology"] + ["median_" + x for x in z_cols]])


  # form concatenated dataframe with all permutations and save it
  allperm_df = pd.concat(dfs,axis=0)
#  allperm_df["abs_mz"] = abs(allperm_df["median_" + args.z_col])
  allperm_df.to_csv("{}{}_{}_allperm{}.tsv".format(outpath,args.dataname,args.num_perms,args.suffix),sep="\t",index=False)

  # separate mz distributions per quantile
  quant_dfs = [allperm_df[allperm_df["n_cells_gene_ont_quant"] == x] for x in range(num_quants)]
  
  for z_col in z_cols:
    pvals = []
  
#    median_val = orig_df.drop_duplicates("cell_gene")[z_col].median()
  
    # try out method from simulation
    # for each row in real data, subset to cell quantile and find frac perm with >= mz val
    for index, row in tqdm(fdr.iterrows()):
      temp = quant_dfs[int(row["n_cells_gene_ont_quant"])]
#      temp.reset_index(inplace=True) 
      # try subsetting to only given gene
#      temp = temp[temp["geneR1A_uniq"] == row["geneR1A_uniq"]]
      if temp.shape[0] > 0:
#        delta = abs(median_val - row["median_" + z_col])
    #    frac = temp[temp["abs_mz"] >= row["abs_mz"]].shape[0]/(temp.shape[0])
        try:
#          frac = temp[(temp["median_" + z_col] <= (median_val - delta)) | (temp["median_" + z_col] >= (median_val + delta))].shape[0]/(temp.shape[0])
          frac = temp[temp["median_" + z_col] < row["median_" + z_col]].shape[0]/(temp.shape[0])

          
        except Exception as e:
          print("ERROR")
          print(row["geneR1A_uniq"])
          print(temp.shape)
          print(row["ontology"])
          print(temp.head())
          temp.to_csv("temp.tsv",sep="\t")
          print("median", median_val)
          print("delta",delta)
          print(z_col)
  
          print()
      else:
        frac = np.nan
      pvals.append(frac)
  
  
    # should they be adjusted separately per quantile?
  #  pvals_adj = multipletests(pvals,0.05, method="fdr_bh")[1]
  
    fdr["cdf_quant_pval_" + z_col] = pvals
    fdr["inv_cdf_quant_pval_" + z_col] = 1 - fdr["cdf_quant_pval_" + z_col] 
  
    # save unadjusted p values
    fdr["quant_pval_" + z_col] = 2*fdr[["inv_cdf_quant_pval_" + z_col,"cdf_quant_pval_" + z_col]].min(axis=1)
  
    # save adjusted p values
    fdr.loc[~fdr["quant_pval_" + z_col].isna(),"quant_pval_adj_" + z_col] = multipletests(fdr.loc[~fdr["quant_pval_" + z_col].isna(),"quant_pval_" + z_col],0.05, method="fdr_bh")[1]
  
  fdr.to_csv("{}{}_fdr_{}{}.tsv".format(outpath,args.dataname,args.num_perms,args.suffix),sep="\t",index=False)
  
main()
