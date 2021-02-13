import argparse
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import pandas as pd
import random
from tqdm import tqdm

def get_args():
  parser = argparse.ArgumentParser(description="calculate fdr based on permutations")
  parser.add_argument("--dataname",help="name of dataset to use")
  parser.add_argument("--suffix",help="suffix to save with")

  parser.add_argument("--num_perms",type=int,help="number of permutations to run for")
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


def sig_genes(df):
  z_col = "scZ"
  df = df.drop_duplicates("cell_gene")
  df = df[~df["ontology_gene"].isna()]

  # calculate the median SpliZ for each ontology
  df["mz"] = df["ontology_gene"].map(df.groupby("ontology_gene")[z_col].median())

  # Find number of cells in this gene + ontology
  df["n_cells_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell"].nunique())

  # only keep ontologies with at least 20 cells
  df = df[df["n_cells_ont"] > 20]
  return df.drop_duplicates("ontology_gene")


def main():

  outpath = "scripts/output/perm_pvals/"

  # load data
  args = get_args()
  num_quants = 5
  usecols = ["ontology","geneR1A_uniq","scZ","cell","cell_gene","n_cells_gene_ont_quant","ontology_gene"]
  columns = ["geneR1A_uniq","tissue","compartment","free_annotation","ontology","scZ","cell","n.g"]
  inpath = "scripts/output/rijk_zscore/"

  # read in z score file
  if args.temp:
    df = pd.read_parquet("{}{}_sym_temp_S_0.1_z_0.0_b_5.pq".format(inpath,args.dataname),columns=columns)
  else:
    df = pd.read_parquet("{}{}_sym_S_0.1_z_0.0_b_5.pq".format(inpath,args.dataname),columns=columns)

  df["cell_gene"] = df["cell"] + df["geneR1A_uniq"]
  
  # find ontologies, tissues, compartments, cell types
  df["ontology"] = df["tissue"] + df["compartment"] + df["free_annotation"]
  df["ontology_gene"] = df["ontology"] + df["geneR1A_uniq"]
  df["numReads_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["n.g"].sum())
  # remove those that lack annotation or gene name
  df = df[(df["free_annotation"] != "") & (~df["geneR1A_uniq"].isin(["","unknown"]))]
  tiss_dict = pd.Series(df.tissue.values,index=df.ontology).to_dict()
  comp_dict = pd.Series(df.compartment.values,index=df.ontology).to_dict()
  freeann_dict = pd.Series(df.free_annotation.values,index=df.ontology).to_dict()
  df["numReads_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["n.g"].sum())
  
  # mark which quantile each gene is in (by cell)
  df["num_cells_gene_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell"].nunique())
  temp = df.drop_duplicates("ontology_gene")
  temp["n_cells_gene_ont_quant"] = pd.qcut(temp['num_cells_gene_ont'], num_quants, labels=False,duplicates="drop") 
  quant_dict = pd.Series(temp.n_cells_gene_ont_quant.values,index=temp.ontology_gene).to_dict()
  df["n_cells_gene_ont_quant"] = df["ontology_gene"].map(quant_dict)

  # save this so we don't have to read it in at each permutation
  orig_df = df.copy()

  # get median z scores for original data
  fdr = sig_genes(df)
  fdr["abs_mz"] = abs(fdr["mz"])

  # perform permutations
  dfs = []
  for num in tqdm(range(args.num_perms)):
  
    print("round",num)

    # randomize ontologies
    df = rand_rijk(orig_df.copy()[usecols],tiss_dict,comp_dict,freeann_dict)
    print("rank_rijk")

    # find median z scores
    df = sig_genes(df)
    print("sig")

    # add relevant info to list
    dfs.append(df[["mz","n_cells_gene_ont_quant","geneR1A_uniq"]])


  # form concatenated dataframe with all permutations and save it
  allperm_df = pd.concat(dfs,axis=0)
  allperm_df["abs_mz"] = abs(allperm_df["mz"])
  allperm_df.to_csv("{}{}_{}_allperm{}.tsv".format(outpath,args.dataname,args.num_perms,args.suffix),sep="\t",index=False)

  # separate mz distributions per quantile
  quant_dfs = [allperm_df[allperm_df["n_cells_gene_ont_quant"] == x] for x in range(num_quants)]
  
  pvals = []

  # for each row in real data, subset to cell quantile and find frac perm with >= mz val
  for index, row in tqdm(fdr.iterrows()):
    temp = quant_dfs[int(row["n_cells_gene_ont_quant"])]
    frac = temp[temp["abs_mz"] >= row["abs_mz"]].shape[0]/(temp.shape[0])
    pvals.append(frac)


  # should they be adjusted separately per quantile?
  pvals_adj = multipletests(pvals,0.05, method="fdr_bh")[1]
  fdr["quant_pval"] = pvals_adj

  
  fdr.to_csv("{}{}_fdr_{}{}.tsv".format(outpath,args.dataname,args.num_perms,args.suffix),sep="\t",index=False)
  
main()
