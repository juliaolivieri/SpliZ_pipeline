import argparse
from collections import defaultdict
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
import numpy as np
import pandas as pd
import random
import seaborn as sns
from scipy import linalg
from scipy import stats
from statsmodels.stats.proportion import proportion_confint
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({'font.size': 14})

outpath = "/scratch/PI/horence/JuliaO/single_cell/SZS_pipeline3/scripts/output/power_simulation/"

def get_args():
  parser = argparse.ArgumentParser(description="power analysis simulation for spliz")
  parser.add_argument("--num_perms",type=int,help="number of permutations to perform for significance test")
  parser.add_argument("--num_trials",type=int,help="number of trials to perform at each read depth")
  parser.add_argument("--num_cells",type=int,help="number of cells in each cell type")

  parser.add_argument("--max_depth",type=int,help="max depth to sample range(1,max_depth + 1)")
  parser.add_argument("--regime",help="Which regime to simulate from")

  parser.add_argument("--psi",type=float,help="percent spliced in to simulate")


  args = parser.parse_args()
  return args

def calc_pval(var_df):

  # calculate the inner sum that's subtracted
  num = 0
  denom = 0
  for index, row in var_df.iterrows():
    num += row["num_cells_ont"]*row["ont_median"]/row["ont_var"]
    denom += row["num_cells_ont"]/row["ont_var"]
  const = num/denom

  # calculate the outer sum
  sum_vals = 0
  for index, row in var_df.iterrows():
    sum_vals += (row["num_cells_ont"]/row["ont_var"])*(row["ont_median"] - const)**2

  # return the chi^2 p value and the chi^2 statistic
  return 1 - stats.chi2.cdf(sum_vals , var_df.shape[0] - 1), sum_vals

def get_var_df(sub_df, z_col):

  sub_df["num_cells_ont"] = sub_df["cell_type"].map(sub_df.groupby("cell_type")["cell"].nunique())
  sub_df["ont_median"] = sub_df["cell_type"].map(sub_df.groupby("cell_type")[z_col].median())
  sub_df["ont_var"] = sub_df["cell_type"].map(sub_df.groupby("cell_type")[z_col].var())

  var_df = sub_df.drop_duplicates("cell_type")[["cell_type","ont_median","num_cells_ont","ont_var"]]
#   display(var_df)

  # helps when the eigenvector just has a 1 in one place, so all values are the same
  var_df = var_df[var_df["ont_var"] > 0]
  return var_df

def variance_adj_pval(df, z_cols, num_perms = 10, alpha = 0.05):
  var_adj = 0.1
  
  adj_var = False
  df = df.drop_duplicates("cell_gene")
#   z_cols = ["scZ","svd_z0","svd_z1","svd_z2"]
  out = {"pval" : [], "geneR1A_uniq" : [], "num_onts" : [],"z_col" : [],"max_abs_median" : []}

  perm_pval = True
  if perm_pval:
    out["perm_pval"] = []

  for gene, sub_df in df.groupby("geneR1A_uniq"):

    for z_col in z_cols:

      var_df = get_var_df(sub_df, z_col)
#       print(z_col)
#       display(var_df)

      if var_df.shape[0] > 1:
        if adj_var:
          var_df["ont_var"] = var_df["ont_var"] + var_adj
        pval, Tn1 = calc_pval(var_df)
        out["pval"].append(pval)
        out["geneR1A_uniq"].append(gene)
        out["num_onts"].append(var_df.shape[0])
        out["z_col"].append(z_col)
        out["max_abs_median"].append((var_df["ont_median"].abs()).max())

        if perm_pval:
          sub_df_perm = sub_df.copy()
          if (pval < 1):
            Tn1_dist = []
  #            for i in range(args.num_perms):
            while len(Tn1_dist) < num_perms:
              sub_df_perm["cell_type"] = np.random.permutation(sub_df_perm["cell_type"])
              var_df = get_var_df(sub_df_perm, z_col)
              if var_df.shape[0] > 1:
                if adj_var:
                  var_df["ont_var"] = var_df["ont_var"] + var_adj
                pval, Tn1_perm = calc_pval(var_df)
                Tn1_dist.append(Tn1_perm)
            out["perm_pval"].append(len([x for x in Tn1_dist if x < Tn1])/num_perms)
          else:
            out["perm_pval"].append(np.nan)
        else:
          out["perm_pval"].append(np.nan)
  out_df = pd.DataFrame.from_dict(out)

  out_df["perm_pval_inv"] = 1 - out_df["perm_pval"]
  out_df["perm_pval2"] = 2*out_df[["perm_pval","perm_pval_inv"]].min(axis=1)

  # only one gene in these situations, so no need for adjustment
  out_df["pval_adj"] = out_df["pval"]
  out_df["perm_pval2_adj"] = out_df["perm_pval2"]
#  out_df["pval_adj"] = multipletests(out_df["pval"],alpha, method="fdr_bh")[1]
#  out_df.loc[~out_df["perm_pval2"].isna(),"perm_pval2_adj"] = multipletests(out_df.loc[~out_df["perm_pval2"].isna(),"perm_pval2"], alpha, method = "fdr_bh")[1]

  # reformat output
  new_out = {"geneR1A_uniq" : [], "num_onts" : []}
  for z_col in z_cols:
    new_out["chi2_pval_adj_" + z_col] = []
    new_out["perm_pval_adj_" + z_col] = []
    new_out["max_abs_median_" + z_col] = []
    new_out["perm_cdf_" + z_col] = []
  for gene, gene_df in out_df.groupby("geneR1A_uniq"):
    new_out["geneR1A_uniq"].append(gene)
    new_out["num_onts"].append(gene_df["num_onts"].iloc[0])
    temp_z_cols = []
    for z_col, z_df in gene_df.groupby("z_col"):
      new_out["chi2_pval_adj_" + z_col].append(z_df["pval_adj"].iloc[0])
      new_out["perm_pval_adj_" + z_col].append(z_df["perm_pval2_adj"].iloc[0])
      new_out["max_abs_median_" + z_col].append(z_df["max_abs_median"].iloc[0])
      new_out["perm_cdf_" + z_col].append(z_df["perm_pval"].iloc[0])
      temp_z_cols.append(z_col)
    for z_col in [x for x in z_cols if x not in temp_z_cols]:
      new_out["chi2_pval_adj_" + z_col].append(np.nan)
      new_out["perm_pval_adj_" + z_col].append(np.nan)
      new_out["max_abs_median_" + z_col].append(np.nan)
      new_out["perm_cdf_" + z_col].append(np.nan)
  new_out_df = pd.DataFrame.from_dict(new_out).sort_values("perm_pval_adj_scZ")
  return new_out_df

def calc_PSI(df, inc_juncs, exc_juncs):
  
  # subset to only exon skipped junctions
  temp = df[df["refName_newR1"].isin(inc_juncs + exc_juncs)]
  cells_with_event = set(temp["cell"].unique()) 

  # find number of reads assigned to these junctions
  temp["cell_sum"] = temp.groupby("cell")["numReads"].transform("sum")

  #subset to only exc_junc (this will remove cells with no exc_junc)
  temp = temp[temp["refName_newR1"].isin(exc_juncs)]
  
  # count number of exc reads in each cell
  temp["exc_cell_sum"] = temp.groupby("cell")["numReads"].transform("sum")
  temp.drop_duplicates("cell")
  
  # calculate psi
  temp["psi"] = 1 - temp["exc_cell_sum"]/temp["cell_sum"] 
  df["psi"] = df["cell"].map(pd.Series(temp["psi"].values,index=temp["cell"]).to_dict())
  
  # for cells with no exon skipped junctions psi is 1
  df["psi"] = df["psi"].fillna(1)

  # if the event isn't present, don't calculate PSI for it
  df.loc[~(df["cell"].isin(cells_with_event)), "psi"] = np.nan
  
  # NOTE: error if there are no reads for this exon skipping event in the cell (will return 1 as psi)
  return df

def perm_pval(df, num_perms = 1000, num_quants = 5, plot = False):


  # z_col = "scZ"
  z_cols = ["psi","scZ","svd_z0","svd_z1","svd_z2"]
  df["n.g"] = df.groupby("cell_gene")["numReads"].transform("sum")
  # find ontologies, tissues, compartments, cell types
  # df["cell_type"] = df["tissue"] + df["compartment"] + df["free_annotation"]
  df["cell_type_gene"] = df["cell_type"] + df["geneR1A_uniq"]
  df["numReads_ont"] = df["cell_type_gene"].map(df.groupby("cell_type_gene")["numReads"].sum())

  # df["numReads_ont"] = df["cell_type_gene"].map(df.groupby("cell_type_gene")["n.g"].sum())

  # mark which quantile each gene is in (by cell)
  df["num_cells_gene_ont"] = df["cell_type_gene"].map(df.groupby("cell_type_gene")["cell"].nunique())
  temp = df.drop_duplicates("cell_type_gene")

  temp["n_cells_gene_ont_quant"] = pd.qcut(temp['num_cells_gene_ont'], num_quants, labels=False,duplicates="drop")
  quant_dict = pd.Series(temp.n_cells_gene_ont_quant.values,index=temp.cell_type_gene).to_dict()
  df["n_cells_gene_ont_quant"] = df["cell_type_gene"].map(quant_dict)

  # save this so we don't have to read it in at each permutation
  orig_df = df.copy()

  # get median z scores for original data
  fdr = sig_genes(df)
#   fdr["abs_mz"] = abs(fdr["mz"])

  usecols = ["cell_type","geneR1A_uniq","cell","cell_gene","n_cells_gene_ont_quant","cell_type_gene"] + z_cols

  # perform permutations
  dfs = []
  for num in range(num_perms):

  #   print("round",num)

    # randomize ontologies
    df = rand_rijk(orig_df.copy()[usecols])

    # find median z scores
    df = sig_genes(df)


    # add relevant info to list
    dfs.append(df[["n_cells_gene_ont_quant","geneR1A_uniq"] + ["median_" + x for x in z_cols]])

  # form concatenated dataframe with all permutations and save it
  allperm_df = pd.concat(dfs,axis=0)

  if allperm_df["n_cells_gene_ont_quant"].nunique() > 1:
    # separate mz distributions per quantile
    quant_dfs = [allperm_df[allperm_df["n_cells_gene_ont_quant"] == x] for x in range(num_quants)]

  for z_col in z_cols:
    pvals = []
    median_val = orig_df.drop_duplicates("cell_gene")[z_col].median()
    if plot:
      plt.hist(allperm_df["median_" + z_col],100)
      plt.axvline(x=median_val,label="median")
    for index, row in fdr.iterrows():
    #   print(index)
      if allperm_df["n_cells_gene_ont_quant"].nunique() > 1:
        temp = quant_dfs[int(row["n_cells_gene_ont_quant"])]
      else:
        temp = allperm_df
      if temp.shape[0] > 0:

        # psi isn't centered at 0; need to re-center; do for other scores as well to make consistent
    #       if z_col in ["psi","center_psi"]:
        delta = abs(median_val - row["median_" + z_col])
        if plot:
          plt.axvline(x=median_val - delta, label=row["cell_type"],linestyle="--")
          plt.axvline(x=median_val + delta, label=row["cell_type"], linestyle="--")

        frac = temp[(temp["median_" + z_col] <= (median_val - delta)) | (temp["median_" + z_col] >= (median_val + delta))].shape[0]/(temp.shape[0])

      else:
        print("temp.shape[0] == 0")
        print("ERROR")
        display(allperm_df["n_cells_gene_ont_quant"].value_counts(dropna=False))

        display(allperm_df)
        frac = np.nan
    #   print(frac)
      pvals.append(frac)
    if plot:
      plt.legend()
      plt.title(z_col)
      plt.show()
    # should they be adjusted separately per quantile?
    pvals_adj = multipletests(pvals,0.05, method="fdr_bh")[1]
    fdr["quant_pval_" + z_col] = pvals_adj
  return fdr

def sig_genes(df):
  z_cols = ["psi","scZ","svd_z0", "svd_z1", "svd_z2"]
  df = df.drop_duplicates("cell_gene")
  df = df[~df["cell_type_gene"].isna()]

  for z in z_cols:
    df["median_" + z] = df["cell_type_gene"].map(df.groupby("cell_type_gene")[z].median())

  # calculate the median SpliZ for each cell_type
#   df["mz"] = df["cell_type_gene"].map(df.groupby("cell_type_gene")[z_col].median())

  # Find number of cells in this gene + cell_type
  df["n_cells_ont"] = df["cell_type_gene"].map(df.groupby("cell_type_gene")["cell"].nunique())

  # only keep ontologies with at least 20 cells
#   df = df[df["n_cells_ont"] > 20]
  return df.drop_duplicates("cell_type_gene")

def rand_rijk(df):

  shuffle_onts = {}

  if df["n_cells_gene_ont_quant"].nunique() > 1:
    # shuffle independently for each ncells quantile
    for quant, quant_df in df.groupby("n_cells_gene_ont_quant"):
#         print("diff quants")
        # shuffle independently for each gene
        for gene, gene_df in quant_df.groupby("geneR1A_uniq"):

          # get all cells+genes and distribution of ontologies
          cell_genes = list(gene_df["cell_gene"].unique())
          ontologies = list(gene_df.drop_duplicates("cell_gene")["cell_type"])

          # if there's only one ontology for this gene in this range don't assign (will be dropped)
          if len(set(ontologies)) > 1:

            # randomize ontologies (occurs in place)
            random.shuffle(ontologies)

            # map back to cell + gene
            so = {cell_gene : ont for cell_gene, ont in zip(cell_genes,ontologies)}

            # update dictionary in place
            shuffle_onts.update(so)

  else:
#       print("one quant")
      quant_df = df
      # shuffle independently for each gene
      for gene, gene_df in quant_df.groupby("geneR1A_uniq"):

        # get all cells+genes and distribution of ontologies
        cell_genes = list(gene_df["cell_gene"].unique())
        ontologies = list(gene_df.drop_duplicates("cell_gene")["cell_type"])

        # if there's only one ontology for this gene in this range don't assign (will be dropped)
        if len(set(ontologies)) > 1:

          # randomize ontologies (occurs in place)
          random.shuffle(ontologies)

          # map back to cell + gene
          so = {cell_gene : ont for cell_gene, ont in zip(cell_genes,ontologies)}

          # update dictionary in place
          shuffle_onts.update(so)
  # map ontologies in dataframe
  df["cell_type"] = df["cell_gene"].map(shuffle_onts)
  df = df[~df["cell_type"].isna()]
  df["cell_type_gene"] = df["cell_type"] + df["geneR1A_uniq"]

  # return shuffled dataframe
  return df

def prepare_df(df, let, rank_by_donor, rev_let, let_dict):

  # create donor identifier
  df["pos{}_group".format(let)] = df["juncPosR1{}".format(let)].astype(str) + df["geneR1A_uniq"]
  df["rank_" + let_dict[let]] = df.groupby("pos{}_group".format(let))["juncPosR1{}".format(rev_let[let])].rank(method="dense")

  # remove consitutive splicing
  df["max_rank"] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["rank_" + let_dict[let]].max())
  df = df[df["max_rank"] > 1]

  if not rank_by_donor:
    df["rank_" + let_dict[let]] = df.groupby("geneR1A_uniq")["juncPosR1B"].rank(method="dense")
  return df

def calc_Sijk(df,let, pinning_S, let_dict):
  # calculate the average rank calculation per gene
  # same as this calculation (this one's slower): df["rank_mean"] = df.groupby("pos{}_group".format(let)).apply(lambda x: (x["numReads"] * x["rank_acc"])/x["numReads"].sum()).reset_index(level=0,drop=True)

  # number of reads with this donor across all cells
  df["sum_reads_group"] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["numReads"].sum())

  df["read_x_" + let_dict[let]] = df["numReads"] * df["rank_" + let_dict[let]]

  # the sum of acceptors for all reads in all cells with this donor
  df["num"] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["read_x_" + let_dict[let]].sum())

  # average acceptor for a read with this donor (donor has one value for this)
  df["rank_mean"]= df["num"] / df["sum_reads_group"]

#   df["rank_mean"] = df.groupby("pos{}_group".format(let)).apply(lambda x: (x["numReads"] * x["rank_acc"])/x["numReads"].sum())
  # sum squared difference in rank for ever read
  df["sq_diff"] = df["numReads"] * (df["rank_" + let_dict[let]] - df["rank_mean"])**2

  # Get the sum of these squared differences for each donor
  df["don_num"] = df["pos{}_group".format(let)].map(df.groupby("pos{}_group".format(let))["sq_diff"].sum())

  # sum of squared differences normalized by total number of reads
  # changed to make it the sample standard deviation (added minus 1)
  df["don_sigma"] = df["don_num"] / (df["sum_reads_group"])

#   df["don_sigma"] = df.groupby("pos_{}_group".format(let)).apply(lambda x: (x["numReads"] * (x["rank_acc"] - x["rank_mean"])**2).sum()/x["numReads"].sum())
  # this is the S_ijk value (difference normalized by sd) - should be normal 0/1
  df["S_ijk_{}".format(let)] = (df["rank_" + let_dict[let]] - df["rank_mean"])/np.sqrt(df["don_sigma"])

  # round outlying S values
  low_quant = df["S_ijk_{}".format(let)].quantile(q=pinning_S)
  high_quant = df["S_ijk_{}".format(let)].quantile(q=1 - pinning_S)
  df["S_ijk_{}_unpinned".format(let)] = df["S_ijk_{}".format(let)]

  df.loc[df["S_ijk_{}".format(let)] < low_quant,"S_ijk_{}".format(let)] = low_quant
  df.loc[df["S_ijk_{}".format(let)] > high_quant,"S_ijk_{}".format(let)] = high_quant



  # correct for those with no variance
  df.loc[df["don_sigma"] == 0, "S_ijk_{}".format(let)] = 0
  df["n_sijk"] = df["numReads"]
  df.loc[df["don_sigma"] == 0,"n_sijk"] = 0
  return df

def normalize_Sijks(df,let):
  # calculate mean of SijkA's per gene
  df["n_s"] = df["numReads"] * df["S_ijk_" + let]
  df["num"] = df["geneR1A_uniq"].map(df.groupby("geneR1A_uniq")["n_s"].sum())
  df["n_gene"] = df["geneR1A_uniq"].map(df.groupby("geneR1A_uniq")["numReads"].sum())
  df["sijk{}_mean".format(let)] = df["num"] / df["n_gene"]

  # calculate standard deviation of SijkA's per gene
  # temp_df["n_g"] = temp_df["cell_gene"].map(temp_df.groupby("cell_gene")["n_sijk"].sum())
  # temp_df["sijkA_mean"] = temp_df["num"] / temp_df["n_gene"]
  # temp_df["S_ijk_A_mean"] = temp_df["geneR1A_uniq"].map(temp_df.groupby("geneR1A_uniq")["S_ijk_A"].mean())
  df["sd_num"] = df["numReads"] * (df["S_ijk_" + let] - df["sijk{}_mean".format(let)])**2
  df["num"] = df["geneR1A_uniq"].map(df.groupby("geneR1A_uniq")["sd_num"].sum())
  df["sijk{}_var".format(let)] = df["num"] / df["n_gene"]
#  df["S_ijk_{}".format(let)] = (df["S_ijk_{}".format(let)] - df["sijk{}_mean".format(let)])/np.sqrt(df["sijk{}_var".format(let))
  return df

def plot_SVD(gene_mat, vh, s, df, name, psi=True):
  sns.heatmap(gene_mat)
  plt.title(name)
  plt.savefig("{}{}_gene_mat.png".format(outpath, name),bbox_inches="tight")
  plt.close()
  
  sns.heatmap(np.transpose(vh[:,:]),annot=True,fmt=".02f",xticklabels=["s{}=\n{:.03f}".format(v1,v2) for v1, v2 in zip(range(gene_mat.shape[1]),s)])
  plt.title(name)
  plt.savefig("{}{}_v.png".format(outpath, name),bbox_inches="tight")
  plt.close()
  
  temp = gene_mat.dot(np.transpose(vh[:,:]))
  temp["scZ"] = temp.index.map(pd.Series(df["scZ"].values,index=df["cell_gene"]).to_dict())
  
  temp = temp.rename(columns={x : "svd_z{}".format(x) for x in range(gene_mat.shape[1])})
  temp = temp[sorted(temp.columns)]
  
  sns.heatmap(temp)
  plt.title(name)
  plt.savefig("{}{}_newzs.png".format(outpath, name),bbox_inches="tight")
  plt.close()
  
  sub_df = df.drop_duplicates("cell")


  cols = ["scZ","svd_z0","svd_z1","svd_z2"]
  if psi:
    cols.append("psi")
  for col in cols:
    try:
      for celltype, ct_df in sub_df.groupby("cell_type"):
        plt.hist(ct_df[col],30,label="{}; median: {:.03f}".format(celltype,ct_df[col].median()),density=True)
      plt.legend(bbox_to_anchor=(1.5,1.05))
      if col.startswith("svd"):
        title_str = "{}\nfrac: {:,.03f}".format(col,ct_df["f" + col[-1]].unique()[0])
      else:
        title_str = col
      plt.title("{}\n{}".format(title_str,name))
      plt.savefig("{}{}_{}_hist.png".format(outpath, name,col),bbox_inches="tight")
      plt.close()
    except Exception as e:
      print(e)

class CellType:
  def __init__(self, label, num_reads, num_cells, junc_tuples, junc_fracs):
    self.num_reads = num_reads
    self.num_cells = num_cells
    self.junc_tuples = junc_tuples
    self.junc_fracs = [y/sum(junc_fracs) for y in junc_fracs]
    self.label = label
    self.sample_cells()
    self.make_df()
    
  def sample_cells(self):
    
#     self.sample = np.random.multinomial(self.num_reads, self.junc_fracs, size=self.num_cells)
    vals = []
    # number of reads for each cell is poisson
    for i in range(self.num_cells):
      vals.append(np.random.multinomial(np.random.poisson(self.num_reads), self.junc_fracs))
    self.sample = np.vstack(vals)
  def make_df(self):
    dfs = []
    for i in range(len(self.junc_tuples)):
      temp = pd.DataFrame.from_dict({"numReads" : self.sample[:,i]})
      temp["juncPosR1A"] = self.junc_tuples[i][0]
      temp["juncPosR1B"] = self.junc_tuples[i][1]
      temp["cell_type"] = self.label
      temp["geneR1A_uniq"] = "gene1"
      temp["cell"] = ["{}_{}".format(self.label,x) for x in range(temp.shape[0])]
      temp["cell_gene"] = temp["cell"] + "_" + temp["geneR1A_uniq"]
      temp["sign"] = 1
      temp["refName_newR1"] = temp["juncPosR1A"].astype(str) + "_" + temp["juncPosR1B"].astype(str)
      dfs.append(temp)
    self.df = pd.concat(dfs,axis=0)
    self.df = self.df[self.df["numReads"] != 0]
    
  def __str__(self):
    return "{}\nnum reads: {}\nnum cells: {}\njunc_tuples: {}\njunc_fracs: {}".format(self.label,self.num_reads,self.num_cells,self.junc_tuples,self.junc_fracs)

def get_pval_df(cell_types, num_perms, num_trials, num_reads, z_cols):
  pvals = {z : {x : [] for x in [c.label for c in cell_types]} for z in z_cols}
  
  for i in range(num_trials):
    if i == 0:
      df = get_rijk_df(cell_types,num_reads, inc_juncs, exc_juncs, plot=True)
    else:
      df = get_rijk_df(cell_types,num_reads, inc_juncs, exc_juncs, plot=False)
    fdr = perm_pval(df, num_perms = num_perms, num_quants = 1)
    for z_col in z_cols:
      for key in pvals[z_col].keys():
  #       display(fdr[fdr["cell_type"] == key]["quant_pval"])
        print("pvals",pvals)
        print(z_col, key)
        try:
          pvals[z_col][key].append(float(fdr[fdr["cell_type"] == key]["quant_pval_" + z_col]))
        except Exception as e:
          pvals[z_col][key].append(np.nan)

#   pval_df = pd.DataFrame.from_dict(pvals)
  return pvals

def get_rijk_df(cell_types, num_reads, inc_juncs, exc_juncs, plot = False):
  
  pinning_S = 0.1
  # draw new sample
  for c in cell_types:
    c.sample_cells()
    c.make_df()
    
  dfs = []
  for c in cell_types:
    dfs.append(c.df)
  df = pd.concat(dfs,axis=0)

  rank_by_donor=True
  rev_let = {"A" : "B", "B" : "A"}
  let_dict = {"A" : "acc", "B" : "don"}
  full_df = df.copy()
  calc_dfs = {}
  for let in ["A","B"]:
    df = full_df
    df = prepare_df(df, let, rank_by_donor, rev_let, let_dict)
    df = calc_Sijk(df,let,pinning_S, let_dict)
    df = normalize_Sijks(df,let)
    
    
#     display(df)
    # remove those with variance == 0
    df = df[df["sijk{}_var".format(let)] != 0]

    df["n.g_" + let] = df["cell_gene"].map(df.groupby("cell_gene")["numReads"].sum())

    df["nSijk" + let] = (df["S_ijk_" + let] - df["sijk{}_mean".format(let)]) / np.sqrt(df["sijk{}_var".format(let)])
    
    df["mult"] = df["numReads"] * df["nSijk" + let]  / np.sqrt(df["n.g_" + let])
    df["z_" + let] = df["sign"] * df["cell_gene"].map(df.groupby("cell_gene")["mult"].sum())
    df["scaled_z_" + let] = df["z_" + let] / np.sqrt(df["n.g_" + let])
    df["cell_gene_junc"] = df["cell_gene"] + df["refName_newR1"]
    calc_dfs[let] = df
    ## test debug by returning
    
  df = calc_dfs["A"].merge(calc_dfs["B"],on="cell_gene_junc",how="outer",suffixes=("","_x"))
  
  for cx in [x for x in df.columns if x.endswith("_x")]:
    c = cx[:-2]
    df.loc[df[c].isna(),c] = df.loc[df[c].isna(),cx]
  df.drop([x for x in df.columns if x.endswith("_x")],inplace=True,axis=1)
  grouped = df.groupby('geneR1A_uniq')
  
#   return df, calc_dfs
  for let in ["A","B"]:
    z_dict = pd.Series(calc_dfs[let]["z_" + let].values,index=calc_dfs[let].cell_gene).to_dict()
    df["z_" + let] = df["cell_gene"].map(z_dict)
    scz_dict = pd.Series(calc_dfs[let]["scaled_z_" + let].values,index=calc_dfs[let].cell_gene).to_dict()
    df["scaled_z_" + let] = df["cell_gene"].map(scz_dict)
#   df["cov"] = df["geneR1A_uniq"].map(grouped.apply(lambda x: x['z_A'].cov(x['z_B'])))

  idx = df[df["z_A"].isna()].index
  df.loc[idx,"z"] = -df.loc[idx,"z_B"]
  df.loc[idx,"scZ"] = -df.loc[idx,"scaled_z_B"]

  idx = df[df["z_B"].isna()].index
  df.loc[idx,"z"] = df.loc[idx,"z_A"]
  df.loc[idx,"scZ"] = df.loc[idx,"scaled_z_A"]


  idx = df[(~df["z_A"].isna()) & (~df["z_B"].isna())].index
  df.loc[idx,"z"] = (df.loc[idx,"z_A"] - df.loc[idx,"z_B"])/np.sqrt(2 )
  df.loc[idx,"scZ"] = (df.loc[idx,"scaled_z_A"] - df.loc[idx,"scaled_z_B"])/np.sqrt(2 )

  df["n.g"] = df.groupby("cell_gene")["numReads"].transform("sum")
  for let in ["A","B"]:
    df["zcontrib" + let] = df["numReads"] * df["nSijk" + let] / np.sqrt(df["n.g"])


  ##################### SVD calculation #####################
  letters = ["A","B"]  
  for let in letters:
    # find number of reads per donor (or acceptor) per cell
    df["cell_gene_pos" + let] = df["cell_gene"] + df["juncPosR1" + let].astype(str)
    df["n.g_pos" + let] = df.groupby("cell_gene_pos" + let)["numReads"].transform("sum")

    # normalize on a donor/acceptor rather than a gene basis
    # TRY OUT NOT SQRT-ING denominator
    df["zcontrib_posnorm" + let] = df["numReads"] * df["nSijk" + let] / df["n.g_pos" + let]
  zcontrib_col = "zcontrib_posnorm"
  
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

  for gene, gene_df in df.groupby("geneR1A_uniq"):

    # get zcontrib matrix
    gene_mats = []
    for let in letters:
      gene_mat = gene_df.drop_duplicates("cell_gene_pos" + let).pivot_table(index="cell_gene",columns="str_juncPosR1{}".format(let),values="summed_zcontrib" + let,fill_value=0)

      gene_mats.append(gene_mat)
    gene_mat = gene_mats[0].merge(gene_mats[1],on="cell_gene")

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

  for i in range(k):
    df["f{}".format(i)] = df["geneR1A_uniq"].map(loads["f{}".format(i)])
    df["svd_z{}".format(i)] = df["cell_gene"].map(zs["svd_z{}".format(i)])
  df["svd_z_sumsq"] = (df[["svd_z{}".format(i) for i in range(k)]]**2).sum(axis=1)

  
  
  
  # calc PSI
  df = calc_PSI(df, inc_juncs, exc_juncs)
  df["center_psi"] = df["psi"] - 1
  if plot:
    name = "{}_{}".format("_".join([c.label for c in cell_types]),num_reads)
    plot_SVD(gene_mat, vh, s, df, name,psi=True)
  # sub_cols = ["cell","geneR1A_uniq","tissue","compartment","free_annotation","ontology","scZ","svd_z_sumsq","n.g_A","n.g_B"] + ["f{}".format(i) for i in range(k)] + ["svd_z{}".format(i) for i in range(k)]
  return df

def main():
  args = get_args()
#  isoform1 = .2
#  isoform2 = .25
#  isoform3 = .25
  frac1 = 0.25
  frac2 = 0.30
  frac3 = 0.45
  const_frac = 0.99
#  isoform4 = .3

  isoform1 = .15
  isoform2 = .2
  isoform3 = .25
  isoform4 = .4

#  isoform1 = .3
#  isoform2 = .1
#  isoform3 = .1
#  isoform4 = .5

  out_dict = {"read_depth" : [], "z_col" : [], "power" : [], "num_trials" : [], "num_sig" : []}
  z_cols = ["scZ","psi","svd_z0","svd_z1","svd_z2"]
  alpha = 0.05
  
  read_vals = range(1,args.max_depth + 1)
  for num_reads in tqdm(read_vals):
    cell_types = []

    if args.regime == "exon_skip_twoends":

      cell_types.append(CellType("inc",num_reads,args.num_cells,[(0,1),(0,2),(3,4),(0,4)],[args.psi/8,3*args.psi/8,args.psi/2,1 - args.psi]))
      cell_types.append(CellType("exc",num_reads,args.num_cells,[(0,1),(0,2),(3,4),(0,4)],[3*(1-args.psi)/8,(1-args.psi)/8,(1-args.psi)/2,args.psi]))

      inc_juncs = ["0_1","3_4"]
      exc_juncs = ["0_4"]

    if args.regime == "exon_skip_twoends_null":

      cell_types.append(CellType("inc",num_reads,args.num_cells,[(0,1),(0,2),(3,4),(0,4)],[args.psi/4,args.psi/4,args.psi/2,1 - args.psi]))
      cell_types.append(CellType("exc",num_reads,args.num_cells,[(0,1),(0,2),(3,4),(0,4)],[args.psi/4,args.psi/4,args.psi/2,1 - args.psi]))

      inc_juncs = ["0_1","3_4"]
      exc_juncs = ["0_4"]

    if args.regime == "exon_skip":

      cell_types.append(CellType("exc",num_reads,args.num_cells,[(0,1),(2,3),(0,3)],[args.psi/2,args.psi/2,1 - args.psi]))
      cell_types.append(CellType("inc",num_reads,args.num_cells,[(0,1),(2,3),(0,3)],[(1-args.psi)/2,(args.psi)/2,args.psi]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
      inc_juncs = ["0_1","2_3"]
      exc_juncs = ["0_3"]

    if args.regime == "exon_skip_null":

      cell_types.append(CellType("exc",num_reads,args.num_cells,[(0,1),(2,3),(0,3)],[args.psi/2,args.psi/2,1 - args.psi]))
      cell_types.append(CellType("inc",num_reads,args.num_cells,[(0,1),(2,3),(0,3)],[args.psi/2,args.psi/2,1 - args.psi]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
      inc_juncs = ["0_1","2_3"]
      exc_juncs = ["0_3"]

    if args.regime == "simple_skip":

      cell_types.append(CellType("exc",num_reads,args.num_cells,[(0,1),(0,3)],[args.psi,1 - args.psi]))
      cell_types.append(CellType("inc",num_reads,args.num_cells,[(0,1),(0,3)],[(1-args.psi),args.psi]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
      inc_juncs = ["0_1","2_3"]
      exc_juncs = ["0_3"]

    if args.regime == "simple_skip_null":

      cell_types.append(CellType("exc",num_reads,args.num_cells,[(0,1),(0,3)],[args.psi,1 - args.psi]))
      cell_types.append(CellType("inc",num_reads,args.num_cells,[(0,1),(0,3)],[args.psi,1 - args.psi]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
      inc_juncs = ["0_1","2_3"]
      exc_juncs = ["0_3"]
    
    if args.regime == "simple_skip_unbalanced":

      cell_types.append(CellType("small_inc",num_reads,20,[(0,1),(0,3)],[args.psi,1 - args.psi]))
      cell_types.append(CellType("large_inc",num_reads,1000,[(0,1),(0,3)],[args.psi, 1 - args.psi]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(0,3)],[0.25,0.25,0.5]))
      inc_juncs = ["0_1","2_3"]
      exc_juncs = ["0_3"]

    elif args.regime == "comp_skip":
      cell_types.append(CellType("r1_dom",num_reads,args.num_cells,[(0,2),(1,2),(1,3)],[args.psi,0,1 - args.psi]))
      cell_types.append(CellType("r3_dom",num_reads,args.num_cells,[(0,2),(1,2),(1,3)],[1 - args.psi,0,args.psi]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,2),(1,2),(1,3)],[1/3,1/3,1/3]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,2),(1,2),(1,3)],[1/3,1/3,1/3]))
      inc_juncs = ["0_2","1_2"]
      exc_juncs = ["1_3"]

    elif args.regime == "comp_skip_null":
      cell_types.append(CellType("r1_dom",num_reads,args.num_cells,[(0,2),(1,2),(1,3)],[args.psi,0,1 - args.psi]))
      cell_types.append(CellType("r3_dom",num_reads,args.num_cells,[(0,2),(1,2),(1,3)],[args.psi,0,1 - args.psi]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,2),(1,2),(1,3)],[1/3,1/3,1/3]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,2),(1,2),(1,3)],[1/3,1/3,1/3]))
      inc_juncs = ["0_2","1_2"]
      exc_juncs = ["1_3"]
    elif args.regime == "casset2":

      cell_types.append(CellType("casset1",num_reads,args.num_cells,[(0,1),(2,5),(0,3),(4,5)],[args.psi/2, args.psi/2, (1 - args.psi)/2, (1-args.psi)/2]))
      cell_types.append(CellType("casset2",num_reads,args.num_cells,[(0,1),(2,5),(0,3),(4,5)],[ (1 - args.psi)/2, (1-args.psi)/2,args.psi/2, args.psi/2]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,5),(0,3),(4,5)],[.25, .25, .25, .25]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,5),(0,3),(4,5)],[.25, .25, .25, .25]))
      inc_juncs = ["0_3","4_5"]
      exc_juncs = ["2_5"]

    elif args.regime == "casset2_null":

      cell_types.append(CellType("casset1",num_reads,args.num_cells,[(0,1),(2,5),(0,3),(4,5)],[args.psi/2, args.psi/2, (1 - args.psi)/2, (1-args.psi)/2]))
      cell_types.append(CellType("casset2",num_reads,args.num_cells,[(0,1),(2,5),(0,3),(4,5)],[args.psi/2, args.psi/2, (1 - args.psi)/2, (1-args.psi)/2]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,5),(0,3),(4,5)],[.25, .25, .25, .25]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,5),(0,3),(4,5)],[.25, .25, .25, .25]))
      inc_juncs = ["0_3","4_5"]
      exc_juncs = ["2_5"]

    elif args.regime == "casset1":

      cell_types.append(CellType("casset1",num_reads,args.num_cells,[(0,1),(2,5),(0,3),(4,5)],[args.psi/2, args.psi/2, (1 - args.psi)/2, (1-args.psi)/2]))
      cell_types.append(CellType("casset2",num_reads,args.num_cells,[(0,1),(2,5),(0,3),(4,5)],[ (1 - args.psi)/2, (1-args.psi)/2,args.psi/2, args.psi/2]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,5),(0,3),(4,5)],[.25, .25, .25, .25]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,5),(0,3),(4,5)],[.25, .25, .25, .25]))
      inc_juncs = ["0_1","2_5"]
      exc_juncs = ["0_3"]

    elif args.regime == "casset1_null":

      cell_types.append(CellType("casset1",num_reads,args.num_cells,[(0,1),(2,5),(0,3),(4,5)],[args.psi/2, args.psi/2, (1 - args.psi)/2, (1-args.psi)/2]))
      cell_types.append(CellType("casset2",num_reads,args.num_cells,[(0,1),(2,5),(0,3),(4,5)],[args.psi/2, args.psi/2, (1 - args.psi)/2, (1-args.psi)/2]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,5),(0,3),(4,5)],[.25, .25, .25, .25]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,5),(0,3),(4,5)],[.25, .25, .25, .25]))
      inc_juncs = ["0_1","2_5"]
      exc_juncs = ["0_3"]

    elif args.regime == "CD47_type_set2":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform1/4 + isoform2/3 + isoform3/2,isoform1/4 + isoform2/3, isoform1/4, isoform1/4,isoform4,isoform3/2,isoform2/3]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform4/4 + isoform3/3 + isoform2/2,isoform4/4 + isoform3/3, isoform4/4, isoform4/4,isoform1,isoform2/2,isoform3/3]))

      inc_juncs = ["2_3","4_5","4_7"]
      exc_juncs = ["0_7","2_7"]

    elif args.regime == "CD47_type_set2_null":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform1/4 + isoform2/3 + isoform3/2,isoform1/4 + isoform2/3, isoform1/4, isoform1/4,isoform4,isoform3/2,isoform2/3]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform1/4 + isoform2/3 + isoform3/2,isoform1/4 + isoform2/3, isoform1/4, isoform1/4,isoform4,isoform3/2,isoform2/3]))

      inc_juncs = ["2_3","4_5","4_7"]
      exc_juncs = ["0_7","2_7"]

    elif args.regime == "CD47_type_set3":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform1/4 + isoform2/3 + isoform3/2,isoform1/4 + isoform2/3, isoform1/4, isoform1/4,isoform4,isoform3/2,isoform2/3]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform4/4 + isoform3/3 + isoform2/2,isoform4/4 + isoform3/3, isoform4/4, isoform4/4,isoform1,isoform2/2,isoform3/3]))

      inc_juncs = ["4_5","6_7"]
      exc_juncs = ["0_7","2_7","4_7"]

    elif args.regime == "CD47_type_set3_null":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform1/4 + isoform2/3 + isoform3/2,isoform1/4 + isoform2/3, isoform1/4, isoform1/4,isoform4,isoform3/2,isoform2/3]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform1/4 + isoform2/3 + isoform3/2,isoform1/4 + isoform2/3, isoform1/4, isoform1/4,isoform4,isoform3/2,isoform2/3]))
      inc_juncs = ["4_5","6_7"]
      exc_juncs = ["0_7","2_7","4_7"]

    elif args.regime == "CD47_type_set1":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform1/4 + isoform2/3 + isoform3/2,isoform1/4 + isoform2/3, isoform1/4, isoform1/4,isoform4,isoform3/2,isoform2/3]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform4/4 + isoform3/3 + isoform2/2,isoform4/4 + isoform3/3, isoform4/4, isoform4/4,isoform1,isoform2/2,isoform3/3]))

      inc_juncs = ["0_1","2_3","2_7"]
      exc_juncs = ["0_7"]

    elif args.regime == "CD47_type_set1_null":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform1/4 + isoform2/3 + isoform3/2,isoform1/4 + isoform2/3, isoform1/4, isoform1/4,isoform4,isoform3/2,isoform2/3]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[isoform1/4 + isoform2/3 + isoform3/2,isoform1/4 + isoform2/3, isoform1/4, isoform1/4,isoform4,isoform3/2,isoform2/3]))


      inc_juncs = ["0_1","2_3","2_7"]
      exc_juncs = ["0_7"]

    elif args.regime == "CD47_type_set1_const":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac1/4)*const_frac + (frac1/3)*(1-const_frac) + frac2/2, (frac1/4)*const_frac + (frac1/3)*(1-const_frac), (frac1/4)*const_frac, (frac1/4)*const_frac, frac3, frac2/2, (frac1/3)*(1-const_frac)]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac3/4)*const_frac + (frac3/3)*(1-const_frac) + frac2/2, (frac3/4)*const_frac + (frac3/3)*(1-const_frac), (frac3/4)*const_frac, (frac3/4)*const_frac, frac1, frac2/2, (frac3/3)*(1-const_frac)]))



      inc_juncs = ["0_1","2_3","2_7"]
      exc_juncs = ["0_7"]

    elif args.regime == "CD47_type_set1_const_null":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac1/4)*const_frac + (frac1/3)*(1-const_frac) + frac2/2, (frac1/4)*const_frac + (frac1/3)*(1-const_frac), (frac1/4)*const_frac, (frac1/4)*const_frac, frac3, frac2/2, (frac1/3)*(1-const_frac)]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac1/4)*const_frac + (frac1/3)*(1-const_frac) + frac2/2, (frac1/4)*const_frac + (frac1/3)*(1-const_frac), (frac1/4)*const_frac, (frac1/4)*const_frac, frac3, frac2/2, (frac1/3)*(1-const_frac)]))


      inc_juncs = ["0_1","2_3","2_7"]
      exc_juncs = ["0_7"]

    elif args.regime == "CD47_type_set2_const":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac1/4)*const_frac + (frac1/3)*(1-const_frac) + frac2/2, (frac1/4)*const_frac + (frac1/3)*(1-const_frac), (frac1/4)*const_frac, (frac1/4)*const_frac, frac3, frac2/2, (frac1/3)*(1-const_frac)]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac3/4)*const_frac + (frac3/3)*(1-const_frac) + frac2/2, (frac3/4)*const_frac + (frac3/3)*(1-const_frac), (frac3/4)*const_frac, (frac3/4)*const_frac, frac1, frac2/2, (frac3/3)*(1-const_frac)]))


      inc_juncs = ["2_3","4_5","4_7"]
      exc_juncs = ["0_7","2_7"]

    elif args.regime == "CD47_type_set2_const_null":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac1/4)*const_frac + (frac1/3)*(1-const_frac) + frac2/2, (frac1/4)*const_frac + (frac1/3)*(1-const_frac), (frac1/4)*const_frac, (frac1/4)*const_frac, frac3, frac2/2, (frac1/3)*(1-const_frac)]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac1/4)*const_frac + (frac1/3)*(1-const_frac) + frac2/2, (frac1/4)*const_frac + (frac1/3)*(1-const_frac), (frac1/4)*const_frac, (frac1/4)*const_frac, frac3, frac2/2, (frac1/3)*(1-const_frac)]))


      inc_juncs = ["2_3","4_5","4_7"]
      exc_juncs = ["0_7","2_7"]

    elif args.regime == "CD47_type_set3_const":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac1/4)*const_frac + (frac1/3)*(1-const_frac) + frac2/2, (frac1/4)*const_frac + (frac1/3)*(1-const_frac), (frac1/4)*const_frac, (frac1/4)*const_frac, frac3, frac2/2, (frac1/3)*(1-const_frac)]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac3/4)*const_frac + (frac3/3)*(1-const_frac) + frac2/2, (frac3/4)*const_frac + (frac3/3)*(1-const_frac), (frac3/4)*const_frac, (frac3/4)*const_frac, frac1, frac2/2, (frac3/3)*(1-const_frac)]))



      inc_juncs = ["4_5","6_7"]
      exc_juncs = ["0_7","2_7","4_7"]

    elif args.regime == "CD47_type_set3_const_null":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac1/4)*const_frac + (frac1/3)*(1-const_frac) + frac2/2, (frac1/4)*const_frac + (frac1/3)*(1-const_frac), (frac1/4)*const_frac, (frac1/4)*const_frac, frac3, frac2/2, (frac1/3)*(1-const_frac)]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(frac1/4)*const_frac + (frac1/3)*(1-const_frac) + frac2/2, (frac1/4)*const_frac + (frac1/3)*(1-const_frac), (frac1/4)*const_frac, (frac1/4)*const_frac, frac3, frac2/2, (frac1/3)*(1-const_frac)]))


      inc_juncs = ["4_5","6_7"]
      exc_juncs = ["0_7","2_7","4_7"]

    elif args.regime == "CD47_type":

      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[1/4 + (1 + args.psi)/3 + (1 + 2*args.psi)/2, 1/4 + (1 + args.psi)/3, 1/4, 1/4,1 + 3*args.psi,(1 + 2*args.psi)/2,(1 + args.psi)/3]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[(1 + 3*args.psi)/4 + (1 + 2*args.psi)/3 + (1 + args.psi)/2, (1 + 3*args.psi)/4 + (1 + 2*args.psi)/3, (1 + 3*args.psi)/4, (1 + 3*args.psi)/4,1,(1 + args.psi)/2,(1 + 2*args.psi)/3]))
      inc_juncs = ["0_1","2_3","2_7"]
      exc_juncs = ["0_7"]
      
    elif args.regime == "CD47_type_null":

      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[1/4 + (1 + args.psi)/3 + (1 + 2*args.psi)/2, 1/4 + (1 + args.psi)/3, 1/4, 1/4,1 + 3*args.psi,(1 + 2*args.psi)/2,(1 + args.psi)/3]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,7),(2,7),(4,7)],[1/4 + (1 + args.psi)/3 + (1 + 2*args.psi)/2, 1/4 + (1 + args.psi)/3, 1/4, 1/4,1 + 3*args.psi,(1 + 2*args.psi)/2,(1 + args.psi)/3]))

      inc_juncs = ["0_1","2_3","2_7"]
      exc_juncs = ["0_7"]


    elif args.regime == "one_donor":

#      cell_types.append(CellType("ascending",num_reads,100,[(0,1),(0,2),(0,3)],[1 + args.psi,1, 1 + 2*args.psi]))
#      cell_types.append(CellType("descending",num_reads,100,[(0,1),(0,2),(0,3)],[1,1 + 2*args.psi,1 + args.psi]))
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(0,2),(0,3)],[1,1 + args.psi, 1 + 2*args.psi]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(0,2),(0,3)],[1 + 2*args.psi,1 + args.psi,1]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(0,2),(0,3)],[1/3,1/3,1/3]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(0,2),(0,3)],[1/3,1/3,1/3]))
      inc_juncs = ["0_2"]
      exc_juncs = ["0_3"]

    elif args.regime == "one_donor_fair":

      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(0,2),(0,3)],[1 + args.psi,1, 1 + 2*args.psi]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(0,2),(0,3)],[1,1 + 2*args.psi,1 + args.psi]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(0,2),(0,3)],[1/3,1/3,1/3]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(0,2),(0,3)],[1/3,1/3,1/3]))
#      inc_juncs = ["0_2"]
#      exc_juncs = ["0_3"]

      inc_juncs = ["0_1"]
      exc_juncs = ["0_2","0_3"]

    elif args.regime == "one_donor_fair_null":

      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(0,2),(0,3)],[1,1 + args.psi, 1 + 2*args.psi]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(0,2),(0,3)],[1,1 + args.psi, 1 + 2*args.psi]))

#      cell_types.append(CellType("ascending",num_reads,100,[(0,1),(0,2),(0,3)],[1 + args.psi,1, 1 + 2*args.psi]))
#      cell_types.append(CellType("descending",num_reads,100,[(0,1),(0,2),(0,3)],[1 + args.psi,1, 1 + 2*args.psi]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(0,2),(0,3)],[1/3,1/3,1/3]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(0,2),(0,3)],[1/3,1/3,1/3]))
      inc_juncs = ["0_1"]
      exc_juncs = ["0_2","0_3"]

    elif args.regime == "one_donor_null":
      cell_types.append(CellType("ascending",num_reads,args.num_cells,[(0,1),(0,2),(0,3)],[1,1 + args.psi,1 + 2*args.psi]))
      cell_types.append(CellType("descending",num_reads,args.num_cells,[(0,1),(0,2),(0,3)],[1,1 + args.psi,1 + 2*args.psi]))

#      cell_types.append(CellType("ascending",num_reads,100,[(0,1),(0,2),(0,3)],[1 + args.psi,1, 1 + 2*args.psi]))
#      cell_types.append(CellType("descending",num_reads,100,[(0,1),(0,2),(0,3)],[1 + args.psi,1, 1 + 2*args.psi]))

#     cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(0,2),(0,3)],[1/3,1/3,1/3]))
#     cell_types.append(CellType("equal",num_reads,1000,[(0,1),(0,2),(0,3)],[1/3,1/3,1/3]))
      inc_juncs = ["0_2"]
      exc_juncs = ["0_3"]

    elif args.regime == "double":
      cell_types.append(CellType("double_inc",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[args.psi/4,args.psi/4,args.psi/4,args.psi/4,(1-args.psi)/2,(1-args.psi)/2]))
      cell_types.append(CellType("double_exc",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[(1-args.psi)/4,(1-args.psi)/4,(1-args.psi)/4,(1-args.psi)/4,args.psi/2,args.psi/2]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[.125,.125,.125,.125,.25,.25]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[.125,.125,.125,.125,.25,.25]))
      inc_juncs = ["0_1","2_3"]
      exc_juncs = ["0_3"]

    elif args.regime == "double2":
      cell_types.append(CellType("double_inc",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[args.psi/4,args.psi/4,args.psi/4,args.psi/4,(1-args.psi)/2,(1-args.psi)/2]))
      cell_types.append(CellType("double_exc",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[(1-args.psi)/4,(1-args.psi)/4,(1-args.psi)/4,(1-args.psi)/4,args.psi/2,args.psi/2]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[.125,.125,.125,.125,.25,.25]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[.125,.125,.125,.125,.25,.25]))
      inc_juncs = ["4_5","6_7"]
      exc_juncs = ["4_7"]
      
    elif args.regime == "double2_null":
      cell_types.append(CellType("double_inc",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[args.psi/4,args.psi/4,args.psi/4,args.psi/4,(1-args.psi)/2,(1-args.psi)/2]))
      cell_types.append(CellType("double_exc",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[args.psi/4,args.psi/4,args.psi/4,args.psi/4,(1-args.psi)/2,(1-args.psi)/2]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[.125,.125,.125,.125,.25,.25]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[.125,.125,.125,.125,.25,.25]))
      inc_juncs = ["4_5","6_7"]
      exc_juncs = ["4_7"]

    elif args.regime == "double_null":
      cell_types.append(CellType("double_inc",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[args.psi/4,args.psi/4,args.psi/4,args.psi/4,(1-args.psi)/2,(1-args.psi)/2]))
      cell_types.append(CellType("double_exc",num_reads,args.num_cells,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[args.psi/4,args.psi/4,args.psi/4,args.psi/4,(1-args.psi)/2,(1-args.psi)/2]))
#      cell_types.append(CellType("equal_small",num_reads,100,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[.125,.125,.125,.125,.25,.25]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(4,5),(6,7),(0,3),(4,7)],[.125,.125,.125,.125,.25,.25]))
      inc_juncs = ["0_1","2_3"]
      exc_juncs = ["0_3"]


    elif args.regime == "double_switch":
      cell_types.append(CellType("double1",num_reads,args.num_cells,[(0,1),(2,3),(4,7),(0,3),(4,5),(6,7)],[args.psi/3,args.psi/3,args.psi/3,(1 - args.psi)/3,(1 - args.psi)/3,(1 - args.psi)/3]))
      cell_types.append(CellType("double2",num_reads,args.num_cells,[(0,1),(2,3),(4,7),(0,3),(4,5),(6,7)],[(1 - args.psi)/3,(1 - args.psi)/3,(1 - args.psi)/3,args.psi/3,args.psi/3,args.psi/3]))
#      cell_types.append(CellType("equal_small",num_reads,1000,[(0,1),(2,3),(4,7),(0,3),(4,5),(6,7)],[1,1,1,1,1,1]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(4,7),(0,3),(4,5),(6,7)],[1,1,1,1,1,1]))
      inc_juncs = ["0_1","2_3"]
      exc_juncs = ["0_3"]

    elif args.regime == "double_switch_null":
      cell_types.append(CellType("double1",num_reads,args.num_cells,[(0,1),(2,3),(4,7),(0,3),(4,5),(6,7)],[args.psi/3,args.psi/3,args.psi/3,(1 - args.psi)/3,(1 - args.psi)/3,(1 - args.psi)/3]))
      cell_types.append(CellType("double2",num_reads,args.num_cells,[(0,1),(2,3),(4,7),(0,3),(4,5),(6,7)],[args.psi/3,args.psi/3,args.psi/3,(1 - args.psi)/3,(1 - args.psi)/3,(1 - args.psi)/3]))
#      cell_types.append(CellType("equal_small",num_reads,1000,[(0,1),(2,3),(4,7),(0,3),(4,5),(6,7)],[1,1,1,1,1,1]))
#      cell_types.append(CellType("equal",num_reads,1000,[(0,1),(2,3),(4,7),(0,3),(4,5),(6,7)],[1,1,1,1,1,1]))
      inc_juncs = ["0_1","2_3"]
      exc_juncs = ["0_3"]


    pvals = {y : [] for y in z_cols}
    for i in range(args.num_trials):
      df = get_rijk_df(cell_types,num_reads, inc_juncs, exc_juncs, plot=False)
      out = variance_adj_pval(df, z_cols, num_perms = args.num_perms)
      for z_col in z_cols:
        if not np.isnan(out["perm_pval_adj_" + z_col][0]):
          pvals[z_col].append(out["perm_pval_adj_" + z_col][0])
    for z_col in z_cols:
      out_dict["read_depth"].append(num_reads)
      out_dict["z_col"].append(z_col)
      out_dict["num_trials"].append(args.num_trials)

      if len(pvals[z_col]) == 0:
        out_dict["power"].append(0)
        out_dict["num_sig"].append(0)
      else:
        out_dict["num_sig"].append(len([x for x in pvals[z_col] if x < alpha]))
        out_dict["power"].append(len([x for x in pvals[z_col] if x < alpha])/args.num_trials)

  out_df = pd.DataFrame.from_dict(out_dict)

  out_df.to_csv("{}{}_{}_{}_{}_{}.tsv".format(outpath,args.regime,args.num_trials,args.num_perms,args.psi,args.num_cells),sep="\t",index=False)
  z_color = {"psi" : u'#ff7f0e', "scZ" : u'#1f77b4', "svd_z0" : u'#2ca02c', "svd_z1" : u'#9467bd', "svd_z2" : u'#8c564b'}
  for z_col, z_df in out_df.groupby("z_col"):
    lower_ci, upper_ci = proportion_confint(z_df["num_sig"], z_df["num_trials"])
    z_df["lower_ci"] = lower_ci
    z_df["upper_ci"] = upper_ci
    
    plt.errorbar(z_df["read_depth"],z_df["power"],yerr=[z_df["power"] - z_df["lower_ci"], z_df["upper_ci"] - z_df["power"]],label=z_col,marker="o",alpha=0.6,color=z_color[z_col])

  plt.axhline(y=0.05, label="y = 0.05",color="k",linestyle="--")
  plt.legend(bbox_to_anchor=(1.05, 1))
  plt.xlabel("mean number of reads")
  plt.ylabel("power")
  plt.ylim([-.025,1.025])
  plt.title("{}\nnum perms: {} num trials: {} psi: {} num cells: {}".format(args.regime,args.num_perms, args.num_trials, args.psi, args.num_cells))
  plt.savefig("{}{}_{}_{}_{}_{}.png".format(outpath,args.regime,args.num_trials,args.num_perms,args.psi,args.num_cells),bbox_inches="tight")
  plt.close()
      
  # pval_df = pd.DataFrame.from_dict(pvals)


main()
