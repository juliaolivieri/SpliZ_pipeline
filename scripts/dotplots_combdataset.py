import argparse
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
import numpy as np
import os
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({'font.size': 14})

## create a color dictionary for donors
#def compartment_colors(compartments):
#    
#    import matplotlib.colors as pltcolors
#    
#    cmap = plt.cm.get_cmap("YlOrRd")
#        
#    compartment_color_dict = {}
#    j=1/len(compartments)
#    for c in compartments:
#        compartment_color_dict[c] = pltcolors.to_hex(cmap(j))
#        j+=1/len(compartments)
#        
#    return compartment_color_dict

# create a color dictionary for donors (our colors)
def compartment_colors(compartments):
    compartment_color_dict = {comp : col for comp, col in zip(compartments,sns.color_palette("deep",len(compartments)))}
    return compartment_color_dict


def get_args():
  parser = argparse.ArgumentParser(description="make dotplots per donor")

  parser.add_argument("--params",help="tab separated tsv file with columns 'gene','let', and 'end'")


  args = parser.parse_args()
  return args

def annotation_plot(gtf, domains, gene, end,outpath,datanames,zoom=True, plot_all=True):
  rev_dict = {"A" : "B","B" : "A"}

  # gene = "CALD1"
  # end = 134928754
  don_df = pd.read_csv("/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/dotplots_combdataset/{}/{}_{}_{}_coords.tsv".format(gene,"_".join(datanames),gene,end),sep="\t")
  don_df = don_df.rename(columns={"rank_acc" : "rank", "rank_don" : "rank"}).astype(int)


  if don_df["juncPosR1A"].nunique() == 1:
    let = "A"
  else:
    let = "B"

  shared_ends = list(don_df["juncPosR1" + rev_dict[let]])

  # for gene in genes:
  gene_gtf = gtf[gtf["gene_id"] == gene]
  gene_gtf = gene_gtf[gene_gtf["feature"].isin(["exon"])]

  legend = True
  name_end = ""

  if don_df["juncPosR1" + rev_dict[let]].nunique() > 1:

    colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf',u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
    someX, someY = 0.5, 0.5
    plt.figure(figsize=(12, 6))
    h = 1
    offset = 1
    currentAxis = plt.gca()
    count = 1
    y_labels = []
    arc_height = 10
    y_ticks = []
    arcs = False
    chromosome = gene_gtf["seqname"].iloc[0]
    gene_min_all = gene_gtf["start"].min()
    gene_max_all = gene_gtf["end"].max()

    legend_elements = []

    gene_domains = domains[(domains[1] == chromosome) & (domains[2] < gene_max_all) & (domains[3] > gene_min_all)]

    # if arcs:
    # count = 1
    if gene_gtf["strand"].iloc[0] == "+":
      asc = True
    else:
      asc = False

  #   plt.text(row["juncPosR1" + rev_dict[let]],gene_gtf["transcript_id"].nunique() + 1,row["rank"],horizontalalignment="center")
    for ind, row in don_df.iterrows():
        if plot_all:
          name_end += "_all"
          plt.text(row["juncPosR1" + rev_dict[let]],gene_gtf["transcript_id"].nunique() + 1,row["rank"],horizontalalignment="center")
          plt.plot([row["juncPosR1" + rev_dict[let]],row["juncPosR1" + rev_dict[let]]],[ -0.5,gene_gtf["transcript_id"].nunique() + 0.5],color="red")
    plt.plot([row["juncPosR1" + let],row["juncPosR1" + let]],[ -0.5,gene_gtf["transcript_id"].nunique() + 0.5],color="blue")


    count = 0
    for transcript, trans_df in gene_gtf.groupby("transcript_id"):
      y_labels.append(transcript)
      y_ticks.append(count * offset)
      gene_min = trans_df["start"].min()
      gene_max = trans_df["end"].max()
      plt.plot([gene_min,gene_max],[offset * count,offset * count],color="k")

      for exons in set([tuple(x) for x in trans_df[["start","end"]].to_numpy()]):
        plot_exon(exons,currentAxis,offset = count * offset)
      i = 0
      for d in set([tuple(x) for x in gene_domains[[2,3,4]].to_numpy()]):
        plot_exon(d[:2],currentAxis,offset = count * offset,alpha = 0.4,color=colors[i],ecolor=None,h=0.5)
        legend_elements.append(Patch(facecolor=colors[i], edgecolor=None,label=d[2],alpha=0.4))
        i += 1

      count += 1

    if arcs:
      plt.ylim([-1,gene_gtf["transcript_id"].nunique() + arc_height + 1])
    else:
      plt.ylim([-1,gene_gtf["transcript_id"].nunique() + 2])

    currentAxis.ticklabel_format(useOffset=False,style="plain")
    if legend:
      currentAxis.legend(handles=legend_elements[:gene_domains.shape[0]],bbox_to_anchor=(1., 1.0))
    plt.yticks(y_ticks,y_labels)
  #   plt.xlim([gene_min_all,gene_max_all])
    if zoom:
      print("shared_ends",shared_ends,"end",end)
      buff = max([abs(int(x) - int(end)) for x in shared_ends])/12
      plt.xlim(min([int(end)] + shared_ends) - buff,max([int(end)] + shared_ends) + buff)
      name_end += "_zoom"
    else:
      plt.xlim([gene_min_all,gene_max_all])

    plt.title("{} {} {} don".format(gene,chromosome,end))
    try:
      plt.savefig("{}{}/{}_{}_{}_don_ann{}.png".format(outpath,gene,"_".join(datanames),gene,end, name_end),bbox_inches = "tight")
      plt.close()
    except Exception as e:
      print(e)

def plot_exon(locs,ax,h=0.25,offset = 0,bin_size=400,alpha=1,color="k",ecolor="k"):
#   print("plotting exon")
  ax.add_patch(Rectangle((locs[0], -h + offset), locs[1] - locs[0], 2*h,edgecolor=ecolor,color=color,alpha=alpha,linewidth=0))
  
def load_gtf(gtf_file,filt_chr):
  gtf = pd.read_csv(gtf_file, names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"], sep="\t")
  if 'NC_000023.11' in gtf["seqname"].unique():
    can_chrom = [x for x in gtf["seqname"].unique() if x.startswith("NC_")]
    name_dict = {x : "chr" + str(int(x.split("_")[1].split(".")[0])) for x in can_chrom}
    name_dict['NC_000023.11'] = "chrX"
    name_dict['NC_000024.10'] = "chrY"
    name_dict['NC_012920.1'] = "chrM"
    gtf["seqname"] = gtf["seqname"].map(name_dict)
    gtf = gtf[~gtf["seqname"].isna()]
    filt_chr = False
  #   gtf = gtf[gtf["feature"] == "exon"]
  try:
    gtf["gene_id"] = gtf["attribute"].str.split("gene_name").str[1].str.split(";").str[0].str.split('"').str[1]
  except:
    gtf["gene_id"] = gtf["attribute"].str.split("gene_id").str[1].str.split(";").str[0].str.split('"').str[1]
  gtf["transcript_id"] = gtf["attribute"].str.split("transcript_id").str[1].str.split(";").str[0].str.split('"').str[1]
#  filt_chr = True
  if filt_chr:
    # don't include immature scaffolds
    chromosomes = gtf["seqname"].unique()
    chromosomes = [x for x in chromosomes if "_" not in x and not x.startswith("KN")]
#    print("chromosomes",chromosomes)
    gtf = gtf[gtf["seqname"].isin(chromosomes)]
  gtf["chr_gene"] = gtf["seqname"] + gtf["gene_id"]
  return gtf


def tissue_colors():
    
    tissue_color_dict = {'Bladder': '#e7969c',
             'Blood': '#d6616b',
             'Bone_Marrow': '#cedb9c',
#              'Fat': '#e7cb94',
#              'Heart': '#637939',
             'Kidney': '#7b4173',
             'Large_Intestine': '#31a354',
             'Lung': '#3182bd',
             'Lymph_Node': '#8c6d31',
             'Muscle': '#e7ba52',
             'Pancreas': '#fd8d3c',
             'Skin': '#ce6dbd',
             'Small_Intestine': '#6baed6',
             'Spleen': '#393b79',
             'Thymus': '#9c9ede',
             'Trachea': '#969696',
             'Vasculature': '#843c39'}
    
    return tissue_color_dict

def test_forward(array):
  
  return [float(x) / 2 for x in array]

def test_backward(array):
  return [x * 2 for x in array]

def dot_plot(don_df, let, let_dict, palette, onts, outpath, gene, don, tiss, dataname, rev_dict):
#  print("in dot")
  print("a")
  don_df["ontology_rank"] = don_df["ontology"] + don_df["rank_" + let_dict[let]].astype(str)
  don_df["rank_count"] = don_df["ontology_rank"].map(don_df.groupby("ontology_rank")["numReads"].sum())
  ont_dict = {o : i for o, i in zip(onts,range(1,don_df["ontology"].nunique() + 1))}
  don_df["ont_num"] = don_df["ontology"].map(ont_dict)
  pdf = don_df.drop_duplicates("ontology_rank")
  pdf["rank_sum"] = pdf["ontology"].map(pdf.groupby("ontology")["rank_count"].sum())
  pdf["frac_rank"] = pdf["rank_count"] / pdf["rank_sum"]
#  pdf["rank_sum"] = pdf["ontology"].map(pdf.groupby("ontology")["scaled_rank"].sum())
#  pdf["scaled_rank"] = pdf["scaled_rank"] / pdf["rank_sum"]
  print("b")
  sns.relplot(x="rank_" + let_dict[let], y="ont_num", size="frac_rank",
              sizes=(10, 400), alpha=.5, palette=palette,hue="compartment",
              height=max(4,pdf["ontology"].nunique()*0.3), data=pdf)

#     onts = []
#     vals = []
#     for key, value in ont_dict.items():
#       onts.append(key)
#       vals.append(value)
  plt.yticks(range(1,don_df["ontology"].nunique() + 1),onts)
  plt.title("{}\n{} {} {} {}".format(dataname,gene,tiss, don, let_dict[rev_dict[let]]))
#   plt.savefig("{}{}_{}_{}_{}_{}_dot.png".format(outpath, gene, don, tiss, dataname, let_dict[rev_dict[let]]),bbox_inches="tight")
#  print("saved")
  plt.close()
  print("c")
  return 0

def plot_df(df, let, cell_lim, outpath, gene, dataname, let_dict, palette, rev_dict, don,tiss, don_df, comp_sort):
#    df["num_cells"] = df["ontology"].map(df.groupby("ontology")["cell"].nunique())
#
#    df = df[df["num_cells"] > cell_lim]

    # calculate ontology-wide values to sort on
    df["ont_rank"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].median())
    df["ont_75"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].quantile(0.75))
    df["ont_25"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].quantile(0.25))

    df["ont_max"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].max())
    df["ont_min"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].min())
    df = df.sort_values(["ont_rank","ont_75","ont_25","ont_max","ont_min"])

    num_cells = list(df.drop_duplicates("ontology")["num_cells"])
    medians = list(df.drop_duplicates("ontology")["ont_rank"])
    if df.shape[0] > 0:
      fig,(ax2) = plt.subplots(1)
      g = sns.boxplot(x="avg_rank", y="ontology",hue="compartment",dodge=False,
                       data=df,
                       orient="h", palette=palette)
#       g = sns.violinplot(x="avg_rank", y="ontology",hue="compartment",dodge=False,
#                        data=df,
#                        orient="h", palette=palette)
      for i,artist in enumerate(ax2.artists):
          # Set the linecolor on the artist to the facecolor, and set the facecolor to None
          col = artist.get_facecolor()
          artist.set_edgecolor(col)
          artist.set_facecolor('None')

          # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
          # Loop over them here, and use the same colour as above
          for j in range(i*6,i*6+6):
              line = ax2.lines[j]
              line.set_color(col)
              line.set_mfc(col)
              line.set_mec(col)

      # Also fix the legend
      for legpatch in ax2.get_legend().get_patches():
          col = legpatch.get_facecolor()
          legpatch.set_edgecolor(col)
          legpatch.set_linewidth(10)
          legpatch.set_facecolor('None')
#         medians = df.groupby("ontology")["avg_rank"].median()
      for i in range(len(num_cells)):
        plt.text(df["avg_rank"].max() + (df["avg_rank"].max() - df["avg_rank"].min())/12,i, num_cells[i])


        plt.scatter([medians[i]],[i],color = "k",s=20,zorder=100)
      plt.title("{}\n{} {} {} {}\nmean: {:0.2f} median: {:0.2f}".format(dataname,gene,tiss, don, let_dict[rev_dict[let]], don_df["avg_rank"].mean(), don_df["avg_rank"].median()))
      plt.legend(bbox_to_anchor=(1.5, 1.05))
#       plt.savefig("{}{}_{}_{}_{}_{}.png".format(outpath, gene, don, tiss, dataname, let_dict[rev_dict[let]]),bbox_inches="tight")
      plt.close()
      return df

def box(df, let, cell_lim, outpath, gene, dataname, let_dict, palette, rev_dict, comp_sort):

  for don, don_df in df.groupby("pos{}_group".format(let)):
    temp = plot_df(don_df.drop_duplicates("pos{}_cell".format(let)), let, cell_lim, outpath, gene, dataname, let_dict, palette, rev_dict, don,"all", don_df, comp_sort)
    if not temp is None:
      if temp["ontology"].nunique() > 0:
        onts = list(temp.sort_values(["ont_rank","ont_75","ont_25","ont_max","ont_min"]).drop_duplicates("ontology")["ontology"].unique())

        onts.reverse()
        response = dot_plot(don_df, let, let_dict, palette,onts, outpath, gene, don, "all", dataname, rev_dict)
    for tiss, tiss_df in don_df.groupby("tissue"):
      temp = plot_df_vioin(tiss_df.drop_duplicates("pos{}_cell".format(let)), let, cell_lim, outpath, gene, dataname, let_dict, palette, rev_dict, don,tiss, don_df)
      if not temp is None:
        if temp["ontology"].nunique() > 0:
          onts = list(temp.sort_values(["ont_rank","ont_75","ont_25","ont_max","ont_min"]).drop_duplicates("ontology")["ontology"].unique())

          onts.reverse()
          response = dot_plot(tiss_df, let, let_dict, palette,onts,outpath, gene, don, tiss, dataname, rev_dict)
          if type(response) != int:
            print("ERROR")
            return response

def main():
  args = get_args()
  outpath = "/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/dotplots_combdataset/"
  plot_domains = True
  if plot_domains:
    domains = pd.read_csv("/oak/stanford/groups/horence/Roozbeh/single_cell_project/utility_files/ucscGenePfam.txt",sep="\t",header=None)
  else:
    domains = pd.DataFrame(columns=[1,2,3,4,5,6,7])
  # load gtf
  # gtf_file = "/scratch/PI/horence/JuliaO/single_cell/STAR_wrapper/gtf_files/hg38.ncbiRefSeq.gtf" 
  # gtf_file = "/oak/stanford/groups/horence/circularRNApipeline_Cluster/index/hg38_genes.gtf"
  gtf_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/utility_files/GRCh38_latest_genomic.gtf"
  gtf = load_gtf(gtf_file,True)
  sig_color = False
  comp_sort = True
#  gtf["gene_id"] = gtf["gene_id"].replace("FYB1","FYB")
#  gtf["gene_id"] = gtf["gene_id"].replace("ATP5F1C","ATP5C1")


#  datanames = ["TS_10x_redo","TSP2_10x_rerun_3prime","TSP1_SS2"]
#  datanames = ["TSP1_10x_nopanc_with_postprocessing","TSP2_10x_3prime_with_postprocessing","TSP1_SS2"]
#  datanames = ["TSP1_10x_nopanc_with_postprocessing","TSP2_10x_3prime_with_postprocessing"]
  datanames = ["HLCA4_P2_10x_with_postprocessing_lung","HLCA4_P3_10x_with_postprocessing_lung"]

  if sig_color: 
    if datanames[0].startswith("HLCA"):
      FDR_col = "FDR_HLCA4_P2_10x_with_postprocessing_lung_HLCA4_P3_10x_with_postprocessing_lung"
    elif datanames[0].startswith("TSP"):
      FDR_col = "FDR_TSP1_10x_nopanc_with_postprocessing_TSP2_10x_3prime_with_postprocessing"

  sig_onts = defaultdict(lambda : set())
  if sig_color:
    for dataname in datanames:
      sig = pd.read_csv("/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/final_FDRs_mz/{}_FDR_S_0.1_z_0.0_b_5.tsv".format(dataname),sep="\t")
      sig = sig[sig[FDR_col] < 0.05] 
      for g, gene_df in sig.groupby("geneR1A_uniq"):
        sig_onts[g] = set.union(sig_onts[g],set(gene_df["ontology"]))
  
  # datanames = ["TS_10x_redo","TSP2_10x_rerun_3prime"]
  df_dict = {}
  gene = "all"
  for dataname in datanames:
  #   dataname = "TS_10x_redo"
    df = pd.read_parquet("/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/rijk_zscore/{}_sym_S_0.1_z_0.0_b_5.pq".format(dataname),columns=["tissue","compartment","free_annotation","geneR1A_uniq","posA_group","posB_group","cell", "juncPosR1A", "juncPosR1B", "numReads","sign", "splice_ann"])
    compartments = sorted([x for x in list(df["compartment"].unique()) if x != None])
  
    palette = compartment_colors(compartments)

#    outpath = "/scratch/PI/horence/JuliaO/single_cell/Differential_Splicing/scripts/output/donor_boxplots/"
    if gene != "all":
      df = df[df["geneR1A_uniq"] == gene]
#      df = df[df["juncPosR1B"].isin([56161387,56160626])]
      df = df[df["juncPosR1B"].isin([108047292])]

    # only include ontology/gene pairs with at least 20 cells
    cell_lim = 10
    df["ontology"] = df["tissue"] + df["compartment"] + df["free_annotation"]
    df["ontology1"] = df["compartment"] + df["tissue"] +  df["free_annotation"]
    df["ontology2"] = df["free_annotation"]

    df["ontology_gene"] = df["ontology"] + df["geneR1A_uniq"]
    df["ontology_don"] = df["ontology"] + df["posA_group"]
    df["ontology_acc"] = df["ontology"] + df["posB_group"]
    df["num_cells"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell"].nunique())
  
    df = df[df["num_cells"] > cell_lim]
    df["num_cell_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell"].nunique())
    df = df[df["num_cell_ont"] > cell_lim]
    df_dict[dataname] = df
#    print("df shape",df.shape)
#  print("keys",df_dict.keys())
  
#  df["ontology"] = df["compartment"] + " " + df["tissue"] + " " + df["free_annotation"]
  let_dict = {"A" : "acc", "B" : "don"}
  rev_dict = {"A" : "B", "B" : "A"}

  param_df = pd.read_csv(args.params,sep="\t")

#  param_df = pd.DataFrame.from_dict({"gene" : ["MYL6"], "let" : ["A"], "end" : ["56160320"]})
#  param_df = pd.DataFrame.from_dict({"gene" : ["COMMD6","COMMD6"], "let" : ["A","B"], "end" : [75530311,75530311]})


  print("param_df",param_df)

  for index, row in tqdm(param_df.iterrows()):
      try:
       
#        print("outpath",outpath + gene + "/")
        gene = row["gene"]
        print("gene",gene)
        if not os.path.exists(outpath + gene + "/"):
            os.makedirs(outpath + gene + "/")
    
        let = row["let"]
        end = row["end"]
    #  gene = "MYL6"
    #  end = 56160320
    #  let = "A"
        gene_dfs = {}
        # palette = {"Endothelial" : u'#1f77b4', "Stromal" : u'#ff7f0e', "Epithelial" : u'#2ca02c', "Immune" : u'#d62728'}
#        palette = compartment_colors()
        
        
        for dataname in datanames:
          
    
          print(dataname)
#          print("end",end,"gene",gene)
#          print(df_dict[dataname][df_dict[dataname]["geneR1A_uniq"] == gene]["juncPosR1" + let].value_counts().head())
#          end2 = df_dict[dataname][df_dict[dataname]["geneR1A_uniq"] == gene]["juncPosR1" + let].value_counts().index[0]
#          print("end2",end2)
#          print("end2 == end", end2 == end)
#          print("int(end2) == end",int(end2) == end)
#          print("end2 == float(end)",end2 == float(end))
#          print("try 1",df_dict[dataname][(df_dict[dataname]["geneR1A_uniq"] == gene) & (df_dict[dataname]["juncPosR1" + let] == end)].shape[0])
#          print("try 2",df_dict[dataname][(df_dict[dataname]["juncPosR1" + let] == end)].shape[0])
   
#          gene_df = df_dict[dataname][(df_dict[dataname]["geneR1A_uniq"] == gene)]
#          gene_df["juncPosR1" + let] = gene_df["juncPosR1" + let].astype(int)
#          print("try 3",gene_df[(gene_df["juncPosR1" + let] == end)].shape[0])
   
          gene_df = df_dict[dataname][(df_dict[dataname]["geneR1A_uniq"] == gene) & (df_dict[dataname]["juncPosR1" + let] == float(end))]
          print("gene df 1",dataname,end,gene_df.shape[0])
            
          gene_df["num_" + let_dict[let]] = gene_df["pos{}_group".format(let)].map(gene_df.groupby("pos{}_group".format(let))["pos{}_group".format(rev_dict[let])].nunique())
          gene_df = gene_df[gene_df["num_" + let_dict[let]] > 1]
          gene_df["pos{}_cell".format(let)] = gene_df["pos{}_group".format(let)] + gene_df["cell"]
          gene_df["rank_{}".format(let_dict[let])] = gene_df.groupby("pos{}_group".format(let))["juncPosR1{}".format(rev_dict[let])].rank(method="dense")
          gene_df["scaled_rank"] = gene_df["rank_" + let_dict[let]] * gene_df["numReads"]
          gene_df["num"] = gene_df["pos{}_cell".format(let)].map(gene_df.groupby("pos{}_cell".format(let))["scaled_rank"].sum())
          gene_df["denom"] = gene_df["pos{}_cell".format(let)].map(gene_df.groupby("pos{}_cell".format(let))["numReads"].sum())
          gene_df["avg_rank"] = gene_df["num"]/gene_df["denom"]
          gene_dfs[dataname] = gene_df
        ont_list = []
        for dataname, gene_df in gene_dfs.items():
          print("A")
  #        print("gene_df",gene_df.shape)
          temp = plot_df(gene_df.drop_duplicates("pos{}_cell".format(let)), let, cell_lim, "", gene, dataname, let_dict, palette, rev_dict, end,"all", gene_df, comp_sort)
          if not temp is None:
            print("B")
            onts = list(temp.sort_values(["ont_rank","ont_75","ont_25","ont_max","ont_min"]).drop_duplicates("ontology")["ontology"].unique())

            onts.reverse()
            ont_list.append(set(onts))
            print("E")
            #   dot_plot(gene_df, let, let_dict, palette,onts, outpath, gene, don, "all", dataname, rev_dict)
        print("len ont list",len(ont_list))
        if len(ont_list) > 0:
          if len(ont_list) == 2:
            print("F")
            if len(set.union(*ont_list)) < 20:
              print("1")
              shared_onts = set.union(*ont_list)
            elif len(set.intersection(*ont_list).union(ont_list[0])) < 20:
              print("2")
              shared_onts = set.intersection(*ont_list).union(ont_list[0])
            elif len(set.intersection(*ont_list).union(ont_list[1])) < 20:
              print("3")
              shared_onts = set.intersection(*ont_list).union(ont_list[1])
            else:
              print("4")
              print("ont_list",ont_list)
              shared_onts = set.intersection(*ont_list)
              print("shared_onts",len(shared_onts),shared_onts)
  
          elif len(ont_list) == 3:
            if len(set.union(*ont_list)) < 20:
              shared_onts = set.union(*ont_list)
            elif len(set.union(*ont_list[:2]).intersection(ont_list[2]).union(set.union(*ont_list[:2]))) < 20:
              shared_onts = set.union(*ont_list[:2]).intersection(ont_list[:3]).union(set.union(*ont_list[:2]))
            elif len(ont_list[1].intersection(ont_list[2]).union(ont_list[0])) < 20:
              shared_onts = ont_list[1].intersection(ont_list[2]).union(ont_list[0])
            elif len(set.union(ont_list[0].intersection(ont_list[1]), ont_list[1].intersection(ont_list[2]), ont_list[0].intersection(ont_list[2]))) < 20:
              shared_onts = set.union(ont_list[0].intersection(ont_list[1]), ont_list[1].intersection(ont_list[2]), ont_list[0].intersection(ont_list[2]))
            else:
              shared_onts = set.intersection(*ont_list)
                
          
          print("G")
          for dataname, gene_df in gene_dfs.items():
            print("H") 
            dot_plot(gene_df[gene_df["ontology"].isin(shared_onts)], let, let_dict, palette,shared_onts, outpath, gene, end, "all", dataname, rev_dict)
          count = 0
          pdfs = []
          shift = 0.2
          print("C") 
          
          for dataname, gene_df in gene_dfs.items():
            print("D")
            print(dataname)
            gene_df["rank_" + let_dict[let]] = gene_df["rank_" + let_dict[let]].rank(method="dense")
            tiss = "all"
            don_df = gene_df[gene_df["ontology"].isin(shared_onts)]
            don_df["ontology_rank"] = don_df["ontology"] + don_df["rank_" + let_dict[let]].astype(str)
            don_df["rank_count"] = don_df["ontology_rank"].map(don_df.groupby("ontology_rank")["numReads"].sum())
            
            ont_dict = {o : i for o, i in zip(onts,range(1,don_df["ontology"].nunique() + 1))}
            
            don_df["ont_num"] = don_df["ontology"].map(ont_dict)
            pdf = don_df.drop_duplicates("ontology_rank")
            pdf["rank_sum"] = pdf["ontology"].map(pdf.groupby("ontology")["rank_count"].sum())
            pdf["frac_rank"] = pdf["rank_count"] / pdf["rank_sum"]
            #  pdf["rank_sum"] = pdf["ontology"].map(pdf.groupby("ontology")["scaled_rank"].sum())
            #  pdf["scaled_rank"] = pdf["scaled_rank"] / pdf["rank_sum"]
            pdf["rank_" + let_dict[let]] = pdf["rank_" + let_dict[let]] + count * shift
            pdf["dataset"] = dataname
            pdfs.append(pdf)
      
            count += 1
          full_pdf = pd.concat(pdfs,axis=0)
          print("created full pdf")
          ax = sns.relplot(x="rank_" + let_dict[let], y="ont_num", size="frac_rank",
                      sizes=(10, 400), alpha=.5, palette=palette,hue="compartment",
                      height=max(4,pdf["ontology"].nunique()*0.3), data=full_pdf)
          plt.yticks(range(1,don_df["ontology"].nunique() + 1),onts)
          # if let == "B":
          #   don = acc
          plt.title("{}\n{} {} {} {}".format(dataname,gene,tiss, end, let_dict[rev_dict[let]]))
          plt.xticks([1 + shift*(len(datanames) - 1)/2,2 + shift*(len(datanames) - 1)/2],[1,2])
          # plt.savefig("{}{}_{}.png".format(outpath,gene,end),bbox_inches="tight")
          plt.close()
          ont_dict = {}
    #         print(full_pdf.groupby("rank_" + let_dict[let])["frac_rank"].sum().sort_values().index[-1])
      
#          for i in range(len(compartments)):
          for ont, ont_df in full_pdf[full_pdf["rank_" + let_dict[let]] == full_pdf.groupby("rank_" + let_dict[let])["frac_rank"].sum().sort_values().index[-1]].groupby("ontology"):
            ont_dict[ont] = ont_df["frac_rank"].sum()
          onts = sorted(ont_dict, key=ont_dict.get)
          print("in dot")
          # onts = shared_onts
          count = 0
          pdfs = []
          shift = 0.2
  #            shared_ends = (set(gene_dfs["TSP1_10x_nopanc_with_postprocessing"]["juncPosR1" +rev_dict[let]]).intersection(set(gene_dfs["TSP2_10x_3prime_with_postprocessing"]["juncPosR1" + rev_dict[let]]))).intersection(set(gene_dfs["TSP1_SS2"]["juncPosR1" + rev_dict[let]]))
          shared_ends = set.intersection(*[set(v["juncPosR1" + rev_dict[let]].unique()) for  v in gene_dfs.values()])
          print("found shared ends") 

          ont_dict = {o : i for o, i in zip(onts,range(1,len(onts) + 1))}
        
          
          for dataname, gene_df in gene_dfs.items():
            print("dataname",2)
          #   gene_df = gene_dfs[datanames[0]]
          #   dataname = datanames[0]
            gene_df = gene_df[gene_df["juncPosR1" + rev_dict[let]].isin(shared_ends)]
          #   gene_df = gene_df[gene_df["rank_" + let_dict[let]].isin(gene_df["rank_" + let_dict[let]].value_counts().index[:2])]
            gene_df["rank_" + let_dict[let]] = gene_df["rank_" + let_dict[let]].rank(method="dense")
            tiss = "all"
            don_df = gene_df[gene_df["ontology"].isin(shared_onts)]
            don_df["ontology_rank"] = don_df["ontology"] + don_df["rank_" + let_dict[let]].astype(str)
            don_df["rank_count"] = don_df["ontology_rank"].map(don_df.groupby("ontology_rank")["numReads"].sum())
          #   ont_dict = {o : i for o, i in zip(onts,range(1,don_df["ontology"].nunique() + 1))}
            don_df["ont_num"] = don_df["ontology"].map(ont_dict)
            pdf = don_df.drop_duplicates("ontology_rank")
            pdf["rank_sum"] = pdf["ontology"].map(pdf.groupby("ontology")["rank_count"].sum())
            pdf["frac_rank"] = pdf["rank_count"] / pdf["rank_sum"]
            #  pdf["rank_sum"] = pdf["ontology"].map(pdf.groupby("ontology")["scaled_rank"].sum())
            #  pdf["scaled_rank"] = pdf["scaled_rank"] / pdf["rank_sum"]
            pdf["rank_" + let_dict[let]] = pdf["rank_" + let_dict[let]] + count * shift
            pdf["dataset"] = dataname
            pdfs.append(pdf)
            count += 1
          factor = 0.5
#          factor = 0.32

          full_pdf = pd.concat(pdfs,axis=0)
          ann_dict = pd.Series(full_pdf.splice_ann.values,index=full_pdf["rank_" + let_dict[let]]).to_dict() 
          if comp_sort:
#            full_pdf = full_pdf.sort_values("ontology1")
            full_pdf["ont_num"] = full_pdf["ontology1"].map({k : v for k, v in zip(sorted(full_pdf["ontology1"].unique()),reversed(range(1,len(onts) +1)))})
  
          coords = full_pdf.drop_duplicates("rank_" + let_dict[let])[["rank_" + let_dict[let],"juncPosR1A","juncPosR1B"]]
          coords = coords[coords["rank_" + let_dict[let]].isin(range(1,int(coords["rank_" + let_dict[let]].max() + 1)))].sort_values("rank_" + let_dict[let])
          coords.to_csv("{}{}_{}_{}_coords.tsv".format(outpath + gene + "/","_".join(datanames),gene,end),sep="\t",index=False)
          annotation_plot(gtf, domains, gene, end,outpath,datanames)
          annotation_plot(gtf, domains, gene, end,outpath,datanames,plot_all=False)

          annotation_plot(gtf, domains, gene, end,outpath,datanames,zoom=False)
          annotation_plot(gtf, domains, gene, end,outpath,datanames,zoom=False,plot_all=False)


         
          print("before relplot") 
          g = sns.relplot(x="rank_" + let_dict[let], y="ont_num", size="frac_rank",
                      sizes=(10, 400), alpha=1, palette=palette,hue="compartment", 
                      height=max(4,pdf["ontology"].nunique()*factor),
                           aspect=(len(datanames)*len(shared_ends)/1.5)/max(4,pdf["ontology"].nunique()*factor), 
                           data=full_pdf)
  #            for tick_label in g.ax.get_xticklabels():
  #              try:
  #                if float(tick_label.get_text()) in ann_dict.keys():
  #                  if ann_dict[int(float(tick_label.get_text()))]:
  #                    tick_label.set_color("blue")
  #                  else:
  #                    tick_label.set_color("red")
  #              except Exception as e:
  #                print(Exception)
          ontology1 = list(full_pdf.drop_duplicates("ontology1")["ontology1"])
          ontology2 = list(full_pdf.drop_duplicates("ontology1")["ontology2"])
          ontology1, ontology2 = (list(t) for t in zip(*sorted(zip(ontology1, ontology2))))
          if comp_sort:
#            print(reverse(range(1,len(onts) + 1)),ontology2)
            plt.yticks(range(1,len(onts) + 1),reversed(ontology2))
          else:
            plt.yticks(range(1,len(onts) + 1),onts)
      #    if let == "B":
      #      don = acc
          plt.title("{}\n{} {} {} {}".format(" ".join(datanames),gene,tiss, end, let_dict[rev_dict[let]]))
          plt.xticks([x + shift*(len(datanames) - 1)/2 for x in range(1,len(shared_ends) + 1)],range(1,len(shared_ends) + 1))#[1 + shift*(len(datanames) - 1)/2,2 + shift*(len(datanames) - 1)/2],[1,2])
          ax = plt.gca()
          
          for i in range(len(ax.get_xticklabels())):
          #   try:
              x = ax.get_xticklabels()[i]
              if float(x.get_text()) in ann_dict.keys():
  #                  print("in dict",x.get_text())
  #                  print("ann dict",ann_dict[int(float(x.get_text()))])
                if ann_dict[int(float(x.get_text()))]:
  #                    print("changed blue")
                  ax.get_xticklabels()[i].set_color('blue') 
                else:
                    ax.get_xticklabels()[i].set_color('red') 

          for i in range(len(ax.get_yticklabels())):
            x = ax.get_yticklabels()[i]
            if x.get_text() in sig_onts[gene]:
              ax.get_yticklabels()[i].set_color('green')  
          # KEEP
          plt.savefig("{}{}_{}_{}.png".format(outpath + gene + "/","_".join(datanames),gene,end),bbox_inches="tight")
          print("saved at","{}{}_{}_{}.png".format(outpath + gene + "/","_".join(datanames),gene,end))
          plt.close()

      except Exception as e:
        print(e)
  
  
#    # genes = ["CAST"]
#    # end = 96675539
#    # let = "B"
#    zoom = False
#    # for gene in genes:
#    gene_gtf = gtf[gtf["gene_id"] == gene]
#    gene_gtf = gene_gtf[gene_gtf["feature"].isin(["exon"])]
#    # gene_df = df[(df["geneR1A_uniq"] == gene)]
#    # gene_df = gene_df.drop_duplicates("refName_newR1")
#    legend = True
#    #   gene_df = gene_df[gene_df["rank_acc"].isin(gene_df["rank_acc"].value_counts().index[:2])]
#    #   gene_df["rank_acc"] = gene_df["rank_acc"].rank(method="dense")
#    # for don, don_df in gene_df.drop_duplicates("refName_newR1").groupby("juncPosR1A"):
#    #   if don_df["juncPosR1B"].nunique() > 1:
#    
#    # colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
#    colors = [ u'#2ca02c',  u'#9467bd',    u'#17becf']
#    
#    someX, someY = 0.5, 0.5
#    plt.figure(figsize=(12, 6))
#    h = 1
#    offset = 1
#    currentAxis = plt.gca()
#    # plt.plot([gene_min,gene_max],[0,0],color="k")
#    count = 1
#    y_labels = []
#    arc_height = 10
#    y_ticks = []
#    arcs = False
#    chromosome = gene_gtf["seqname"].iloc[0]
#    gene_min_all = gene_gtf["start"].min()
#    gene_max_all = gene_gtf["end"].max()
#    
#    legend_elements = []
#    
#    gene_domains = domains[(domains[1] == chromosome) & (domains[2] < gene_max_all) & (domains[3] > gene_min_all)]
#    
#    # if arcs:
#    # count = 1
#    if gene_gtf["strand"].iloc[0] == "+":
#      asc = True
#    else:
#      asc = False
#      
#    plt.plot([end,end],[ -0.5,gene_gtf["transcript_id"].nunique() + 0.5],color="blue")
#    
#    shared_ends = list(shared_ends)
#    shared_ends.sort()
#    count = 0
#    for shared_end in shared_ends:
#      count += 1
#      plt.plot([shared_end,shared_end],[ -0.5,gene_gtf["transcript_id"].nunique() + 0.5],color="red")
#      plt.text(shared_end,gene_gtf["transcript_id"].nunique() + 1,count,horizontalalignment="center")
#        
#        
#    count = 0
#    for transcript, trans_df in gene_gtf.groupby("transcript_id"):
#      y_labels.append(transcript)
#      y_ticks.append(count * offset)
#      gene_min = trans_df["start"].min()
#      gene_max = trans_df["end"].max()
#      plt.plot([gene_min,gene_max],[offset * count,offset * count],color="k")
#    
#    #         print(set([tuple(x) for x in trans_df[["start","end"]].to_numpy()]))
#      for exons in set([tuple(x) for x in trans_df[["start","end"]].to_numpy()]):
#        plot_exon(exons,currentAxis,offset = count * offset)
#      i = 0
#      for d in set([tuple(x) for x in gene_domains[[2,3,4]].to_numpy()]):
#        plot_exon(d[:2],currentAxis,offset = count * offset,alpha = 0.4,color=colors[i],ecolor=None,h=0.5)
#        legend_elements.append(Patch(facecolor=colors[i], edgecolor=None,label=d[2],alpha=0.4))
#        i += 1
#    
#      count += 1
#    
#    if arcs:
#      plt.ylim([-1,gene_gtf["transcript_id"].nunique() + arc_height + 1])
#    else:
#      plt.ylim([-1,gene_gtf["transcript_id"].nunique() + 2])
#    
#    currentAxis.ticklabel_format(useOffset=False,style="plain")
#    if legend:
#      currentAxis.legend(handles=legend_elements,bbox_to_anchor=(1., 1.0))
#    plt.yticks(y_ticks,y_labels)
#    
#    if zoom:
#      buff = max([abs(x - end) for x in shared_ends])/12
#      plt.xlim(min([end] + shared_ends) - buff,max([end] + shared_ends) + buff)
#    else:
#      plt.xlim([gene_min_all,gene_max_all])
#    #       plt.xlim([7802500,gene_max_all])
#    # matplotlib.scale.FuncScale(currentAxis, (test_forward,test_backward))
#    # currentAxis.set_xscale(functions=(test_forward,test_backward))
#    plt.title("{} {} {} ".format(gene,chromosome,end))
#    #   plt.title("{} {} ({})".format(gene,reordered_bins[i],i))
#    try:
#      # KEEP
#      plt.savefig("{}{}_{}_{}_don_ann.png".format(outpath + gene + "/",gene,end, dataname),bbox_inches = "tight")
#      plt.close()
#    except Exception as e:
#      print(e)
#    #         print(don,"error")
    
main()
