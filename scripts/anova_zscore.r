library(data.table)

args = commandArgs(trailingOnly=TRUE)


#datanames = c("TS_10x_redo","TSP1_SS2","TSP2_10x_rerun_3prime","TSP2_SS2")
#datanames = c("TS_10x_redo","TSP1_SS2","TSP2_10x_rerun_3prime")
datanames <- args[2:length(args)]

#datanames = c("TS_10x_redo_bestrefseq","TSP1_SS2_bestrefseq","TSP2_10x_rerun_3prime_bestrefseq")

inpath <- 'scripts/output/significant_genes/'
outpath <- 'scripts/output/anova_zscore/'
#suff <- '_unfilt'
suff <- args[1]
out_suff <- suff
z_col <- "scZ"
resid <- FALSE
if (resid) {
  out_suff <- paste0(out_suff,"_residual")
}
for (dataname in datanames) {
  
  # read in data
  df <- read.table(file = paste0(inpath,dataname,'-',z_col,'_allp',suff,'.tsv'), sep = '\t', header = TRUE)
  
  # require at least 2 compartments, 2 tissues, and 5 ontologies per gene 
  setDT(df)[, num_tiss := uniqueN(tissue),by = geneR1A_uniq]
  setDT(df)[, num_comp := uniqueN(compartment), by = geneR1A_uniq]
  setDT(df)[, num_ont := uniqueN(ontology), by = geneR1A_uniq]
  df <- df[which(df$num_comp > 1 & df$num_ont > 4) ]
  if (length(unique(df$tissue)) > 1) {
    df <- df[which(df$num_tiss > 1) ]

  }
  
  # make empty output dataframe
#  freeanns <- unique(df$condensed_freeann)
  compartments = unique(df$compartment)
  tissues <- unique(df$tissue)
  genes <- unique(df$geneR1A_uniq)
  
#  out <- data.frame(matrix(ncol = length(freeanns) + length(compartments) + length(tissues) + 1, nrow = length(genes)))
#  colnames(out) <- c(c("gene"),as.character(freeanns), as.character(compartments), as.character(tissues))
  out <- data.frame(matrix(ncol =  1, nrow = length(genes)))
  colnames(out) <- c("gene")
  out$gene <- genes
  out_unweight <- data.frame(out)

  out_coeff <- data.frame(matrix(ncol =  1, nrow = length(genes)))
  colnames(out_coeff) <- c("gene")
  out_coeff$gene <- genes
  out_coeff_unweight <- data.frame(out_coeff)
 
  # add one hot encoding for free annotations
#  for (freeann in freeanns) {
#    df[,freeann] <- df$condensed_freeann == freeann
#  }
  
  # add one hot encoding for tissue
  for (tissue in tissues) {
    df[,tissue] <- df$tissue == tissue
  }
  
  # add one hot encoding for compartment
  for (compartment in compartments) {
    df[,compartment] <- df$compartment == compartment
  }
  
  # get only columns necessary for linear model
#  cols <- c(compartments,tissues,c("is_stem","mz"))
  cols <- c(compartments,tissues,"mz")

  
  # loop over all genes
  for (i in 1:length(genes)) {
    
    gene <- out[i,"gene"]
    temp <- subset(df,subset=geneR1A_uniq == gene)
    
    # remove variables that would overdefine the model (only one entry for a value)
    temp_cols <- c("mz")
#    if (length(which(table(temp$is_stem) > 1)) == 2) {
#      temp_cols <- c(temp_cols,"is_stem")
#    }
#    temp_cols <- c(temp_cols,names(which(table(temp$condensed_freeann) > 1 )))
    temp_cols <- c(temp_cols,names(which(table(temp$compartment) > 1 )))

    if (resid) {
      tiss_cols <- c("mz",names(which(table(temp$tissue) > 1 )))
    } else {
      temp_cols <- c(temp_cols,names(which(table(temp$tissue) > 1 )))
    }

    if (resid) {
      # fit weighted linear model using only tissues as predictors
      model1 <- lm(mz~. -1,weights = temp$n_cells_ont,data=temp[,..tiss_cols])
      temp[["mz"]] <- resid(model1)
    }
    
    # fit weighted linear model on tissue, compartment, and all condensed free annotations
    model2 <- lm(mz~. -1,weights = temp$n_cells_ont,data=temp[,..temp_cols])
    anova_out2 <- anova(model2)  

    # remove weighting
    model3 <- lm(mz~. -1,data=temp[,..temp_cols])
    anova_out3 <- anova(model3)  


    # parse output p values
    for (row in rownames(anova_out2)) {
      if (row != "Residuals"){
        
        out[i,row] <- anova_out2[row,"Pr(>F)"]
      }
    }
    # save linear model coefficients
    for (name in names(model2$coefficients)) {
      out_coeff[i,name] <- model2$coefficients[name]
    
      
    }

    # parse output p values not weighted
    for (row in rownames(anova_out3)) {
      if (row != "Residuals"){
        
        out_unweight[i,row] <- anova_out3[row,"Pr(>F)"]
      }
    }
    # save linear model coefficients
    for (name in names(model3$coefficients)) {
      out_coeff_unweight[i,name] <- model3$coefficients[name]
    
      
    }

    
  }
  
  write.table(out,file=paste0(outpath,dataname,"_",z_col,'_anova_out',out_suff,'.tsv'),row.names=FALSE,sep="\t",quote=FALSE)
  write.table(out_coeff,file=paste0(outpath,dataname,"_",z_col,'_anova_coeff',out_suff,'.tsv'),row.names=FALSE,sep="\t",quote=FALSE)
  write.table(out_unweight,file=paste0(outpath,dataname,"_",z_col,'_anova_out_unweight',out_suff,'.tsv'),row.names=FALSE,sep="\t",quote=FALSE)
  write.table(out_coeff_unweight,file=paste0(outpath,dataname,"_",z_col,'_anova_coeff_unweight',out_suff,'.tsv'),row.names=FALSE,sep="\t",quote=FALSE)

}
