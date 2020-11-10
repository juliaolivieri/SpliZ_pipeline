library(data.table)
#dataname <- "TS_10x_redo"
args = commandArgs(trailingOnly=TRUE)


inpath <- 'scripts/output/significant_genes/'
outpath <- 'scripts/output/anova_factor/'
#suff <- '_unfilt'
suff <- args[1]
dataname <- args[2]

out_suff <- suff
z_col <- "scZ"

#  read in data
df <- read.table(file = paste0(inpath,dataname,'-',z_col,'_allp',suff,'.tsv'), sep = '\t', header = TRUE)

# require at least 2 compartments, 2 tissues, and 5 ontologies per gene 
setDT(df)[, num_tiss := uniqueN(tissue),by = geneR1A_uniq]
setDT(df)[, num_comp := uniqueN(compartment), by = geneR1A_uniq]
setDT(df)[, num_ont := uniqueN(ontology), by = geneR1A_uniq]
df <- df[which(df$num_comp > 1 & df$num_ont > 4) ]
if (length(unique(df$tissue)) > 1) {
  df <- df[which(df$num_tiss > 1) ]

}

compartments = unique(df$compartment)
tissues <- unique(df$tissue)
genes <- unique(df$geneR1A_uniq)
  
out <- data.frame(matrix(ncol =  1, nrow = length(genes)))
colnames(out) <- c("gene")
out$gene <- genes
out_unweight <- data.frame(out)

out_coeff <- data.frame(matrix(ncol = 1, nrow = length(genes)))

colnames(out_coeff) <- c("gene")
out_coeff$gene <- genes
out_coeff_unweight <- data.frame(out_coeff)

for (i in 1:length(genes)) {
  gene <- out[i,"gene"]
  temp <- subset(df,subset=geneR1A_uniq == gene)
  temp <- temp[temp$compartment %in% names(which(table(temp$compartment) > 1)),] 
  if (length(tissues) > 1){
    temp <- temp[temp$tissue %in% names(which(table(temp$tissue) > 1)),] 
  }
  if (length(unique(temp$compartment)) > 1) {
    if ( (length(unique(temp$tissue)) > 1) | (length(tissues) == 1)) {

      if (length(tissues) > 1) {
        model <- lm(mz ~ tissue + compartment - 1, weights = temp$n_cells_ont,data=temp)   # fit weighted linear model on tissue, compartment, and all condensed free annotations

      } else {
        model <- lm(mz ~ compartment - 1, weights = temp$n_cells_ont,data=temp)   # fit weighted linear model on tissue, compartment, and all condensed free annotations

      }
      if (length(tissues) > 1) {
        model2 <- lm(mz ~ tissue + compartment - 1,data=temp)   # fit weighted linear model on tissue, compartment, and all condensed free annotations

      } else {
        model2 <- lm(mz ~ compartment - 1, data=temp)   # fit weighted linear model on tissue, compartment, and all condensed free annotations

      }

      anova_out <- anova(model)
      anova_out2 <- anova(model2)

      for (row in rownames(anova_out)) {
        if (row != "Residuals"){

          out[i,row] <- anova_out[row,"Pr(>F)"]

        }
        for (name in names(model$coefficients)) {
          out_coeff[i,name] <- model$coefficients[name]
        }
    }
      for (row in rownames(anova_out2)) {
        if (row != "Residuals"){

          out_unweight[i,row] <- anova_out2[row,"Pr(>F)"]

        }
        for (name in names(model2$coefficients)) {
          out_coeff_unweight[i,name] <- model2$coefficients[name]
        }
    }
    }
}
}  
write.table(out,file=paste0(outpath,dataname,"_",z_col,'_anova_out',out_suff,'.tsv'),row.names=FALSE,sep="\t",quote=FALSE)
write.table(out_coeff,file=paste0(outpath,dataname,"_",z_col,'_anova_coeff',out_suff,'.tsv'),row.names=FALSE,sep="\t",quote=FALSE)

write.table(out_unweight,file=paste0(outpath,dataname,"_",z_col,'_anova_out_unweight',out_suff,'.tsv'),row.names=FALSE,sep="\t",quote=FALSE)
write.table(out_coeff_unweight,file=paste0(outpath,dataname,"_",z_col,'_anova_coeff_unweight',out_suff,'.tsv'),row.names=FALSE,sep="\t",quote=FALSE)

