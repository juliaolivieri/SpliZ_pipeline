# Written By : Roozbeh Dehghannasiri (rdehghan@stanford.edu)
# This scripts annotates the the most variable splice sites of the genes across celltypes found based on the eigenvectors
# It determines whether splice sites are annotated exon boundaries, or known to be involved in AS, are in an UTR, and are conserved in mouse and lemur genomes (i.e., also found to be a top splice site in mouse/lemur genomes).

library(data.table)
library(ggplot2)

############  input files #########################
known_splice_sites_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/SICILIAN_references/human/ncbi_refseq/hg38_refseq_known_splice_sites.txt"
top_coords_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/processed_files/zscores/PCA/coordinates_to_plot_normdonor/TSP2_10x.tsv"
top_coords_second_evec_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/processed_files/zscores/PCA/coordinates_to_plot_normdonor/second_evec_TSP2_10x.tsv"
top_coords_third_evec_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/processed_files/zscores/PCA/coordinates_to_plot_normdonor/third_evec_TSP2_10x.tsv"
zscores_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/processed_files/zscores/zscore_files/TSP2_10x_rerun_with_postprocessing_3prime_cellann_sym_SVD_normdonor_S_0.1_z_0.0_b_5.tsv"
UTR_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/utility_files/biomart_hg38_UTR.txt"
exon_annotations_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/processed_files/splicing_calls/TSP2_10x_with_postprocessing_exon_annotation.txt"
mouse_liftover_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/annotate_donors/TSP2_10x_top_coords_mouse_liftover_mapping.txt"
lemur_liftover_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/annotate_donors/TSP2_10x_top_coords_lemur_liftover_mapping.txt"
lemur_top_coords_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/processed_files/zscores/PCA/coordinates_to_plot_normdonor/Lemur_10x_Antoine.tsv"
mouse_top_coords_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/processed_files/zscores/PCA/coordinates_to_plot_normdonor/Tabula_muris_senis_P1_10x.tsv"
gtf_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/SICILIAN_references/human/ncbi_refseq/hg38_refseq_UTR_coords.txt"
output_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/TS_biology_paper/outputs/TSP2_10x_top_splice_sites_annotation.txt"
##################################################

##### read in input files ######################
known_splice_sites = fread(known_splice_sites_file)
UTR = fread(UTR_file)
eigen_values = fread(zscores_file,select=c("geneR1A_uniq","f0","f1","f2"),sep="\t",header=TRUE)
top_coords = fread(top_coords_file)
top_coords_second_evec = fread(top_coords_second_evec_file)
top_coords_third_evec = fread(top_coords_third_evec_file)
exon_annotations = fread(exon_annotations_file)
mouse_liftover = fread(mouse_liftover_file,header=TRUE)
lemur_liftover = fread(lemur_liftover_file,header=TRUE)
lemur_top_coords = fread(lemur_top_coords_file)
mouse_top_coords = fread(mouse_top_coords_file)
gtf = fread(gtf_file)
################################################

### merging top coords from all 3 eigen vectors ######
top_coords = merge(top_coords,gtf[!duplicated(gene_name),list(gene_name,chr,strand)],all.x=TRUE,all.y=FALSE,by.x="gene",by.y="gene_name")
top_coords = top_coords[!is.na(chr)]
top_coords[,evec:=1]
a = ave(top_coords$gene, top_coords$gene, FUN=seq_along)
top_coords$site_rank=a
top_coords_second_evec = merge(top_coords_second_evec,gtf[!duplicated(gene_name),list(gene_name,chr,strand)],all.x=TRUE,all.y=FALSE,by.x="gene",by.y="gene_name")
top_coords_second_evec = top_coords_second_evec[!is.na(chr)]
top_coords_second_evec[,evec:=2]
a = ave(top_coords_second_evec$gene, top_coords_second_evec$gene, FUN=seq_along)
top_coords_second_evec$site_rank=a
top_coords_third_evec = merge(top_coords_third_evec,gtf[!duplicated(gene_name),list(gene_name,chr,strand)],all.x=TRUE,all.y=FALSE,by.x="gene",by.y="gene_name")
top_coords_third_evec = top_coords_third_evec[!is.na(chr)]
top_coords_third_evec[,evec:=3]
a = ave(top_coords_third_evec$gene, top_coords_third_evec$gene, FUN=seq_along) ## add the tank of each site in the given eigenvector
top_coords_third_evec$site_rank=a
top_coords = rbind(top_coords,top_coords_second_evec,top_coords_third_evec)
######################################################


# now we want to do the same analysis by looking at the splice sites that are involved in at least two distinct junctions
known_splice_sites[,chr_V2:=paste(V1,V2,sep=":"),by=paste(V1,V2)]
known_splice_sites[,chr_V3:=paste(V1,V3,sep=":"),by=paste(V1,V3)]
known_splice_sites[,num_uniq_V2_for_V3:=length(unique(chr_V2)),by = chr_V3] # number of partner splice cites for each V2
known_splice_sites[,num_uniq_V3_for_V2:=length(unique(chr_V3)),by = chr_V2] # number of partner splice cites for each V3
alt_v2 = known_splice_sites[num_uniq_V3_for_V2>1]$V2 # the V2 coordinates that have more than one splice site partner
alt_v3 = known_splice_sites[num_uniq_V2_for_V3>1]$V3 # the V3 coordinates that have more than one splice site partner
total = c(alt_v2, alt_v3, alt_v2-1, alt_v3-1, alt_v2+1, alt_v3+1) # I concatenate all coordinates with their +-1 counterparts
top_coords[,AS_annot:=0]
top_coords[end%in%total,AS_annot:=1]

## Also if the driving splice sites are annotated exons
exon_annotations_1 = exon_annotations[,list(juncPosR1A,exon_annR1A)]
exon_annotations_2 = exon_annotations[,list(juncPosR1B,exon_annR1B)]
names(exon_annotations_1) = c("pos","exon_ann")
names(exon_annotations_2) = c("pos","exon_ann")
exon_annotations = rbind(exon_annotations_1,exon_annotations_2)
top_coords = merge(top_coords,exon_annotations[!duplicated(pos)],all.x=TRUE,all.y=FALSE,by.x="end",by.y="pos")


## prepare coords for liftover
top_coords[,end_1:=end+1,by=end]
top_coords[,for_liftover:=paste(chr,":",end,"-",end_1,sep="")]
#write.table(top_coords,"/oak/stanford/groups/horence/Roozbeh/single_cell_project/annotate_donors/TSP2_10x_top_coords_for_liftover.txt",sep="\t",row.names=FALSE,quote=FALSE)


## now we want to compare with lemur and mouse liftovers and top coordinates
a = ave(mouse_top_coords$gene, mouse_top_coords$gene, FUN=seq_along)
mouse_top_coords$site_rank_mouse = a 
a = ave(lemur_top_coords$gene, lemur_top_coords$gene, FUN=seq_along)
lemur_top_coords$site_rank_lemur = a
top_coords = merge(top_coords,lemur_liftover[!duplicated(human)],all.x=TRUE,all.y=FALSE,by.x="for_liftover",by.y="human")
top_coords = merge(top_coords,mouse_liftover[!duplicated(human)],all.x=TRUE,all.y=FALSE,by.x="for_liftover",by.y="human")
top_coords[,lemur_coord:=strsplit(lemur_pos,split="[:-]")[[1]][2],by=1:nrow(top_coords)]
top_coords[,mouse_coord:=strsplit(mouse_pos,split="[:-]")[[1]][2],by=1:nrow(top_coords)]
top_coords$lemur_coord = as.numeric(top_coords$lemur_coord)
top_coords$mouse_coord = as.numeric(top_coords$mouse_coord)
top_coords[,lemur_ever_sig:=0]
top_coords[,mouse_ever_sig:=0]
mouse_top_coords$end=as.numeric(mouse_top_coords$end)
lemur_top_coords$end=as.numeric(lemur_top_coords$end)
top_coords = merge(top_coords,mouse_top_coords[,list(end,site_rank_mouse)],all.x=TRUE,all.y=FALSE,by.x="mouse_coord",by.y="end")
top_coords = merge(top_coords,lemur_top_coords[,list(end,site_rank_lemur)],all.x=TRUE,all.y=FALSE,by.x="lemur_coord",by.y="end")


## sometimes it is possible that just by adding one to the coordinates they become conserved so I check that one as well
top_coords[,lemur_coord_1:=lemur_coord+1]
top_coords[,mouse_coord_1:=mouse_coord+1]
top_coords = merge(top_coords,mouse_top_coords[,list(end,site_rank_mouse)],all.x=TRUE,all.y=FALSE,by.x="mouse_coord_1",by.y="end")
top_coords = merge(top_coords,lemur_top_coords[,list(end,site_rank_lemur)],all.x=TRUE,all.y=FALSE,by.x="lemur_coord_1",by.y="end")
top_coords[is.na(site_rank_mouse.x) & !is.na(site_rank_mouse.y),site_rank_mouse:=site_rank_mouse.y]
top_coords[!is.na(site_rank_mouse.x) & is.na(site_rank_mouse.y),site_rank_mouse:=site_rank_mouse.x]
top_coords[is.na(site_rank_lemur.x) & !is.na(site_rank_lemur.y),site_rank_lemur:=site_rank_lemur.y]
top_coords[!is.na(site_rank_lemur.x) & is.na(site_rank_lemur.y),site_rank_lemur:=site_rank_lemur.x]
top_coords[!is.na(site_rank_lemur),lemur_ever_sig:=1]
top_coords[!is.na(site_rank_mouse),mouse_ever_sig:=1]
top_coords[,c("site_rank_lemur.x","site_rank_lemur.y","site_rank_mouse.x","site_rank_mouse.y"):=NULL]
top_coords[,end_1:=NULL]; top_coords[,mouse_pos:=NULL]; top_coords[,lemur_pos:=NULL]
top_coords[,for_liftover:=NULL]; top_coords[,lemur_coord_1:=NULL]; top_coords[,mouse_coord_1:=NULL] 

## see if the splice site is in a UTR region
for (counter in 1:nrow(top_coords)){
  UTR_gene = UTR[`Gene name`==top_coords[counter]$gene]
  p = top_coords$end[counter]
  is.5UTR <-  any(UTR_gene[(!is.na(`5' UTR start`))]$`5' UTR start` <= p & UTR_gene[(!is.na(`5' UTR start`))]$`5' UTR end` >= p)
  is.3UTR <-  any(UTR_gene[(!is.na(`3' UTR start`))]$`3' UTR start` <= p & UTR_gene[(!is.na(`3' UTR start`))]$`3' UTR end` >= p)
  top_coords$is.5UTR[counter]=is.5UTR
  top_coords$is.3UTR[counter]=is.3UTR
}

top_coords = merge(top_coords,eigen_values[!duplicated(geneR1A_uniq)],all.x=TRUE,all.y=FALSE,by.x="gene",by.y="geneR1A_uniq")
top_coords[,ever_sig:=0]
top_coords[lemur_ever_sig+mouse_ever_sig>0,ever_sig:=1]


# since for some of the splice sites in TSP2 10x the annotation status for some of the exons were wrong I fix them manually here.  
top_coords[,chr_end:=paste(chr,end)]
annotated_coords = c("chrX 77902643","chrX 77899593","chrX 41342653","chrX 41343351","chrX 41343737","chrX 72181828","chrX 48902042","chrX 23675106","chrX 23682396","chrX 23683705","chrX 48577109","chrX 154400464","chrX 154399941","chrX 72275544","chrX 53403566","chrX 53414870","chrX 153796552","chrX 153797725","chrX 47585331","chrX 47585543","chrX 47657572","chrX 47657782","chrX 48902517","chrX 53412253","chr4 121681671","chr2 177219314","chr22 37956842","chrX 48902732","chr20 50267663","chr1 24648859","chrX 47657183","chrX 41341484")
unannotated_coords = c("chrX 12976968","chrX 12976947","chrX 72272705")
top_coords[chr_end%in%annotated_coords,exon_ann:=TRUE]
top_coords[chr_end%in%unannotated_coords,exon_ann:=FALSE]
top_coords[,chr_end:=NULL]

## collect driving splice sites based on their eigenvalues
driving_top_coords = rbind(top_coords[evec==1],top_coords[(f0<0.7) & (evec==2)& (site_rank==1)])
driving_top_coords = rbind(driving_top_coords,top_coords[(f0<0.7) & (f2>0.2)& (evec==3) &(site_rank==1)])


# ## plot summary of the driving splice sites
 counts = top_coords[!is.na(mouse_coord)&!is.na(lemur_coord)&(evec==1),.N,by=paste(site_rank,AS_annot,ever_sig)]
 counts[,site_rank:=strsplit(paste,split=" ")[[1]][1],by=paste]
 counts[,AS_annot:=strsplit(paste,split=" ")[[1]][2],by=paste]
 counts[,ever_sig:=strsplit(paste,split=" ")[[1]][3],by=paste]
 ggplot(counts, aes(x = ever_sig, y = N, fill =AS_annot)) + geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ site_rank) + theme_bw()  +scale_fill_manual(values = c("lightblue","tomato","darkolivegreen","blueviolet"))
# 
# 
 counts = top_coords[!is.na(mouse_coord)&!is.na(lemur_coord)&(evec==2),.N,by=paste(site_rank,AS_annot,ever_sig)]
 counts[,site_rank:=strsplit(paste,split=" ")[[1]][1],by=paste]
 counts[,AS_annot:=strsplit(paste,split=" ")[[1]][2],by=paste]
 counts[,ever_sig:=strsplit(paste,split=" ")[[1]][3],by=paste]
 ggplot(counts, aes(x = ever_sig, y = N, fill =AS_annot)) + geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ site_rank) + theme_bw()  +scale_fill_manual(values = c("lightblue","tomato","darkolivegreen","blueviolet"))
 
 
 counts = top_coords[!is.na(mouse_coord)&!is.na(lemur_coord)&(evec==3),.N,by=paste(site_rank,AS_annot,ever_sig)]
 counts[,site_rank:=strsplit(paste,split=" ")[[1]][1],by=paste]
 counts[,AS_annot:=strsplit(paste,split=" ")[[1]][2],by=paste]
 counts[,ever_sig:=strsplit(paste,split=" ")[[1]][3],by=paste]
 ggplot(counts, aes(x = ever_sig, y = N, fill =AS_annot)) + geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ site_rank) + theme_bw()  +scale_fill_manual(values = c("lightblue","tomato","darkolivegreen","blueviolet"))


write.table(driving_top_coords,output_file,sep="\t",row.names=FALSE,quote=FALSE)
