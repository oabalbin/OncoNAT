# TODO: Add comment
# 
# Author: alebalbin
# This script produces the results presented in Figure 2 and 3 of the 
# Balbin et al The landscape of antisense gene expression in human cancers
# Genome Research manuscript
###############################################################################

################################################################################
### code chunk number 1: SETUP AND SUB-ROUTINES
################################################################################

source(file.path(source_code_dir,'sas_general_functions.R'))
source(file.path(source_code_dir,'functions_tocompute_correlation_between_saspairs.R'))
source(file.path(source_code_dir,'plot_correlation_btw_overlapping_pairs.R'))

results_dir_fig2<-paste(results_dir,'figure2',sep='/')
###############################################################################
### code chunk number 2: Precalculated Data structures to load
###############################################################################
# CORRELATIONS COMPUTED WITHIN COHORTS
# FOR KNOWN OVERLAPPING GENES
CORR_LOCI_K2OV_BC_FILE<-paste(data_dir,"known_overlaps_correlation_table_trans.R",sep="/")
CORR_LOCI_K2OV_MEANS_BC_FILE<-paste(data_dir,
		"known_overlaps_correlation_table_using_meansbycohort_trans.R",sep="/")
# FOR ALL LOCI
CORR_ALLLOCI_BC_FILE<-paste(data_dir,"all_overlaps_correlation_table_trans.R",sep="/")
# RANDOM MATRICES BY COHORT
RANDOM_BYCOHORT_FILE<-paste(data_dir,"random_correlationbycohort_table_trans.R",sep='/')
RANDOM_BYCOHORT_MEANS_FILE<-paste(data_dir,
		"random_overlaps_correlation_table_using_meansbycohort_trans.R",sep='/')

# CORRELATIONS COMPUTED ACCROSS COHORTS
# FOR KNOWN OVERLAPPING GENES
CORR_LOCI_K2OV_AC_FILE<-paste(data_dir,"known_overlaps_correlation_table_accross_cohorts.R",sep='/')
# FOR ALL LOCI
CORR_LOCI_K2OV_AC_FILE<-paste(data_dir,"all_correlations_acohort_trans_table.R",sep='/')
# RANDOM MATRICES BY COHORT
RANDOM_ACOHORT_FILE<-paste(data_dir,"random_correlation_table_trans.R",sep='/')
# TABLE OF CISNAT PAIRS CALLED BY THE OPS AND NASTISEQ METHODS
FULL_UNIQNASTISEQ_LOCI_FILE <-  paste(data_dir,"all_cisNATpairs_OPS_nastiseq_loci.tsv",sep='/')

### Load the Data structures

# Load correlation table for known pairs by cohort
load(file.path(data_dir,CORR_LOCI_K2OV_BC_FILE))
CORR_LOCI_K2OV_BC<-correlation_tables; rm(correlation_tables)
# Load correlation table for random pairs by cohort
load(file.path(data_dir,RANDOM_BYCOHORT_FILE));
RANDOM_BYCOHORT<-NOT_NEIGHBORS_RANDOMCOR_BYCOHORT[[1]]; rm(NOT_NEIGHBORS_RANDOMCOR_BYCOHORT)
# Select only genes in nastiseq loci
corr_table<-CORR_LOCI_K2OV_BC$corr_table


# Load the full list of unique nastiseq loci

full_uniqnastiseq_loci<-as.vector(read.table(FULL_UNIQNASTISEQ_LOCI_FILE,sep='\t',header=FALSE,skip=1)[,2])
uniq_full_uniqnastiseq_loci <- unique(full_uniqnastiseq_loci)
uniq_akto_cohort_corr_table<- corr_table[which(as.vector(corr_table$gene_id) %in% uniq_full_uniqnastiseq_loci),]


###############################################################################
### code chunk number 3: Plots Figures 2 A, B, C
###############################################################################
if(!file.exists(results_dir_fig2)){
	system(paste("mkdir",results_dir_fig2,sep=" "))
}
setwd(results_dir_fig2)

# Plot Figure 2a. Summary correlation plot
plot_correlations_allcohorts(uniq_akto_cohort_corr_table,RANDOM_BYCOHORT,EXONS=FALSE,"correlation_plot_alltrans_bycohort.pdf")
# Plot Figure 2b,c for each pairs type
plot_correlation_btw_pairs_trans(uniq_akto_cohort_corr_table, RANDOM_BYCOHORT)
# Make supplementary table with the statistics and pvalues for the differences in correlations
corr_summary_tables <- summarize_correlation_btw_pairs_exons(uniq_akto_cohort_corr_table, RANDOM_BYCOHORT,EXONS=FALSE)
write.table(corr_summary_tables[[1]],file='correlation_summary_table.tsv',row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
write.table(corr_summary_tables[[2]],file='correlation_pairstype_comparison_table.tsv',row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
plot_correlations_allcohorts(corrtables_ovgenes_using_means,randomcorrtables_ovgenes_using_means[[1]],
		EXONS=FALSE,"correlation_plot_alltrans_usingmeans.pdf")


###############################################################################
### code chunk number 4: Density plots for Figures 2E, 3A, 3E
###############################################################################

### Plot all scatterplots for all pairs

for(i in seq(1,nrow(uniq_akto_cohort_corr_table))){
	
	gene_name=as.vector(ov_trans_corrtable[i,"gene_left_name"]);
	gene_name_as=as.vector(ov_trans_corrtable[i,"gene_right_name"]);
	gene_symbol=as.vector(ov_trans_corrtable[i,"gene_name"]);
	
	print(c(gene,gene_name_as))
	plot_sense_asense_smooth(sense_norm, sense_norm,
			as.data.frame(sense_norm_mean_cohort),as.data.frame(sense_norm_mean_cohort),
			gene_name,cohort_list_names, cohort_name="",scaled=FALSE,log=TRUE,gene_symbol,gene_name_as);
}

###############################################################################
###
###############################################################################

###############################################################################
### code chunk number 5: Supplementary statistics: 
### For the number of pairs that belong to each structural category
###############################################################################
count_sasgenes_by_structuralcategory <- function(ovft){
	intronic<-which(ovft$overlap_type=="INTRONIC")
	EMB<- setdiff(which(ovft$trans_overlap_type=="EMB"),intronic)
	H2H<- setdiff(which(ovft$trans_overlap_type=="H2H"),intronic)
	T2T<- setdiff(which(ovft$trans_overlap_type=="T2T"),intronic)
	# Pairs type
	p2p<-c("protein_coding&protein_coding")
	p2a<-c("protein_coding&antisense","antisense&protein_coding")
	p2l<-c("protein_coding&lincRNA","lincRNA&protein_coding")
	p2ncRNA<-c("protein_coding&processed_transcript","processed_transcript&protein_coding",
			"pseudogene&protein_coding","protein_coding&pseudogene","pseudogene&protein_coding",
			"protein_coding&miRNA","miRNA&protein_coding","protein_coding&sense_intronic","sense_intronic&protein_coding")
	
	a2ncRNA<-c("antisense&antisense","pseudogene&antisense","antisense&pseudogene", "lincRNA&antisense",
			"processed_transcript&antisense","antisense&processed_transcript","antisense&lincRNA",
			"antisense&miRNA","antisense&sense_intronic","3prime_overlapping_ncrna&antisense",
			"antisense&3prime_overlapping_ncrna", "antisense&misc_RNA","antisense&non_coding",                             
			"antisense&polymorphic_pseudogene","antisense&rRNA","antisense&sense_overlapping",
			"antisense&snoRNA","antisense&snRNA")                                
	l2l<-c("lincRNA&lincRNA")
	other <- c(p2p,p2a,p2l,p2ncRNA,a2ncRNA,l2l)
	
	p2p_ind <- which(ovft$gene_biotype %in% p2p)
	p2a_ind <- which(ovft$gene_biotype %in% p2a)
	p2l_ind <- which(ovft$gene_biotype %in% p2l)
	p2ncRNA_ind <- which(ovft$gene_biotype %in% p2ncRNA)
	a2ncRNA_ind <- which(ovft$gene_biotype %in% a2ncRNA)
	l2l_ind <- which(ovft$gene_biotype %in% l2l)
	other_ind <- which(!ovft$gene_biotype %in% other)
	
	m<-as.data.frame(matrix(0,nrow=4,ncol=9))
	rownames(m)<-c("H2H","T2T","EMB","INTRONIC")
	colnames(m)<-c("Number_pairs","median_length_overlap",
			"p2a","p2p","p2l","p2ncRNA","a2ncRNA","l2l","other")
	m["H2H","Number_pairs"]<-length(H2H)
	m["T2T","Number_pairs"]<-length(T2T)
	m["EMB","Number_pairs"]<-length(EMB)
	m["INTRONIC","Number_pairs"]<-length(intronic)
	m["H2H","median_length_overlap"]<-median(ovft$overlap_length[H2H])
	m["T2T","median_length_overlap"]<-median(ovft$overlap_length[T2T])
	m["EMB","median_length_overlap"]<-median(ovft$overlap_length[EMB])
	m["INTRONIC","median_length_overlap"]<-median(ovft$overlap_length[intronic])
	
	# Counting the types of pairs
	m["H2H",c("p2a","p2p","p2l","p2ncRNA","a2ncRNA","l2l","other")] <- c(length(intersect(H2H,p2a_ind)),
			length(intersect(H2H,p2p_ind)),length(intersect(H2H,p2l_ind)),length(intersect(H2H,p2ncRNA_ind)),
			length(intersect(H2H,a2ncRNA_ind)),length(intersect(H2H,l2l_ind)), length(intersect(H2H,other_ind)))
	
	m["T2T",c("p2a","p2p","p2l","p2ncRNA","a2ncRNA","l2l","other")] <- c(length(intersect(T2T,p2a_ind)),
			length(intersect(T2T,p2p_ind)),length(intersect(T2T,p2l_ind)),length(intersect(T2T,p2ncRNA_ind)),
			length(intersect(T2T,a2ncRNA_ind)),length(intersect(T2T,l2l_ind)), length(intersect(T2T,other_ind)))
	
	m["EMB",c("p2a","p2p","p2l","p2ncRNA","a2ncRNA","l2l","other")] <- c(length(intersect(EMB,p2a_ind)),
			length(intersect(EMB,p2p_ind)),length(intersect(EMB,p2l_ind)),length(intersect(EMB,p2ncRNA_ind)),
			length(intersect(EMB,a2ncRNA_ind)),length(intersect(EMB,l2l_ind)), length(intersect(EMB,other_ind)))
	
	m["INTRONIC",c("p2a","p2p","p2l","p2ncRNA","a2ncRNA","l2l","other")] <- c(length(intersect(intronic,p2a_ind)),
			length(intersect(intronic,p2p_ind)),length(intersect(intronic,p2l_ind)),length(intersect(intronic,p2ncRNA_ind)),
			length(intersect(intronic,a2ncRNA_ind)),length(intersect(intronic,l2l_ind)), length(intersect(intronic,other_ind)))
	
	colnames(m)<-c("Number_pairs","median_length_overlap",
			"protein2antisense","protein2protein","protein2lincRNA","prontein2ncRNA","antisense2ncRNA","lincRNA2lincRNA","Other")
	
	return(m)
	
}

uniq_akto_feature_ov<-feature_overlaps_additionals[which(as.vector(feature_overlaps_additionals$gene_id) %in% uniq_full_uniqnastiseq_loci),]
summary_table_overlaps<-count_sasgenes_by_structuralcategory(uniq_akto_feature_ov)

write.table(summary_table_overlaps,file=paste(cohort_name,"summary_table_overlapping_pairs","tsv",sep="."),
		row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)

################################################################################
#### Done
################################################################################






