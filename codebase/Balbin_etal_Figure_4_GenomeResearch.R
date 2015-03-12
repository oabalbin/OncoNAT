# TODO: Add comment
# 
# Author: alebalbin
# This script produces the results presented in Figure 4 of the 
# Balbin et al The landscape of antisense gene expression in human cancers
# Genome Research manuscript
###############################################################################
### code chunk number 1: SETUP AND SUB-ROUTINES
###############################################################################
source(paste(source_code_dir,'sas_general_functions.R',sep=''))

results_dir_fig4 <- paste(results_dir,"figure4",sep="/")
if(!file.exists(results_dir_fig4)){
	system(paste("mkdir",results_dir_fig4,sep=" "))
}
###############################################################################
### code chunk number 2: Load sense_norm and asense_norm data structures
###############################################################################
load(TRANSCRIPT_MATRICES_ROBJ)

###############################################################################
### code chunk number 3: Load annotation of tumor suppresors and oncogenes
### Annotate transcripts with cancer_related_information
###############################################################################


tumor_supressors <- read.table(TUMORSUP_FILE,header=TRUE,skip=1,sep='\t')
oncogenes <- read.table(ONOCOGENES_FILE,header=TRUE,skip=1,sep='\t')

cancer_related_genesymbol <- unique(union(oncogenes$GeneSymbol,tumor_supressors$GeneSymbol))
cancer_related_genes<-rownames(geneAnnot_ori)[which(geneAnnot_ori$gene_name %in% cancer_related_genesymbol)]

###############################################################################
### code chunk number 3: Find genes with OPS ratio greater than pe_i_th per sample
###############################################################################

load(CISNAT_TRAINING_DATASET_FILE)
tetha_threshold=10 #0
pe_modified_norm_data<-pe_modified_calculation2(sense_norm,asense_norm,negpairs,den_pad=0,min_avgcov=100)

asloci_pei_sd1 <-compute_genes_with_asratio_above_pei_by_sample(sense_norm,asense_norm,
		cohortAnnot,pe_modified_norm_data,
		bycohort=TRUE,
		min_samples=1,
		min_num_sds=1,
		min_avgcov=tetha_threshold,
		count_by_cohort=TRUE)

sample_thr=0.3
OPS_genes_alltissues<-which(rowSums(asloci_pei_sd1$asloci_bysample) >= sample_thr*ncol(sense_norm))

###############################################################################
### code chunk number 4: Read in precomputed NASTI-seq and OPS results
### for each cohort and determine ubiquitous, tissue enriched and tissue 
### specific antisense loci
###############################################################################

cohort_names<-c("luad","lusc","lucl","benign","brca","prca","panc","meningioma") #"ovarian",
cancer_names<-c("luad","lusc","lucl","brca","prca","panc","meningioma")
cnastiseq_list<-c()
common_natgenes<-c()
all_natgenes<-matrix(0,nrow=nrow(sense_norm),ncol=length(cohort_names))
rownames(all_natgenes)<-rownames(sense_norm); colnames(all_natgenes)<-cohort_names;
all_POSgenes<-matrix(0,nrow=nrow(sense_norm),ncol=length(cohort_names))
rownames(all_POSgenes)<-rownames(sense_norm); colnames(all_POSgenes)<-cohort_names;
all_POSgenes_20percent<-matrix(0,nrow=nrow(sense_norm),ncol=length(cohort_names))
rownames(all_POSgenes_20percent)<-rownames(sense_norm); colnames(all_POSgenes_20percent)<-cohort_names;


annoterror<-matrix(0,nrow=nrow(sense_norm),ncol=1);
rownames(annoterror)<-rownames(sense_norm); colnames(annoterror)<-c("ANNOT_ERROR");

all_natgenes_ASscore<-matrix(0,nrow=nrow(sense_norm),ncol=length(cohort_names))
rownames(all_natgenes_ASscore)<-rownames(sense_norm); colnames(all_natgenes_ASscore)<-cohort_names;
first=TRUE

### Load the precomputed reults for the OPS and NASTI-seq calculations
file_basename_final<-'NASTIseq_fullresults_trans_mincount5_rebuttal'

sample_20thr=0.20
		
for(cname in cohort_names){
	cohort<-cohort_list_names[[cname]]
	
	filename_tmp<-paste(file_basename_final,cname,"R",sep=".")
	setwd(data_dir_OPSnastiseq)
	# Loads the nastiseq results
	load(filename_tmp)
	
	print(cname)
	geneAnnot_natcand<-nastiseq_results$geneAnnot_natcand
	
	cg<-rownames(geneAnnot_natcand)[which(geneAnnot_natcand$OPS_common==1)]
	#cg<-rownames(geneAnnot_natcand)
	ops_20per<-which(rowSums(asloci_pei_sd1$asloci_bysample[,cohort]) >= sample_20thr*length(cohort))
	
	all_natgenes[cg,cname]<-1
	all_natgenes_ASscore[cg,cname]<-nastiseq_results$ASscore[cg]
	annoterror[cg,"ANNOT_ERROR"]<- (annoterror[cg,"ANNOT_ERROR"] + nastiseq_results$geneAnnot_natcand[cg,"ANNOT_ERROR"])
	all_POSgenes[cg,cname]<-geneAnnot_natcand[cg,"OPS_common"]
	
	all_POSgenes_20percent[ops_20per,cname]<-1
	
	print("Annotation error...")
	print(length(which(nastiseq_results$geneAnnot_natcand[cg,"ANNOT_ERROR"]==1)))
	
	if(first){
		common_natgenes <- cg
		first=FALSE
	}else{
		common_natgenes <- intersect(common_natgenes,cg)
	}	
}

###############################################################################
### code chunk number 5: Prepare the structure for plotting heatmap in figure 4
###############################################################################
#
geneAnnot_ori$cancer_related<-0;
geneAnnot_ori[cancer_related_genes,"cancer_related"]<-1
geneAnnot_ori$maxASscore <- NA;
geneAnnot_ori[rownames(all_natgenes_ASscore),"maxASscore"] <-apply(all_natgenes_ASscore,1, max,na.rm=TRUE)
#
common_natgenes<-unique(common_natgenes)
all_natgenes<-as.data.frame(all_natgenes);
all_natgenes_ASscore<-as.data.frame(all_natgenes_ASscore)
all_natgenes$numcohorts<-apply(all_natgenes[,cancer_names],1,sum)
all_natgenes$annoterror<-as.vector(annoterror)
#
all_POSgenes<-as.data.frame(all_POSgenes)
all_POSgenes$OPSnumcohorts<-apply(all_POSgenes[,cancer_names],1,sum)
all_POSgenes_20percent<-as.data.frame(all_POSgenes_20percent)
all_POSgenes_20percent$OPSnumcohorts<-apply(all_POSgenes_20percent[,cancer_names],1,sum)
all_POSgenes$OPS_20percent<-all_POSgenes_20percent$OPSnumcohorts
# Join matrices
all_natgenes <- cbind(all_natgenes,all_POSgenes,geneAnnot_ori)
#
o<-order(all_natgenes$numcohorts,
		all_natgenes$benign,
		all_natgenes$luad,
		all_natgenes$lusc,
		all_natgenes$brca,
		all_natgenes$prca,
		all_natgenes$panc,
		#all_natgenes$ovarian,
		all_natgenes$meningioma,
		all_natgenes$lucl,
		decreasing=TRUE);
all_natgenes<-all_natgenes[o,];
all_natgenes_ori <- all_natgenes

###############################################################################
### code chunk number 6: Save objects to an R structur for future load
###############################################################################
setwd(data_dir)
save(list(all_natgenes,all_natgenes_ASscore,all_POSgenes),
		file=paste(data_dir,'all_natgenes_objects.R'))

###############################################################################
### code chunk number 7: Select antisense loci present in at least one cohort
###############################################################################

all_natgenes<-all_natgenes[all_natgenes$numcohorts!=0 & all_natgenes$annoterror==0 & all_natgenes$OPSnumcohorts !=0 ,]

###############################################################################
### code chunk number 8: Find ubiquitous, tissue enriched and tissue specific 
### antisense loci
###############################################################################
### Total number of cancer tissue cohorts
totalnumcohorts=length(cancer_names)

# Tissue specific antisense loci
all_natgenes_unique<-all_natgenes[all_natgenes$benign<=1 &
				all_natgenes$numcohorts==1 & all_natgenes$OPS_20percent==1,]

# Ubiquitous
ubiquitous<-intersect(which(all_natgenes$benign<=1),
		which(all_natgenes$numcohorts==totalnumcohorts))
potential_ubiquitous<-all_natgenes[all_natgenes$benign<=1 &
				all_natgenes$numcohorts==1 & all_natgenes$OPS_20percent==totalnumcohorts,]

# Select ubiquitous loci
all_natgenes_common<-rbind(all_natgenes[ubiquitous,],potential_ubiquitous)

# Determine tissue enriched loci.
tissue_enriched <- setdiff(intersect(which(all_natgenes$benign<=1),
				which(all_natgenes$numcohorts>=2)),ubiquitous)
# Add loci that look like tissue specific but they are found to be antisense loci
# in at least 20% of cohort samples but less than 30% and therefore appear as
# expressed in only one tissue type. They are erroneously included in the tissue specific set.

potential_enriched<-all_natgenes[all_natgenes$benign<=1 &
				all_natgenes$numcohorts==1 & all_natgenes$OPS_20percent > 1,]

potential_enriched<-potential_enriched[setdiff(rownames(potential_enriched),
				rownames(potential_ubiquitous)),]
# Select all tissue enriched loci 
all_natgenes_tissue_enriched<-rbind(all_natgenes[tissue_enriched,],potential_enriched)
all_natgenes_tissue_enriched<-all_natgenes_tissue_enriched[order(all_natgenes_tissue_enriched$numcohorts,decreasing=FALSE),]

#
all_natgenes_protcods<-all_natgenes[all_natgenes$gene_biotype=="protein_coding",]

###############################################################################
### code chunk number 9: Aggregate the number of antisense loci in each category
###############################################################################

heatmap_genes_unique<-intersect(rownames(all_natgenes_unique),rownames(all_natgenes_protcods))
heatmap_genes_tissueenriched<-intersect(rownames(all_natgenes_tissue_enriched),rownames(all_natgenes_protcods))
heatmap_genes_ubiquitous<-intersect(rownames(all_natgenes_common),rownames(all_natgenes_protcods))
#
heatmap_genes<-c(heatmap_genes_unique,heatmap_genes_tissueenriched,heatmap_genes_ubiquitous)

figure4_summarytable<-matrix(NA,nrow=3,ncol=3)
rownames(figure4_summarytable)<-c("tissue_specific","tissue_enriched","tissue_ubiquitous")
colnames(figure4_summarytable)<-c("total_number_loci","protcods_loci","cancer_related")
figure4_summarytable["tissue_specific","total_number_loci"]<-nrow(all_natgenes_unique);
figure4_summarytable["tissue_specific","protcods_loci"]<-length(heatmap_genes_unique);
figure4_summarytable["tissue_specific","cancer_related"]<-nrow(all_natgenes_unique[all_natgenes_unique$cancer_related==1,]);
figure4_summarytable["tissue_enriched","total_number_loci"]<-nrow(all_natgenes_tissue_enriched);
figure4_summarytable["tissue_enriched","protcods_loci"]<-length(heatmap_genes_tissueenriched);
figure4_summarytable["tissue_enriched","cancer_related"]<-nrow(all_natgenes_tissue_enriched[all_natgenes_tissue_enriched$cancer_related==1,]);
figure4_summarytable["tissue_ubiquitous","total_number_loci"]<-nrow(all_natgenes_common);
figure4_summarytable["tissue_ubiquitous","protcods_loci"]<-length(heatmap_genes_ubiquitous);
figure4_summarytable["tissue_ubiquitous","cancer_related"]<-nrow(all_natgenes_common[all_natgenes_common$cancer_related==1,]);

#total_number_loci protcods_loci cancer_related
#tissue_specific                 894           653             44
#tissue_enriched                7942          4798            282
#tissue_ubiquitous              4828          3273            178

###############################################################################
### code chunk number 10: Plot Figure 4
### Heat Maps for tissue_specific, tissue_enriched and tissue_ubiquitous 
### antisense loci
###############################################################################
setwd(results_dir_fig4)
library("gplots")
colsorder<-c("benign","luad","lusc","brca","prca","panc","meningioma","lucl")
colsorder<-c("luad","lusc","brca","prca","panc","meningioma","lucl")

heatmap_genes<-c(heatmap_genes_unique,heatmap_genes_tissueenriched,heatmap_genes_ubiquitous)
asnorm_meanbycohort<-compute_mean_by_cohort(asense_norm,cohort_list_names)

heatmap_matrix_asense_tissue_specific<-asnorm_meanbycohort[heatmap_genes_unique,toupper(colsorder)]
heatmap_matrix_asense_tissue_enriched<-asnorm_meanbycohort[heatmap_genes_tissueenriched,toupper(colsorder)]
heatmap_matrix_asense_ubiquitous<-asnorm_meanbycohort[heatmap_genes_ubiquitous,toupper(colsorder)]

pdf('figure4_tissue_specific_protcodsloci.pdf')
gn<-character(nrow(all_natgenes_unique))
cr<-which(all_natgenes_unique$cancer_related==1)
gn[cr]<-as.vector(all_natgenes_unique$gene_name[cr])

heatmap.2(as.matrix(heatmap_matrix_asense_tissue_specific),
		Rowv=NA,Colv=NA,
		col = rev(heat.colors(256)),#length(palette.breaks)-1)), #256
		#breaks=palette.breaks,
		scale="row",
		dendrogram="none",
		key=TRUE,
		labRow=gn,#geneAnnot_ori[heatmap_genes_unique,"gene_name"],
		symkey=FALSE, density.info="none", trace="none", cexRow=0.2) #,cexCol=0.1
dev.off()

pdf('figure4_tissue_enriched_protcodsloci.pdf')
gn<-character(nrow(all_natgenes_tissue_enriched))
cr<-which(all_natgenes_tissue_enriched$cancer_related==1 &
						all_natgenes_tissue_enriched$maxASscore>1000)
gn[cr]<-as.vector(all_natgenes_tissue_enriched$gene_name[cr])

heatmap.2(as.matrix(heatmap_matrix_asense_tissue_enriched),
		Rowv=NA,Colv=NA,
		col = rev(heat.colors(256)),#length(palette.breaks)-1)), #256
		#breaks=palette.breaks,
		scale="row",
		dendrogram="none",
		key=TRUE,
		labRow=gn,#geneAnnot_ori[heatmap_genes_tissueenriched,"gene_name"],
		symkey=FALSE, density.info="none", trace="none", cexRow=0.2) #,cexCol=0.1
dev.off()

pdf('figure4_tissue_ubiquitous_protcodsloci.pdf')
gn<-character(nrow(all_natgenes_common))
cr<-which(all_natgenes_common$cancer_related==1 &
				all_natgenes_common$maxASscore>1000 )
gn[cr]<-as.vector(all_natgenes_common$gene_name[cr])
heatmap.2(as.matrix(heatmap_matrix_asense_ubiquitous),
		Rowv=NA,Colv=NA,
		col = rev(heat.colors(256)),#length(palette.breaks)-1)), #256
		#breaks=palette.breaks,
		scale="row",
		dendrogram="none",
		key=TRUE,
		labRow=gn,#geneAnnot_ori[heatmap_genes_ubiquitous,"gene_name"],
		symkey=FALSE, density.info="none", trace="none", cexRow=0.2) #,cexCol=0.1
dev.off()

###############################################################################
### code chunk number 11: Plot full heat map of antisense expression
### dendogram on the tissues and the samples.
### This plot is not included in the manuscript figures.
###############################################################################


# HeatMaps:
patientcolors <- rep("white",ncol(asense_norm_heatmap))
patientcolors[which(cohortAnnot$tissuetype=="BRCA")]<-"yellow"
patientcolors[which(cohortAnnot$tissuetype=="LUAD")]<-"blue"
patientcolors[which(cohortAnnot$tissuetype=="LUSC")]<-"red"
patientcolors[which(cohortAnnot$tissuetype=="LUCL")]<-"pink"
patientcolors[which(cohortAnnot$tissuetype=="PRCA")]<-"cyan"
patientcolors[which(cohortAnnot$tissuetype=="PANC")]<-"green"
patientcolors[which(cohortAnnot$tissuetype=="OVARIAN")]<-"purple"
patientcolors[which(cohortAnnot$tissuetype=="MENINGIOMA")]<-"orange"

asense_norm_heatmap<-asense_norm[heatmap_genes,]
#asense_norm_heatmap<-asnorm_meanbycohort[heatmap_genes,]
x<-log2(asense_norm_heatmap+1)
xn<-t(scale(t(x),center=TRUE, scale=TRUE))
finite_values=as.vector(which(is.finite(rowMeans(xn))))
xn <- xn[finite_values,]
major_tissue_cohorts <- c("BRCA","LUAD","LUSC","LUCL","PRCA","PANC","MENINGIOMA","OV")
cols<-which(cohortAnnot$tissuetype %in% major_tissue_cohorts)


pdf('asense_norm_tissue_heatmap.pdf')
palette.breaks <- seq(-4, 4, 0.01)
heatmap.2(as.matrix(xn[,cols]),
		Rowv=TRUE,Colv=TRUE,
		#col = rev(heat.colors(256)),#length(palette.breaks)-1)), #256
		col=greenred(length(palette.breaks)-1), 
		breaks=palette.breaks,
		scale="row",
		dendrogram="both",
		key=TRUE,
		ColSideColors=patientcolors[cols],
		labRow=geneAnnot_ori[heatmap_genes,"gene_name"],
		symkey=FALSE, density.info="none", trace="none", cexRow=0.8) #,cexCol=0.1
dev.off()


###############################################################################
### code chunk number 12: Write the tables to a file
###############################################################################

all_heatmap_matrix_asense <- as.data.frame(asnorm_meanbycohort[heatmap_genes,toupper(colsorder)])
all_heatmap_matrix_asense$meanASscore<-rowMeans(all_natgenes_ASscore[heatmap_genes,],na.rm=TRUE)
all_heatmap_matrix_asense <- cbind(all_heatmap_matrix_asense,geneAnnot_ori[heatmap_genes,])
#
all_heatmap_matrix_asense$tissue_specific<-0; 
all_heatmap_matrix_asense[heatmap_genes_unique,"tissue_specific"]<-1
all_heatmap_matrix_asense$tissue_enriched<-0; 
all_heatmap_matrix_asense[heatmap_genes_tissueenriched,"tissue_enriched"]<-1
all_heatmap_matrix_asense$tissue_ubiquitous<-0; 
all_heatmap_matrix_asense[heatmap_genes_ubiquitous,"tissue_ubiquitous"]<-1

all_heatmap_matrix_asense<-cbind(all_heatmap_matrix_asense,all_natgenes[rownames(all_heatmap_matrix_asense),])


o<-order(all_heatmap_matrix_asense$tissue_specific,
		all_heatmap_matrix_asense$tissue_enriched,
		all_heatmap_matrix_asense$tissue_ubiquitous,
		all_heatmap_matrix_asense$luad,
		all_heatmap_matrix_asense$lusc,
		all_heatmap_matrix_asense$brca,
		all_heatmap_matrix_asense$prca,
		all_heatmap_matrix_asense$panc,
		all_heatmap_matrix_asense$meningioma,
		all_heatmap_matrix_asense$lucl,
		all_heatmap_matrix_asense$maxASscore,
		all_heatmap_matrix_asense$cancer_related,
		decreasing=TRUE);

all_heatmap_matrix_asense<-all_heatmap_matrix_asense[o,]

write.table(all_heatmap_matrix_asense,file="all_heatmap_matrix_asense_figure4.tsv",
		row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)

###############################################################################
### Done
###############################################################################




