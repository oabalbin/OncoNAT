# TODO: Add comment
# 
# Author: alebalbin
###############################################################################
# TODO: Add comment
# 
# Author: alebalbin
###############################################################################

################################################################################
### code chunk number 1: SETUP AND SUB-ROUTINES
################################################################################
# Load functions for the simulations for the rebuttal.
source(file.path(source_code_dir,'sas_general_functions.R'))

###############################################################################
### code chunk number 2: Load precomputed data structures
###############################################################################
#TRANSCRIPT_MATRICES_ROBJ object loads the basic data structures for the analyses.
#sense_cds,asense_cds, 
#sense_norm,asense_norm,
#geneAnnot_ori, cohortAnnot,
#cohort_list_names,
#protocolerror=pe,
#cohortLibAnnot,

load(TRANSCRIPT_MATRICES_ROBJ)
feature_overlaps <- read.table(OVERLAP_ANNOTATION_FILE,header=TRUE,sep='\t')


###############################################################################
### code chunk number 3: Plots Figure 1B
###############################################################################
#Smooth plot of sense/antisense expression for protein coding genes
#Not involved in protein&protein coding gene pairs
#geneAnnot <- geneAnnot_ori[reliable_transcripts,]
library(RColorBrewer)
library(KernSmooth)
library(MASS)
geneAnnot <- geneAnnot_ori
cis_pairs <- c("convergent","divergent")
inpairs<-unique(as.vector(as.matrix(feature_overlaps[intersect(which(
												feature_overlaps$gene_biotype=="protein_coding&protein_coding"),
										which(!as.vector(feature_overlaps$overlap_type %in% cis_pairs))),
								c("gene_left","gene_right")])))
pc<-rownames(geneAnnot)[which(geneAnnot$gene_biotype=="protein_coding")]
pcd <-setdiff(pc,inpairs)
spcdm<-as.vector(as.matrix(sense_norm[pcd,]))
apcdm<-as.vector(as.matrix(asense_norm[pcd,]))

# Create density plot for genes with counts != 0 in both strands
z <- kde2d(log10(apcdm[apcdm>0 & spcdm>0 ]+1), log10(spcdm[apcdm>0 & spcdm>0]+1), n=100)

#Plot only the genes that correspond to protein coding genes and do not overlap
#Other protein coding genes
g = 9
my.cols <- rev(brewer.pal(g, "RdYlBu"))
#my.cols <- brewer.pal(g, "OrRd")
mycols<-colorRampPalette(c("white", "red"))(g)


results_dir_fig1 <- paste(results_dir,"figure1",sep="/")
if(!file.exists(results_dir_fig1)){
	system(paste("mkdir",results_dir_fig1,sep=" "))
}
# Change directory to output results
setwd(results_dir_fig1)

pdf("smooth_antisense_vs_sense_overNOToverlappingPrortcoding2.pdf")
percent_outliers=0.0	
smoothScatter(log10(apcdm+1),
		log10(spcdm+1),
		nbin=200,
		nrpoints = percent_outliers*length(spcdm),
		ylab="Sense expression (log10(norm(count)))",xlab="Antisense expresion (log10(norm(count)))",
		cex=0.15,bg="black",pch=21, col="green1"
#colramp=colorRampPalette(mycols)
)
contour(z, drawlabels=FALSE, nlevels=g, col=my.cols, add=TRUE, lwd = 2)
dev.off()


###############################################################################
### code chunk number 4: Find genes with antisense ratio above pe_i_th per
### sample
###############################################################################
###
### Script starts here
#1. Compute the reliable transcripts for each cohort.
#2. Generate figure 1 B, C, D and Supplementary figure
###############################################################################

print("Computing ssRNASeq protocol error rate")
all_genesInpairs<-unique(as.vector(as.matrix(feature_overlaps[,c("gene_left","gene_right")])))
all_genesNOTinpairs<-setdiff(rownames(geneAnnot_ori),all_genesInpairs)
# Restrict to protein coding genes to stablish the negative distribution
negpairsAnnot<-geneAnnot_ori[all_genesNOTinpairs,]
negpairsAnnot<-negpairsAnnot[which(negpairsAnnot$gene_biotype=="protein_coding"),]
negpairs <- unique(rownames(negpairsAnnot))

# Get reliable transcripts
#expressed_loci_persample<-get_reliable_transcripts_matrix(sense_norm,asense_norm,tetha_threshold)

pe_modified_norm_data<-pe_modified_calculation2(sense_norm,asense_norm,negpairs,den_pad=0,min_avgcov=100)

# Count transcripts that are expressed above pe_i_th for each sample
# Counting done per sample
asloci_pei_sd1_countbysample <-compute_genes_with_asratio_above_pei_by_sample(sense_norm,asense_norm,
		cohortAnnot,pe_modified_norm_data,
		bycohort=TRUE,
		min_samples=1,
		min_num_sds=1,
		min_avgcov=10,
		count_by_cohort=FALSE)

# Counting done per cohort
asloci_pei_sd1_countbycohort <-compute_genes_with_asratio_above_pei_by_sample(sense_norm,asense_norm,
		cohortAnnot,pe_modified_norm_data,
		bycohort=TRUE,
		min_samples=1,
		min_num_sds=1,
		min_avgcov=10,
		count_by_cohort=TRUE)

###############################################################################
### code chunk number 5: Plots Figures 1 C, D 
###############################################################################


major_cohorts <- c("BRCA","LUAD","LUSC","LUCL","PRCA","PANC","OVARIAN","MENINGIOMA")
this_rows<-seq(8,nrow(asloci_pei_sd1_countbycohort$summary_bycohort)-1)
barplot_table<-asloci_pei_sd1_countbycohort$summary_bycohort
barplot_table<-asloci_pei_sd1$summary_bycohort
this_rows<-rownames(barplot_table)[this_rows]

### This plot corresponds to Figure 1 C
pdf("barplot_OPSratio_fullcohorts.pdf")
barplot(barplot_table[this_rows,"ALL"]/barplot_table["Total_cohort_numloci_mincov_>=1_sample","ALL"]*100,
		ylab=c("Percentage of loci consistently expressing the opposite strand"),
		xlab=c("Percentage of total samples"),
		beside=TRUE,las=2,
		col="gray",
		border="black")
dev.off()

### This plot corresponds to Figure 1D
pdf("barplot_OPSratio_by_tissue.2.pdf")
barplot(barplot_table[this_rows,major_cohorts]/barplot_table["Total_cohort_numloci_mincov_>=1_sample",
				major_cohorts]*100,
		ylab=c("Percentage of loci consistently expressing the opposite strand"),
		xlab=c("Percentage of total samples"),
		beside=TRUE,las=2,
		col="gray",
		border="white")
dev.off()

###############################################################################
### code chunk number 6: Plots Supplementary Figure 6
###############################################################################
## This generates all panels of Supplementary Figure 6
## Barplots for each tissue type, supplementary table to figure 1D
for(cname in major_cohorts){
	pdf(paste("barplot_OPSratio_by_tissue",cname,"pdf",sep="."))
	barplot(barplot_table[this_rows,cname]/barplot_table["Total_cohort_numloci_mincov_>=1_sample",cname]*100,
			ylab=c("Percentage of loci consistently expressing the opposite strand"),
			xlab=c("Percentage of total samples"),
			beside=TRUE,las=2,
			col="gray",
			border="black")
	dev.off()
}


###############################################################################
### Done
###############################################################################




