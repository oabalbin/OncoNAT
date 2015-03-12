# TODO: Add comment
# 
# Author: alebalbin
# This script produces the tables of the OncoNAT catalogue provided as a resource in
# Balbin et al The landscape of antisense gene expression in human cancers
# Genome Research manuscript
###############################################################################

###############################################################################
### code chunk number 1: Load required R objects

source(paste(source_code_dir,'sas_general_functions.R',sep=''))
source(paste(source_code_dir,'functions_tocompute_correlation_between_saspairs.R',sep=''))

###############################################################################

load(TRANSCRIPT_MATRICES_ROBJ)

CORR_LOCI_K2OV_AC_FILE_TABLE<-paste(data_dir,'known_overlaps_correlation_table_accross_cohorts.corrtable.tsv',sep="/")
ALL_NATGENES_OBJECTS <- paste(data_dir,'all_natgenes_objects.R',sep="/")
ALL_NATGENES_OBJECTS <- paste(data_dir,'all_natgenes_objects.R')
load(ALL_NATGENES_OBJECTS)
all_natgenes<-all_natgenes_objects[[1]]
all_natgenes_ASscore<-all_natgenes_objects[[2]]
#all_natgenes<-all_natgenes[all_natgenes$numcohorts!=0 & all_natgenes$annoterror==0 & all_natgenes$OPSnumcohorts !=0 ,]

load(TRANSCRIPT_MATRICES_ROBJ)
total_correlations_table <- read.table(CORR_LOCI_K2OV_AC_FILE_TABLE,header=TRUE,sep="\t")

results_dir_onconat <- paste(results_dir,'onconat',sep='/')
if(!file.exists(results_dir_onconat)){
	system(paste("mkdir",results_dir_onconat,sep=" "))
}

###############################################################################
### code chunk number 2: SELECT FEATURES THAT MAP TO EITHER TUMOR SUPRESSOR OR ONCOGENES
###############################################################################
tumor_supressors <- read.table(TUMORSUP_FILE,header=TRUE,skip=1,sep='\t')
oncogenes <- read.table(ONOCOGENES_FILE,header=TRUE,skip=1,sep='\t')
# split those gene names in three non-overlapping groups
genes_classified_asboth<-intersect(tumor_supressors$GeneSymbol,oncogenes$GeneSymbol)
tumor_supressors<-setdiff(tumor_supressors$GeneSymbol,genes_classified_asboth)
oncogenes<-setdiff(oncogenes$GeneSymbol,genes_classified_asboth)
oncots<-genes_classified_asboth

###############################################################################
### code chunk number 3:
### Question 1: How many oncogenes, and tumor suppressor have an annotated
### antisense partner.
### Get the full annotation
###############################################################################

feature_overlaps_additionals <- read.table(OVERLAP_ANNOTATION_FILE,header=TRUE,sep='\t')

feature_overlaps_additionals$overlap_type <- feature_overlaps_additionals$exon_overlap_type
feature_overlaps_additionals$overlap_length <- feature_overlaps_additionals$exon_overlaping_length
rmpairs <- which(feature_overlaps_additionals$exon_overlap_type %in% c("convergent","divergent","downstream","upstream"))
newfeature_sasoverlaps<-feature_overlaps_additionals[-rmpairs,]
#
tumor_supressors_ensl<-rownames(geneAnnot_ori)[which(geneAnnot_ori$gene_name %in% tumor_supressors)]
oncogenes_ensl<-rownames(geneAnnot_ori)[which(geneAnnot_ori$gene_name %in% oncogenes)]
oncotusup_ensl<-rownames(geneAnnot_ori)[which(geneAnnot_ori$gene_name %in% oncots)]
all_other_genes<-setdiff(rownames(geneAnnot_ori),unique(c(tumor_supressors_ensl,oncogenes_ensl,oncotusup_ensl)))
all_other_genes_protcod<-intersect(all_other_genes,rownames(geneAnnot_ori)[which(geneAnnot_ori$gene_biotype == "protein_coding")])
actionable<-c(tumor_supressors_ensl,oncogenes_ensl,oncotusup_ensl)

#Determine how many have overlap
tumor_supressors_table<-intersect(tumor_supressors_ensl,unique(as.vector(as.matrix(newfeature_sasoverlaps[,c("gene_left","gene_right")]))))
oncogenes_table<-intersect(oncogenes_ensl,unique(as.vector(as.matrix(newfeature_sasoverlaps[,c("gene_left","gene_right")]))))
oncotusup_table<-intersect(oncotusup_ensl,unique(as.vector(as.matrix(newfeature_sasoverlaps[,c("gene_left","gene_right")]))))
all_other_genes_table<-intersect(all_other_genes,unique(as.vector(as.matrix(newfeature_sasoverlaps[,c("gene_left","gene_right")]))))
all_other_genes_protcod_table<-intersect(all_other_genes_protcod,unique(as.vector(as.matrix(newfeature_sasoverlaps[,c("gene_left","gene_right")]))))

#Actionable genes with not overlap
actionable_table<-unique(c(tumor_supressors_table,oncogenes_table,oncotusup_table))
actionable_notoverlap<-setdiff(actionable,actionable_table)
#geneAnnot_actionable_notoverlap<-geneAnnot_ori[actionable_notoverlap,]

# Actionable genes that intersect with nastiseq genes.
setwd(results_dir_onconat)
# Select loci called by nastiseq and the OPS ratio method
all_natgenes<-all_natgenes[all_natgenes$numcohorts!=0 & all_natgenes$annoterror==0 & all_natgenes$OPSnumcohorts !=0 ,]

actionable_notoverlap_natgene<-intersect(rownames(all_natgenes),actionable_notoverlap)
geneAnnot_actionable_notoverlap<-geneAnnot_ori[actionable_notoverlap_natgene,]
geneAnnot_actionable_notoverlap$ASscore<-all_natgenes[rownames(geneAnnot_actionable_notoverlap),"maxASscore"]
geneAnnot_actionable_notoverlap <- cbind(geneAnnot_actionable_notoverlap,all_natgenes[rownames(geneAnnot_actionable_notoverlap),])

#Determine the not overlapping features that are very close
notoverlapfeature_sasoverlaps<-feature_overlaps_additionals[which(feature_overlaps_additionals$exon_overlap_type %in% c("convergent","divergent")),]
notoverlapfeature_sasoverlaps<-notoverlapfeature_sasoverlaps[notoverlapfeature_sasoverlaps$overlap_length<=500,]
actionable_notoverlap_natgene_closeNeighbor <-intersect(actionable_notoverlap_natgene,
		unique(as.vector(as.matrix(notoverlapfeature_sasoverlaps[,c("gene_left","gene_right")]))))

# Number of actionable genes called by nastiseq and the OPS ratio method

allactionable_genes_nastiseq<-intersect(actionable,rownames(all_natgenes))
# Actionable genes with overlap and nastiseq.
allactionable_genes_nastiseq_table<-intersect(allactionable_genes_nastiseq,unique(as.vector(as.matrix(newfeature_sasoverlaps[,c("gene_left","gene_right")]))))


#Comparing Tumor Suppressors and oncogenes to the rest of protein coding genes
a<-length(all_other_genes_protcod_table)
b<-length(all_other_genes_protcod)-a
c<-length(tumor_supressors_table)
d<-length(tumor_supressors_ensl)-c
fishertable <-matrix(c(a,c,b,d),ncol=2,byrow=T)
rownames(fishertable)<-c("OV","NOTOV")
colnames(fishertable)<-c("ALL","ACT")
prop.test(fishertable)
fisher.test(fishertable,alternative="two.sided")
# p-value = 0.005272


#Comparing oncogenes  to the rest of protein coding genes
a<-length(all_other_genes_protcod_table)
b<-length(all_other_genes_protcod)-a
c<-length(oncogenes_table)
d<-length(oncogenes_ensl)-c
fishertable <-matrix(c(a,c,b,d),ncol=2,byrow=T)
rownames(fishertable)<-c("OV","NOTOV")
colnames(fishertable)<-c("ALL","ACT")
prop.test(fishertable)
fisher.test(fishertable,alternative="two.sided")

#Comparing tumor sup and oncogenes 
a<-length(tumor_supressors_table)
b<-length(tumor_supressors_ensl)-a
c<-length(oncogenes_table)
d<-length(oncogenes_ensl)-c
fishertable <-matrix(c(a,c,b,d),ncol=2,byrow=T)
rownames(fishertable)<-c("OV","NOTOV")
colnames(fishertable)<-c("ALL","ACT")
prop.test(fishertable)
fisher.test(fishertable,alternative="two.sided")

#Comparing genes considered dual  to the rest of protein coding genes
a<-length(all_other_genes_protcod_table)
b<-length(all_other_genes_protcod)-a
c<-length(oncotusup_table)
d<-length(oncotusup_ensl)-c
fishertable <-matrix(c(a,c,b,d),ncol=2,byrow=T)
rownames(fishertable)<-c("OV","NOTOV")
colnames(fishertable)<-c("ALL","ACT")
prop.test(fishertable)
fisher.test(fishertable,alternative="two.sided")

###############################################################################
### code chunk number 4:
### Question 2: is there any difference when looking at
### the correlation between the sense and antisense genes
### for tumor suppressors. 
### A: There were not correlations computed for overlapping pairs of tumor
### suppressor genes. Well there are for couple of them in the 
### category of dual genes. 
### The main reason is that antisense expression in tumor suppressor genes is overall
### lower than what was observed for oncogenes. mean=30 normalized counts, vs 168 respectively.
###############################################################################

#Now looking at the correlations.
setwd(results_dir_onconat)
total_correlations_table <- read.table(CORR_LOCI_K2OV_AC_FILE_TABLE,header=TRUE,sep="\t")
total_correlations_table_only_sas<-total_correlations_table[which(!as.vector(total_correlations_table$overlap_type) %in%
						c("convergent","divergent","downstream","upstream")),]

chrfun<-function(x){return(paste("chr",x[1],":",as.numeric(x[2]),"-",as.numeric(x[3]),sep=""))}
namefun<-function(x){return(paste(x[1],x[2],x[3],sep="_"))}
total_correlations_table_only_sas$chr_reg<-as.vector(apply(total_correlations_table_only_sas[,c("chr","start","end")],1,chrfun))
total_correlations_table_only_sas$cisNATid<-as.vector(apply(total_correlations_table_only_sas[,c("gene_left","gene_right","overlap_type")],1,namefun))
# From the cisnat script I got the all_natgenes_ASscore table
corrtable_sas_ASscore <- as.matrix(cbind(all_natgenes[as.vector(total_correlations_table_only_sas$gene_left),"maxASscore"],
				all_natgenes[as.vector(total_correlations_table_only_sas$gene_right),"maxASscore"]))
#
total_correlations_table_only_sas$ASscore<-as.vector(apply(corrtable_sas_ASscore,1,max,na.rm=TRUE))
# From the cisnat script I got the all_natgenes table
natind<-union(which(as.vector(total_correlations_table_only_sas$gene_left) %in% rownames(all_natgenes)),
				which(as.vector(total_correlations_table_only_sas$gene_right) %in% rownames(all_natgenes)))

# Label gene pairs in which one at least one gene was identified as asloci by nastiseq
total_correlations_table_only_sas$natgene<-0
total_correlations_table_only_sas$natgene[natind]<-1

# Label gene pairs in which  at least one gene was identified as a cancer speficic gene
all_natgenes_notbenign<-all_natgenes[all_natgenes$benign==0,]
canind<-union(which(as.vector(total_correlations_table_only_sas$gene_left) %in% rownames(all_natgenes_notbenign)),
		which(as.vector(total_correlations_table_only_sas$gene_right) %in% rownames(all_natgenes_notbenign)))

# Label gene pairs in which  at least one gene was identified as a cancer speficic gene
total_correlations_table_only_sas$natgene_canspecific<-0
total_correlations_table_only_sas$natgene_canspecific[canind]<-1

#from the cpg islands script I got the cpgi_mapped_table
###
cpgi_mapped_table <- read.table(CPGI_MAPPED_TABLE_FILE,header=TRUE,sep="\t");

cpgind<-intersect(rownames(cpgi_mapped_table),rownames(total_correlations_table_only_sas))
total_correlations_table_only_sas$cpgisland<-0
total_correlations_table_only_sas[cpgind,"cpgisland"]<-1

#Extract tumor supressors and oncogenes
corr_th=0.0
actgenes_ov<-get_oncogene_tumorsup_annot_corrtable2(total_correlations_table_only_sas,oncogenes_table, tumor_supressors_table,corr_th, absv=TRUE)
actgenes_ov<-featureov_eliminate_redundancy(actgenes_ov,"trans_overlap_type")
o<-order(actgenes_ov$natgene,actgenes_ov$srho,actgenes_ov$cpgisland,actgenes_ov$trans_overlap_type,
		actgenes_ov$transcript_biotype,actgenes_ov$ASscore,decreasing=TRUE)
actgenes_ov<-actgenes_ov[o,]
#Select oncogenes and tumor sup
onco_corrtable<-actgenes_ov[actgenes_ov$oncogene==1,]
tsup_corrtable<-actgenes_ov[actgenes_ov$tumsup==1,]
#
#poscorr_th=0.3; negcorr_th=0.1
corr_th=0.0
actgenes_ov<-get_oncogene_tumorsup_annot_corrtable2(total_correlations_table_only_sas,oncogenes_table, oncotusup_table,corr_th, absv=TRUE)
actgenes_ov<-featureov_eliminate_redundancy(actgenes_ov,"trans_overlap_type")
o<-order(actgenes_ov$natgene,actgenes_ov$srho,actgenes_ov$cpgisland,actgenes_ov$trans_overlap_type,
		actgenes_ov$transcript_biotype,actgenes_ov$ASscore,decreasing=TRUE)
actgenes_ov<-actgenes_ov[o,]
#onco_corrtable<-actgenes_ov[actgenes_ov$oncogene==1,]
#oncotusup_corrtable<-actgenes_ov[actgenes_ov$tumsup==1,]
oncotusup_corrtable <-actgenes_ov

# All actioable genes with correlation and annotated.
allactionable_corrtable<-rbind(onco_corrtable,tsup_corrtable,oncotusup_corrtable)
allactionable_corrtable<-featureov_eliminate_redundancy(allactionable_corrtable,"trans_overlap_type")
allactionable_corrtable_nastiseq<-allactionable_corrtable[which(allactionable_corrtable$natgene==1),]
#allactionable_genes_nastiseq_table2<-intersect(allactionable_genes_nastiseq,unique(as.vector(as.matrix(allactionable_corrtable_nastiseq[,c("gene_left","gene_right")]))))

#collect pair for the ones that have a neighbor annotated
actgenes_ov<-get_oncogene_tumorsup_annot_corrtable2(total_correlations_table,actionable_notoverlap_natgene_closeNeighbor, 
		actionable_notoverlap_natgene_closeNeighbor,corr_th, absv=TRUE)
actgenes_ov<-featureov_eliminate_redundancy(actgenes_ov,"trans_overlap_type")
#o<-order(actgenes_ov$natgene,actgenes_ov$srho,actgenes_ov$cpgisland,actgenes_ov$trans_overlap_type,
#		actgenes_ov$transcript_biotype,actgenes_ov$ASscore,decreasing=TRUE)
#actgenes_ov<-actgenes_ov[o,]
#onco_corrtable<-actgenes_ov[actgenes_ov$oncogene==1,]
#notoverlap_corrtable<-actgenes_ov[intersect(which(actgenes_ov$tumsup==1),which(actgenes_ov$overlaping_length<=500)),]
notoverlap_corrtable<-actgenes_ov[actgenes_ov$exon_overlaping_length<=500,]
notoverlap_corrtable_farneighbors<-actgenes_ov[actgenes_ov$exon_overlaping_length>500,]

notoverlap_corrtable$cisNATid<-as.vector(apply(notoverlap_corrtable[,c("gene_left","gene_right","overlap_type")],1,namefun))
notoverlap_corrtable_farneighbors$cisNATid<-as.vector(apply(notoverlap_corrtable_farneighbors[,c("gene_left","gene_right","overlap_type")],1,namefun))

###############################################################################
### code chunk number 5:
### Write the tables for the OncoNAT catalogue
###############################################################################


cols2print<-c("cisNATid","chr_reg","transcript_biotype","gene_name","exon_overlaping_length","trans_overlap_type",
		"overlap_type","srho","spvalue_adjusted","spvalue_fdr","cohort","ASscore","natgene","cpgisland","oncogene","tumsup","natgene_canspecific")

cols2print2<-c("cisNATid","chr","start","end","transcript_biotype","gene_name","exon_overlaping_length","trans_overlap_type",
		"overlap_type","srho","spvalue_adjusted","spvalue_fdr","cohort","oncogene","tumsup")

### Main OncoNAT DB tables
### EMB H2H NOV T2T
o<-order(allactionable_corrtable_nastiseq$trans_overlap_type,
		allactionable_corrtable_nastiseq$srho,decreasing=TRUE)
allactionable_corrtable_nastiseq<-allactionable_corrtable_nastiseq[o,]
write.table(allactionable_corrtable_nastiseq[allactionable_corrtable_nastiseq$trans_overlap_type=="H2H",
				cols2print],file="DataS11_OncoNAT.HTH.tsv",
		row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
write.table(allactionable_corrtable_nastiseq[allactionable_corrtable_nastiseq$trans_overlap_type=="T2T",
				cols2print],file="DataS12_OncoNAT.TTT.tsv",
		row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
write.table(allactionable_corrtable_nastiseq[allactionable_corrtable_nastiseq$trans_overlap_type=="EMB",
				cols2print],file="DataS13_OncoNAT.EMB.tsv",
		row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
### Misannotated and genes with not overlapping transcripts in Emsembl v69.
write.table(notoverlap_corrtable[,cols2print2],file="DataS14_OncoNAT.MISANNOTATED.tsv",
		row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)

write.table(geneAnnot_actionable_notoverlap,file="DataS15_OncoNAT.NotOverlap.tsv",
		row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)

###############################################################################
### Done
###############################################################################
