# TODO: Add comment
# 
# Author: alebalbin
# This script produces the results presented in Figure 6 of the 
# Balbin et al The landscape of antisense gene expression in human cancers
# Genome Research manuscript

###############################################################################
### code chunk number 1: Specific functions of this script 
###############################################################################

calculate_diffexp_genes_DESeq <-function(smat,condition,controlgroup_name, querygroup_name,smallest_group,strand){
	library(DESeq)
	min_num_samples=7 # recommended from the DESeq viggnette
	sense_cds<-newCountDataSet(smat,condition)
	###################################################
	### code chunk number : estimateSizeFactors
	### and dispersion
	###################################################
	print("Starting DESeq size factor and dispersion estimation")
	sense_cds = estimateSizeFactors( sense_cds )
	if((length(condition) < min_num_samples) || (length(smallest_group) < min_num_samples)){
		mysharingMode="maximum"
	}else{mysharingMode="gene-est-only"}
	# calculate dispersion
	print(paste("Sharing mode, ",mysharingMode,sep=" "))
	sense_cds = estimateDispersions( sense_cds,method="per-condition",sharingMode=mysharingMode,fitType="local")
	# compute differential expression
	print("Computing DESeq nbinomTest calculation")
	res = nbinomTest( sense_cds, controlgroup_name, querygroup_name )
	pdf("DESeq_MAplot.pdf")
	plotMA(res)
	dev.off()
	pdf("DESeq_histogram.pdf")
	hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
	dev.off()
	save(res,file=paste("nBimTestRes",controlgroup_name, querygroup_name,strand,"R",sep="."))
	return(res)
}

select_diffgenes2<-function(res,pvalth){
	resSig = res[ res$padj < pvalth, ]
	return(resSig[ order( -resSig$foldChange, -resSig$baseMean ), ])
}

#evaluate_assignificant_pair <- function(x,y,pvalth){
#	if( x$padj < pvalth & y$padj < pvalth){}
#}

select_diffgenes<-function(res,gene_pairs,pvalth){
	glsig<-c();grsig<-c()
	rownames(res)<-res$id
	gl_names<-as.vector(gene_pairs[,1]); gr_names<-as.vector(gene_pairs[,2]);
	
	for(i in seq(1,nrow(gene_pairs))){
		gl<-gl_names[i];gr<-gr_names[i];
		if( (gl %in% rownames(res) & gr %in% rownames(res)) & (!is.na(res[gl,"padj"]) & !is.na(res[gr,"padj"]) ) 
				& (res[gl,"padj"] < pvalth & res[gr,"padj"] < pvalth)){
			glsig<-c(glsig,gl);grsig<-c(grsig,gr)
		}
	}
	res_grsig<-res[grsig,]; colnames(res_grsig) <- paste(colnames(res_grsig),"1",sep=".")
	resSig<-cbind(res[glsig,],res_grsig)
	save(res,file=paste("resSignificantPairs",controlgroup_name, querygroup_name,"R",sep="."))
	return(resSig);
}

select_diffgenes_twomat<-function(res,ares,gl_names,pvalth){
	glsig<-c();grsig<-c()
	rownames(res)<-res$id
	#gl_names<-as.vector(gene_pairs[,1]); gr_names<-as.vector(gene_pairs[,2]);	
	withpval<-intersect(which(!is.na(res[,"padj"])),which(!is.na(ares[,"padj"])))
	withpval_th<-intersect(which(res[,"padj"] < pvalth),which(ares[,"padj"] < pvalth))
	rsig<-intersect(withpval,withpval_th)
	
	res_grsig<-ares[rsig,]; colnames(res_grsig) <- paste(colnames(res_grsig),"1",sep=".")
	resSig<-cbind(res[rsig,],res_grsig)
	save(res,file=paste("resSignificantPairs.bothstrands",controlgroup_name, querygroup_name,"R",sep="."))
	return(resSig);
}



select_diffgenes_bylfc<-function(res,lfcth){
	resSig1 = res[ intersect(which(res$log2FoldChange >= lfcth),which(res$log2FoldChange.1 >= lfcth)),]
	resSig2 = res[ intersect(which(res$log2FoldChange <= -1*lfcth), which(res$log2FoldChange.1 <= -1*lfcth)),]
	resSig3 = res[ intersect(which(res$log2FoldChange >= lfcth),which(res$log2FoldChange.1 <= -1*lfcth)),]
	resSig4 = res[ intersect(which(res$log2FoldChange <= -1*lfcth), which(res$log2FoldChange.1 >= lfcth)),]
	
	return( list(consistent=rbind(resSig1,resSig2), inconsistent=rbind(resSig3,resSig4) ))
}


plot_scatter_lfc<-function(res,resSig,consistent, inconsistent,name,pdffile=FALSE){

	a<-res$log2FoldChange
	s<-res$log2FoldChange.1
	as<-resSig$log2FoldChange
	ss<-resSig$log2FoldChange.1
	amci<-consistent$log2FoldChange
	smci<-consistent$log2FoldChange.1
	amciv<-inconsistent$log2FoldChange
	smciv<-inconsistent$log2FoldChange.1
	if(pdffile==TRUE){
		pdf(paste('scaterplot.mean_sasLFC.3',name,'pdf',sep="."))
	}else{
		svg(filename=paste('scaterplot.mean_sasLFC',name,'svg',sep='.'),bg = "white")
	}
	plot(a,s,
			xlab="Mean antisense log fold change Tumor vs Normal",
			ylab="Mean sense log fold change Tumor vs Normal",
			cex=0.3,col="gray"
			,xlim=c(-8,8),ylim=c(-8,8))
	abline(v=0,h=0,col="gray")
# red candidates are loci for which the nature of the loci changes
# from sense to antisense or viceversa. And the lfc between sense and
# antisense strands is significant.
#	points(as,ss,
#			xlab="Mean antisense log fold change Tumor vs Normal",
#			ylab="Mean sense log fold change Tumor vs Normal",col="coral1",cex=0.3)

	points(amci,smci,
			xlab="Mean antisense log fold change Tumor vs Normal",
			ylab="Mean sense log fold change Tumor vs Normal",col="violetred",cex=0.6,pch=16)#coral1
	points(amciv,smciv,
			xlab="Mean antisense log fold change Tumor vs Normal",
			ylab="Mean sense log fold change Tumor vs Normal",col="forestgreen",cex=0.6,pch=16)
	
	dev.off()
	pdf(paste('scaterplot.mean_sasLFC.frame',name,'pdf',sep="."))
	plot(1,1,
			xlab="Mean antisense log fold change Tumor vs Normal",
			ylab="Mean sense log fold change Tumor vs Normal",
			cex=0.3,col="gray",
			xlim=c(-8,8),ylim=c(-10,10))#xlim=c(-10,10),ylim=c(-10,10))
	abline(v=0,h=0,col="gray")
	dev.off()
	
	
	
}


split_pair_name<-function(ov_exons,sep="_"){
	results <- matrix("NA",nrow=length(ov_exons),ncol=2)
	for(i in seq(1:length(ov_exons))){
		results[i,]<-unlist(strsplit(ov_exons[i],sep))[c(1,2)]
	}
	return(results)
}

add_gene_names <- function(resSig){
	gene_pairs<-split_pair_name(as.vector(resSig$gene_pair),"_"); 
	colnames(gene_pairs)<-c("gene_left_name","gene_right_name")
	resSig<-cbind(resSig,gene_pairs)
	gene_pairs<-split_pair_name(as.vector(resSig$gene_name),"_"); 
	colnames(gene_pairs)<-c("gene_left_hugo","gene_right_hugo")
	resSig<-cbind(resSig,gene_pairs)
	return(resSig)
}



order_overlapping_features <- function(res){
	biotype<-as.vector(res$transcript_biotype);
	gene_left_name<-as.vector(res$gene_left_name)
	gene_right_name<-as.vector(res$gene_right_name)
	new_gene_left_name<-gene_left_name
	new_gene_right_name<-gene_right_name
	
	for(i in seq(1,length(biotype))){
		b<-biotype[i]
		bp<-unlist(strsplit(b,split="&"));
		gl<-bp[1];gr<-bp[2]
		gl1<-gene_left_name[i];
		gr1<-gene_right_name[i];
		print(c(gl,gr))
		print(c(gl1,gr1))
		
		if((gl=="protein_coding") && (gr!="protein_coding")){
			new_gene_left_name[i]<-gl1
			new_gene_right_name[i]<-gr1
			
		}else if((gl!="protein_coding") && (gr=="protein_coding")){
			new_gene_left_name[i]<-gr1
			new_gene_right_name[i]<-gl1
		}else{
			new_gene_left_name[i]<-gl1
			new_gene_right_name[i]<-gr1
		}
	}
	res$new_gene_left_name<-new_gene_left_name
	res$new_gene_right_name<-new_gene_right_name
	
	return(res)
}


################################################################################
### code chunk number 2: SETUP AND SUB-ROUTINES
################################################################################

root_dir <- paste(base_dir,'sascompendia_codebase',sep='/')
data_dir <-paste(root_dir,'data',sep='')
data_dir_OPSnastiseq<-paste(data_dir,'nastiseq_final_R_objects',sep='/')
references_dir<-paste(root_dir,'reference_files/',sep='')
source_code_dir<-paste(root_dir,'codebase',sep='/')
results_dir<-paste(root_dir,'results',sep='/')

source(paste(source_code_dir,'sas_general_functions.R',sep=''))
source(paste(source_code_dir,'functions_tocompute_correlation_between_saspairs.R',sep=''))

results_dir_fig6 <- paste(results_dir,'figure6',sep='/')
if(!file.exists(results_dir_fig6)){
	system(paste("mkdir",results_dir_fig6,sep=" "))
}

###############################################################################
### code chunk number 3: Data structures to load
###############################################################################

trans_DESeq <- 'cohort.matrices.transcripts.DESeq.R'

# Overlapping features
overlap_annotation_filename<-paste(references_dir,"GRCh37.69.saspairs_trans.tsv",sep="")
feature_overlaps <- read.table(overlap_annotation_filename,header=TRUE,sep='\t')
feature_overlaps$overlap_type<-feature_overlaps$exon_overlap_type
feature_overlaps$gene_left_name<-feature_overlaps$gene_left; feature_overlaps$gene_right_name<-feature_overlaps$gene_right;	
gene_pairs<-split_exonpair_name(as.vector(feature_overlaps$gene_biotype),"&"); 
colnames(gene_pairs)<-c("gene_left_biotype","gene_right_biotype")
feature_overlaps_ori <- cbind(feature_overlaps, gene_pairs)
# Find the unique features. Eliminate pairs that involve the same to genes
# and the same type of overlap, in order to avoid calculating the correlation
# twice.
rmpairs <- which(as.vector(feature_overlaps_ori$exon_overlap_type) %in%  c("convergent","divergent","downstream","upstream"))	# upstream and downstream genes represent almost duplicated sets.
feature_overlaps <- feature_overlaps_ori[-rmpairs,]
feature_overlaps <- featureov_eliminate_redundancy(feature_overlaps)


library(DESeq)
load(TRANSCRIPT_MATRICES_ROBJ) # returns sense_cds, asense_cds
sense_counts<-counts(sense_cds,normalized=FALSE);
asense_counts<-counts(asense_cds,normalized=FALSE);

###############################################################################
### code chunk number 4: Determine concordant/discordant loci
###############################################################################

cohort_list_names[["luad_benign"]]<-c("pt_lung_A65N","pt_lung_A63N","pt_lung_A35N")
cohort_list_names[["luad_match"]]<-c("pt_lung_A65","pt_lung_A63","pt_lung_A35")
cohort_list_names[["lusc_benign"]]<-c("pt_lung_SCC07","pt_lung_SCC01N","pt_lung_SCC03N")

comparisons_table<-matrix(c("luad","luad_benign"),ncol=2,byrow=TRUE)

###############################################################################
### code chunk number 5:
### QUANTIFICATION OF DIFFERENTIALLY EXPRESSED TRANSCRIPTS
### SELECTION OF PAIRS CONSISTENTLY/INCONSISTENTLY REGULATED
###############################################################################

for(i in seq(1,nrow(comparisons_table))){

	querygroup_name<-comparisons_table[i,1]
	controlgroup_name<-comparisons_table[i,2]
	print(c(querygroup_name,controlgroup_name))
	querygroup <- cohort_list_names[[querygroup_name]];
	controlgroup<- cohort_list_names[[controlgroup_name]];
	
	sppairs <- c(querygroup,controlgroup)
	print(sppairs)
	
	min_count=1
	theta_q=min_count*length(querygroup);theta_c=min_count*length(controlgroup);
	rs_q = which(rowSums ( sense_counts[,querygroup] ) >= theta_q)
	rs_c = which(rowSums ( sense_counts[,querygroup] ) >= theta_c)
	use = union(rs_q,rs_c)
		
	sense_expmat <- sense_counts[use,sppairs]
	antisense_expmat <- asense_counts[use,sppairs]
	condition<-factor(c(rep(querygroup_name,length(querygroup)),rep(controlgroup_name,length(controlgroup))))
	pvalth=0.1
	results_dir_fig6_sub<-paste(results_dir_fig6,"/figure6_DESeq_correction",querygroup_name,"_",controlgroup_name,"_","pvalth",pvalth,sep="");
	
	if(!file.exists(results_dir_fig6_sub)){
		system(paste("mkdir",results_dir_fig6_sub,sep=" "))		
	}	
	setwd(results_dir_fig6_sub)
	# Calculate binomial test for each gene using DESeq
	smallest_group<-c()
	if(length(querygroup) >= length(controlgroup)){smallest_group<-controlgroup}else{smallest_group<-querygroup}
	res<-calculate_diffexp_genes_DESeq(sense_expmat,condition,controlgroup_name, querygroup_name,smallest_group,"forward")
	ares<-calculate_diffexp_genes_DESeq(antisense_expmat,condition,controlgroup_name, querygroup_name,smallest_group,"reverse")
	
	# Determine significant genes
	gene_pairs<-feature_overlaps[,c("gene_left","gene_right")]
	resSig<-select_diffgenes(res,gene_pairs,pvalth)
	
	#Annotation of the differentially expressed genes
	resSig$gene_pair <- as.vector(apply(resSig[,c(1,9)],1,paste,collapse="_"))
	annot<-feature_overlaps[which(feature_overlaps$gene_id %in% as.vector(resSig$gene_pair)),
			c("gene_name","transcript_biotype","trans_overlap_type","overlap_type","chr","start","end")]
	resSig <-cbind(resSig,annot)
	write.table(resSig,file=paste("resSigannot",querygroup_name,controlgroup_name,"tsv",sep="."),
			row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)
}

###############################################################################
### code chunk number 6: DETERMINE CONCORDANT AND DISCORDANT GENES
###############################################################################

#POST PROCESSING
#IF THESE ANALYSIS ARE DONE AS POST-PROCESSING
#load("nBimTestRes.luad_benign.luad_match.R")
#load("nBimTestRes.luad_benign.luad_match.forward.R")
#load("nBimTestRes.luad_benign.luad_match.reverse.R")

load("nBimTestRes.luad_benign.luad.reverse.R")
ares<-res
load("nBimTestRes.luad_benign.luad.forward.R")


gene_pairs<-feature_overlaps[,c("gene_left","gene_right")]
resSig<-select_diffgenes(res,gene_pairs,pvalth)

#Annotation of the differentially expressed genes
resSig$gene_pair <- as.vector(apply(resSig[,c(1,9)],1,paste,collapse="_"))
annot<-feature_overlaps[which(feature_overlaps$gene_id %in% as.vector(resSig$gene_pair)),
		c("gene_name","transcript_biotype","trans_overlap_type","overlap_type","chr","start","end")]
resSig <-cbind(resSig,annot)


lfcth=1.0
resSig2filter<-select_diffgenes_bylfc(resSig,lfcth)
consist<-resSig2filter$consistent;inconsist<-resSig2filter$inconsistent;

# Annotate all genes in results to see which ones have pairs in common
res2<-select_diffgenes(res,gene_pairs,pvalth=1)
res2$gene_pair <- as.vector(apply(res2[,c(1,9)],1,paste,collapse="_"))
annot<-feature_overlaps[which(feature_overlaps$gene_id %in% as.vector(res2$gene_pair)),
		c("gene_name","trans_overlap_type","transcript_biotype","overlap_type","chr","start","end")]
res2 <-cbind(res2,annot)

###############################################################################
### code chunk number 7:
### DO SCATTER PLOT OF CONSISTENT AND INCOSISTENTLY REGULATED
### TRANSCRIPTS ACROSS THE COMPARISON
###############################################################################

plot_scatter_lfc(res2,resSig,consist, inconsist,'pairs')
plot_scatter_lfc(res2,resSig,consist, inconsist,'pairs',TRUE)

# Write the final matrices for enrichment analysis.
resSig <- add_gene_names(resSig)
consist <- add_gene_names(consist)
inconsist <- add_gene_names(inconsist)
res2<-add_gene_names(res2)
querygroup_name<-"luad";controlgroup_name<-"luad_benign"
#querygroup_name<-"lusc";controlgroup_name<-"lusc_benign"
write.table(res2,file=paste("res_annot_alloverlapping",querygroup_name,controlgroup_name,"tsv",sep="."),row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)
write.table(resSig,file=paste("resSigannot2",querygroup_name,controlgroup_name,"tsv",sep="."),row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)
write.table(consist,file=paste("resSigannot_consistent",querygroup_name,controlgroup_name,"tsv",sep="."),row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)
write.table(inconsist,file=paste("resSigannot_inconsistent",querygroup_name,controlgroup_name,"tsv",sep="."),row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)


sense_norm<-counts(sense_cds,normalized=TRUE);
asense_norm<-counts(asense_cds,normalized=TRUE);
sense_expmat <- sense_norm[use,sppairs]
antisense_expmat <- asense_norm[use,sppairs]

sensedelta <- log2(sense_expmat+1) - log2(res$baseMeanA+1)
#deltasheatmap_tmp<- cbind(sensedelta[sigchangers1,querygroup],antisensedelta[sigchangers1,querygroup])
consist<-order_overlapping_features(consist)
deltasheatmap_tmp<- cbind(sensedelta[consist$new_gene_left_name,querygroup],sensedelta[consist$new_gene_right_name,querygroup])

#o3<-order(rowMeans(sensedelta[sigchangers1,querygroup],na.rm=TRUE),decreasing=TRUE)#
o3<-order(rowMeans(deltasheatmap_tmp,na.rm=TRUE),decreasing=TRUE)
o3<-order(rowMeans(sensedelta[consist$new_gene_left_name,querygroup],na.rm=TRUE),decreasing=TRUE)
deltasheatmap<-deltasheatmap_tmp[o3,]
heatmap_figure_fn1<-"DESeq_consistent_genes.pdf"

ramp<-colorRamp(c("blue", "white","red"))
MIN<-min(deltasheatmap,na.rm=TRUE); MAX<-max(deltasheatmap,na.rm=TRUE)
LIM<-round(min(abs(MIN),abs(MAX)),0)
LIM<-2
#palette.breaks <- seq(MIN, MAX, 0.1)
palette.breaks <- seq(-1*LIM, LIM, 0.1)
palette.breaks <- c(MIN,palette.breaks,MAX)
breaks=1/(length(palette.breaks)-2)
color.palette <- rgb( ramp(seq(0, 1, breaks)), max = 255)
#sppairs2<-c("pt_lung_A65N","pt_lung_A63N","pt_lung_A35N","pt_lung_A65","pt_lung_A63","pt_lung_A35")
#gnames<-as.vector(geneAnnot[rownames(deltasheatmap),"gene_name"])
gnames<-as.vector(consist$gene_name)
pdf(heatmap_figure_fn1)
image(t(as.matrix(deltasheatmap[order(seq(1:nrow(deltasheatmap)),decreasing=TRUE),])), 
		col= color.palette,breaks = palette.breaks,axes=FALSE)# 
axis(2,at=seq(0,1,length=length(gnames)),labels=gnames[order(seq(1:nrow(deltasheatmap)),decreasing=TRUE)],las=2,cex.axis=0.3)
dev.off()


# Change the order of protein_coding/antisense
inconsist2<-order_overlapping_features(inconsist)
inconsist<-inconsist2
deltasheatmap_tmp<- cbind(sensedelta[inconsist$new_gene_left_name,querygroup],sensedelta[inconsist$new_gene_right_name,querygroup])
o3<-order(rowMeans(sensedelta[inconsist$new_gene_left_name,querygroup],na.rm=TRUE),decreasing=TRUE)
deltasheatmap<-deltasheatmap_tmp[o3,]
heatmap_figure_fn1<-"DESeq_inconsistent_genes.pdf"

ramp<-colorRamp(c("blue", "white","red"))
MIN<-min(deltasheatmap,na.rm=TRUE); MAX<-max(deltasheatmap,na.rm=TRUE)
LIM<-round(min(abs(MIN),abs(MAX)),0)
LIM<-2
#palette.breaks <- seq(MIN, MAX, 0.1)
palette.breaks <- seq(-1*LIM, LIM, 0.1)
palette.breaks <- c(MIN,palette.breaks,MAX)
breaks=1/(length(palette.breaks)-2)
color.palette <- rgb( ramp(seq(0, 1, breaks)), max = 255)
#sppairs2<-c("pt_lung_A65N","pt_lung_A63N","pt_lung_A35N","pt_lung_A65","pt_lung_A63","pt_lung_A35")
#gnames<-as.vector(geneAnnot[rownames(deltasheatmap),"gene_name"])
gnames<-as.vector(inconsist$gene_name)
pdf(heatmap_figure_fn1)
image(t(as.matrix(deltasheatmap[order(seq(1:nrow(deltasheatmap)),decreasing=TRUE),])), 
		col= color.palette,breaks = palette.breaks,axes=FALSE)# 
axis(2,at=seq(0,1,length=length(gnames)),labels=gnames[order(seq(1:nrow(deltasheatmap)),decreasing=TRUE)],las=2,cex.axis=0.3)
dev.off()


###############################################################################
### code chunk number 8:
### SECOND ANALYSIS USING THE EXPRESSIONF OVER FORWARD AND REVERSE STRAND OF A 
### PARTICULAR LOCUS
###############################################################################

lfcth=1.0
resSigSAS<-select_diffgenes_twomat(res,ares,as.vector(res$id),pvalth)

#resSigSAS$gene_pair <- as.vector(apply(resSigSAS[,c(1,9)],1,paste,collapse="_"))
#annot<-feature_overlaps[which(feature_overlaps$gene_id %in% as.vector(resSigSAS$gene_pair)),
#		c("gene_name","trans_overlap_type","transcript_biotype","overlap_type","chr","start","end")]
#resSigSAS <-cbind(resSigSAS,annot)

resSigSAS$new_gene_left_name<-resSigSAS$id;
resSigSAS$new_gene_right_name<-resSigSAS$id.1;
resSigSAS2filter<-select_diffgenes_bylfc(resSigSAS,lfcth)
sasconsist<-resSigSAS2filter$consistent;sasinconsist<-resSigSAS2filter$inconsistent;
#
plot_scatter_lfc(res2,resSigSAS2filter,sasconsist, sasinconsist,'sas')
plot_scatter_lfc(res2,resSigSAS2filter,sasconsist, sasinconsist,'sas',TRUE)

write.table(resSigSAS,file=paste("resSigannot2_sas",querygroup_name,controlgroup_name,"tsv",sep="."),row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)
write.table(sasconsist,file=paste("resSigannot_consistent_sas",querygroup_name,controlgroup_name,"tsv",sep="."),row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)
write.table(sasinconsist,file=paste("resSigannot_inconsistent_sas",querygroup_name,controlgroup_name,"tsv",sep="."),row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)


sense_norm<-counts(sense_cds,normalized=TRUE);
asense_norm<-counts(asense_cds,normalized=TRUE);
sense_expmat <- sense_norm[use,sppairs]
antisense_expmat <- asense_norm[use,sppairs]

sensedelta <- log2(sense_expmat+1) - log2(res$baseMeanA+1)
asensedelta <- log2(antisense_expmat+1) - log2(ares$baseMeanA+1)
#deltasheatmap_tmp<- cbind(sensedelta[sigchangers1,querygroup],antisensedelta[sigchangers1,querygroup])
#consist<-order_overlapping_features(consist)
consist<-sasconsist
deltasheatmap_tmp<- cbind(sensedelta[consist$new_gene_left_name,querygroup],asensedelta[consist$new_gene_right_name,querygroup])

#o3<-order(rowMeans(sensedelta[sigchangers1,querygroup],na.rm=TRUE),decreasing=TRUE)#
o3<-order(rowMeans(deltasheatmap_tmp,na.rm=TRUE),decreasing=TRUE)
o3<-order(rowMeans(sensedelta[consist$new_gene_left_name,querygroup],na.rm=TRUE),decreasing=TRUE)
deltasheatmap<-deltasheatmap_tmp[o3,]
heatmap_figure_fn1<-"DESeq_consistent_genes.sas.pdf"

ramp<-colorRamp(c("blue", "white","red"))
MIN<-min(deltasheatmap,na.rm=TRUE); MAX<-max(deltasheatmap,na.rm=TRUE)
LIM<-round(min(abs(MIN),abs(MAX)),0)
LIM<-2
#palette.breaks <- seq(MIN, MAX, 0.1)
palette.breaks <- seq(-1*LIM, LIM, 0.1)
palette.breaks <- c(MIN,palette.breaks,MAX)
breaks=1/(length(palette.breaks)-2)
color.palette <- rgb( ramp(seq(0, 1, breaks)), max = 255)
#sppairs2<-c("pt_lung_A65N","pt_lung_A63N","pt_lung_A35N","pt_lung_A65","pt_lung_A63","pt_lung_A35")
#gnames<-as.vector(geneAnnot[rownames(deltasheatmap),"gene_name"])
gnames<-as.vector(geneAnnot_ori[consist$new_gene_left_name,"gene_name"])
pdf(heatmap_figure_fn1)
image(t(as.matrix(deltasheatmap[order(seq(1:nrow(deltasheatmap)),decreasing=TRUE),])), 
		col= color.palette,breaks = palette.breaks,axes=FALSE)# 
axis(2,at=seq(0,1,length=length(gnames)),labels=gnames[order(seq(1:nrow(deltasheatmap)),decreasing=TRUE)],las=2,cex.axis=0.3)
dev.off()


# Change the order of protein_coding/antisense
#inconsist2<-order_overlapping_features(inconsist)
inconsist<-sasinconsist
deltasheatmap_tmp<- cbind(sensedelta[inconsist$new_gene_left_name,querygroup],asensedelta[inconsist$new_gene_right_name,querygroup])
o3<-order(rowMeans(sensedelta[inconsist$new_gene_left_name,querygroup],na.rm=TRUE),decreasing=TRUE)
deltasheatmap<-deltasheatmap_tmp[o3,]
heatmap_figure_fn1<-"DESeq_inconsistent_genes.sas.pdf"

ramp<-colorRamp(c("blue", "white","red"))
MIN<-min(deltasheatmap,na.rm=TRUE); MAX<-max(deltasheatmap,na.rm=TRUE)
LIM<-round(min(abs(MIN),abs(MAX)),0)
LIM<-2
#palette.breaks <- seq(MIN, MAX, 0.1)
palette.breaks <- seq(-1*LIM, LIM, 0.1)
palette.breaks <- c(MIN,palette.breaks,MAX)
breaks=1/(length(palette.breaks)-2)
color.palette <- rgb( ramp(seq(0, 1, breaks)), max = 255)
#sppairs2<-c("pt_lung_A65N","pt_lung_A63N","pt_lung_A35N","pt_lung_A65","pt_lung_A63","pt_lung_A35")
#gnames<-as.vector(geneAnnot[rownames(deltasheatmap),"gene_name"])
gnames<-as.vector(geneAnnot_ori[inconsist$new_gene_left_name,"gene_name"])
pdf(heatmap_figure_fn1)
image(t(as.matrix(deltasheatmap[order(seq(1:nrow(deltasheatmap)),decreasing=TRUE),])), 
		col= color.palette,breaks = palette.breaks,axes=FALSE)# 
axis(2,at=seq(0,1,length=length(gnames)),labels=gnames[order(seq(1:nrow(deltasheatmap)),decreasing=TRUE)],las=2,cex.axis=0.3)
dev.off()


fishertable <-matrix(c(6807,1944,2152,125,133,65,60,14,37),ncol=3,byrow=T)
rownames(fishertable)<-c("ALL","CONST","INCO")
colnames(fishertable)<-c("EMB","HTH","TTT")
fisher.test(fishertable,simulate.p.value = TRUE, B = 1e5)

fishertable <-matrix(c(6807,1944,125,133),ncol=2,byrow=T)
rownames(fishertable)<-c("ALL","CONST")
colnames(fishertable)<-c("EMB","HTH")
fisher.test(fishertable,alternative="two.sided")

fishertable <-matrix(c(6807,1944,2152,171,131,84,74,11,64),ncol=3,byrow=T)
rownames(fishertable)<-c("ALL","CONST","INCO")
colnames(fishertable)<-c("EMB","HTH","TTT")
fisher.test(fishertable,simulate.p.value = TRUE, B = 1e5)

fishertable <-matrix(c(6807,1944,171,131),ncol=2,byrow=T)
rownames(fishertable)<-c("ALL","INCO")
colnames(fishertable)<-c("EMB","HTH")
fisher.test(fishertable,alternative="two.sided")

fishertable <-matrix(c(6807,1944,74,11),ncol=2,byrow=T)
rownames(fishertable)<-c("ALL","INCO")
colnames(fishertable)<-c("EMB","HTH")
fisher.test(fishertable,alternative="two.sided")

##### LOOK FOR ACTIONABLE GENES
tumor_supressors <- read.table(TUMORSUP_FILE,header=TRUE,skip=1,sep='\t')
oncogenes <- read.table(ONOCOGENES_FILE,header=TRUE,skip=1,sep='\t')
# split those gene names in three non-overlapping groups
genes_classified_asboth<-intersect(tumor_supressors$GeneSymbol,oncogenes$GeneSymbol)
tumor_supressors<-setdiff(tumor_supressors$GeneSymbol,genes_classified_asboth)
oncogenes<-setdiff(oncogenes$GeneSymbol,genes_classified_asboth)
oncots<-genes_classified_asboth
#
tumor_supressors_ensl<-rownames(geneAnnot_ori)[which(geneAnnot_ori$gene_name %in% tumor_supressors)]
oncogenes_ensl<-rownames(geneAnnot_ori)[which(geneAnnot_ori$gene_name %in% oncogenes)]
oncotusup_ensl<-rownames(geneAnnot_ori)[which(geneAnnot_ori$gene_name %in% oncots)]
all_other_genes<-setdiff(rownames(geneAnnot_ori),unique(c(tumor_supressors_ensl,oncogenes_ensl,oncotusup_ensl)))
all_other_genes_protcod<-intersect(all_other_genes,rownames(geneAnnot_ori)[which(geneAnnot_ori$gene_biotype == "protein_coding")])
actionable<-c(tumor_supressors_ensl,oncogenes_ensl,oncotusup_ensl)

consist_actionable<-intersect(consist$id,actionable)
inconsist_actionable<-intersect(inconsist$id,actionable)

consist_actionable<-cbind(geneAnnot_ori[consist_actionable,],consist[consist_actionable,])
inconsist_actionable<-cbind(geneAnnot_ori[inconsist_actionable,],inconsist[inconsist_actionable,])

write.table(consist_actionable,file=paste("consist_actionable",querygroup_name,controlgroup_name,"tsv",sep="."),
		row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)

write.table(inconsist_actionable,file=paste("inconsist_actionable",querygroup_name,controlgroup_name,"tsv",sep="."),
		row.names=TRUE,col.names=NA,sep='\t',eol='\n',quote=FALSE)

