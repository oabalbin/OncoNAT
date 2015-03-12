# TODO: Add comment
# 
# Author: alebalbin
# /var/folders/kr/krqck8AdEIKtDx3cghx5kk+++TI/-Tmp-//RtmpsBHMh3/downloaded_packages
# Aug/2013
# This functio will load the matrices of read counts or either full transcript
# or exons
# It will return a DESeq count dataset with normalized counts, estimated size factors
# and estimated dispersion.
###############################################################################
#Libraries
load_data<-function(data_dir,expr_matrix_for_filename,
	expr_matrix_rev_filename,EXONS,
	number_lines_transcripts,number_lines_exons,
	header_lines,annotation_lines,gene_annot_cols,row_names_col_transcripts,
	row_names_col_exons,positive_reads_counts,negative_reads_counts,feature_overlaps
		){
	library(DESeq)
	#data_dir <-paste(root_dir,'runs/aug_gene_cov',sep='')
#	number_lines_transcripts=50971
#	number_lines_exons=660716
#	header_lines=6
#	annotation_lines=2#6
#	gene_annot_cols=8
#	row_names_col_transcripts=2
#	row_names_col_exons=3
#	positive_reads_counts=4
#	negative_reads_counts=5
	
	# Matrix of genes
	#cohort_for_transcripts_filename <- 'cohort.posstrand.cov.matrix.08_2013.transcripts.tsv'
	#cohort_rev_transcripts_filename <- 'cohort.negstrand.cov.matrix.08_2013.transcripts.tsv'
	# Matrix of exons
	#cohort_for_exons_filename <- 'cohort.posstrand.cov.matrix.08_2013.exons.tsv'
	#cohort_rev_exosn_filename <- 'cohort.negstrand.cov.matrix.08_2013.exons.tsv'
	
	
	###################################################
	### code chunk number : load count data matrices
	###################################################
	setwd(data_dir)
	min_pe=0.2
	# Load positivestrand counts
	if(EXONS==TRUE){
		cohort_for_filename<-expr_matrix_for_filename
		cohort_rev_filename<-expr_matrix_rev_filename
		number_lines<-number_lines_exons
		row_names_col<-row_names_col_exons
		deseq_filename<-"cohort.matrices.exons.DESeq.2014.R"
	}else{
		cohort_for_filename<-expr_matrix_for_filename
		cohort_rev_filename<-expr_matrix_rev_filename
		number_lines<-number_lines_transcripts
		row_names_col<-row_names_col_transcripts
		deseq_filename<-"cohort.matrices.transcripts.DESeq.2014.R"
	}
	print("Loading positive strand expression matrix")
	
	mf<-load_counts_data(cohort_for_filename, number_lines,header_lines,
			gene_annot_cols,row_names_col,annotation_lines)
		
	print("Done loading positive strand expression matrix")
	# Load negativestrand counts
	print("Loading negative strand expression matrix")
	mr<-load_counts_data(cohort_rev_filename, number_lines,header_lines,
			gene_annot_cols,row_names_col,annotation_lines)
	print("Done loading negative strand expression matrix")
	
	# Select the sense and antisense expression matrices
	print("Creating sense and antisense matrices")
	sample_names<-colnames(mf$mat)
	loci_for <- which(mf$row_annot[,"strand"]=="+");loci_rev <- which(mr$row_annot[,"strand"]=="-");
	smat<-rbind(mf$mat[loci_for,sample_names],mr$mat[loci_rev,sample_names]);
	amat<-rbind(mr$mat[loci_for,sample_names],mf$mat[loci_rev,sample_names]);
	geneAnnot_ori <- rbind(mf$row_annot[loci_for,],mr$row_annot[loci_rev,])
	print("Done loading the expression matrices")
	
	###################################################
	### code chunk number: calculate the p.e using 
	### genes that do not belong to pairs
	###################################################
	print("Computing ssRNASeq protocol error rate")
	all_genesInpairs<-unique(as.vector(as.matrix(feature_overlaps[,c("gene_left","gene_right")])))
	all_genesNOTinpairs<-setdiff(rownames(geneAnnot_ori),all_genesInpairs)
	# Restrict to protein coding genes to stablish the negative distribution
	negpairsAnnot<-geneAnnot_ori[all_genesNOTinpairs,]
	negpairsAnnot<-negpairsAnnot[which(negpairsAnnot$gene_biotype=="protein_coding"),]
	negpairs <- unique(rownames(negpairsAnnot))
	smatsum <- apply(smat[negpairs,], 2, sum)
	asmatsum <- apply(amat[negpairs,], 2, sum)
		
	pe <- asmatsum/(smatsum + asmatsum)
	write.table(pe,file="strand_specific_error_matrix.tsv",row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
	###################################################
	### code chunk number: calculate the p.e using 
	### for each genes that do not belong to pairs (instead of sample level)
	### generate the boxplots protcolerror_genelevel...
	###################################################
	
#	gsmatsum <- apply(sense_norm_trasn[negpairs,], 1, sum)
#	gasmatsum <- apply(asense_norm_trans[negpairs,], 1, sum)
#	pe_genelevel <- gasmatsum/(gsmatsum + gasmatsum)
#	
#	gtotalcov<-smatsum+asmatsum+1
#	hcmat<-cbind(pe_genelevel,totalcov)
#	hcmat<-hcmat[order(hcmat[,2]),]
#	hcmat<-hcmat[which(-is.na(hcmat[,1])),]
#	bins<-seq(1,nrow(hcmat),200)
#	repbins<-sort(c(rep(bins[-length(bins)],200),rep(nrow(hcmat),(nrow(hcmat)-200*(length(bins)-1)))))
#	hcmatd<-data.frame(expbin=as.factor(repbins),pegenelevel=hcmat[,1],readcounts=hcmat[,2])
#	meanexp<-aggregate(hcmatd$readcounts,list(hcmatd$expbin),mean)
#	
#	pdf("protocolerror_genelevel_vs_total.boxplot.pdf")
#	boxplot(hcmat[,1] ~ as.factor(repbins),las=2,names=round(log10(meanexp[,2]),2),xlab="Bin expression (log10(norm(count)))",
#			ylab="Fraction of reads mapping to the unintended (opposite) strand")
#	dev.off()
#	pdf("protocolerror_genelevel_vs_total.boxplot.subset.pdf")
#	boxplot(hcmat[,1] ~ as.factor(repbins),las=2,names=round(log10(meanexp[,2]),2),
#			ylim=c(0,0.20),xlab="Bin expression (log10(norm(count)))",ylab="Fraction of reads mapping to the unintended (opposite) strand")
#	dev.off()
#	
	
	
	#Samples for which the strand specificity is very low and the pe is quite high.
    sp_blacklisted<-which(mf$col_annot["annotation8_disease_type",]=="BENIGN_BODYMAP")
	sp_blacklisted<-union(sp_blacklisted,which(pe>=min_pe))
	
	# Adjust the smat, and amat
	smat<-smat[,-sp_blacklisted]
	amat<-amat[,-sp_blacklisted]
	mf$col_annot <-mf$col_annot[,-sp_blacklisted]
	pe<-pe[-sp_blacklisted]
	cohortAnnot<-data.frame(
			sample_names=colnames(smat),
			sampleid=as.vector(as.matrix(mf$col_annot["annotation1_sample_id",])),
			sampletype=as.vector(as.matrix(mf$col_annot["annotation2_sample_type",])),
			nucleotidetype=as.vector(as.matrix(mf$col_annot["annotation3_nucleotide_type",])),
			mctpcohort=as.vector(as.matrix(mf$col_annot["annotation4_cohort",])),
			disease=as.vector(as.matrix(mf$col_annot["annotation5_disease",])),
			cancerprogression=as.vector(as.matrix(mf$col_annot["annotation6_cancer_progression",])),
			tumorcontent=as.vector(as.matrix(mf$col_annot["annotation7_tumor_content",])),			
			tissuetype=as.vector(as.matrix(mf$col_annot["annotation8_disease_type",]))
			)
	cohortLibAnnot<-mf$library_annot#hda
	rm(mr); rm(mf); gc();		
	# tmp: Change once the matrix annotation is corrected
	condition<-cohortAnnot$tissuetype#as.factor(c(as.vector(as.matrix(mf$col_annot["annotation1_tissue_type",])),"LUAD","LUAD"))
	
	###################################################
	### code chunk number : Create a count matrix for DESeq
	# The idea is to use DESeq to normalize the data and avoid the rpkm normalization.
	# DESeq could also be used to determine differential expressed genes or exons.
	###################################################
	print("Starting DESeq data set")
	sense_cds<-newCountDataSet(smat,condition)
	asense_cds<-newCountDataSet(amat,condition)
	# Remove duplicate matrices and garbage memory collection
	rm(smat);rm(amat);gc();
	
	
	###################################################
	### code chunk number : estimateSizeFactors
	### and dispersion
	###################################################
	print("Starting DESeq size factor and dispersion estimation")
	sense_cds = estimateSizeFactors( sense_cds )
	asense_cds = estimateSizeFactors( asense_cds )
	sense_cds = estimateDispersions( sense_cds,method="per-condition",sharingMode="gene-est-only",fitType="local")
	asense_cds = estimateDispersions( asense_cds,method="per-condition",sharingMode="gene-est-only",fitType="local")
	
	# Save the DESeq representation
	#save(sense_cds,asense_cds, geneAnnot_ori, cohortAnnot, file = deseq_filename)
	
	
	# Define the samples that belong to each cohort.
	luad<-which(cohortAnnot$tissuetype=="LUAD"); lusc<-which(cohortAnnot$tissuetype=="LUSC");
	brca<-which(cohortAnnot$tissuetype=="BRCA"); prca<-which(cohortAnnot$tissuetype=="PRCA");
	panc<-which(cohortAnnot$tissuetype=="PANC"); ov<-which(cohortAnnot$tissuetype=="OVARIAN");
	lucl<-which(cohortAnnot$tissuetype=="LUCL"); benign<-which(cohortAnnot$tissuetype=="BENIGN" & cohortAnnot$sampletype=="tissue")
	meningioma<-which(cohortAnnot$tissuetype=="MENINGIOMA");#benign_bm<-which(cohortAnnot$tissuetype=="BENIGN_BODYMAP");
	cholangio<-which(cohortAnnot$tissuetype=="CHOLANGIO");sarcoma<-which(cohortAnnot$tissuetype=="SARCOMA");
	
	# Define the cohort that you want to work with.
	cohort_list_names <- list(luad=as.vector(cohortAnnot$sample_names)[luad],
			lusc=as.vector(cohortAnnot$sample_names)[lusc],
			brca=as.vector(cohortAnnot$sample_names)[brca],
			prca=as.vector(cohortAnnot$sample_names)[prca],
			lucl=as.vector(cohortAnnot$sample_names)[lucl],
			benign=as.vector(cohortAnnot$sample_names)[benign],
			meningioma=as.vector(cohortAnnot$sample_names)[meningioma],
			ov=as.vector(cohortAnnot$sample_names)[ov],
			panc=as.vector(cohortAnnot$sample_names)[panc],
			cholangio=as.vector(cohortAnnot$sample_names)[cholangio],
			sarcoma=as.vector(cohortAnnot$sample_names)[sarcoma]
	)
#	cohort_list <- list(luad=luad,lusc=lusc,lucl=lucl)
#	cohort_list_names <- list(luad=as.vector(cohortAnnot$colname)[luad],
#			lusc=as.vector(cohortAnnot$colname)[lusc],
#			lucl=as.vector(cohortAnnot$colname)[lucl]
#	)
	# This function should return the sense_cds, and asense_cds
	sense_norm <- counts(sense_cds,normalized=TRUE)
	asense_norm<- counts(asense_cds,normalized=TRUE)

	save(sense_cds,asense_cds, 
			sense_norm,asense_norm,
			geneAnnot_ori, cohortAnnot,
			cohort_list_names,
			protocolerror=pe,
			cohortLibAnnot,
			file = deseq_filename)
		
	return(list(sense_cds=sense_cds,asense_cds=asense_cds,
					msense=sense_norm,masense=asense_norm,
					geneAnnot=geneAnnot_ori,cohortAnnot=cohortAnnot,
					cohort_list_names=cohort_list_names,
					protocolerror=pe,
					cohortLibAnnot=cohortLibAnnot
					))

}