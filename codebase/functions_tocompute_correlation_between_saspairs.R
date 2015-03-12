# TODO: Add comment
# 
# Author: alebalbin
# Aug 13/2013
###############################################################################

perform_vector_correlation <- function(mat1,mat2,pval_th,cmethod,
		gene_annot_cols,log,all){
	#cmethod= spearman,pearson, kendall
	pad=1
	m1<-mat1[,gene_annot_cols:ncol(mat1)]
	m2<-mat2[,gene_annot_cols:ncol(mat2)]
	results<-matrix(NA,nrow=nrow(m1),ncol=4)
	
	for(i in seq(1,nrow(m1))){
		#print(i)
		if(log==TRUE){
			print(i)
			x<-log10(as.numeric(as.vector(m1[i,]))+pad)
			y<-log10(as.numeric(as.vector(m2[i,]))+pad)
			
		}else{
			x<-as.numeric(as.vector(m1[i,]))
			y<-as.numeric(as.vector(m2[i,]))
		}
		
		t<-cor.test(x,y,alternative ="two.sided",
				method=cmethod,
				conf.level = 0.95,
				exact=FALSE
		)
		results[i,1]<-i
		results[i,2]<-t$estimate;results[i,3]<-t$p.value;results[i,4]<-t$statistic;
		
		#if(t$p.value <= pval_th){
		#	results<-rbind(results,c(i,t$estimate,t$p.value,t$statistic))
		#	#print(c(i,t$estimate,t$p.value,t$statistic))
		#}
	}
	colnames(results)<-c("gene_ind","rho","pvalue","statistic")
	results <- as.data.frame(results)
	results$pval_adjusted <- p.adjust(results$pvalue, method = "hochberg")
	results$fdr <- p.adjust(results$pvalue, method = "fdr")
	results <- cbind(mat1[,1:(gene_annot_cols-1)],results)
	#o<-order(results$rho,decreasing=TRUE)
	#results <- results[o,]
	
	return(results)
}


perform_allvector_correlation <- function(m1,m2,
		pval_th,
		gene_annot_cols,log,
		add_row_labels=TRUE){
	#cmethod= spearman,pearson, kendall
	library(robust)
	library(entropy)
	
	pad=1
	results<-matrix(NA,nrow=nrow(m1),ncol=13)
	MIN_NUM_POINTS=10#5
	irecord=0
	epsilon=1e3
	
	for(i in seq(1,nrow(m1))){
		gene_left<-rownames(m1)[i];gene_right<-rownames(m2)[i];
		#print(c(gene_left,gene_right))
		if(log==TRUE){
			x<-log10(as.numeric(as.vector(m1[i,]))+pad); y<-log10(as.numeric(as.vector(m2[i,]))+pad)
		}else{
			xo<-as.numeric(as.vector(m1[i,])); yo<-as.numeric(as.vector(m2[i,]));
			x <-scale(xo,center=TRUE,scale=TRUE); y <-scale(yo,center=TRUE,scale=TRUE);
		}
		results[i,1]<-gene_left;results[i,2]<-gene_right;
		if((length(intersect(which(!is.na(x)),which(!is.na(y))))> MIN_NUM_POINTS)&(
					length(intersect(which(xo!=0),which(yo!=0))) > MIN_NUM_POINTS)){
			p<-cor.test(x,y,alternative ="two.sided",
					method='pearson',
					conf.level = 0.95,
					exact=FALSE)
			s<-cor.test(x,y,alternative ="two.sided",
					method='spearman',
					conf.level = 0.95,
					exact=FALSE)
			
			results[i,3]<-p$estimate;results[i,4]<-p$p.value;results[i,5]<-p$statistic;
			results[i,6]<-s$estimate;results[i,7]<-s$p.value;results[i,8]<-s$statistic;
			# # INCLUDE HERE THE MUTUAL INFORMATION 
			xmax=max(as.numeric(as.vector(m1[i,]))+1,na.rm=TRUE);
			ymax=max(as.numeric(as.vector(m2[i,]))+1,na.rm=TRUE);
			ybin=10;
			xbin=ybin*ceiling(xmax/ymax)
			freqs2d = discretize2d( as.numeric(as.vector(m1[i,])), as.numeric(as.vector(m2[i,])), 
					xbin, ybin, r1=range(as.numeric(as.vector(m1[i,]))), r2=range(as.numeric(as.vector(m2[i,]))) )
			
			results[i,12]<-mi.plugin(freqs2d)
			results[i,13]<-0.5*chi2indep.plugin(freqs2d)			
			#print(results[i,])
			#Compute Robust correlation
			error.flag=FALSE
			t <- tryCatch(covRob(cbind(x,y),corr=TRUE,distance=TRUE), error=function(e){
						warning(e);error.flag<<-TRUE})
			if (error.flag){
				#print(paste("there was an error in covRob for genes ",gene_left,gene_right,sep=" "))
				results[i,9]<-NA;results[i,10]<-NA;results[i,11]<-NA;
			}else{
				results[i,9]<-t$cov[2];results[i,10]<-t$center[1];results[i,11]<-t$center[2];
			}
		}else{
			results[i,3]<-NA;results[i,4]<-NA;results[i,5]<-NA;
			results[i,6]<-NA;results[i,7]<-NA;results[i,8]<-NA;
			results[i,9]<-NA;results[i,10]<-NA;results[i,11]<-NA;
			results[i,12]<-NA;results[i,13]<-NA;
		}
	}
	colnames(results)<-c("gene_left_corr","gene_right_corr","prho","ppvalue","pstatistic","srho","spvalue","sstatistic","rrho","center1","center2","mutualinfo","chi2ind")
	results <- as.data.frame(results)
	results$spvalue_adjusted <- p.adjust(as.numeric(as.vector(results$spvalue)), method = "hochberg")
	results$spvalue_fdr <- p.adjust(as.numeric(as.vector(results$spvalue)), method = "fdr")
	if(add_row_labels==TRUE){
		results <- cbind(gene_annot_cols,results)
	}
	#
	results <- as.data.frame(results)
	
	#o<-order(results$rho,decreasing=TRUE)
	#results <- results[o,]
	print("DONE COMPUTING CORRELATION BETWEEN PAIRS")
	
	return(results)
}

perform_allvector_correlation_simple <- function(m1,m2,
		pval_th,
		gene_annot_cols,log,
		add_row_labels=TRUE){
	#cmethod= spearman,pearson, kendall
	library(robust)
	library(entropy)
	
	pad=1
	results<-matrix(NA,nrow=nrow(m1),ncol=8)
	MIN_NUM_POINTS=10#5
	irecord=0
	epsilon=1e3
	
	for(i in seq(1,nrow(m1))){
		gene_left<-rownames(m1)[i];gene_right<-rownames(m2)[i];
		if(log==TRUE){
			x<-log10(as.numeric(as.vector(m1[i,]))+pad); y<-log10(as.numeric(as.vector(m2[i,]))+pad)
		}else{
			xo<-as.numeric(as.vector(m1[i,])); yo<-as.numeric(as.vector(m2[i,]));
			x <-scale(xo,center=TRUE,scale=TRUE); y <-scale(yo,center=TRUE,scale=TRUE);
		}
		results[i,1]<-gene_left;results[i,2]<-gene_right;
		if((length(intersect(which(!is.na(x)),which(!is.na(y))))> MIN_NUM_POINTS)&(
					length(intersect(which(xo!=0),which(yo!=0))) > MIN_NUM_POINTS)){
			p<-cor.test(x,y,alternative ="two.sided",
					method='pearson',
					conf.level = 0.95,
					exact=FALSE)
			s<-cor.test(x,y,alternative ="two.sided",
					method='spearman',
					conf.level = 0.95,
					exact=FALSE)
			
			results[i,3]<-p$estimate;results[i,4]<-p$p.value;results[i,5]<-p$statistic;
			results[i,6]<-s$estimate;results[i,7]<-s$p.value;results[i,8]<-s$statistic;
		}else{
			results[i,3]<-NA;results[i,4]<-NA;results[i,5]<-NA;
			results[i,6]<-NA;results[i,7]<-NA;results[i,8]<-NA;
		}
	}
	colnames(results)<-c("gene_left_corr","gene_right_corr","prho","ppvalue","pstatistic","srho","spvalue","sstatistic")
	results <- as.data.frame(results)
	results$spvalue_adjusted <- p.adjust(as.numeric(as.vector(results$spvalue)), method = "hochberg")
	results$spvalue_fdr <- p.adjust(as.numeric(as.vector(results$spvalue)), method = "fdr")
	if(add_row_labels==TRUE){
		results <- cbind(gene_annot_cols,results)
	}
	#
	results <- as.data.frame(results)
	
	#o<-order(results$rho,decreasing=TRUE)
	#results <- results[o,]
	print("DONE COMPUTING CORRELATION BETWEEN PAIRS")
	
	return(results)
}

perform_vector_robustcorrelation <- function(mat1,mat2,pval_th,cmethod,
		gene_annot_cols,log){
	#cmethod= spearman,pearson, kendall
	pad=1
	m1<-mat1[,gene_annot_cols:ncol(mat1)]
	m2<-mat2[,gene_annot_cols:ncol(mat2)]
	results<-matrix(NA,nrow=nrow(m1),ncol=4)
	
	for(i in seq(1,nrow(m1))){
		#print(i)
		if(log==TRUE){
			print(i)
			x<-log10(as.numeric(as.vector(m1[i,]))+pad)
			y<-log10(as.numeric(as.vector(m2[i,]))+pad)
			print(x)
			print(y)
		}else{
			x<-as.numeric(as.vector(m1[i,]))
			y<-as.numeric(as.vector(m2[i,]))
		}
		
		t<-covRob(cbind(x,y),corr=TRUE,distance=TRUE)
		
		results[i,1]<-i
		results[i,2]<-t$cov[2];results[i,3]<-t$center[1];results[i,4]<-t$center[2];
		
		#if(t$p.value <= pval_th){
		#	results<-rbind(results,c(i,t$estimate,t$p.value,t$statistic))
		#	#print(c(i,t$estimate,t$p.value,t$statistic))
		#}
	}
	colnames(results)<-c("gene_ind","rho","center1","center2")
	results <- as.data.frame(results)
	#results$pval_adjusted <- p.adjust(results$pvalue, method = "hochberg")
	#results$fdr <- p.adjust(results$pvalue, method = "fdr")
	results <- cbind(mat1[,1:(gene_annot_cols-1)],results)
	#o<-order(results$rho,decreasing=TRUE)
	#results <- results[o,]
	
	return(results)
}


compute_sense_antisense_correlation <- function(groupA_for,groupA_rev,
		gene_annot_cols,
		log_data,pval_th,
		add_row_labels
){
	
	################################################################
	### INPUT
	### 1. expression matrices of sense/antisense
	### 2. gene annotation, gene_annotation_cols, 
	### 3. Tissue type or an already subset of the data.
	### 4. loci_for, loci_rev
	### 5. Need whether to use log data to perform correlations. 
	#specially important for pearson correlation
	### 6. pval_th: To determine what correlations are considered
	#  as significant correlations.
	### 7. A name for the correlation plot
	################################################################
	
	# Needs: reliable_transcripts, gene_annotation_for_reliable_prot
	#log_data=TRUE;pval_th=0.001
	#groupA_for,groupA_rev
	srho_th1=0.9; srho_th2=0.7; srho_th3=0.5;srho_th4=0.3
	srho_th5=-0.3; srho_th6=-0.5; srho_th7=-0.7; srho_th8=-0.9
	#CORRELATION_PLOT='correlation_plot_sense_antisense.pdf'
	### COMPUTE THE CORRELATION BETWEN SENSE AND ANTISENSE EXPRESSION
	sas_pair_all_linear_corr_feature <- perform_allvector_correlation(groupA_for,groupA_rev,pval_th,gene_annot_cols,log_data,add_row_labels)
	### 
	#sas_pair_all_linear_corr_feature <- rbind(sas_pair_all_linear_corr[loci_rev,],sas_pair_all_linear_corr[loci_for,])
	### TABULATE ALL CORRELATIONS OBSERVED BY THE CORRELATION VALUE
	### REPORT THIS TABLE
	rhona<-which(!is.na(sas_pair_all_linear_corr_feature$srho))
	corr_counts<-c(length(intersect(which(sas_pair_all_linear_corr_feature$srho > srho_th1),rhona)),
			length(intersect(which(sas_pair_all_linear_corr_feature$srho > srho_th2),rhona)),
			length(intersect(which(sas_pair_all_linear_corr_feature$srho > srho_th3),rhona)),
			length(intersect(which(sas_pair_all_linear_corr_feature$srho > srho_th4),rhona)),
			length(intersect(which(sas_pair_all_linear_corr_feature$srho < srho_th5),rhona)),
			length(intersect(which(sas_pair_all_linear_corr_feature$srho < srho_th6),rhona)),
			length(intersect(which(sas_pair_all_linear_corr_feature$srho < srho_th7),rhona)),
			length(intersect(which(sas_pair_all_linear_corr_feature$srho < srho_th8),rhona))
	)
	corr_counts_table <- rbind(corr_counts,(corr_counts/length(rhona))*100)
	corr_counts_table<-as.data.frame(corr_counts_table)
	
	# RETURN CORR_COUNTS_TABLE, AND sas_pair_all_linear_corr_feature
	# PRINT THE HISTOGRAM FOR DATA EXPLORATION.
	print("DONE COUNTING HOW MANY PAIRS FALL IN BIN X OF CORRELATION")
	return(list(sas_pair_all_linear_corr_feature,corr_counts_table))
	
}



get_reliable_transcripts <- function(sense_norm,asense_norm,theta){
	sense_rel <- rowMeans(sense_norm,na.rm=TRUE);
	asense_rel <- rowMeans(asense_norm,na.rm=TRUE);
	reliable_transcripts <- union(which(sense_rel >=theta),which(asense_rel >=theta));
	return(reliable_transcripts)
}


compute_pairs_correlation_in_cohorts <- function(groupA_for,groupA_rev,
		log_data,pval_th,theta,
		cohort_list,
		gene_annot_cols,
		mode_overlap,
		add_row_labels){
	
	################################################################
	### INPUT
	### 1. expression matrices of sense/antisense
	### 2. gene annotation, gene_annotation_cols, 
	### 3. List of Tissue types on which to perform correlation
	### 4. loci_for, loci_rev
	### 5. Need whether to use log data to perform correlations. 
	#specially important for pearson correlation
	### 6. pval_th: To determine what correlations are considered
	#  as significant correlations.
	### 7. A name for the correlation plot
	### 8. List of tissues
	################################################################
	total_correlations<-c()
	total_corr_counts<-c()
	for(cohort_name in names(cohort_list)){
		print(cohort_name)
		cohort=cohort_list[cohort_name][[1]]
		#
		
		# Get the reliable transcripts: transcripts with expression
		# > theta reads on average for the cohort.
		reliable_trans<-get_reliable_transcripts(groupA_for[,cohort],groupA_rev[,cohort],theta)
		if(length(reliable_trans)==0){
			print(paste(cohort_name, "DID NOT HAVE ANY RELIABLE TRANSCRIPTS FOR THIS TYPE OF OVERLAPPING",mode_overlap,sep=" "))
			next
		}else{
			print(paste(cohort_name, "HAVE",length(reliable_trans) ,"RELIABLE TRANSCRIPTS FOR THIS TYPE OF OVERLAPPING",mode_overlap,sep=" "))
		}
		sas_corr <- compute_sense_antisense_correlation(groupA_for[reliable_trans,cohort],groupA_rev[reliable_trans,cohort],
				gene_annot_cols[reliable_trans,],
				log_data,pval_th,
				add_row_labels)
		sas_correlations <- sas_corr[[1]]; corr_counts_table <-sas_corr[[2]];
		sas_correlations$cohort<-cohort_name;corr_counts_table$cohort$cohort_name;
		
		total_correlations <- rbind(total_correlations, sas_correlations)
		total_corr_counts <- rbind(total_corr_counts, corr_counts_table)
	}
	print(paste("DONE WITH PAIRS CORRELATION FOR",cohort_name,sep=' '))
	
	return(list(total_correlations,total_corr_counts))
}


corrtable_eliminate_redundancy<-function(corr_table){
	ids<-c();seenids<-c()
	for(i in seq(1,dim(corr_table)[1])){
		id<-paste(as.vector(as.matrix(corr_table[i,c("name","overlap_type","cohort")])),collapse=";");
		if(!id %in% seenids){seenids<-c(seenids,id);
		ids<-c(ids,i)}
	}
	unique_corr_table <- matrix(NA,nrow=length(ids),ncol=dim(corr_table)[2]);
	unique_corr_table <- corr_table[ids,];
	return(unique_corr_table)
}

featureov_eliminate_redundancy<-function(corr_table,thiscol=""){
	ids<-c();seenids<-c()
	for(i in seq(1,dim(corr_table)[1])){
		if(thiscol!=""){
			id<-paste(as.vector(as.matrix(corr_table[i,c("name",thiscol)])),collapse=";");	
		}
		else if("exon_overlap_type" %in% colnames(corr_table)){
			id<-paste(as.vector(as.matrix(corr_table[i,c("name","exon_overlap_type")])),collapse=";");	
		}else{
			id<-paste(as.vector(as.matrix(corr_table[i,c("name","overlap_type")])),collapse=";");
		}
		
		if(!id %in% seenids){seenids<-c(seenids,id);
			ids<-c(ids,i)}
	}
	unique_corr_table <- matrix(NA,nrow=length(ids),ncol=dim(corr_table)[2]);
	unique_corr_table <- corr_table[ids,];
	return(unique_corr_table)
}
table_eliminate_redundancy<-function(corr_table,features){
	ids<-c();seenids<-c()
	for(i in seq(1,dim(corr_table)[1])){
		id<-paste(as.vector(as.matrix(corr_table[i,features])),collapse=";");
		if(!id %in% seenids){seenids<-c(seenids,id);
			ids<-c(ids,i)}
	}
	unique_corr_table <- matrix(NA,nrow=length(ids),ncol=dim(corr_table)[2]);
	unique_corr_table <- corr_table[ids,];
	return(unique_corr_table)
}


sub_matrix_slicer<-function(groupA,genepairs,EXONS=FALSE){
	if(EXONS==TRUE){
		gids_left<-as.vector(genepairs$gene_left);gids_right<-as.vector(genepairs$gene_right);
		
	}else{
		gids_left<-as.vector(genepairs$gene_left_name);gids_right<-as.vector(genepairs$gene_right_name);
	}
	m1<-groupA[gids_left,]
	m2<-groupA[gids_right,]	
	return(list(m1,m2))
}


filter_genepairs <- function(genespairs,reliable_transcripts,EXONS){
	if(EXONS==TRUE){
		gids_left<-as.vector(genespairs$gene_left);gids_right<-as.vector(genespairs$gene_right);
		
	}else{
		gids_left<-as.vector(genespairs$gene_left_name);gids_right<-as.vector(genespairs$gene_right_name);
	}
	newind<-c()
	for(i in seq(1,length(gids_left))){
		if(length(union(intersect(gids_left[i],reliable_transcripts),intersect(gids_right[i],reliable_transcripts)))==2){
			newind<-c(newind,i)
		}
	}
	return(genespairs[newind,])
}


compute_correlation_between_genepairs<-function(expmat,feature_overlaps_annotation,
		cohort_list_names,EXONS=TRUE,
		log_data=FALSE,pval_th=0.05,theta,
		gene_annot_cols=1,
		output_dir,
		data_dir,corrtable_filename){
	#expmat[expmat==0]<-NA
	###
	add_row_labels=TRUE
	# DETERMINE OVERLAPPING TYPES TO BE CONSIDERED FOR THE CORRELATION ANALYSIS
	if("exon_overlap_type" %in% colnames(feature_overlaps_annotation)){
		overlaping_types <- as.vector(unique(feature_overlaps_annotation$exon_overlap_type)) 
	}else{
		overlaping_types <- as.vector(unique(feature_overlaps_annotation$overlap_type)) 
	}
	overlaping_types <- overlaping_types[!is.na(overlaping_types)]
	# FOR KNOWN PAIRS
	total_correlations_table <- c()
	total_corr_counts_table <- c()
	# Create the correlation matrix
	for(ov in overlaping_types){
		if("exon_overlap_type" %in% colnames(feature_overlaps_annotation)){
			guse<-which(feature_overlaps_annotation$exon_overlap_type==ov)
		}else{
			guse<-which(feature_overlaps_annotation$overlap_type==ov)
		}
		
		genespairs <- feature_overlaps_annotation[guse,];
		#genespairs <- filter_genepairs(genespairs,reliable_transcripts,EXONS)
		
		print(paste("Slicing the matrix for overlaping type ",ov," pairs = ",length(guse),sep=""))
		newmat<-sub_matrix_slicer(expmat,genespairs,EXONS)
		
		correlation_table <- compute_pairs_correlation_in_cohorts(newmat[[1]],newmat[[2]],
				log_data, pval_th,theta,
				cohort_list_names,genespairs,
				ov,
				add_row_labels)
		total_correlations_table <- rbind(total_correlations_table, correlation_table[[1]])
		total_corr_counts_table <- rbind(total_corr_counts_table, correlation_table[[2]])	
	}
	
	correlation_tables <- list(corr_table=total_correlations_table,corr_count=total_corr_counts_table)
	setwd(data_dir)
	save(correlation_tables,
			file = paste(corrtable_filename,".R",sep=""))
	
	setwd(output_dir)
	print("Writing the table of total correlations")
	write.table(total_correlations_table,file=paste(corrtable_filename,".corrtable.tsv",sep=""),row.names=TRUE, col.names=NA,sep='\t')
	write.table(total_corr_counts_table,file=paste(corrtable_filename,".counts.tsv",sep=""),row.names=TRUE, col.names=NA,sep='\t')
	
	
	return(correlation_tables)
}

compute_random_correlations <- function(expmat,pval_th,gene_annot_cols,log_data,
		data_dir,corrtable_filename){
#	# COMPUTE THE CORRELATION BETWEEN RANDOM PAIRS
	if(nrow(expmat)>1e5){random_pairs<-1e5}else{random_pairs<-nrow(expmat)}
	EXPMAT_RANDOM1 <- expmat[sample.int(random_pairs),]
	EXPMAT_RANDOM2 <- expmat[sample.int(random_pairs),]
	print(paste("Computing ",nrow(EXPMAT_RANDOM1),"random correlations",sep=""))
	NOT_NEIGHBORS_RANDOMCORR <- perform_allvector_correlation(EXPMAT_RANDOM1,EXPMAT_RANDOM2,pval_th,gene_annot_cols,log_data,add_row_labels=FALSE)
#	#NOT_NEIGHBORS_RANDOMCORR2 <- perform_allvector_correlation(EXPMAT_RANDOM1,EXPMAT_RANDOM,pval_th,"spearman",gene_annot_cols,log_data)
	print("Writing random correlations")
	write.table(NOT_NEIGHBORS_RANDOMCORR,file='total_random_correlations_table.tsv',row.names=TRUE, col.names=NA,sep='\t')
	setwd(data_dir)
	save(NOT_NEIGHBORS_RANDOMCORR,
			file = corrtable_filename)
	
	return(NOT_NEIGHBORS_RANDOMCORR)
}


compute_random_correlations_by_cohort <- function(expmat,
		pval_th,gene_annot_cols,log_data,
		theta,cohort_list_names,add_row_labels,
		data_dir,corrtable_filename){
	
	if(nrow(expmat)>1e5){random_pairs<-1e5}else{random_pairs<-nrow(expmat)}
	EXPMAT_RANDOM1 <- expmat[sample.int(random_pairs),]
	EXPMAT_RANDOM2 <- expmat[sample.int(random_pairs),]

	print(paste("Computing ",nrow(EXPMAT_RANDOM1)," random correlations",sep=""))
	NOT_NEIGHBORS_RANDOMCOR_BYCOHORT <- compute_pairs_correlation_in_cohorts(EXPMAT_RANDOM1,
			EXPMAT_RANDOM2,
			log_data, pval_th,theta,
			cohort_list_names,gene_annot_cols,
			"ALL",
			add_row_labels)
	print("Writing random correlations")
	setwd(data_dir)
	write.table(NOT_NEIGHBORS_RANDOMCOR_BYCOHORT[[1]],file=paste(corrtable_filename,'.corrbycohort.tsv',sep=""),row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
	write.table(NOT_NEIGHBORS_RANDOMCOR_BYCOHORT[[2]],file=paste(corrtable_filename,'.countsbycohort.tsv',sep=""),row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
	
	save(NOT_NEIGHBORS_RANDOMCOR_BYCOHORT,
			file = paste(corrtable_filename,'.corrbycohort.R',sep=""))
	return(NOT_NEIGHBORS_RANDOMCOR_BYCOHORT) # It is a list
}
	
	

compute_all_correlations <- function(sense_norm,asense_norm,pval_th,gene_annot_cols,log_data,
		data_dir,corrtable_filename,theta){
	reliable<-get_reliable_transcripts(sense_norm,asense_norm,theta)
	
	print(paste("Computing ",length(reliable)," correlations",sep=""))
	cohort<-colnames(sense_norm)
	ALL_CORR_ACOHORT <- perform_allvector_correlation(sense_norm[reliable,cohort],asense_norm[reliable,cohort],pval_th,gene_annot_cols,
			log_data,add_row_labels=FALSE)
	
	print("Writing random correlations")
	filename_table<-paste(corrtable_filename,"tsv",sep=".")
	filename_RObj<-paste(corrtable_filename,"R",sep=".")
	
	write.table(ALL_CORR_ACOHORT,file=filename_table,row.names=TRUE, col.names=NA,sep='\t',)
	setwd(data_dir)
	save(ALL_CORR_ACOHORT,reliable,
			file = filename_RObj)
	return(list(corr_table=ALL_CORR_ACOHORT,used_exons=reliable))
}

compute_correlations_among_overlaping_genes <- function(sense_norm,asense_norm,pval_th,gene_annot_cols,log_data,
		data_dir,corrtable_filename,theta){
	reliable<-get_reliable_transcripts(sense_norm,asense_norm,theta)
	
	print(paste("Computing ",length(reliable)," correlations",sep=""))
	cohort<-colnames(sense_norm)
	corr_overlapping_genes <- perform_allvector_correlation(sense_norm[reliable,cohort],asense_norm[reliable,cohort],pval_th,gene_annot_cols,log_data,add_row_labels=FALSE)
	
	print("Writing correlations for overlapping genes")
	filename_table<-paste(corrtable_filename,"tsv",sep=".")
	filename_RObj<-paste(corrtable_filename,"R",sep=".")
	
	write.table(corr_overlapping_genes,file=filename_table,row.names=TRUE, col.names=NA,sep='\t',)
	setwd(data_dir)
	save(corr_overlapping_genes,reliable,
			file = filename_RObj)
	return(list(corr_table=corr_overlapping_genes,used_exons=reliable))
}
