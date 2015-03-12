# TODO: Add comment
# 
# Author: alebalbin
###############################################################################


pe_original_calculation<-function(smat,amat,negpairs){
	smatsum <- apply(smat[negpairs,], 2, sum);
	asmatsum <- apply(amat[negpairs,], 2, sum);
	pe <- asmatsum/(smatsum + asmatsum);
	return(pe)
}

pe_modified_calculation<-function(smat,amat,negpairs,
		den_pad=0,min_avgcov=100){
	#Modified to be at gene level
	ncols=ncol(smat);
	# Calculate the reverse comparison
	ms<-rowMeans(smat[negpairs,],na.rm=TRUE)
	negpairs<-negpairs[which(ms>=min_avgcov)]
	
	
	total <- (smat[negpairs,] + amat[negpairs,])+den_pad;#+1)
	asratio <- amat[negpairs,]/total;
	# Compute the mean protocol error rate
	# for each sample
	pei <- colMeans(asratio, na.rm=TRUE);
	sd_pei<-apply(asratio, 2, sd,na.rm=TRUE);
	return(list(pe=pei,sd_pe=sd_pei));
}


pe_modified_calculation2<-function(smat,amat,negpairs,
		den_pad=0,min_avgcov=100){
	#Modified to be at gene level
	ncols=ncol(smat);
	# Calculate the reverse comparison
	ms<-rowMeans(smat[negpairs,],na.rm=TRUE)
	#negpairs<-negpairs[which(ms>=min_avgcov)]
	pei<-matrix(NA, nrow=3,ncol = ncol(smat))
	colnames(pei) <- colnames(smat)	
	rownames(pei) <- c('pei','sd_pei','ngenes_used')	
	
	for (sp in colnames(smat)){
		ng_sp <- negpairs[which(smat[negpairs,sp] > min_avgcov)];
		total <- smat[ng_sp,sp] + amat[ng_sp,sp] +den_pad
		asratio <- amat[ng_sp,sp]/total;
		pei['pei',sp] <- mean(asratio,na.rm=TRUE) 
		pei['sd_pei',sp] <- sd(asratio,na.rm=TRUE)
		pei['ngenes_used',sp] <- length(ng_sp)
	}
	# Compute the mean protocol error rate
	# for each sample
	return(list(pe=pei['pei',],sd_pe=pei['sd_pei',],num_genes_tested=pei['ngenes_used',]));
}



compute_genes_with_asratio_above_pei_by_sample <- function(smat,amat,cohortAnnot,
		protocolerror_bysample,
		bycohort=TRUE,min_samples=1,
		min_num_sds=2,
		min_avgcov=100,
		count_by_cohort=TRUE){
	# Compute the genes with antisense expression > than protocol error
	cohort_names<-unique(as.vector(cohortAnnot$tissuetype))
	#cohort_names<-names(cohort_list_names)
	lociForcorr_by_chort<-matrix(0,nrow=nrow(smat),ncol=length(cohort_names)+1);
	lociForcorr_by_sample <- matrix(0,nrow=nrow(smat),ncol=ncol(smat));
	rownames(lociForcorr_by_sample) <- rownames(smat);
	colnames(lociForcorr_by_sample) <- colnames(smat);
	
	lociCounts_by_sample<-matrix(0,nrow=5,ncol=ncol(smat));	
	colnames(lociCounts_by_sample) <- colnames(smat);
	rownames(lociCounts_by_sample) <- c('numloci_mincov','numloci_minasratio','numloci_mincov_asratio','pe_threshold','min_senseavgcov');
	
	lociMinCov_by_sample <- matrix(0,nrow=nrow(smat),ncol=ncol(smat));
	rownames(lociMinCov_by_sample) <- rownames(smat);
	colnames(lociMinCov_by_sample) <- colnames(smat);
	
	
	#Mean protocolerror by sample
	pei<-protocolerror_bysample$pe
	sd_pei<-protocolerror_bysample$sd_pe
	
	#Determine the protocolerror_th for each sample as 
	# pei + min_num_sds * sd_pei
	protocolerror <- pei + min_num_sds*sd_pei
	
	#Compute asratio for each sample
	tumor_cov <- compute_sas_ratio(smat,amat,den_pad=0)
	sense_ratio <- tumor_cov[[1]]
	asratio <- tumor_cov[[2]]
	total_tumor_cov <- tumor_cov[[3]]
	rm(tumor_cov)
	
	# Compare each gene asratio to the protocol error in that sample
	for(i in seq(1,length(protocolerror))){
		loci_with_mincov<- which(total_tumor_cov[,i] >= min_avgcov);
		loci_with_minasratio<- which(asratio[,i] > protocolerror[i]);
		ind<-intersect(loci_with_mincov,loci_with_minasratio)
		#
		lociCounts_by_sample['numloci_mincov',i]<-length(loci_with_mincov);
		lociCounts_by_sample['numloci_minasratio',i]<-length(loci_with_minasratio);
		lociCounts_by_sample['numloci_mincov_asratio',i]<-length(ind);
		lociCounts_by_sample['pe_threshold',i] <- protocolerror[i];
		lociCounts_by_sample['min_senseavgcov',i] <- min_avgcov;
		lociForcorr_by_sample[ind,i]<-1;
		lociMinCov_by_sample[loci_with_mincov,i]<-1
		
	}
	
	# Calculate the number of genes that pass pe at the cohort
	# level
	cohort_names2<-c(cohort_names,"ALL")
	pe_by_cohort<-matrix(NA,nrow=7,ncol=length(cohort_names2));	
	
	
	for(j in seq(1,length(cohort_names2))){
		#cohort_samples <- cohort_list_names[[cohort_names[j]]];
		cohort_samples <- as.vector(cohortAnnot$sample_name[cohortAnnot$tissuetype==cohort_names2[j]])
		if(cohort_names2[j] == "BENIGN"){
			cohort_samples <- as.vector(cohortAnnot$sample_name[cohortAnnot$tissuetype==cohort_names2[j] & cohortAnnot$sampletype=="tissue"])
		}else if(cohort_names2[j] == "ALL"){
			cohort_samples<-as.vector(cohortAnnot$sample_name)
		}
		
		#print(cohort_samples)
		if(length(cohort_samples) >= 3){
			min_samples=nsamples<-round(quantile(seq(1,length(cohort_samples)),0.05),0)
			
			inc<-which(rowSums(lociForcorr_by_sample[,cohort_samples])>=min_samples) # seen in at least 2 samples
			lociForcorr_by_chort[inc,j]<-1;
			pe_by_cohort[1,j]<-mean(protocolerror[cohort_samples])
			pe_by_cohort[2,j]<-sd(protocolerror[cohort_samples])
			pe_by_cohort[3,j]<-length(cohort_samples)
			# Mean fraction of transcripts for which asratio > pei
			# At the cohort level
			mean_numgenes2use<-mean(colSums(lociForcorr_by_sample[,cohort_samples],na.rm=TRUE))/mean(lociCounts_by_sample['numloci_mincov',cohort_samples]) #nrow(lociForcorr_by_sample))	
			pe_by_cohort[4,j]<-mean_numgenes2use;
			pe_by_cohort[5,j]<-mean(lociCounts_by_sample['numloci_mincov',cohort_samples])
			pe_by_cohort[6,j]<-mean(lociCounts_by_sample['numloci_mincov_asratio',cohort_samples])
			
			pe_by_cohort[7,j]<-length(which(rowSums(lociMinCov_by_sample[,cohort_samples])>=1))
			
		}
	}
	
	#Quantify the number of transcripts with asratio > pe_th
	#in at >= percentage of samples in the cohort.
	spquant<-c(0.01,0.05,seq(0.1,1,0.05))
	pe_by_cohort2<-matrix(NA,nrow=length(spquant)+1,ncol=length(cohort_names2));	
	for(j in seq(1,length(cohort_names2))){
		#cohort_samples <- cohort_list_names[[cohort_names[j]]];
		cohort_samples <- as.vector(cohortAnnot$sample_name[cohortAnnot$tissuetype==cohort_names2[j]])
		if(cohort_names2[j] == "BENIGN"){
			cohort_samples <- as.vector(cohortAnnot$sample_name[cohortAnnot$tissuetype==cohort_names2[j] & cohortAnnot$sampletype=="tissue"])
		}else if(cohort_names2[j] == "ALL"){
			cohort_samples<-as.vector(cohortAnnot$sample_name)
		}		
		#print(cohort_samples)
		if(length(cohort_samples) >= 3){
			#nsamples<-floor(c(1,2,3,summary(seq(1,length(cohort_samples)))[2:6]))
			nsamples<-round(quantile(seq(1,length(cohort_samples)),spquant),0)
			for(k in seq(1,length(nsamples))){
				min_samples <- nsamples[k]
				inc<-which(rowSums(lociForcorr_by_sample[,cohort_samples])>=min_samples)
				
				lociForcorr_by_sample2 <- matrix(0,nrow=nrow(lociForcorr_by_sample),ncol=length(cohort_samples));
				colnames(lociForcorr_by_sample2)<-cohort_samples
				for(sp in cohort_samples){
					ones<-which(lociForcorr_by_sample[,sp]==1)
					incones<-intersect(inc,ones)
					lociForcorr_by_sample2[incones,sp]<-1
				}
				
				if(count_by_cohort){
					#Do the counting at the cohort level.
					numloci_mincov_asratio<-which(rowSums(lociForcorr_by_sample[,cohort_samples])>=1)
					pe_by_cohort2[k,j]<-length(inc)#/length(numloci_mincov_asratio) #nrow(lociForcorr_by_sample)
				}else{
					#Do the counting at the sample level
					pe_by_cohort2[k,j]<- mean(colSums(lociForcorr_by_sample2,na.rm=TRUE))
				}
			}
			pe_by_cohort2[nrow(pe_by_cohort2),j]<-length(numloci_mincov_asratio)
			
		}
	}
	
	rownames(lociForcorr_by_chort) <- rownames(smat);
	colnames(lociForcorr_by_chort) <- toupper(cohort_names2);
	colnames(pe_by_cohort)<-toupper(cohort_names2);
	rownames(pe_by_cohort)<-c("cohort_pe","cohort_pe_sd","cohort_size","cohort_mean_percentage_GenesAbovePE",
			"cohort_numloci_mincov","cohort_numloci_mincov_asratio","Total_cohort_numloci_mincov_>=1_sample");
	rownames(pe_by_cohort2)<-c(0.01,0.05,seq(0.1,1,0.05),"cohort_numloci_mincov_asratio_>=1_sample");
	colnames(pe_by_cohort2)<-toupper(cohort_names2);
	pe_by_cohort_all<-rbind(pe_by_cohort,pe_by_cohort2);
	
	return(list(asloci_bycohort=lociForcorr_by_chort,
					summary_bycohort=pe_by_cohort_all,
					asloci_bysample=lociForcorr_by_sample,
					summary_by_sample=lociCounts_by_sample));
}



simulation_effect_sas_threshold <- function(lociCounts_by_sample,
		smat, asmat, 
		min_avgcov=100,
		min_sp_th=0.05,
		number_simulations=100){
	
	cohort_samples<-colnames(smat)
	#simulation <- vector("numeric",length=number_simulations);
	
	all_loci <- rownames(smat);
	#min_nsamples<-round(ncol(lociForcorr_by_sample)*min_sp_th,0)
	#reliable_transcripts <- get_reliable_transcripts(smat,asmat,theta)
	
	tumor_cov <- compute_sas_ratio(smat,asmat,den_pad=0)
	sense_ratio <- tumor_cov[[1]]
	asratio <- tumor_cov[[2]]
	total_tumor_cov <- tumor_cov[[3]]
	rm(tumor_cov)
	
	#Quantify the number of transcripts with asratio > pe_th
	#in at >= percentage of samples in the cohort.		
	spquant<-c(0.01,0.05,seq(0.1,1,0.05))
	simulation<-matrix(NA,nrow=number_simulations,ncol=length(spquant));
	nsamples<-round(quantile(seq(1,length(cohort_samples)),spquant),0)
	
	for (i in seq(1,number_simulations)){
		print(number_simulations - i)
		# Initialize the matrix that contains the loci > pe_i for each simulation
		lociForcorr_by_sample <- matrix(0,nrow=nrow(smat),ncol=ncol(smat));
		rownames(lociForcorr_by_sample) <- rownames(smat);
		colnames(lociForcorr_by_sample) <- colnames(smat);
		
		for(sample_name in colnames(lociCounts_by_sample)){
			pe_i <- lociCounts_by_sample["pe_threshold",sample_name];
			num_loci2select <- lociCounts_by_sample["numloci_minasratio",sample_name];
			
			# select randomly the number of transcripts
			rand_loci <- sample(all_loci,num_loci2select,replace=FALSE);
			loci_with_minasratio<- all_loci[which(asratio[rand_loci,sample_name] > pe_i)];		
			# Control for the minimun coverage
			loci_with_mincov<- all_loci[which(total_tumor_cov[rand_loci,sample_name] >= min_avgcov)];
			# Interect loci with minimum coverage in the sample and OPS>pe_i in the sample 
			ind_rand_loci<-intersect(loci_with_mincov,loci_with_minasratio)
			
			#Store in the matrix
			lociForcorr_by_sample[ind_rand_loci,sample_name]<-1;
			
		}
		
		#Quantify the number of transcripts with asratio > pe_th
		#in at >= percentage of samples in the cohort.		
		loci_present_nsamples<-rowSums(lociForcorr_by_sample)
		for(k in seq(1,length(nsamples))){
			min_samples <- nsamples[k]
			ASloci_above_noise<-which(loci_present_nsamples >= min_samples)
			
			simulation[i,k]<-length(ASloci_above_noise)#/nrow(lociForcorr_by_sample)
			
		}		
	}
	colnames(simulation)<-c(0.01,0.05,seq(0.1,1,0.05));
	return(simulation)
}


compute_asexp_for_genes_below200reads <- function(smat,amat,
		min_avgcov=200){
	#
	#
	#Compute asratio for each sample
	tumor_cov <- compute_sas_ratio(smat,amat,den_pad=0)
	sense_ratio <- tumor_cov[[1]]
	asratio <- tumor_cov[[2]]
	total_tumor_cov <- tumor_cov[[3]]
	rm(tumor_cov)
	
	lociCounts_by_sample<-matrix(0,nrow=6,ncol=ncol(smat));
	
	colnames(lociCounts_by_sample) <- colnames(smat);
	rownames(lociCounts_by_sample) <- c('numloci_above150asreads','numloci_above100asreads','numloci_above50asreads','numloci_above10asreads',
			'numloci_above2asreads','numloci_0asreads');
	
	for(i in seq(1,ncol(smat))){
		loci_with_mincov<- which(total_tumor_cov[,i] < min_avgcov);
		#
		lociCounts_by_sample['numloci_above150asreads',i]<-length(which(amat[loci_with_mincov,i] > 150));
		lociCounts_by_sample['numloci_above100asreads',i]<-length(which(amat[loci_with_mincov,i] > 100));
		lociCounts_by_sample['numloci_above50asreads',i]<-length(which(amat[loci_with_mincov,i] > 50)); 
		lociCounts_by_sample['numloci_above10asreads',i] <- length(which(amat[loci_with_mincov,i] > 10));
		lociCounts_by_sample['numloci_above2asreads',i] <- length(which(amat[loci_with_mincov,i] > 2));
		lociCounts_by_sample['numloci_0asreads',i]<-length(which(amat[loci_with_mincov,i] > 0));
		
	}
	return(lociCounts_by_sample)
}


loci_min_coverate<-function(x,min_avgcov){
	loci_with_mincov<- length(which(x < min_avgcov));
	return(loci_with_mincov)
}
