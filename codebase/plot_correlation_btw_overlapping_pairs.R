# TODO: Add comment
# 
# Author: alebalbin
###############################################################################
# Aug/2013
# Plotting all tissues in a density plot.
###############################################################################
#EXONS
#plot_correlation_btw_pairs_exons(correlation_tables$corr_table, NOT_NEIGHBORS_RANDOMCORR)
#plot_correlation_btw_pairs_exons(correlation_tables$corr_table, random_correlations_bycohort)
#corr_summary_tables <- summarize_correlation_btw_pairs_exons(correlation_tables$corr_table, NOT_NEIGHBORS_RANDOMCORR,EXONS)
#corr_summary_tables <- summarize_correlation_btw_pairs_exons(correlation_tables$corr_table, random_correlations_bycohort,EXONS)
#write.table(corr_summary_tables[[1]],file='correlation_summary_table.tsv',row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
#write.table(corr_summary_tables[[2]],file='correlation_pairstype_comparison_table.tsv',row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)

#TRANSCRIPTS
#plot_correlation_btw_pairs_trans(correlation_tables$corr_table, random_correlations_bycohort)
#corr_summary_tables <- summarize_correlation_btw_pairs_exons(correlation_tables$corr_table, random_correlations_bycohort,EXONS)
#write.table(corr_summary_tables[[1]],file='correlation_summary_table.tsv',row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
#write.table(corr_summary_tables[[2]],file='correlation_pairstype_comparison_table.tsv',row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)
#plot_correlations_allcohorts(correlation_tables$corr_table,random_correlations_bycohort,EXONS=TRUE,"correlation_plot_alltrans_bycohort.pdf")
#plot_correlations_allcohorts(all_corr_acohort[[1]],random_correlations_acohort,EXONS=TRUE,"correlation_plot_alltrans_across_cohort.pdf")

#TO PLOT PARTICULAR EXONS/TRANSCRIPTS WITH ABS(CORRELATION)> TH.
#plot_correlated_pairs(correlation_tables$corr_table, sense_norm, corr_th=-0.5,cohort_list_names,absv=FALSE,EXONS=TRUE)
#plot_correlated_pairs(correlation_tables$corr_table, sense_norm, corr_th=0.8,cohort_list_names,absv=FALSE,EXONS=TRUE)
# To plot the correlation for all exons regardless of overlaping
# Results were save to /exds/users/oabalbin/projects/sascompendia/version2/analysis/loci_proportions/all_exons
#plot_correlated_sense_asense(all_correlation_table[[1]], sense_norm, asense_norm, corr_th=-0.5,cohort_list_names2,absv=FALSE,EXONS=TRUE,scaled=FALSE)

plot_correlations_allcohorts<-function(corr_table,rand_corr_table,
		EXONS=TRUE,figname){
	notna<-which(!is.na(as.numeric(as.vector(corr_table$srho))))
	srho<-as.numeric(as.vector(corr_table$srho))
	# For the random matrix
	notna_random<-which(!is.na(as.numeric(as.vector(rand_corr_table$srho))))
	srho_random<-as.numeric(as.vector(rand_corr_table$srho))
	srhorandom<-srho_random[notna_random];
	if(length(notna) < length(srhorandom)){srhorandom<-srhorandom[sample.int(length(notna))]}
	XMIN=-1;XMAX=1;YMIN=0;YMAX=2.0;
	
	pdf(figname)
	plot(density(srho[notna]),xlim=c(XMIN,XMAX),ylim=c(YMIN,YMAX),
			col="black",lwd=3,main=paste("median",median(srho[notna]),sep="="))
	lines(density(srhorandom),xlim=c(XMIN,XMAX),
			col="gray",lty=2,main=paste("median",median(srho[notna]),sep="="))
	abline(v=median(srhorandom),col="gray",lty=2)
	abline(v=median(srho[notna]),col="black",lty=2)
	
	dev.off()	
}

summarize_correlation_btw_pairs_exons<-function(total_correlations_table,NOT_NEIGHBORS_RANDOMCORR,
		EXONS=TRUE){

	#if(!EXONS){}
#	conv<-which(total_correlations_table$overlap_type=="convergent")
#	diverg<-which(total_correlations_table$overlap_type=="divergent")
#	ups<-which(total_correlations_table$overlap_type=="upstream")
#	downs<-which(total_correlations_table$overlap_type=="downstream")
	MIN_NUM_OBS=10
	SM<-matrix(NA,nrow=13,ncol=6)
	notna<-which(!is.na(as.numeric(as.vector(total_correlations_table$srho))))
	srho<-as.numeric(as.vector(total_correlations_table$srho))
	# For the random matrix
	notna_random<-which(!is.na(as.numeric(as.vector(NOT_NEIGHBORS_RANDOMCORR$srho))))
	srho_random<-as.numeric(as.vector(NOT_NEIGHBORS_RANDOMCORR$srho))
	srhorandom<-srho_random[notna_random];srhorandom<-srhorandom[sample.int(length(notna))]
	
	overlap_types=c("3UTR_3UTR","3UTR_EXON","EXON_3UTR",
			"5UTR_5UTR","5UTR_EXON","EXON_5UTR",
			"EXON_EXON","SINGLE_EXON_EXON",
			"convergent","divergent","upstream","downstream","random")
	overlap_types_ind<-list()
	for(i in seq(1,length(overlap_types))){
		ov<-overlap_types[i]
		if(ov=="random"){
			c(SM[i,1],SM[i,2],SM[i,3],SM[i,4],SM[i,5],SM[i,6]) := summary(srhorandom)
			overlap_types_ind[[ov]]<-seq(1,length(srhorandom))
			next};
		ind<-intersect(notna,which(total_correlations_table$overlap_type==ov));
		overlap_types_ind[[ov]]<-ind
		if(length(ind>0)){
			c(SM[i,1],SM[i,2],SM[i,3],SM[i,4],SM[i,5],SM[i,6]) := summary(srho[ind])
		}
	}
	rownames(SM)<-overlap_types;colnames(SM)<-c("MIN","Q25","Q50","AVG","Q75","MAX")
	# COMPUTE T.TEST BETWEEN DIFFERENT TYPES OF OVERLAP
	TT<-matrix(NA, nrow=length(overlap_types),ncol=length(overlap_types))
	rownames(TT)<-overlap_types;colnames(TT)<-overlap_types;
	overlap_types_comp<-list()
	for(ov1 in names(overlap_types_ind)){
		for(ov2 in names(overlap_types_ind)){
			comp_name<-paste(ov1,"_vs_",ov2,sep="")
			ind1<-overlap_types_ind[[ov1]];ind2<-overlap_types_ind[[ov2]];
			if(length(ind1)>MIN_NUM_OBS & length(ind2)>MIN_NUM_OBS){
				if(ov1=="random"){
					tt<-t.test(srhorandom[ind1],srho[ind2],alternative="greater");
					TT[ov1,ov2]=pval=tt[[3]];
					overlap_types_comp[[comp_name]]<-list(pval=tt[[3]],n1=length(ind1),n2=length(ind2));
				}else if(ov2=="random"){
					tt<-t.test(srho[ind1],srhorandom[ind2],alternative="greater");
					TT[ov1,ov2]=pval=tt[[3]];
					overlap_types_comp[[comp_name]]<-list(pval=tt[[3]],n1=length(ind1),n2=length(ind2));
				}else if (ov1!="random" & ov2!="random"){
				tt<-t.test(srho[ind1],srho[ind2],alternative="greater");
				TT[ov1,ov2]=pval=tt[[3]];
				overlap_types_comp[[comp_name]]<-list(pval=tt[[3]],n1=length(ind1),n2=length(ind2));
			}else{
				overlap_types_comp[[comp_name]]<-list(pval=NA,n1=length(ind1),n2=length(ind2));
			}
		}
	  }
   }
	return(list(sumtable=SM,comptable=TT, comp_table_list=overlap_types_comp))
 }

plot_correlation_btw_pairs_exons<-function(total_correlations_table,
		NOT_NEIGHBORS_RANDOMCORR){
	
	# PLOT PARAMETERS
	XMIN=-1;XMAX=1
	YMIN=0;YMAX=2.0
	##################################################################
	notna<-which(!is.na(as.numeric(as.vector(total_correlations_table$srho))))
	
	utr3s<-which(total_correlations_table$overlap_type=="3UTR_3UTR")
	utr3s_exon<-which(total_correlations_table$overlap_type=="3UTR_EXON")
	exon_utr3s<-which(total_correlations_table$overlap_type=="EXON_3UTR")
	utr5s<-which(total_correlations_table$overlap_type=="5UTR_5UTR")
	utr5s_exon<-which(total_correlations_table$overlap_type=="5UTR_EXON")
	exon_utr5s<-which(total_correlations_table$overlap_type=="EXON_5UTR")

	exon_exon<-which(total_correlations_table$overlap_type=="EXON_EXON")
	single_exon<-which(total_correlations_table$overlap_type=="SINGLE_EXON_EXON")
		
	luadc<-which(total_correlations_table$cohort=="luad")
	luscc<-which(total_correlations_table$cohort=="lusc")
	brcac<-which(total_correlations_table$cohort=="brca")
	prcac<-which(total_correlations_table$cohort=="prca")
	luclc<-which(total_correlations_table$cohort=="lucl")
	srho<-as.numeric(as.vector(total_correlations_table$srho))
	
	# Random not neighboring genes
	notna_random<-which(!is.na(as.numeric(as.vector(NOT_NEIGHBORS_RANDOMCORR$srho))))
	srho_random<-as.numeric(as.vector(NOT_NEIGHBORS_RANDOMCORR$srho))
	
	################################################################
	### PLOT THE OVERLAPING GENES

	pdf('UTRs3_alltissues.pdf')
	plot(density(srho[intersect(notna,intersect(utr3s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main="Overlapping 3UTR_3UTR genes",xlab="")
	lines(density(srho[intersect(notna,intersect(utr3s,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(utr3s,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(utr3s,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	lines(density(srho[intersect(notna,intersect(utr3s,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('3UTR_EXON_alltissues.pdf')
	plot(density(srho[intersect(notna,intersect(utr3s_exon,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main="Overlapping 3UTR_EXON genes",xlab="")
	lines(density(srho[intersect(notna,intersect(utr3s_exon,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(utr3s_exon,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(utr3s_exon,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(utr3s_exon,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('EXON_3UTR_alltissues.pdf')
	plot(density(srho[intersect(notna,intersect(exon_utr3s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main="Overlapping EXON_3UTR genes",xlab="")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))	
	dev.off()
	
	pdf('UTR5s_alltissues.pdf')
	plot(density(srho[intersect(notna,intersect(utr5s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main="Overlapping 5UTR_5UTR genes",xlab="")
	lines(density(srho[intersect(notna,intersect(utr5s,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(utr5s,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(utr5s,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(utr5s,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('5UTR_EXON_alltissues.pdf')
	plot(density(srho[intersect(notna,intersect(utr5s_exon,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main="Overlapping 5UTR_EXON genes",xlab="")
	lines(density(srho[intersect(notna,intersect(utr5s_exon,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(utr5s_exon,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(utr5s_exon,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(utr5s_exon,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('EXON_5UTR_alltissues.pdf')
	plot(density(srho[intersect(notna,intersect(exon_utr5s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main="Overlapping 5UTR_EXON genes",xlab="")
	lines(density(srho[intersect(notna,intersect(exon_utr5s,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(exon_utr5s,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(exon_utr5s,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(exon_utr5s,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('EXON_EXON_alltissues.pdf')
	plot(density(srho[intersect(notna,intersect(exon_exon,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon",main="Overlapping EXON_EXON genes",xlab="")
	lines(density(srho[intersect(notna,intersect(exon_exon,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(exon_exon,luadc))]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,intersect(exon_exon,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(exon_exon,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()

	pdf('SINGLE_EXON_EXON_alltissues.pdf')
	plot(density(srho[intersect(notna,intersect(single_exon,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main="Overlapping SINGLE_EXON_EXON genes",xlab="")
	lines(density(srho[intersect(notna,intersect(single_exon,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(single_exon,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(single_exon,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(single_exon,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	################################################################
	### PLOT ONLY THE GENE RELATIONSHIPS WITHOUT REGARDLESS OF TISSUE
	pdf('gene_relations_alltissues_with_null.pdf')
	srhorandom<-srho_random[notna_random];
	if(length(notna) < length(srhorandom)){srhorandom<-srhorandom[sample.int(length(notna))]}
	
	plot(density(srho_random[notna_random]),xlim=c(XMIN,XMAX),ylim=c(YMIN,YMAX),
			col="gray",lty=2,main="")
    #pdf('gene_relations_alltissues_with_null.pdf')
	lines(density(srho[intersect(notna,utr5s)]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,utr3s)]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,utr3s_exon)]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,exon_utr3s)]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,utr5s)]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,utr5s_exon)]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,exon_utr5s)]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,exon_exon)]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,single_exon)]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	abline(v=median(srhorandom),col="gray",lty=2)
	legend("topleft",c("neighbors not overlap","3UTR-overlap","5UTR-overlap","Nested","Random"),lty=c(1,1,1,1,2),
			col=c("gray","blue","salmon","forestgreen","black"))
	dev.off()
}

plot_correlation_btw_pairs_trans<-function(total_correlations_table,
		NOT_NEIGHBORS_RANDOMCORR){
	
	# PLOT PARAMETERS
	XMIN=-1;XMAX=1
	YMIN=0;YMAX=2
	
	##################################################################
	notna<-which(!is.na(as.numeric(as.vector(total_correlations_table$srho))))
	
	# Plotting all tissues in a density plot.
	notna<-which(!is.na(as.numeric(total_correlations_table$srho)))
	conv<-which(total_correlations_table$overlap_type=="convergent")
	diverg<-which(total_correlations_table$overlap_type=="divergent")
	ups<-which(total_correlations_table$overlap_type=="upstream")
	downs<-which(total_correlations_table$overlap_type=="downstream")
	
	utr3s<-which(total_correlations_table$overlap_type=="3UTR_3UTR")
	utr3s_exon<-which(total_correlations_table$overlap_type=="3UTR_EXON")
	exon_utr3s<-which(total_correlations_table$overlap_type=="EXON_3UTR")
	#T2T <- unique(c(utr3s,utr3s_exon,utr3s_exon))
	
	utr5s<-which(total_correlations_table$overlap_type=="5UTR_5UTR")
	utr5s_exon<-which(total_correlations_table$overlap_type=="5UTR_EXON")
	exon_utr5s<-which(total_correlations_table$overlap_type=="EXON_5UTR")
	#H2H<-unique(c(utr5s,utr5s_exon,exon_utr5s))
	
	exon_exon<-which(total_correlations_table$overlap_type=="EXON_EXON")
	single_exon<-which(total_correlations_table$overlap_type=="SINGLE_EXON_EXON")
	#EMB<- unique(c(exon_exon,single_exon))
	#print(c(length(EMB),length(H2H),length(T2T)))
	
	
	#Plot the relationships by general types. H2H,T2T, EMB
	intronic<-which(total_correlations_table$overlap_type=="INTRONIC")
	EMB<- setdiff(which(total_correlations_table$trans_overlap_type=="EMB"),intronic)
	H2H<- setdiff(which(total_correlations_table$trans_overlap_type=="H2H"),intronic)
	T2T<- setdiff(which(total_correlations_table$trans_overlap_type=="T2T"),intronic)
	
	print(c(length(EMB),length(H2H),length(T2T)))
	
	luadc<-which(total_correlations_table$cohort=="luad")
	luscc<-which(total_correlations_table$cohort=="lusc")
	brcac<-which(total_correlations_table$cohort=="brca")
	prcac<-which(total_correlations_table$cohort=="prca")
	luclc<-which(total_correlations_table$cohort=="lucl")
	srho<-as.numeric(as.vector(total_correlations_table$srho))
	
	# Random not neighboring genes
	notna_random<-which(!is.na(as.numeric(as.vector(NOT_NEIGHBORS_RANDOMCORR$srho))))
	srho_random<-as.numeric(as.vector(NOT_NEIGHBORS_RANDOMCORR$srho))
	
	################################################################
	### PLOT THE OVERLAPING GENES. GENERAL ORIENTATIONS
	pdf('gene_GRAL_relations_alltissues_with_null.pdf')
	srhorandom<-srho_random[notna_random];
	if(length(notna) < length(srhorandom)){srhorandom<-srhorandom[sample.int(length(notna))]}
	
	plot(density(srho_random[notna_random]),xlim=c(XMIN,XMAX),ylim=c(YMIN,YMAX),
			col="gray",lty=2,main="")
	#pdf('gene_relations_alltissues_with_null.pdf')
	lines(density(srho[intersect(notna,H2H)]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,T2T)]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,EMB)]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intronic)]),xlim=c(XMIN,XMAX),
			col="purple")
	abline(v=median(srhorandom),col="gray",lty=2)
	legend("topleft",c("Random","Tail-to-Tail","Head-to_Head","Embeded","Intronic"),lty=c(2,1,1,1,1),
			col=c("gray","blue","salmon","forestgreen","purple"))
	dev.off()
	
	# TAIL TO TAIL
	pdf('T2T_alltissues.pdf')
	med<-median(srho[intersect(notna,T2T)])
	plot(density(srho[intersect(notna,intersect(utr3s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping tail to tail genes, ",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(T2T,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(T2T,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(T2T,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	lines(density(srho[intersect(notna,intersect(T2T,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	abline(v=median(srho[intersect(notna,utr3s)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	# HEAD TO HEAD
	pdf('H2H_alltissues.pdf')
	med<-median(srho[intersect(notna,H2H)])
	plot(density(srho[intersect(notna,intersect(H2H,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping head to head genes, ",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(H2H,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(H2H,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(H2H,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	lines(density(srho[intersect(notna,intersect(H2H,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	abline(v=median(srho[intersect(notna,H2H)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()

	# HEAD TO HEAD
	pdf('EMB_alltissues.pdf')
	med<-median(srho[intersect(notna,EMB)])
	plot(density(srho[intersect(notna,intersect(EMB,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping embeded genes, ",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(EMB,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(EMB,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(EMB,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	lines(density(srho[intersect(notna,intersect(EMB,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	abline(v=median(srho[intersect(notna,EMB)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	# HEAD TO HEAD
	pdf('INTRONIC_alltissues.pdf')
	med<-median(srho[intersect(notna,intronic)])
	plot(density(srho[intersect(notna,intersect(intronic,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping embeded genes, ",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(intronic,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(intronic,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(intronic,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	lines(density(srho[intersect(notna,intersect(intronic,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	abline(v=median(srho[intersect(notna,intronic)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	
	################################################################
	### PLOT THE OVERLAPING GENES BY FINE TYPE OF OVERLAP. 
	### WHAT ARE THE EXONS INCLUDED
	
	pdf('UTRs3_alltissues.pdf')
	med<-median(srho[intersect(notna,utr3s)])
	plot(density(srho[intersect(notna,intersect(utr3s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping 3UTR_3UTR genes, ",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(utr3s,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(utr3s,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(utr3s,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	lines(density(srho[intersect(notna,intersect(utr3s,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	abline(v=median(srho[intersect(notna,utr3s)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('3UTR_EXON_alltissues.pdf')
	med<-median(srho[intersect(notna,utr3s_exon)])
	plot(density(srho[intersect(notna,intersect(utr3s_exon,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping 3UTR_EXON genes",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(utr3s_exon,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(utr3s_exon,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(utr3s_exon,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(utr3s_exon,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	abline(v=median(srho[intersect(notna,utr3s_exon)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('EXON_3UTR_alltissues.pdf')
	med<-median(srho[intersect(notna,exon_utr3s)])
	plot(density(srho[intersect(notna,intersect(exon_utr3s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping EXON_3UTR genes",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,intersect(exon_utr3s,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	abline(v=median(srho[intersect(notna,exon_utr3s)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))	
	dev.off()
	
	pdf('UTR5s_alltissues.pdf')
	med<-median(srho[intersect(notna,utr5s)])
	plot(density(srho[intersect(notna,intersect(utr5s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping 5UTR_5UTR genes",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(utr5s,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(utr5s,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(utr5s,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(utr5s,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	abline(v=median(srho[intersect(notna,utr5s)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('5UTR_EXON_alltissues.pdf')
	med<-median(srho[intersect(notna,utr5s_exon)])
	plot(density(srho[intersect(notna,intersect(utr5s_exon,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping 5UTR_EXON genes",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(utr5s_exon,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(utr5s_exon,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(utr5s_exon,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(utr5s_exon,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	abline(v=median(srho[intersect(notna,utr5s_exon)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('EXON_5UTR_alltissues.pdf')
	med<-median(srho[intersect(notna,exon_utr5s)])
	plot(density(srho[intersect(notna,intersect(exon_utr5s,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping 5UTR_EXON genes",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(exon_utr5s,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(exon_utr5s,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(exon_utr5s,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(exon_utr5s,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	abline(v=median(srho[intersect(notna,exon_utr5s)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('EXON_EXON_alltissues.pdf')
	med<-median(srho[intersect(notna,exon_exon)])
	plot(density(srho[intersect(notna,intersect(exon_exon,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon",main=paste("Overlapping EXON_EXON genes",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(exon_exon,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(exon_exon,luadc))]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,intersect(exon_exon,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(exon_exon,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	abline(v=median(srho[intersect(notna,exon_exon)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	pdf('SINGLE_EXON_EXON_alltissues.pdf')
	med<-median(srho[intersect(notna,single_exon)])
	plot(density(srho[intersect(notna,intersect(single_exon,luadc))]),xlim=c(XMIN,XMAX),
			col="blue",main=paste("Overlapping SINGLE_EXON_EXON genes",med,sep=""),xlab="")
	lines(density(srho[intersect(notna,intersect(single_exon,prcac))]),xlim=c(XMIN,XMAX),
			col="orange")
	lines(density(srho[intersect(notna,intersect(single_exon,luscc))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,intersect(single_exon,brcac))]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,intersect(single_exon,luclc))]),xlim=c(XMIN,XMAX),
			col="gray")
	abline(v=median(srho[intersect(notna,single_exon)]),col="gray",lty=2)
	legend("topleft",c("luad","lusc","lucl","brca","prca"),lty=c(1,1,1,1,1),
			col=c("blue","forestgreen","gray","salmon","orange"))
	dev.off()
	
	################################################################
	### PLOT ONLY THE GENE RELATIONSHIPS WITHOUT REGARDLESS OF TISSUE
	srhorandom<-srho_random[notna_random];
	if(length(notna) < length(srhorandom)){srhorandom<-srhorandom[sample.int(length(notna))]}
	#print(length(srhorandom))
	pdf('gene_relations_alltissues_with_null.pdf')
	plot(density(srhorandom),xlim=c(XMIN,XMAX),ylim=c(YMIN,YMAX),
			col="gray",lty=5,main="")
	#Plot down and upstream genes
	lines(density(srho[intersect(notna,conv)]),xlim=c(XMIN,XMAX),
			col="gray")#,lty=3)
	lines(density(srho[intersect(notna,diverg)]),xlim=c(XMIN,XMAX),
			col="gray")#,lty=4)
	#lines(density(srho[intersect(notna,ups)]),xlim=c(XMIN,XMAX),
	#		col="gray")
	#lines(density(srho[intersect(notna,downs)]),xlim=c(XMIN,XMAX),
	#		col="gray")
	#Plot overlaping genes
	lines(density(srho[intersect(notna,utr5s)]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,utr3s)]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,utr3s_exon)]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,exon_utr3s)]),xlim=c(XMIN,XMAX),
			col="blue")
	lines(density(srho[intersect(notna,utr5s)]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,utr5s_exon)]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,exon_utr5s)]),xlim=c(XMIN,XMAX),
			col="salmon")
	lines(density(srho[intersect(notna,exon_exon)]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	lines(density(srho[intersect(notna,single_exon)]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	abline(v=median(srhorandom),col="gray",lty=2)
	legend("topleft",c("neighbors not overlap","3UTR-overlap","5UTR-overlap","Nested","Random"),lty=c(1,1,1,1,2),
			col=c("gray","blue","salmon","forestgreen","black"))
	dev.off()
	
	notna<-which(!is.na(as.numeric(total_correlations_table$srho)))
	nonoverlap<-unique(c(conv,diverg,ups,downs))
	non_utr5 <- unique(c(utr3s,utr3s_exon,exon_utr3s,exon_exon,single_exon))
	
	prot2as<-which(total_correlations_table$gene_biotype=="protein_coding&antisense")
	as2prot<-which(total_correlations_table$gene_biotype=="antisense&protein_coding")
	prot2prot<-which(total_correlations_table$gene_biotype=="protein_coding&protein_coding")
	
	# PLOTTING ONLY CORRELATIONS BETWEEN PROTEIN CODING GENES AND ANTISENSE GENES
	pdf('all_protcod_asense_correlations.pdf')
	plot(density(srho[intersect(notna,setdiff(setdiff(prot2prot,nonoverlap),non_utr5))]),xlim=c(XMIN,XMAX),
			col="black",main="protcod&antisense, protcod&protcod all overlpaings")	
	lines(density(srho[intersect(notna,setdiff(setdiff(prot2as,nonoverlap),non_utr5))]),xlim=c(XMIN,XMAX),
			col="blue")	
	lines(density(srho[intersect(notna,setdiff(setdiff(as2prot,nonoverlap),non_utr5))]),xlim=c(XMIN,XMAX),
			col="forestgreen")
	abline(v=median(srho_random[notna_random]),col="black",lty=2)
	legend("topleft",c("protcod-protcod", "protcod-asense","asense vs protcod"),
			lty=c(1,1,1),
			col=c("black","blue","forestgreen"))
	dev.off()
	
}

#####
## Repeat the plot for 3UTR Correlation
#####

#notna<-which(!is.na(as.numeric(total_correlations_table$srho)))
#nonoverlap<-unique(c(conv,diverg,ups,downs))
#non_utr3 <- unique(c(utr5s,utr5s_exon,exon_utr5s,exon_exon,single_exon))
#
#prot2as<-which(total_correlations_table$mode_gene_biotype=="protein_coding&antisense")
#as2prot<-which(total_correlations_table$mode_gene_biotype=="antisense&protein_coding")
#prot2prot<-which(total_correlations_table$mode_gene_biotype=="protein_coding&protein_coding")
#
## PLOTTING ONLY CORRELATIONS BETWEEN PROTEIN CODING GENES AND ANTISENSE GENES
#pdf('ALL_protcod_asense_correlations.pdf')
#plot(density(srho[intersect(notna,setdiff(prot2prot,nonoverlap))]),xlim=c(XMIN,XMAX),
#		col="black",main="protcod&antisense, protcod&protcod all overlpaings")
#lines(density(srho[intersect(notna,setdiff(prot2as,nonoverlap))]),xlim=c(XMIN,XMAX),
#		col="blue")
#lines(density(srho[intersect(notna,setdiff(as2prot,nonoverlap))]),xlim=c(XMIN,XMAX),
#		col="forestgreen")
#abline(v=median(srho_random[notna_random]),col="black",lty=2)
#legend("topleft",c("protcod-protcod", "protcod-asense","asense vs protcod"),
#		lty=c(1,1,1),
#		col=c("black","blue","forestgreen"))
#dev.off()
#
#pdf('UTR3_protcod_asense_correlations.pdf')
#plot(density(srho[intersect(notna,setdiff(setdiff(prot2as,nonoverlap), non_utr3))]),xlim=c(XMIN,XMAX),
#		col="blue",lwd=3,main="protcod&antisense, protcod&protcod only 3UTR")
#lines(density(srho[intersect(notna,setdiff(setdiff(prot2prot,nonoverlap),non_utr3))]),xlim=c(XMIN,XMAX),
#		col="black")
#lines(density(srho[intersect(notna,setdiff(setdiff(as2prot,nonoverlap),non_utr3))]),xlim=c(XMIN,XMAX),
#		col="forestgreen")
#abline(v=median(srho_random[notna_random]),col="black",lty=2)
#legend("topleft",c("protcod-protcod", "protcod3UTR-asense","asense3UTR vs protcod"),
#		lty=c(1,1,1),
#		col=c("black","blue","forestgreen"))
#dev.off()
#
#### SELECT THE 3UTR-EXON PROTEIN&ANTISENSE CANDIDATES
#UTR3_cand<-total_correlations_table[intersect(notna,setdiff(setdiff(prot2as,nonoverlap), non_utr3)),]
#UTR3_cand <- UTR3_cand[which(as.numeric(UTR3_cand$srho) > median(as.numeric(UTR3_cand$srho))),]
#write.table(UTR3_cand,file="UTR3_protasense.tsv",row.names=TRUE,col.names=TRUE,sep="\t")

