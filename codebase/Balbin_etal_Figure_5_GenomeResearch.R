# TODO: Add comment
# 
# Author: alebalbin
# This script produces the results presented in Figure 4 of the 
# Balbin et al The landscape of antisense gene expression in human cancers
# Genome Research manuscript
###############################################################################
###############################################################################
### This scripts depend on the objects created with the scripts in 
### Balbin_etal_Figure_4_GenomeResearch.R e.g. all_natgenes
###############################################################################

load(paste(data_dir,'all_natgenes_objects.R',sep="/"))
all_natgenes<-all_natgenes_objects[[1]]
all_natgenes<-all_natgenes[all_natgenes$numcohorts!=0 & all_natgenes$annoterror==0 & all_natgenes$OPSnumcohorts !=0 ,]

###############################################################################
### code chunk number 1: Set up results directories
###############################################################################

results_dir_fig5 <- paste(results_dir,'figure5',sep='/')
if(!file.exists(results_dir_fig5)){
	system(paste("mkdir",results_dir_fig5,sep=" "))
}

###############################################################################
### code chunk number 2: Determine potential lung cancer speficific genes
###############################################################################
cohort_names<-c("luad","lusc","lucl","benign","brca","prca","panc","meningioma") #"ovarian",
cancer_names<-c("luad","lusc","lucl","brca","prca","panc","meningioma")
totalnumcohorts=length(cancer_names)
		
all_natgenes_tissuelung<- rowSums(all_natgenes[,c("luad","lusc")],na.rm=TRUE);
all_natgenes_lung<- rowSums(all_natgenes[,c("luad","lusc","lucl")],na.rm=TRUE);
all_natgenes_othertissues <- rowSums(all_natgenes[,c("brca","prca","panc","meningioma")],na.rm=TRUE);

all_natgenes_luad_lusc <- all_natgenes[all_natgenes$benign==0 & (all_natgenes$luad==1 || all_natgenes$lusc==1),]
all_natgenes_lung_tissuespecific <- all_natgenes[all_natgenes$benign==0 & all_natgenes_tissuelung == 2 & all_natgenes_othertissues!=0,]
all_natgenes_lung_canspecific <- all_natgenes[all_natgenes$benign==0 & all_natgenes_tissuelung == 2 & all_natgenes_lung ==3 & all_natgenes_othertissues!=0,]
all_natgenes_unique_lung <- all_natgenes[all_natgenes$benign==0 & all_natgenes_tissuelung == 2 & all_natgenes_othertissues==0,]


figure5_summarytable<-matrix(NA,nrow=4,ncol=3)
rownames(figure5_summarytable)<-c("luad_lusc_specific","luad_and_lusc","luad_lusc_lucl","lung_unique")
colnames(figure5_summarytable)<-c("total_number_loci","protcods_loci","cancer_realted_loci")
figure5_summarytable["luad_lusc_specific","total_number_loci"]<-nrow(all_natgenes_luad_lusc);
figure5_summarytable["luad_lusc_specific","cancer_realted_loci"]<-nrow(all_natgenes_luad_lusc[all_natgenes_luad_lusc$cancer_related==1,]);
figure5_summarytable["luad_lusc_specific","protcods_loci"]<-nrow(all_natgenes_luad_lusc[all_natgenes_luad_lusc$gene_biotype=="protein_coding",]);


figure5_summarytable["luad_and_lusc","total_number_loci"]<-nrow(all_natgenes_lung_tissuespecific);
figure5_summarytable["luad_and_lusc","cancer_realted_loci"]<-nrow(all_natgenes_lung_tissuespecific[all_natgenes_lung_tissuespecific$cancer_related==1,]);
figure5_summarytable["luad_and_lusc","protcods_loci"]<-nrow(all_natgenes_lung_tissuespecific[all_natgenes_lung_tissuespecific$gene_biotype=="protein_coding",]);

figure5_summarytable["luad_lusc_lucl","total_number_loci"]<-nrow(all_natgenes_lung_canspecific);
figure5_summarytable["luad_lusc_lucl","cancer_realted_loci"]<-nrow(all_natgenes_lung_canspecific[all_natgenes_lung_canspecific$cancer_related==1,]);
figure5_summarytable["luad_lusc_lucl","protcods_loci"]<-nrow(all_natgenes_lung_canspecific[all_natgenes_lung_canspecific$gene_biotype=="protein_coding",]);

figure5_summarytable["lung_unique","total_number_loci"]<-nrow(all_natgenes_unique_lung);
figure5_summarytable["lung_unique","cancer_realted_loci"]<-nrow(all_natgenes_unique_lung[all_natgenes_unique_lung$cancer_related==1,]);
figure5_summarytable["lung_unique","protcods_loci"]<-nrow(all_natgenes_unique_lung[all_natgenes_unique_lung$gene_biotype=="protein_coding",]);

#total_number_loci protcods_loci cancer_realted_loci
#luad_lusc_specific              6286          3437                 201
#luad_and_lusc                   1772           916                  44
#luad_lusc_lucl                  1226           664                  32
#lung_unique                      154            65                   2

###############################################################################
### code chunk number 3: Plot Figure 5
### Binary matrix
### lung cancer tissue specific heat maps
###############################################################################

setwd(results_dir_fig5)
library("gplots")

colsorder<-c("benign","luad","lusc","brca","prca","panc","meningioma","lucl")
heatmap_matrix_asense <- all_natgenes_lung_tissuespecific[all_natgenes_lung_tissuespecific$gene_biotype=="protein_coding",]
o<-order(heatmap_matrix_asense$numcohorts[heatmap_matrix_asense$gene_biotype=="protein_coding"],decreasing=FALSE)
heatmap_matrix_asense<-heatmap_matrix_asense[o,]

#o<-order(heatmap_matrix_asense$numcohorts,decreasing=FALSE)
heatmap_matrix_asense[heatmap_matrix_asense==1]<-2
heatmap_matrix_asense[heatmap_matrix_asense==0]<-1
color.palette <- c('wheat','red')
palette.breaks <- c(0,1,2)

pdf("Figure5_all_natgenes_lung_tissuespecific.pdf")

gn<-character(nrow(heatmap_matrix_asense))
cr<-which(heatmap_matrix_asense$cancer_related==2 &
				heatmap_matrix_asense$maxASscore>1000 )
gn[cr]<-as.vector(heatmap_matrix_asense$gene_name[cr])

heatmap.2(as.matrix(heatmap_matrix_asense[,colsorder]),
		Rowv=NA,Colv=NA,
		col = color.palette,#,rev(heat.colors(256)),#length(palette.breaks)-1)), #256
		breaks=palette.breaks,#palette.breaks,
		scale="none",
		dendrogram="none",
		key=TRUE,
		labRow=gn,#geneAnnot_ori[rownames(heatmap_matrix_asense),"gene_name"],
		symkey=FALSE, density.info="none", trace="none", cexRow=0.8) #,cexCol=0.1
dev.off()

# All heat map genes
all_natgenes_lung<-rbind(all_natgenes_unique_lung,all_natgenes_lung_tissuespecific)

colsorder<-c("benign","luad","lusc","brca","prca","panc","meningioma","lucl")
heatmap_matrix_asense <- all_natgenes_lung[all_natgenes_lung$gene_biotype=="protein_coding",]
o<-order(heatmap_matrix_asense$numcohorts[heatmap_matrix_asense$gene_biotype=="protein_coding"],decreasing=FALSE)
heatmap_matrix_asense<-heatmap_matrix_asense[o,]

#o<-order(heatmap_matrix_asense$numcohorts,decreasing=FALSE)
heatmap_matrix_asense[heatmap_matrix_asense==1]<-2
heatmap_matrix_asense[heatmap_matrix_asense==0]<-1
color.palette <- c('wheat','red')
palette.breaks <- c(0,1,2)

pdf("Figure5_all_natgenes_lung.pdf")

gn<-character(nrow(heatmap_matrix_asense))
cr<-which(heatmap_matrix_asense$cancer_related==2 &
				heatmap_matrix_asense$maxASscore>1000 )
gn[cr]<-as.vector(heatmap_matrix_asense$gene_name[cr])

heatmap.2(as.matrix(heatmap_matrix_asense[,colsorder]),
		Rowv=NA,Colv=NA,
		col = color.palette,#,rev(heat.colors(256)),#length(palette.breaks)-1)), #256
		breaks=palette.breaks,#palette.breaks,
		scale="none",
		dendrogram="none",
		key=TRUE,
		labRow=gn,#geneAnnot_ori[rownames(heatmap_matrix_asense),"gene_name"],
		symkey=FALSE, density.info="none", trace="none", cexRow=0.8) #,cexCol=0.1
dev.off()


### lung cancer unique antisense loci
colsorder<-c("benign","luad","lusc","brca","prca","panc","meningioma","lucl")
heatmap_matrix_asense <- all_natgenes_unique_lung[all_natgenes_unique_lung$gene_biotype=="protein_coding",]
o<-order(heatmap_matrix_asense$numcohorts[heatmap_matrix_asense$gene_biotype=="protein_coding"],decreasing=FALSE)
heatmap_matrix_asense<-heatmap_matrix_asense[o,]

heatmap_matrix_asense[heatmap_matrix_asense==1]<-2
heatmap_matrix_asense[heatmap_matrix_asense==0]<-1
color.palette <- c('wheat','red')
palette.breaks <- c(0,1,2)

pdf("Figure5_all_natgenes_unique_lung.pdf")

gn<-character(nrow(heatmap_matrix_asense))
cr<-which(heatmap_matrix_asense$cancer_related==2 &
				heatmap_matrix_asense$maxASscore>1000 )
gn[cr]<-as.vector(heatmap_matrix_asense$gene_name[cr])

heatmap.2(as.matrix(heatmap_matrix_asense[,colsorder]),
		Rowv=NA,Colv=NA,
		col = color.palette,#,rev(heat.colors(256)),#length(palette.breaks)-1)), #256
		breaks=palette.breaks,#palette.breaks,
		scale="none",
		dendrogram="none",
		key=TRUE,
		labRow=as.vector(heatmap_matrix_asense$gene_name), #gn,#geneAnnot_ori[rownames(heatmap_matrix_asense),"gene_name"],
		symkey=FALSE, density.info="none", trace="none", cexRow=0.8) #,cexCol=0.1
dev.off()

### lung cancer tissue specific heat map
### genes of interest
colsorder<-c("benign","luad","lusc","brca","prca","panc","meningioma","lucl")
heatmapgenesind<-intersect(which(all_natgenes_lung_tissuespecific$gene_biotype=="protein_coding"),
		which(all_natgenes_lung_tissuespecific$cancer_related==1))

heatmap_matrix_asense <- all_natgenes_lung_tissuespecific[heatmapgenesind,]
o<-order(heatmap_matrix_asense$numcohorts[heatmap_matrix_asense$gene_biotype=="protein_coding"],decreasing=FALSE)
heatmap_matrix_asense<-heatmap_matrix_asense[o,]


#o<-order(heatmap_matrix_asense$numcohorts,decreasing=FALSE)
heatmap_matrix_asense[heatmap_matrix_asense==1]<-2
heatmap_matrix_asense[heatmap_matrix_asense==0]<-1
color.palette <- c('wheat','red')
palette.breaks <- c(0,1,2)

pdf("Figure5_cancer_natgenes_lung_tissuespecific.pdf")

gn<-character(nrow(heatmap_matrix_asense))
cr<-which(heatmap_matrix_asense$cancer_related==2 &
				heatmap_matrix_asense$maxASscore>1000 )
gn[cr]<-as.vector(heatmap_matrix_asense$gene_name[cr])

heatmap.2(as.matrix(heatmap_matrix_asense[,colsorder]),
		Rowv=NA,Colv=NA,
		col = color.palette,#,rev(heat.colors(256)),#length(palette.breaks)-1)), #256
		breaks=palette.breaks,#palette.breaks,
		scale="none",
		dendrogram="none",
		key=TRUE,
		labRow=gn,#geneAnnot_ori[rownames(heatmap_matrix_asense),"gene_name"],
		symkey=FALSE, density.info="none", trace="none", cexRow=0.8) #,cexCol=0.1
dev.off()

write.table(all_natgenes_lung_tissuespecific,file="all_natgenes_lung_tissuespecific.tsv",
		row.names=TRUE, col.names=NA,sep='\t',quote=FALSE)

###############################################################################
### Done
###############################################################################

