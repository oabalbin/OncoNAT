# TODO: Add comment
# 
# Author: alebalbin
# This scripts sets up the paths and required file names for all analysis 
# scripts, which reproduce the results presented in 
# Balbin et al The landscape of antisense gene expression in human cancers
# Genome Research manuscript
#
###############################################################################

###############################################################################
### Set up of the data folders
### gunzip sascompendia_codebase
###############################################################################

base_dir<-''#Replace with the base path in your computer
root_dir <- paste(base_dir,'sascompendia_codebase',sep='/')
data_dir <-paste(root_dir,'data',sep='')
data_dir_OPSnastiseq<-paste(data_dir,'nastiseq_final_R_objects',sep='/')
references_dir<-paste(root_dir,'reference_files/',sep='')
source_code_dir<-paste(root_dir,'codebase',sep='/')
results_dir<-paste(root_dir,'results',sep='/')

###############################################################################
### Load general functions for analyses
###############################################################################
# Source the functions for this pipeline
source(paste(source_code_dir,'sas_general_functions.R',sep=''))
# Check that all matrices have the same size and the same genes in same rows.
# Load the matrices
#source(paste(source_code_dir,'sas_load_count_matrices_forDESeq.R',sep=''))

###############################################################################
### Load pre-calculated tables and objects
###############################################################################
TRANSCRIPT_MATRICES_ROBJ <- paste(data_dir,'cohort.matrices.transcripts.DESeq.R',sep="/")
OVERLAP_ANNOTATION_FILE<-paste(references_dir,"GRCh37.69.saspairs_trans.annotation.tsv",sep="/")
TUMORSUP_FILE <- paste(references_dir,"Cancer_Genes_Tumor_Suppresors.txt")
ONOCOGENES_FILE <- paste(references_dir,"Cancer_Genes_Oncogenes.txt")


###############################################################################
### Figure 1
### Antisense expression is pervasive across the human transcriptome.
###############################################################################

#Figure 1 scripts will calculate compute the number of loci with antisense ratio
#greater than the confidence interval for the protocol error rate: pei_th.

#The main data structures that are produce are: asloci_pei_sd1_countbysample, 
#asloci_pei_sd1_countbycohort.

#Run: figure_1_GenomeResearch.R
source(paste(source_code_dir,'Balbin_etal_Figure_1_GenomeResearch.R',sep=''))
###############################################################################
# 2. Determine the number of bonafide antisense loci combining the information
#obtained in step 1 and the prediction of NASTI-seq
#Run: compute_nastiseq scores_allcohorts.R
#This script will call compute_nastiseq_loci_allcohorts.R
#In order to speed computation time, 
#the main data structures produced for these scripts are provided as 
#precomputed R objects in the ~/data/nastiseq_final_R_objects/

###############################################################################

###############################################################################
### Figure 2
### Correlation of cis-NAT pair expression and regulation by bidirectional 
### promoters.
### Figure 3
### Antisense regulation of cognate sense genes.
###############################################################################

### Data structures to load
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
# Note find this file name and create a relative path
FULL_UNIQNASTISEQ_LOCI_FILE <-  paste(data_dir,"all_cisNATpairs_OPS_nastiseq_loci.tsv",sep='/')

source(paste(source_code_dir,'Balbin_etal_Figure_2_3_GenomeResearch.R',sep=''))
###############################################################################
### Figure 4
### Expression of antisense loci across cancer subtypes/
###############################################################################

### Data structures to load
# LOAD CISNAT TRAINING DATASET
CISNAT_TRAINING_DATASET_FILE <-pasta(data_dir,'cis_natpairs_training_dataset.R',sep="/")

# Directory with the R objects with the matrices 
# with the antisense loci nominated by both the 
# OPS and NASTIseq methods for each cohort

# data_dir_OPSnastiseq
# R Object base name. There is an R object for each cohort.
file_basename_final<-'NASTIseq_fullresults_trans_mincount5_rebuttal'

source(paste(source_code_dir,'Balbin_etal_Figure_4_GenomeResearch.R',sep=''))

###############################################################################
### Figure 5
### Cancer specific antisense loci. 
###############################################################################

### Data structures to load

source(paste(source_code_dir,'Balbin_etal_Figure_5_GenomeResearch.R',sep=''))
###############################################################################
### Create main OncoNAT tables
### 
###############################################################################
### Data structures to load
CORR_LOCI_K2OV_AC_FILE_TABLE<-paste(data_dir,
		'known_overlaps_correlation_table_accross_cohorts.corrtable.tsv',sep="/")
ALL_NATGENES_OBJECTS <- paste(data_dir,'all_natgenes_objects.R',sep="/")
CPGI_MAPPED_TABLE_FILE<-paste(references_dir,'genes.GRCh37.69.saspairs_trans.cpgislands_count.tsv')

source(paste(source_code_dir,'OncoNAT_tables.R',sep=''))

###############################################################################
### Figure 6
### Antisense loci dysregulation in cancer
###############################################################################

source(paste(source_code_dir,'Balbin_etal_Figure_6_GenomeResearch.R',sep=''))

