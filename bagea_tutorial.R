#########################################################
## BAGEA tutorial 
##
## the following is a full example of the bagea workflow
## with reduced data size to increase speed and decrease
## memory footprint.
#########################################################
#########################################################
require(bagea)
require(data.table)
require(Matrix)
#########################################################
## We start by installing gtex data into BAGEA_DATA for the tutorial.
##
## (make sure the function 'install_external_data' has been
## successfully called beforehand.)
#########################################################
install_gtex("https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/all_snp_gene_associations/Cells_Transformed_fibroblasts.allpairs.txt.gz")


#########################################################
## Next, we download directed annotation bed files for the tutorial.
##
## The following code downloads 6 directed  annotation bed files
## for the tutorial and places it in the folder ${BAGEA_PATH}/Data/tutorial_annotation_files/
#########################################################

filenames2download=c("SydhHuvecCfosUcd.bed","SydhK562Cfos.bed","UtaHuvecPol2.bed","UtaK562Pol2.bed","UtaK562Ctcf.bed","UtaHuvecCtcf.bed")

download_example_data=function(filenames2download){
	BAGEA_PATH=Sys.getenv("BAGEA_PATH")
	if(BAGEA_PATH==""){
		stop("BAGEA path has to be set.")
	}
	bagea_data_path=paste(BAGEA_PATH,"/Data",sep="")
	if(!dir.exists(bagea_data_path)){
		stop("${BAGEA_PATH}/Data folder doesn't exist. Make sure you run bagea::install_external_data() first.")
	}
	tutorial_data_path=paste(BAGEA_PATH,"/Data/tutorial_annotation_files/",sep="")
	if(!dir.exists(tutorial_data_path)){
		system(paste("mkdir",tutorial_data_path))
	}
	cmdstr=paste("wget -P ",tutorial_data_path," https://s3-us-west-1.amazonaws.com/bagea-data/bagea_data_freeze/tutorial_annotation_files/",sep="")
	cmdstrs=paste(cmdstr,filenames2download,sep="")
	print("downloading tutorial files..")
	lapply(cmdstrs,function(x){
		print(x)
		system(x)
	})
}

download_example_data(filenames2download)

BAGEA_PATH=Sys.getenv("BAGEA_PATH")
tutorial_data_path=paste(BAGEA_PATH,"/Data/tutorial_annotation_files/",sep="")


#########################################################
## Next, we prepare a `bagea_annotation_mat` object.
##
## RANGE_AROUND_TSS is set lower for the tutorial 
## to increased speed 
## and  decrease memory footprint.
#########################################################

my_annotation_mapping=prepare_annotation_mapping_mat(RANGE_AROUND_TSS = 2e+04, TSS_CUTOFFS = c(250, 500, 1000, 2000, 5000, 10000),SHIFT_BED_FILES_PATHS =tutorial_data_path, SHIFT_BED_FILES=filenames2download)


########################################################
## Next, we have to prepare SHIFT_METAINFO_TBL
##
## This is a data frame (or data table from the package data.table) 
## That groups the directed annotations acording to 'meta'-annotations 
## such as celltype (cellt) and assay type (assayt)
#########################################################

annotation_id=sub("\\.bed$","",filenames2download)
cellt=c("Huvec","K562","Huvec","K562","Huvec","K562")
assayt=c("Cfos","Cfos","Pol2","Pol2","Ctcf","Ctcf")

myshift_metainfo_tbl=data.table(annotation_id,cellt,assayt)


########################################################
## Next, prepare the  bagea_observed_data object. 
##
## For the tutorial we restrict processed genes to chromosome 21 for training and to chromosome 22 for testing.
#########################################################
pval_cutoff=1E-10
myn=300 ## The sample size for the fibroblast data is taken from the gtex home-page

summary_stats_path=paste(BAGEA_PATH,"/Data/gtex_summary/Cells_Transformed_fibroblasts_chr",sep="")

my_observed_data_train=prepare_observed_1KG_data(BAGEA_ANNOTATION_MAT=my_annotation_mapping,NSAMPLES=myn,VARIANCE_PERCENTAGE2KEEP = 95, RANGE_AROUND_TSS=20000,SUMMARY_BASE_STAT_PATH=summary_stats_path,NCORES=3,SINGLESNP_MINPVALUE_CUTOFF=pval_cutoff,CHR=c(1:15),MAF_CUTOFF=0.05,SHIFT_METAINFO_TBL=myshift_metainfo_tbl)

my_observed_data_test=prepare_observed_1KG_data(BAGEA_ANNOTATION_MAT=my_annotation_mapping,NSAMPLES=myn,VARIANCE_PERCENTAGE2KEEP = 95, RANGE_AROUND_TSS=20000,SUMMARY_BASE_STAT_PATH=summary_stats_path,NCORES=3,SINGLESNP_MINPVALUE_CUTOFF=pval_cutoff,CHR=c(16:22),MAF_CUTOFF=0.05,SHIFT_METAINFO_TBL=myshift_metainfo_tbl)
########################################################
## Next, we prepare the preare  bagea_hyperparameter_list object.
##
## We use the default settings.
#########################################################
### extract all undirected annotation names to set them
nu_names=colnames(my_annotation_mapping$annotation_mat)

my_hyperparameter_list=set_hyperparameter_list(p=5,chi1=1,chi2=0.0003,zeta1=1,zeta2=1E2,nu_names=nu_names)
########################################################
## Next, we run the bagea updates
##
## For the tutorial we use only 3 cores and run 50 rounds.
#########################################################
my_results=run_bagea(gene_dat_list=my_observed_data_train,hyperparameter_list=my_hyperparameter_list,calc_L=FALSE,ncores=3,nrounds=50)

########################################################
## Next, print the estimated expectation of omega:
#########################################################
print(my_results$Eu_omega_list[[1]])

########################################################
## Next, print the estimated expectation of upsilon cellt:
#########################################################
print(my_results$Eu_upsilon_list[[1]])
print(my_results$D_names_list[[1]])

########################################################
## Next, we calculate mse_dir for the genes in the test set:
#########################################################
mse_dir_tbl=predict_directed_fast(observed_data=my_observed_data_test,bagea_result=my_results)
########################################################
## Next, we test whether msed_dir in the test set is significantly lower than 1. 
## We take the top quartile of genes (w.r.t. predictor magnitude).
## Then we count how many mse_dir values are below 1.
## We then test via the binomial distribution.
#########################################################
mse_dir_tbl_top_quartile=mse_dir_tbl[S_approx>quantile(S_approx,0.75),]
smallerthan1=mse_dir_tbl_top_quartile[,mse_dir_approx<1]

pvalue=1-pbinom(q=(sum(smallerthan1)-1),size=length(smallerthan1),0.5)
print(pvalue)