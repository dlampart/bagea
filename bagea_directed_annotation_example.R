#########################################################
## BAGEA full example 
##
## the following is a full example of the bagea workflow
## This takes substantial resources. It assumes 15 cores
## and substantial memory.
## Do not run this on your laptop!! 
##
## make sure the function 'install_external_data' has been
## successfully called beforehand.
## This is the case if you ran 'bagea_tutorial.R'
##
##
## We are running bagea for the gtex fibroblast data.
## bagea is using all TF annotations derived from expecto.
## 
#########################################################
#########################################################

unloadNamespace("bagea")
require(bagea)
require(data.table)
require(Matrix)
#########################################################
NCORES=15
#########################################################
SAVE_PATH=""### change this to were you want the output of bagea saved to.
if(SAVE_PATH==""){
	stop("set SAVE_PATH before running this script.")
}

BAGEA_PATH=Sys.getenv("BAGEA_PATH")


my_annotation_mapping=prepare_annotation_mapping_mat(RANGE_AROUND_TSS = 2E5, TSS_CUTOFFS = c(250, 500, 1000, 2000, 5000, 10000, 20000,50000))

########################################################
## Next, we load the met-annotations and restrict the 
## directed annotations to TF annotations.
## We then prepare the SHIFT_METAINFO_TBL to call bagea
## in lasso mode.
##
## For the tutorial we use only 15 cores and run 300 rounds.
#########################################################
full_myshift_metainfo_tbl=load_metaannotation_tbl()
tokeep=full_myshift_metainfo_tbl[,`Assay type`=="TF"]
my_annotation_mapping[[3]]=my_annotation_mapping[[3]][,tokeep]


shift_metainfo_tbl=unique(full_myshift_metainfo_tbl[`Assay type`=="TF",list(annotation_id=annotation_id,cluster=annotation_id)])


myn=300 ## The sample size for the fibroblast data is taken from the gtex home-page
summary_stats_path=paste(BAGEA_PATH,"/Data/gtex_summary/Cells_Transformed_fibroblasts_chr",sep="")

########################################################
## Next, we run the bagea updates
##
## For the tutorial we use only 15 cores and run 300 rounds.
#########################################################
my_observed_data_train=prepare_observed_1KG_data(BAGEA_ANNOTATION_MAT=my_annotation_mapping,NSAMPLES=myn,VARIANCE_PERCENTAGE2KEEP = 95, RANGE_AROUND_TSS=150000,SUMMARY_BASE_STAT_PATH=summary_stats_path,NCORES=3,SINGLESNP_MINPVALUE_CUTOFF=1E-7,CHR=c(1:22),MAF_CUTOFF=0.05,SHIFT_METAINFO_TBL=shift_metainfo_tbl)

nu_names=colnames(my_annotation_mapping$annotation_mat)
myc=c(1,rep(0.3,length(nu_names)))
my_hyperparameter_list=set_hyperparameter_list(nu_names=nu_names,myc=myc)
########################################################
## Next, we run the bagea updates
##
## For the tutorial we use only 15 cores and run 300 rounds.
#########################################################
my_results=run_bagea(gene_dat_list=my_observed_data_train,hyperparameter_list=my_hyperparameter_list,calc_L=FALSE,ncores=NCORES,nrounds=300)
print(my_results$Eu_omega_list[[1]])

########################################################
## Finally, we print omega value fits that substantially
## deviate from 0.
#########################################################
omega_mat=my_results$Eu_omega_list[[1]][,1]
print("example done")
print("omega values deviating substantially from 0:")
print(omega_mat[abs(omega_mat)>0.001])


########################################################
## here, we remove the gene-wise for saving.
########################################################
my_results$gene_dat_list=NULL
save(my_results,file=SAVE_PATH)

