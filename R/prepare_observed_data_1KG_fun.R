require(parallel)
require(Matrix)

add_genotypes2bagea_annotation_mat=function(bagea_annotation_mat,genotype_mat_scaled,ADD_SNP_FREQS2ANNOT){
	genotype_matcher=match(bagea_annotation_mat$snp_tss_map[,snpid],rownames(genotype_mat_scaled))
	genotype_mat_expanded=genotype_mat_scaled[genotype_matcher,]
	if(ADD_SNP_FREQS2ANNOT){
		snpfreq_annot_mat=get_snpfreq_annot(genotype_mat,maf_ranges_lower=c(0,0.05,0.1,0.2,0.35),maf_ranges_upper=c(0.05,0.1,0.2,0.35,0.5))
		snpfreq_annot_mat_expanded=snpfreq_annot_mat[genotype_matcher,]
		mycolnames=c(colnames(bagea_annotation_mat$annotation_mat),colnames(snpfreq_annot_mat_expanded))
		bagea_annotation_mat$annotation_mat=cbind(annotation_mat_chr,snpfreq_annot_mat_expanded)
		colnames(bagea_annotation_mat$annotation_mat)=mycolnames
	}
	bagea_annotation_mat$genotype_mat_expanded=genotype_mat_expanded
	class(bagea_annotation_mat)=c("bagea_annotation_mat_wgenotypes",class(bagea_annotation_mat))
	return(bagea_annotation_mat)
}

subset_bagea_annotation_mat_onsnpnames=function(bagea_annotation_mat,snpnames,symols2filteron=NULL){
	### remove SNPS that are not present in all files # genotype_mat_scaled gets expanded later. 
	if(!is.null(symols2filteron)){
		snpsongenes=bagea_annotation_mat$snp_tss_map[is.element(symbol,symols2filteron),snpid]
		snpnames=intersect(snpnames,snpsongenes)
	}
	snpsintersected=intersect_order_preserving(bagea_annotation_mat$snp_tss_map[,snpid],snpnames)
	isel=bagea_annotation_mat$snp_tss_map[,is.element(snpid,snpsintersected)]
	bagea_annotation_mat$snp_tss_map=bagea_annotation_mat$snp_tss_map[isel,]
	bagea_annotation_mat$annotation_mat=bagea_annotation_mat$annotation_mat[isel,]
	if(shift_mat_is_available(bagea_annotation_mat)){
		bagea_annotation_mat$annotation_shift_mat=bagea_annotation_mat$annotation_shift_mat[isel,]
	}
	return(bagea_annotation_mat)
}

add_summary_stats2bagea_annotation_mat=function(bagea_annotation_mat,summary_stats_tbl){
	##### sort summary_stats_tbl into snp_tss_map_chr order
	## (both genes and snps)
	###########################
	bagea_annotation_mat$snp_tss_map[,fused:=paste(snpid,symbol)]
	summary_stats_tbl[,fused:=paste(snps,gene)]
	myfused=intersect_order_preserving(bagea_annotation_mat$snp_tss_map[,fused],summary_stats_tbl[,fused])
	###########################
	isel=bagea_annotation_mat$snp_tss_map[,is.element(fused,myfused)]
	bagea_annotation_mat$snp_tss_map=bagea_annotation_mat$snp_tss_map[isel,]
	bagea_annotation_mat$annotation_mat=bagea_annotation_mat$annotation_mat[isel,]

	if(shift_mat_is_available(bagea_annotation_mat)){
		bagea_annotation_mat$annotation_shift_mat=bagea_annotation_mat$annotation_shift_mat[isel,]
	}
	summary_stats_tbl=summary_stats_tbl[match(bagea_annotation_mat$snp_tss_map[,fused],fused),]
	###########################
	if(sum(summary_stats_tbl[,fused]!=bagea_annotation_mat$snp_tss_map[,fused])!=0){stop("errs23emi2")}
	bagea_annotation_mat$summary_stats_tbl=summary_stats_tbl
	class(bagea_annotation_mat)=c("bagea_annotation_mat_wsumstats",class(bagea_annotation_mat))
	return(bagea_annotation_mat)
}



subset_annotation_mat_rows2range=function(bagea_annotation_mat,RANGE_AROUND_TSS,gwas_summary_stats_mode,symbols2filteron=NULL){
	if(gwas_summary_stats_mode){
		snp_tss_map_chr_id=bagea_annotation_mat$snp_tss_map[,snpid==snpid]
	}else{
		if(is.null(symbols2filteron)){
			snp_tss_map_chr_id=bagea_annotation_mat$snp_tss_map[,abs(snp_pos-tss_pos)<RANGE_AROUND_TSS]
		}else{
			snp_tss_map_chr_id=bagea_annotation_mat$snp_tss_map[,abs(snp_pos-tss_pos)<RANGE_AROUND_TSS & is.element(symbol,symbols2filteron)]
		}
	}
	bagea_annotation_mat$snp_tss_map=bagea_annotation_mat$snp_tss_map[snp_tss_map_chr_id,]
	bagea_annotation_mat$annotation_mat=bagea_annotation_mat$annotation_mat[snp_tss_map_chr_id,]
	if(shift_mat_is_available(bagea_annotation_mat)){
		bagea_annotation_mat$annotation_shift_mat=bagea_annotation_mat$annotation_shift_mat[snp_tss_map_chr_id,]
	}
	return(bagea_annotation_mat)
}


subset_annotation_mat_rows2chr=function(bagea_annotation_mat,chr){
	curstr=paste("chr",chr,sep="")
	snps_on_chr_inds=which(bagea_annotation_mat$snp_tss_map[,chrom==curstr])
	bagea_annotation_mat$snp_tss_map=bagea_annotation_mat$snp_tss_map[snps_on_chr_inds,]
	bagea_annotation_mat$annotation_mat=bagea_annotation_mat$annotation_mat[snps_on_chr_inds,]
	if(shift_mat_is_available(bagea_annotation_mat)){
		bagea_annotation_mat$annotation_shift_mat=bagea_annotation_mat$annotation_shift_mat[snps_on_chr_inds,]
	}
	return(bagea_annotation_mat)
}

shift_mat_is_available=function(bagea_annotation_mat){
	!is.null(bagea_annotation_mat[["annotation_shift_mat"]])
}

do_simulation=function(param_list){
	!is.null(param_list)
}

load_all_genotype_mats=function(chr,bagea_annotation_mat,gwas_summary_stats_mode,SNPS2KEEP,set){
	if(gwas_summary_stats_mode){
		snp_ids_in_range=bagea_annotation_mat$snp_tss_map[,snpid]
	}else{
		snp_ids_in_range=bagea_annotation_mat$snp_tss_map[abs(snp_pos-tss_pos)<set$RANGE_AROUND_TSS,snpid]
	}
	genotype_mats=lapply(c(1:length(set$GENOTYPE_BASE_FILE_PATHS)),function(i){
		load_genotype_mat(chr,GENOTYPE_BASE_FILE_PATHS=set$GENOTYPE_BASE_FILE_PATHS[i],GENOTYPE_FILE_PATHS_GZ=set$GENOTYPE_FILE_PATHS_GZ[i],snps2keep=snp_ids_in_range)
	})
	if(!is.null(SNPS2KEEP)){
		snps2keep_chr=intersect(SNPS2KEEP,rownames(genotype_mats[[1]]))
	}else{
		snps2keep_chr=NULL
	}
	genotype_mats=intersect_snps_fully(genotype_mats,snps2keep_chr)
	return(genotype_mats)
}

prepare_data_per_chr=function(chr,
	bagea_annotation_mat,	
	simulated_params,
	SNPS2KEEP,
	set
	){
	# restrict annotation data to single chr
	bagea_annotation_mat=subset_annotation_mat_rows2chr(bagea_annotation_mat,chr)
	###########################
	# load summary stats.
	###########################
	summary_stat_path=paste(set$SUMMARY_BASE_STAT_PATH,chr,".txt",sep="")
	summary_stats_tbl=fread(summary_stat_path,header=TRUE)
	###########################
	# END: load summary stats.
	###########################
	###########################
	# check if files are GWAS files and process accordingly.
	###########################
	gwas_summary_stats_mode=is_gwas_summary_stats_tbl(summary_stats_tbl)
	if(gwas_summary_stats_mode){
		summary_stats_tbl=preprocess_gwas_summary_stat_tbl(summary_stats_tbl,bagea_annotation_mat$snp_tss_map)
	}
	###########################
	# END: check if files are GWAS files and process accordingly.
	###########################
	chrstr=paste("chr",chr,sep="")
#	browser()
	genotype_mat=load_all_genotype_mats(chr,bagea_annotation_mat,gwas_summary_stats_mode,SNPS2KEEP,set=set)[[1]]
	###TODO: fix this
	my_variance_transform=FALSE
	all_snps_below_mafcutoff=get_lowmaf_snps(genotype_mat,set$MAF_CUTOFF,variance_transform=my_variance_transform)
	genotype_mat=remove_lowmaf_snps(genotype_mat,set$MAF_CUTOFF,variance_transform=my_variance_transform)
	###########################
	# END: remove low MAF snps.
	###########################
	###########################
	# scale genotype mat
	###########################
	genotype_mat_scaled=scale_genotype_mat_1KG(genotype_mat,set$NSAMPLES)
	rm(genotype_mat)
	###########################
	# END: scale genotype  mats
	###########################
	###########################
	# restrict annotation data to relevant range
	###########################
	bagea_annotation_mat=subset_annotation_mat_rows2range(bagea_annotation_mat,set$RANGE_AROUND_TSS,gwas_summary_stats_mode)
	###########################
	# restrict annotation data to snps on genotypes.
	###########################
	bagea_annotation_mat=subset_bagea_annotation_mat_onsnpnames(bagea_annotation_mat=bagea_annotation_mat,snpnames=rownames(genotype_mat_scaled))

	bagea_annotation_mat=add_summary_stats2bagea_annotation_mat(bagea_annotation_mat,summary_stats_tbl)
	###########################
	###########################
	##### expand genotype matrix into 
	###########################
	bagea_annotation_mat=add_genotypes2bagea_annotation_mat(bagea_annotation_mat,genotype_mat_scaled,set$ADD_SNP_FREQS2ANNOT)
	rm(genotype_mat_scaled)
	do_simulation=do_simulation(simulated_params)
	###########################
	# END: map genotype mat to annotation mat
	###########################	
	###########################
	### get gene list to process
	genes_to_process=unique(bagea_annotation_mat$snp_tss_map[,symbol])
	glen=length(genes_to_process)
	print("prepare observed gene list")
	print(paste(glen," to process",sep=""))
	all_Obs_list=lapply(c(1:glen),process_gene_1KG,
		genes_to_process=genes_to_process,
		bagea_annotation_mat=bagea_annotation_mat,
		simulated_params=simulated_params,
		NSAMPLES=set$NSAMPLES,
		SINGLESNP_MINPVALUE_CUTOFF=set$SINGLESNP_MINPVALUE_CUTOFF,
		VARIANCE_PERCENTAGE2KEEP=set$VARIANCE_PERCENTAGE2KEEP,
		ANNOTATION_SHIFTWEIGHT_LIST=set$ANNOTATION_SHIFTWEIGHT_LIST
	)
	names(all_Obs_list)=genes_to_process
	if(set$SINGLESNP_MINPVALUE_CUTOFF<1){
		all_Obs_list=filter_on_minpvalue(all_Obs_list,set$SINGLESNP_MINPVALUE_CUTOFF)
	}
	attr(all_Obs_list, "all_snps_below_mafcutoff")=all_snps_below_mafcutoff
	return(all_Obs_list)	
}

process_gene_1KG=function(i,
	genes_to_process,
	bagea_annotation_mat,
	simulated_params,
	NSAMPLES,
	SINGLESNP_MINPVALUE_CUTOFF,
	VARIANCE_PERCENTAGE2KEEP,
	ANNOTATION_SHIFTWEIGHT_LIST
	){
	geneid=genes_to_process[i]
	do_simulation=do_simulation(simulated_params)
	if(i %% 50 == 1){
		print(paste(i, " of ",length(genes_to_process)," processed",sep=""))
	}
	chr=sub("chr","",bagea_annotation_mat$snp_tss_map[1,chrom])
	chrnr=as.numeric(chr)
	if(!is.na(chrnr)){
		chr=chrnr
	}
	expanded_ids=which(bagea_annotation_mat$snp_tss_map[,symbol==geneid])
	X=bagea_annotation_mat$genotype_mat_expanded[expanded_ids,,drop=F]
	mapping_mat=bagea_annotation_mat$annotation_mat[expanded_ids,,drop=F]
	if(shift_mat_is_available(bagea_annotation_mat)){
		mapping_shift_mat=bagea_annotation_mat$annotation_shift_mat[expanded_ids,,drop=F]
	}
	ytX=t(bagea_annotation_mat$summary_stats_tbl[gene==geneid,statistic]*sqrt(NSAMPLES))
	if(SINGLESNP_MINPVALUE_CUTOFF<1){
		pval_cutoff=pnorm(-max(abs(ytX))/sqrt(NSAMPLES))*2
		if(pval_cutoff>SINGLESNP_MINPVALUE_CUTOFF){
			return(NULL)
		}
	}
	yty=NSAMPLES
	ns=NSAMPLES
	eigens=rep(0,ns)
	Xsvd=my_caught_svd(X)
	dl=length(Xsvd$d)
	eigens[c(1:dl)]=Xsvd$d^2
	eigens[eigens<0]=0
	cumeigens=cumsum(eigens)
	cumeigens_scaled=cumeigens/cumeigens[length(cumeigens)]
	selection_indices=c(1:min(which(cumeigens_scaled>=VARIANCE_PERCENTAGE2KEEP/100)))
	if(length(selection_indices)>1 && sum(selection_indices)>1 ){
		XtX_sqrt=Xsvd$u[,selection_indices,drop=F]%*%diag(Xsvd$d[selection_indices,drop=F])
	}else{
		XtX_sqrt=Xsvd$u[,selection_indices,drop=F]*Xsvd$d[selection_indices]
	}
	out=list()
	rownames(XtX_sqrt)=rownames(X)
	if(do_simulation){
		if(is.null(simulated_params[["p"]])){
			genewise_param_list=get_genewise_param_list(mapping_mat,mapping_shift_mat=NULL,hyperparameter_list=simulated_params)
		}else{
			genewise_param_list=get_genewise_param_list(mapping_mat,mapping_shift_mat,hyperparameter_list=simulated_params,annotation_nu_names=ANNOTATION_SHIFTWEIGHT_LIST)
		}
		trunc_var_scaled=1-sum(eigens[selection_indices,drop=F])/sum(eigens)
		z=simulate_z(XtX_sqrt,genewise_param_list,my_n=ns,trunc_var_scaled=trunc_var_scaled)
		ytX=t(z)*sqrt(ns)
	}
	out[["chrom"]]=chr
	out[["XtX_sqrt"]]=XtX_sqrt
	out[["ytX"]]=ytX
	out["yty"]=yty
	out[["n"]]=ns
	out[["eig_values"]]=eigens
	out[["mapping_mat"]]=mapping_mat
	if(shift_mat_is_available(bagea_annotation_mat)){
		out[["mapping_shift_mat"]]=mapping_shift_mat
	}
	if(do_simulation){
		out[["genewise_param_list"]]=genewise_param_list
	}
	return(out)
}


filter_on_minpvalue=function(all_Obs_list,SINGLESNP_MINPVALUE_CUTOFF){
	null_inds=sapply(all_Obs_list,function(x){is.null(x)})
	print(paste(sum(null_inds),"genes filtered due to no SNP having p-value below",SINGLESNP_MINPVALUE_CUTOFF))
	all_Obs_list=all_Obs_list[!null_inds]
	print(paste(length(all_Obs_list),"genes left."))
	return(all_Obs_list)
}

is_gwas_summary_stats_tbl=function(summary_stats_tbl){
	all_g=unique(summary_stats_tbl[,gene])
	if(length(all_g)==1 && all_g[1]=="."){
		gwas_summary_stats_mode=TRUE
	}else{
		gwas_summary_stats_mode=FALSE
	}
	return(gwas_summary_stats_mode)
}

scale_genotype_mat_1KG=function(genotype_mat,NSAMPLES){
	print("scale genotype matrix to mean 0 sd and scale to variance it would have had in the full dataset ..")
	genotype_mat_scaled=t(genotype_mat)
	genotype_mat_scaled=my_scale(genotype_mat_scaled,use_minus1=FALSE)
	genotype_mat_scaled=t(genotype_mat_scaled)
	genotype_mat_scaled=genotype_mat_scaled*sqrt(NSAMPLES/dim(genotype_mat)[2])
	return(genotype_mat_scaled)
}

preprocess_gwas_summary_stat_tbl=function(summary_stats_tbl,snp_tss_map_chr){
	blockids=unique(snp_tss_map_chr[,symbol])
	if(sum(!grepl("^block\\d+$",blockids))>0){
		stop("err[aacp,:in gwas_summary_stats_mode annotation mapping has to have format block\\d+")
	}
	summary_stats_tbl_sub=summary_stats_tbl[is.element(snps,snp_tss_map_chr[,snpid]),]
	mymatcher=match(summary_stats_tbl_sub[,snps],snp_tss_map_chr[,snpid])
	myblock_vec=snp_tss_map_chr[mymatcher,symbol]
	summary_stats_tbl_sub[,gene:=myblock_vec]
	summary_stats_tbl=summary_stats_tbl_sub
	return(summary_stats_tbl)
}

simulate_params_from_hyperparams=function(HYPER_PARAM_LIST,D_mat,ANNOTATION_LIST,ANNOTATION_SHIFT_LIST,ANNOTATION_SHIFTWEIGHT_LIST){
	HYPER_PARAM_LIST=add_annotation_params(HYPER_PARAM_LIST,ANNOTATION_LIST)
	HYPER_PARAM_LIST=get_global_param_list(HYPER_PARAM_LIST,D_mat=D_mat)
	HYPER_PARAM_LIST=add_annotation_shift_params(HYPER_PARAM_LIST,ANNOTATION_SHIFT_LIST,D_mat=D_mat)
	if(!is.null(HYPER_PARAM_LIST$p) && sum(class(HYPER_PARAM_LIST)=="bagea_hyperparameter_bool_list")==0){
		HYPER_PARAM_LIST=add_annotation_nu_params(HYPER_PARAM_LIST,ANNOTATION_SHIFTWEIGHT_LIST)		
	}
	return(HYPER_PARAM_LIST)
}

subset_annotation_mat_cols=function(BAGEA_ANNOTATION_MAT,ANNOTATION_LIST,ANNOTATION_SHIFT_LIST,ANNOTATION_SHIFTWEIGHT_LIST){
	annot_list=BAGEA_ANNOTATION_MAT
	snp_tss_map=annot_list[["snp_tss_map"]]
	annotation_mat=annot_list[["annotation_mat"]]
	annotation_mat_full=annot_list[["annotation_mat"]]
	if(sum(names(annot_list)=="annotation_shift_mat")==1){
		shift_mat_available=TRUE
	}else{
		shift_mat_available=FALSE
	}
	if(shift_mat_available){
		annotation_shift_mat=annot_list[["annotation_shift_mat"]]
	}
	annotation_not_found=setdiff(ANNOTATION_LIST,colnames(annotation_mat))
	if(length(annotation_not_found)!=0){
		stopstr=paste("Not all in ANNOTATION_LIST found in BAGEA_ANNOTATION_MAT. first one missing is ",annotation_not_found[1],sep="")
		stop(stopstr)
	}
	shift_annotation_not_found=setdiff(ANNOTATION_SHIFT_LIST,colnames(annotation_shift_mat))
	if(length(shift_annotation_not_found)!=0){
		stopstr=paste("errpomasda:Not all in ANNOTATION_SHIFT_LIST found in BAGEA_ANNOTATION_MAT. first one missing is ",shift_annotation_not_found[1],sep="")
		stop(stopstr)
	}
	annotations_to_keep=intersect(colnames(annotation_mat),ANNOTATION_LIST)
	annotations_shiftweight_to_keep=intersect(colnames(annotation_mat),ANNOTATION_SHIFTWEIGHT_LIST)
	annotations_to_keep=union(annotations_to_keep,annotations_shiftweight_to_keep)
	shift_annotations_to_keep=intersect(colnames(annotation_shift_mat),ANNOTATION_SHIFT_LIST)
	print(paste("trying to find ",length(ANNOTATION_LIST)," annotations",sep=""))
	print(paste("found ",length(annotations_to_keep)," annotations",sep=""))
	annotation_mat=annotation_mat[,is.element(colnames(annotation_mat),annotations_to_keep)]
	print(paste("trying to find ",length(ANNOTATION_SHIFT_LIST)," annotations",sep=""))
	print(paste("found ",length(shift_annotations_to_keep)," annotations",sep=""))
	annotation_mat=annotation_mat[,is.element(colnames(annotation_mat),annotations_to_keep)]
	annotation_shift_mat=annotation_shift_mat[,is.element(colnames(annotation_shift_mat),shift_annotations_to_keep)]
	bagea_annotation_mat=BAGEA_ANNOTATION_MAT
	bagea_annotation_mat[["annotation_mat"]]=annotation_mat
	bagea_annotation_mat[["annotation_shift_mat"]]=annotation_shift_mat
	return(bagea_annotation_mat)
}


###########################
# prepare D_mat from SHIFT_METAINFO_TBL
###########################
prepare_D_mat=function(SHIFT_METAINFO_TBL){
	if(!is.null(SHIFT_METAINFO_TBL)){
		myfun=function(x){
			unq=unique(x)
			mymatch=match(x,unq)
		}
		D_mat=apply(SHIFT_METAINFO_TBL,2,myfun)
		myfun=function(x){
			unq=unique(x)		
		}
		D_names_list=lapply(SHIFT_METAINFO_TBL,myfun)
		names(D_names_list)=names(SHIFT_METAINFO_TBL)
		colnames(D_mat)=names(SHIFT_METAINFO_TBL)
	}else{
		D_mat=NULL
		D_names_list=NULL
	}
	outl=list()
	outl[["D_mat"]]=D_mat
	outl[["D_names_list"]]=D_names_list
	return(outl)
}

setup_D_mat_with_checks=function(SHIFT_METAINFO_TBL,bagea_annotation_mat){
	shift_mat_colnames=colnames(bagea_annotation_mat[["annotation_shift_mat"]])
	if(!is.data.frame(SHIFT_METAINFO_TBL) && !is.null(SHIFT_METAINFO_TBL)){
		stop("SHIFT_METAINFO_TBL has to be either NULL or a data.frame")
	}
	if(is.null(SHIFT_METAINFO_TBL)){
		D_mat_outl=list()
		D_mat=matrix(1,length(bagea_annotation_mat),1)
		colnames(D_mat)="all"
		D_names_list=list()
		D_names_list[["all"]]="all"
		D_mat_outl[["D_mat"]]=D_mat
		D_mat_outl[["D_names_list"]]=D_names_list
		return(D_mat_outl)
	}
	if(!is.element("annotation_id",colnames(SHIFT_METAINFO_TBL))){
		stop("SHIFT_METAINFO_TBL has to contain a column named annotation_id.")
	}	
	shift_metainfo_tbl=as.data.table(SHIFT_METAINFO_TBL)
	annot_ids=shift_metainfo_tbl[,annotation_id]
	if(length(unique(annot_ids))!=length(annot_ids)){
		stop("SHIFT_METAINFO_TBL annotation_id column has to be unique")
	}
	missing_annots=setdiff(shift_mat_colnames,annot_ids)
	if(length(missing_annots)!=0){
		stopstr=paste("SHIFT_METAINFO_TBL annotation_id column does not contain all necessary undirected annotations. first one missing: ",missing_annots[1],sep="")
		stop(stopstr)
	}
	mymatcher=match(shift_mat_colnames,annot_ids)
	if(sum(is.na(mymatcher))>0){
		stop("SHIFT_METAINFO_TBL doesnt match shift_mat_colnames.")
	}
	shift_metainfo_tbl_sorted=shift_metainfo_tbl[mymatcher,]
	shift_metainfo_tbl_sorted[,annotation_id:=NULL]
	D_mat_outl=prepare_D_mat(shift_metainfo_tbl_sorted)
	return(D_mat_outl)
}

## aggregating settings into single object.
get_settingsout_list_obs1KG=function(BAGEA_ANNOTATION_MAT,ANNOTATION_LIST,ANNOTATION_SHIFT_LIST,GENOTYPE_BASE_FILE_PATHS,GENOTYPE_FILE_PATHS_GZ,SUMMARY_BASE_STAT_PATH,NSAMPLES,VARIANCE_PERCENTAGE2KEEP,RANGE_AROUND_TSS,MAF_CUTOFF,HYPER_PARAM_LIST,CHR,NCORES,SNPS2KEEP,ADD_SNP_FREQS2ANNOT,SINGLESNP_MINPVALUE_CUTOFF,ANNOTATION_SHIFTWEIGHT_LIST,SHIFT_METAINFO_TBL){
		settings=list()
	settings[["GENOTYPE_BASE_FILE_PATHS"]]=GENOTYPE_BASE_FILE_PATHS
	settings[["GENOTYPE_FILE_PATHS_GZ"]]=GENOTYPE_FILE_PATHS_GZ
	settings[["SUMMARY_BASE_STAT_PATH"]]=SUMMARY_BASE_STAT_PATH
	settings[["BAGEA_ANNOTATION_MAT"]]=BAGEA_ANNOTATION_MAT$settings
	settings[["ANNOTATION_LIST"]]=ANNOTATION_LIST
	settings[["MAF_CUTOFF"]]=MAF_CUTOFF
	settings[["RANGE_AROUND_TSS"]]=RANGE_AROUND_TSS
	settings[["VARIANCE_PERCENTAGE2KEEP"]]=VARIANCE_PERCENTAGE2KEEP
	settings[["CHR"]]=CHR
#	settings[["KEEP_X_y"]]=KEEP_X_y
	settings[["NCORES"]]=NCORES
	settings[["n_snps2keep"]]=length(SNPS2KEEP)
	settings[["ADD_SNP_FREQS2ANNOT"]]=ADD_SNP_FREQS2ANNOT
	settings[["SINGLESNP_MINPVALUE_CUTOFF"]]=SINGLESNP_MINPVALUE_CUTOFF
	settings[["ANNOTATION_SHIFTWEIGHT_LIST"]]=ANNOTATION_SHIFTWEIGHT_LIST
	settings[["SHIFT_METAINFO_TBL"]]=SHIFT_METAINFO_TBL
	settings[["NSAMPLES"]]=NSAMPLES
	settings[["ANNOTATION_SHIFT_LIST"]]=ANNOTATION_SHIFT_LIST
	return(settings)
}

## no processing here, just checking of input.
check_settings_list_obs1KG=function(BAGEA_ANNOTATION_MAT,SNPS2KEEP,
	set){
	if(is.null(set$HYPER_PARAM_LIST) && is.null(set$SUMMARY_BASE_STAT_PATH)){
		stop("neither HYPER_PARAM_LIST for simulation supplied nor an summary stat file path given.")
	}
	if(!is.null(set$HYPER_PARAM_LIST) && !is.null(set$SUMMARY_BASE_STAT_PATH)){
		stop("both HYPER_PARAM_LIST for simulation  and an SUMMARY_BASE_STAT_PATH file path is set. set one to NULL")
	}
	if(!is.null(set$HYPER_PARAM_LIST) && sum(class(set$HYPER_PARAM_LIST)=="bagea_hyperparameter_list")){
		stop("HYPER_PARAM_LIST is set but is not a bagea_hyperparameter_list object.")
	}
	if(sum(class(BAGEA_ANNOTATION_MAT)=="bagea_annotation_mat")==0){
		stop("BAGEA_ANNOTATION_MAT argument is not a bagea_annotation_mat object.")
	}
	if(!is.vector(set$GENOTYPE_BASE_FILE_PATHS) || !is.character(set$GENOTYPE_BASE_FILE_PATHS)){
			stop("GENOTYPE_BASE_FILE_PATHS argument not a character vector")
	}
	if(!is.character(set$ANNOTATION_LIST) || !is.vector(set$ANNOTATION_LIST)){
		stop("ANNOTATION_LIST argument is not a character vector")
	}
	if(is.null(colnames(BAGEA_ANNOTATION_MAT[["annotation_shift_mat"]]))){
		stop("BAGEA_ANNOTATION_MAT$annotation_shift_mat has no colnames")
	}
	if(is.null(colnames(BAGEA_ANNOTATION_MAT[["annotation_mat"]]))){
		stop("BAGEA_ANNOTATION_MAT$annotation_mat has no colnames")
	}
	if(!is.character(set$ANNOTATION_SHIFT_LIST) || !is.vector(set$ANNOTATION_SHIFT_LIST)){
		stop("ANNOTATION_SHIFT_LIST argument is not a character vector")
	}
	### check GENOTYPE_FILE_PATHS_GZ
	if(!is.vector(set$GENOTYPE_FILE_PATHS_GZ) || !is.logical(set$GENOTYPE_FILE_PATHS_GZ)){
		stop("GENOTYPE_FILE_PATHS_GZ has to be a single boolean")
	}
	if(!is.null(set$SUMMARY_BASE_STAT_PATH)){
		if(length(set$SUMMARY_BASE_STAT_PATH)!=length(set$GENOTYPE_BASE_FILE_PATHS)){
			stop("not the same number of SUMMARY_BASE_STAT_PATH and GENOTYPE_BASE_FILE_PATHS supplied")
		}
	}
	if(length(set$GENOTYPE_FILE_PATHS_GZ)!=length(set$GENOTYPE_BASE_FILE_PATHS)){
		stop("not the same number of GENOTYPE_FILE_PATHS_GZ and GENOTYPE_BASE_FILE_PATHS supplied")
	}
	if(length(set$GENOTYPE_FILE_PATHS_GZ)!=1){
		stop("1KG setting should only have one element in GENOTYPE_FILE_PATHS_GZ.")
	}
	### check KEEP_X_y
	if(!is.null(SNPS2KEEP) && ( !is.vector(SNPS2KEEP) || !is.character(SNPS2KEEP))){
		stop("SNPS2KEEP has to be null or a character vector")
	}
	if(!is.scalar(set$ADD_SNP_FREQS2ANNOT) || !is.logical(set$ADD_SNP_FREQS2ANNOT)){
		stop("ADD_SNP_FREQS2ANNOT has to be a single boolean")
	}
	check_posinteger(set$NCORES)
	### check GENOTYPE_BASE_FILE_PATHS
	if(!is.null(set$SUMMARY_BASE_STAT_PATH)){
		if(!is.vector(set$SUMMARY_BASE_STAT_PATH) || !is.character(set$SUMMARY_BASE_STAT_PATH)){
			stop("SUMMARY_BASE_STAT_PATH argument not a character vector")
		}
	}
	print("checking genotype file existence ..")
	lapply(c(1:length(set$GENOTYPE_FILE_PATHS_GZ)),function(j){
		print(paste("batch ",j,sep=""))
		lapply(set$CHR,function(i){	
			gzstr=(".gz")[set$GENOTYPE_FILE_PATHS_GZ[j]]
			filepath=paste(set$GENOTYPE_BASE_FILE_PATHS[j],i,".txt",gzstr,sep="")
			if(!file.exists(filepath)){
				aa=paste(filepath," file does not exists.",sep="")
				stop(aa)
			}
		})
	})
	### check VARIANCE_PERCENTAGE2KEEP
	if(!is.numeric(set$VARIANCE_PERCENTAGE2KEEP)){
		stop("VARIANCE_PERCENTAGE2KEEP argument not a numeric.")
	}
	if(!is.scalar(set$VARIANCE_PERCENTAGE2KEEP)){
		stop("VARIANCE_PERCENTAGE2KEEP should be single scalar.")
	}
	if(!is.numeric(set$NSAMPLES)){
		stop("NSAMPLES argument not a numeric.")
	}
	if(!is.scalar(set$NSAMPLES)){
		stop("NSAMPLES should be single scalar.")
	}
	if(set$NSAMPLES<=0){
		stop("NSAMPLES should be larger than 0")
	}

	if(set$VARIANCE_PERCENTAGE2KEEP<95 || set$VARIANCE_PERCENTAGE2KEEP>100){
		stop("VARIANCE_PERCENTAGE2KEEP not within range")
	}
	### check RANGE_AROUND_TSS
	if(!is.numeric(set$RANGE_AROUND_TSS)){
		stop("RANGE_AROUND_TSS argument not a numeric.")
	}
	if(!is.scalar(set$RANGE_AROUND_TSS)){
		stop("RANGE_AROUND_TSS should be single scalar.")
	}
	if(set$RANGE_AROUND_TSS<10000  || BAGEA_ANNOTATION_MAT$settings$RANGE_AROUND_TSS<set$RANGE_AROUND_TSS){
		stop("RANGE_AROUND_TSS has to be between 10'000 and BAGEA_ANNOTATION_MAT$settings$RANGE_AROUND_TSS")
	}
	### check MAF_CUTOFF
	if(!is.numeric(set$MAF_CUTOFF)){
		stop("MAF_CUTOFF argument not a numeric.")
	}
	if(!is.scalar(set$MAF_CUTOFF)){
		stop("MAF_CUTOFF should be single scalar.")
	}
	if(set$MAF_CUTOFF<0 || set$MAF_CUTOFF>=0.5){
		stop("MAF_CUTOFF not within range")
	}
	if(!is.numeric(set$CHR) || !is.vector(set$CHR)){
		stop("CHR not set to numeric.")
	}
	if(length(unique(set$CHR))!=length(set$CHR) || length(intersect(c(1:22),set$CHR))!=length(set$CHR)){
		stop("CHR has to be a unique vector of elements between 1 and 22.")
	}
	### check SINGLESNP_MINPVALUE_CUTOFF
	if(!is.numeric(set$SINGLESNP_MINPVALUE_CUTOFF)){
		stop("SINGLESNP_MINPVALUE_CUTOFF argument not a numeric.")
	}
	if(!is.scalar(set$SINGLESNP_MINPVALUE_CUTOFF)){
		stop("SINGLESNP_MINPVALUE_CUTOFF should be single scalar.")
	}
	if(set$SINGLESNP_MINPVALUE_CUTOFF<0 || set$SINGLESNP_MINPVALUE_CUTOFF>1){
		stop("SINGLESNP_MINPVALUE_CUTOFF not within range")
	}
}


intersect_order_preserving=function(x,y){
	x1=x[is.element(x,y)]
	return(x1)
}

my_caught_svd=function(X){
	Xsvd = tryCatch({
			    Xsvd=svd(X)
			}, error = function(e) {
			    Xsvd_tmp=my_caught_svd(X/10)
			    #### my svd crashed here. a simple division by 10 solved the issue
			    Xsvd_tmp$d=Xsvd_tmp$d*10
			    return(Xsvd_tmp)
			})
}

#' prepare observed data for bagea run using external LD (e.g. 1KG). 
#'
#' Takes in external reference genotype data (for instance 1KG) and eqtl summary statistics data and combines them with a \code{bagea_annotation_mat} object to produce  a \code{bagea_observed_data} object that serves as input to \code{run_bagea}.
#' Genotype data has to be supplied indirectly by specifying the file paths of the chromosome-wise dosage files. The same goes for summary statistics data. See Details below for a description of the procedure. 
#' @return
#' An R-object of class bagea_observed_data; an R-list contatining 4 elements. if the phenotype was simulated, the list contains 6 elements (see below).
#' \describe{\item{\code{obs_list}}{A named list of length n where n is the number of genes processed. Each element contains the data for a particular gene that is needed to run \emph{bagea}. If the phenotype was simulated, it also contains the unobserved random variables.}
#' \item{\code{settings}}{A list of parameter settings with which prepare_observed_data was run.(of bagea_annotation_mat object, only the settings-field is stored and of the SNPS2KEEP list only the length).}
#' \item{\code{simulation_unobs_list}}{Only produced if phenotype was simulated, The list contains all the unobserved random variables that are not gene specific. Additionally, it conatins all the hyperparameters.}
#' \item{\code{removed_snps}}{A list of snps that were removed on each chromosome due to low maf (see argument MAF_CUTOFF).}
#' \item{\code{D_mat}}{A matrix that maps the index of the meta-annotation group indices to the indices of the directed annotations.}
#' \item{\code{D_names_list}}{A list of names that contains to the meta-annotation names.}}
#' @param BAGEA_ANNOTATION_MAT
#' An object of class bagea_annotation_mat produced by \code{\link{prepare_annotation_mapping_mat}}.
#' @param ANNOTATION_LIST
#' Character vector or NULL. If not NULL, has to be a selection of column names found in \code{BAGEA_ANNOTATION_MAT$annotation_mat}. Specifies which subset of genome annotations are used. If unspecified, all available annotations are used.
#' @param ANNOTATION_SHIFT_LIST
#' Character vector or NULL. If not NULL, has to be a selection of column names found in \code{BAGEA_ANNOTATION_MAT$annotation_shift_mat}. Specifies which subset of directed genome annotations are used (i.e. which annotations can be used to build up the V-matrix). If unspecified, all available annotations are used.
#' @param ANNOTATION_SHIFTWEIGHT_LIST
#' Character vector or NULL. If not NULL, has to be selection of column names found in \code{BAGEA_ANNOTATION_MAT$annotation_mat}. Specifies which subset of undirected genome annotations are used for weighting the SHIFT component (i.e. which annotations can be used to build up the F-matrix). If unspecified, all available annotations are used.
#' @param GENOTYPE_BASE_FILE_PATHS
#' Character vector. i-th element gives filepath to a directory of dosage files for batch i. See below for details on formatting the path as well as on how to format the dosage files.
#' @param GENOTYPE_FILE_PATHS_GZ
#' Boolean vector.  Set i-th element true if the genotype files for batch i are gzipped. See below for details.
#' @param SUMMARY_BASE_STAT_PATH 
#' String. Gives part of the filepath to summary statistics data.table in tab separated format. Three columns are expected. The first column should be named \code{gene} and should contain gene hgnc symbols. The second column named \code{snps} should contain snpids. The third column should be named \code{statistics} and should contain z-scores. The paths are assembled as \code{\{SUMMARY_BASE_STAT_PATH\}\{chr\}.txt,} where \{chr\} is a number between 1 to 22.
#' @param NSAMPLES
#' number of samples in the data set that were used to calculate the summary statistics.
#' @param VARIANCE_PERCENTAGE2KEEP
#' Scalar between BAGEA 95 and 100. BAGEA allows to approximate the observed LD matrix with a lower rank approximation, by removing eigen vectors that do not contribute substantially to the LD matrix.  VARIANCE_PERCENTAGE2KEEP sets the amount of variance that the kept eigenvectors have to capture.
#' @param RANGE_AROUND_TSS
#' scalar between 10'000 and 500'000. 
#' Window around transcription start site for which snps are considered.
#' @param MAF_CUTOFF
#' Scalar. SNPs with observed maf below or at MAF_CUTOFF will be ignored. has to be between 0 and 0.5. If multiple genotype batchtes are combined, the cutoff is transformed into a variance cutoff assuming hardy-weinberge equilibrium and applied to the combined dataset after mean removal of each batch.
#' @param CHR: 
#' A vector of unique integer numbers between 1 and 22. Gives of chromosome to process. Default is c(1:22) 
#' @param NCORES
#' Positive integer scalar. Number of cores to be used when loading the data.
#' @param SNPS2KEEP
#' Character vector or NULL. if not NULL all SNPS not in SNPS2KEEP will be removed.
#' @param ADD_SNP_FREQS2ANNOT
#' Boolean. if TRUE, SNP frequency will be added to the genomic annotations.
#' @param SINGLESNP_MINPVALUE_CUTOFF
#' A scalar set between 0 and 1. If smaller than 1, genes where no single SNP p-value reaches below this threshold is removed. 
#' @param SHIFT_METAINFO_TBL
#' A data frame. Defines groups shift annotations for modeling group-wise priors. The following has to hold: one column has to be named annotation_id and contain all elements in ANNOTATION_SHIFT_LIST uniquely. Every other column must contain names one particular grouping structure for directed genome annotations are used (for instance which cell type the particular annotation is mapped to).
#' @details Genotype data has to be supplied indirectly by specifying the file paths chromosome wise dosage files. The file paths are specified via 2 arguments.
#' GENOTYPE_BASE_FILE_PATHS and GENOTYPE_FILE_PATHS_GZ, where the i-th element of this vector specifies the paths for batch i.
#' Assuming that only one batch are analysed, then the paths are assembled as \code{\{GENOTYPE_BASE_FILE_PATHS\}\{chr\}.txt,} where \{chr\} is a number between 1 to 22 and \code{GENOTYPE_FILE_PATHS_GZ[i]} specifies whether files in batch i are gzipped.
#' For instance, if GENOTYPE_BASE_FILE_PATHS is \code{~/dosages/mygenotypes_chr,} GENOTYPE_FILE_PATHS_GZ is set to FALSE, then the function looks for 22 files of the format \code{~/dosages/mygenotypes_chr1.txt} to \code{~/dosages/mygenotypes_chr22.txt.} If GENOTYPE_FILE_PATHS_GZ is set to TRUE, the function will look for \code{~/dosages/mygenotypes_chr22.txt.gz} etc instead.
#' Likewise for the chromosome-wise summary statistics, the paths are assembled as \code{\{SUMMARY_BASE_STAT_PATH\}\{chr\}.txt,} where \{chr\} is a number between 1 to 22.
#' All dosage files need to be numeric tab-separated tables (no NAs allowed), except the first column and row. The first row gives the column names, which are sample ids except for the first. The first column should be named 'snpid ' and should contain snp ids in rs numbers. All other column should contain a individuals' SNP dosages.
#' If multiple batches are combined in the analysis,then GENOTYPE_BASE_FILE_PATHS, GENOTYPE_FILE_PATHS_GZ, should be given as vectors, where element i in the vector is the respective parameter setting for the i-th batch and formatting follows the same rules as for one batch. 
#' @export
prepare_observed_1KG_data=function(
	BAGEA_ANNOTATION_MAT=NULL,
	GENOTYPE_BASE_FILE_PATHS=NULL,
	GENOTYPE_FILE_PATHS_GZ=rep(TRUE,length(GENOTYPE_BASE_FILE_PATHS)),
	SUMMARY_BASE_STAT_PATH=NULL,
	NSAMPLES=NULL,
	VARIANCE_PERCENTAGE2KEEP=100,
	RANGE_AROUND_TSS=2e5,
	ANNOTATION_LIST=colnames(BAGEA_ANNOTATION_MAT$annotation_mat),
	ANNOTATION_SHIFT_LIST=colnames(BAGEA_ANNOTATION_MAT$annotation_shift_mat),
	ANNOTATION_SHIFTWEIGHT_LIST=colnames(BAGEA_ANNOTATION_MAT$annotation_mat),
	MAF_CUTOFF=0.05,
	CHR=c(1:22),
	NCORES=1,
	SNPS2KEEP=NULL,
	ADD_SNP_FREQS2ANNOT=FALSE,
	SINGLESNP_MINPVALUE_CUTOFF=1,
	SHIFT_METAINFO_TBL
	){
	###########################
	### set params
	###########################
	HYPER_PARAM_LIST=NULL
	if(is.null(GENOTYPE_BASE_FILE_PATHS)){
		BAGEA_PATH=Sys.getenv("BAGEA_PATH")
		if(BAGEA_PATH==""){
			stop("GENOTYPE_BASE_FILE_PATHS set to NULL but BAGEA_PATH not set.")		
		}
		GENOTYPE_BASE_FILE_PATHS=paste(BAGEA_PATH,"/Data/1KG/1KG_eur_chr",sep="")

	}
	settings=get_settingsout_list_obs1KG(
		BAGEA_ANNOTATION_MAT=BAGEA_ANNOTATION_MAT,
		ANNOTATION_LIST=ANNOTATION_LIST,
		ANNOTATION_SHIFT_LIST=ANNOTATION_SHIFT_LIST,
		GENOTYPE_BASE_FILE_PATHS=GENOTYPE_BASE_FILE_PATHS,
		GENOTYPE_FILE_PATHS_GZ=GENOTYPE_FILE_PATHS_GZ,
		SUMMARY_BASE_STAT_PATH=SUMMARY_BASE_STAT_PATH,
		NSAMPLES=NSAMPLES,
		VARIANCE_PERCENTAGE2KEEP=VARIANCE_PERCENTAGE2KEEP,
		RANGE_AROUND_TSS=RANGE_AROUND_TSS,
		MAF_CUTOFF=MAF_CUTOFF,
		HYPER_PARAM_LIST=HYPER_PARAM_LIST,
		CHR=CHR,
		NCORES=NCORES,
		SNPS2KEEP=SNPS2KEEP,
		ADD_SNP_FREQS2ANNOT=ADD_SNP_FREQS2ANNOT,
		SINGLESNP_MINPVALUE_CUTOFF=SINGLESNP_MINPVALUE_CUTOFF,
		ANNOTATION_SHIFTWEIGHT_LIST=ANNOTATION_SHIFTWEIGHT_LIST,
		SHIFT_METAINFO_TBL=SHIFT_METAINFO_TBL)
	###########################
	### check parameters
	###########################
	check_settings_list_obs1KG(
		BAGEA_ANNOTATION_MAT=BAGEA_ANNOTATION_MAT,
		SNPS2KEEP=SNPS2KEEP,
		set=settings
	)
	###########################
	### END: check parameters
	###########################
	###########################
	### END: set params
	###########################
	###########################
	# prepare annotation data
	###########################
	bagea_annotation_mat=subset_annotation_mat_cols(BAGEA_ANNOTATION_MAT,ANNOTATION_LIST,ANNOTATION_SHIFT_LIST,ANNOTATION_SHIFTWEIGHT_LIST)
	rm(BAGEA_ANNOTATION_MAT)
	###########################
	# prepare D_mat from SHIFT_METAINFO_TBL
	###########################	
	#subset_metainfo_tbl
	D_mat_outl=setup_D_mat_with_checks(SHIFT_METAINFO_TBL,bagea_annotation_mat)
	D_mat=D_mat_outl[["D_mat"]]
	D_names_list=D_mat_outl[["D_names_list"]]
	rm(D_mat_outl)
	###########################
	# END prepare D_mat from SHIFT_METAINFO_TBL
	###########################
	###########################
	### simulate global params
	###########################
	if(do_simulation(HYPER_PARAM_LIST)){
		simulated_params=simulate_params_from_hyperparams(HYPER_PARAM_LIST,D_mat,ANNOTATION_LIST,ANNOTATION_SHIFT_LIST,ANNOTATION_SHIFTWEIGHT_LIST)
	}else{
		simulated_params=NULL
	}
	###########################
	### END: simulate global params
	###########################
	###########################
	### load data
	###########################
	###########################
	# END: expression data
	###########################
	if(length(CHR)!=1){
		if(NCORES!=1){
			tt=mclapply(CHR, prepare_data_per_chr,
				bagea_annotation_mat=bagea_annotation_mat,
				SNPS2KEEP=SNPS2KEEP,
				simulated_params=simulated_params,
				set=settings,
				mc.cores=NCORES,	
				mc.preschedule=F)
		}else{
			tt=lapply(CHR, prepare_data_per_chr,
				bagea_annotation_mat=bagea_annotation_mat,	
				SNPS2KEEP=SNPS2KEEP,
				simulated_params=simulated_params,
				set=settings
				)
		}
		all_removed_snps_list=lapply(c(1:length(CHR)),function(i){attr(tt[[i]], "all_snps_below_mafcutoff")})
	}else{
		tt=lapply(c(CHR),prepare_data_per_chr,
				bagea_annotation_mat=bagea_annotation_mat,
				SNPS2KEEP=SNPS2KEEP,
				simulated_params=simulated_params,
				set=settings				
				)
		all_removed_snps_list=lapply(c(1:22),function(i){NULL})
		all_removed_snps_list[[CHR]]=attr(tt[[1]],"all_snps_below_mafcutoff")
	}
	obs_list=unlist(tt,recursive=FALSE)
	obs_out=list()
	obs_out[["obs_list"]]=obs_list
	obs_out[["settings"]]=settings
	obs_out[["simulation_unobs_list"]]=simulated_params
	obs_out[["removed_snps"]]=all_removed_snps_list
	obs_out[["D_mat"]]=D_mat
	obs_out[["D_names_list"]]=D_names_list
	class(obs_out)=c("bagea_observed_data",class(obs_out))
	return(obs_out)
}
