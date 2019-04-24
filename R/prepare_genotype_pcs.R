get_params_for_global_cov=function(x){
	xtx=t(x)%*%x
	xs=colSums(x)
	len=dim(x)[1]
	outl=list()
	outl[["xtx"]]=xtx
	outl[["xs"]]=xs
	outl[["len"]]=len
	return(outl)
}

assemble_cov=function(xtx,xs,len){
	xm=xs/len
	ret=xtx/len
	ret=ret-(xm%*%t(xm))
	return(ret)
}


#' prepare genotype pcs to control population structure. 
#'
#' Uses genotype data to compute batch-wise eigenvectors which can be used to control population structure when running \code{prepare_observed_data}. Genotype data has to be supplied indirectly by specifying the file paths of the chromosome-wise dosage files. See Details below for a description of the 
#' @return
#' An R-object of class bagea_genotype_pcs_data; an R-list contatining 2 elements.
#' genotype_pcs:
#' A list of matrices that contains the relevant PCs per batch in the columns. 
#' settings:
#' A list of parameter settings with which prepare_genotype_pcs_data was run
#' @param GENOTYPE_BASE_FILE_PATHS character vector. i-th element gives filepath to a directory of dosage files for batch i. See below for details on formatting the path as well as on how to format the dosage files.
#' @param GENOTYPE_FILE_PATHS_GZ boolean vector.  set i-th element true if the genotype files for batch i are gzipped. See below for details.
#' @param MAF_CUTOFF SNPs with observed maf below or at MAF_CUTOFF will be ignored. has to be between 0 and 0.5. The cutoff is transformed into a variance cutoff assuming hardy-weinberg equilibrium and applied to the combined dataset after mean removal of each batch.  
#' @param CHR A scalar between 0 and 22. Gives of chromosome to process. If set to 0, all chromosomes  are processed.
#' @param NCORES positive integer scalar. Number of cores to be used when loading the data.
#' @param SNPS2KEEP character vector or NULL. if not NULL all SNPS not in SNPS2KEEP will be removed.
#' @details genotype data has to be supplied indirectly by specifying the file paths chromosome wise dosage files. The file paths are specified via 2 arguments.
#' GENOTYPE_BASE_FILE_PATHS and GENOTYPE_FILE_PATHS_GZ, where the i-th element of this vector specifies the paths for batch i.
#' Assuming that only one batch are analysed, then the paths are assembled as \preformatted{\{GENOTYPE_BASE_FILE_PATHS\}\{chr\}.txt,} where \{chr\} is a number between 1 to 22 and \code{GENOTYPE_FILE_PATHS_GZ[i]} specifies whether files in batch i are gzipped.
#' For instance, if GENOTYPE_BASE_FILE_PATHS is \preformatted{~/dosages/mygenotypes_chr,} GENOTYPE_FILE_PATHS_GZ is set to FALSE, then the function looks for 22 files of the format \preformatted{~/dosages/mygenotypes_chr1.txt} to \preformatted{~/dosages/mygenotypes_chr22.txt.} If GENOTYPE_FILE_PATHS_GZ is set to TRUE, the function will look for \preformatted{~/dosages/mygenotypes_chr22.txt.gz} etc instead.
#' All dosage files need to be numeric tab-separated tables (no NAs allowed), except the first column and row. The first row gives the column names, which are sample ids except for the first. The first column should be named 'snpid ' and should contain snp ids in rs numbers. All other column should contain a individuals' SNP dosages.
#' If multiple batches are combined in the analysis,then GENOTYPE_BASE_FILE_PATHS, GENOTYPE_FILE_PATHS_GZ,EXPRESSION_FILE_PATHS should be given as vectors, where element i in the vector is the respective parameter setting for the i-th batch and formatting follows the same rules as for one batch. 
#' @export
prepare_genotype_pcs_data=function(GENOTYPE_BASE_FILE_PATHS=NULL,GENOTYPE_FILE_PATHS_GZ=rep(FALSE,length(GENOTYPE_BASE_FILE_PATHS)),MAF_CUTOFF=0.05,CHR=0,NCORES=1,SNPS2KEEP=NULL){
	settings=list()
	settings[["GENOTYPE_BASE_FILE_PATHS"]]=GENOTYPE_BASE_FILE_PATHS
	settings[["GENOTYPE_FILE_PATHS_GZ"]]=GENOTYPE_FILE_PATHS_GZ
	settings[["MAF_CUTOFF"]]=MAF_CUTOFF
	settings[["CHR"]]=CHR
	settings[["NCORES"]]=NCORES
	settings[["n_snps2keep"]]=length(SNPS2KEEP)
	if(length(CHR)==1 && CHR==0){
		chrstorun=c(1:22)
	}else{
		chrstorun=CHR
	}
	print("checking genotype file existence ..")
	lapply(c(1:length(GENOTYPE_FILE_PATHS_GZ)),function(j){
		print(paste("batch ",j,sep=""))
		lapply(chrstorun,function(i){	
			gzstr=(".gz")[GENOTYPE_FILE_PATHS_GZ[j]]
			filepath=paste(GENOTYPE_BASE_FILE_PATHS[j],i,".txt",gzstr,sep="")
			if(!file.exists(filepath)){
				aa=paste(filepath," file does not exists.",sep="")
				stop(aa)
			}
		})
	})
	###########################
	prepare_data_per_chr=function(chr){
		###########################
		# restrict annotation data to single chr
		###########################
		chrstr=paste("chr",chr,sep="")
		genotype_mats=lapply(c(1:length(GENOTYPE_BASE_FILE_PATHS)),function(i){
			load_genotype_mat(chr,GENOTYPE_BASE_FILE_PATHS=GENOTYPE_BASE_FILE_PATHS[i],GENOTYPE_FILE_PATHS_GZ=GENOTYPE_FILE_PATHS_GZ[i],snps2keep=NULL)
		})
		if(!is.null(SNPS2KEEP)){
			snps2keep_chr=intersect(SNPS2KEEP,rownames(genotype_mats[[1]]))
		}else{
			snps2keep_chr=NULL
		}
		genotype_mats=intersect_snps_fully(genotype_mats,snps2keep_chr)
		###########################
		# restrict expression mat to single chr
		###########################
		###########################
		# END: order samples for genotypes.
		###########################
		dd=lapply(genotype_mats,rowSums)
		dd=do.call("cbind",dd)
		dd2=rowSums(dd)
		nsamples=sum(sapply(genotype_mats,function(x){dim(x)[2]}))
		dd=dd2/(nsamples*2)
#		genotype_mat=combine_genotype_mats(genotype_mats)
		all_snps_below_mafcutoff_inds=(dd<MAF_CUTOFF | dd>(1-MAF_CUTOFF))
		all_snps_below_mafcutoff=rownames(genotype_mats[[1]])[all_snps_below_mafcutoff_inds]
		###########################
		# ADD SNPs that are homozygous in one sample
		###########################			
		dd=lapply(genotype_mats,function(x){
			y=x[,1]%*%t(rep(1,dim(x)[2]))
			rowSums((x-y)^2)==0
		})
		dd=do.call("cbind",dd)
		dd2=rowSums(dd)>0
		all_snps_below_mafcutoff=union(all_snps_below_mafcutoff,rownames(genotype_mats[[1]])[dd2])
		###########################
		# remove low maf snps
		###########################		
		genotype_mats=lapply(genotype_mats,function(x){
			snpnames=rownames(x)
			tokeep=!is.element(snpnames,all_snps_below_mafcutoff)
			return(x[tokeep,,drop=FALSE])
		})
		###########################
		# END: remove low maf snps
		###########################
		###########################
		# scale genotype mats
		###########################		
		genotype_mats=lapply(genotype_mats,function(x){
			x=t(x)
			x=my_scale(x)
			x=t(x)
			genotype_mats=t(my_scale(t(x)))
		})
		###########################
		# END scale genotype mats
		###########################	
		myout=lapply(genotype_mats,get_params_for_global_cov)
		return(myout)
	}
	if(NCORES!=1){
		tt=mclapply(chrstorun,prepare_data_per_chr,mc.cores=NCORES,mc.preschedule=F)
	}else{
		tt=lapply(chrstorun,prepare_data_per_chr)
	}
	names(tt)=chrstorun
	all_covs=lapply(c(1:length(GENOTYPE_BASE_FILE_PATHS)),function(i){
		nchrs=length(tt)
		for(chrname in names(tt)){
			if(chrname==names(tt)[1]){
				xtx=tt[[chrname]][[i]][["xtx"]]
				len=tt[[chrname]][[i]][["len"]]
				xs=tt[[chrname]][[i]][["xs"]]
			}else{
				xtx=xtx+tt[[chrname]][[i]][["xtx"]]
				len=len+tt[[chrname]][[i]][["len"]]
				xs=xs+tt[[chrname]][[i]][["xs"]]
			}
		}		
		mycov=assemble_cov(xtx,xs,len)
		return(mycov)
	})	
	###########################
	# END: get covs genotype mats
	###########################
	###########################
	# get pcs mats
	###########################
	vector_mats=lapply(all_covs,function(x){
		mypcs=eigen(x)$vectors
		rownames(mypcs)=rownames(all_covs)
		return(mypcs)
	})
	###########################
	# END: get pcs mats
	###########################
	obs_out=list()
	obs_out[["genotype_pcs"]]=vector_mats
	obs_out[["settings"]]=settings
	class(obs_out)=c("bagea_genotype_pcs_data",class(obs_out))
	return(obs_out)
}