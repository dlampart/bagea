assemble_data_mats=function(expression_mats,genotype_mats,bagea_annotation_mat,do_simulation){
	expression_mats_scaled=list()
	outl=list()
	if(!do_simulation){
		for(i in c(1:length(genotype_mats))){
			genotype_mat=genotype_mats[[i]]
			expression_mat=expression_mats[[i]]
			expression_mat_chr=subset_expression_mat_chr(expression_mat=expression_mat,bagea_annotation_mat_chr=bagea_annotation_mat)
			###########################
			# END: restrict expression mat to single chr
			###########################
			###########################
			# order samples for genotypes.
			###########################
			print("order samples")
			intersected_names=intersect(colnames(genotype_mat),colnames(expression_mat_chr))
			print(paste("number of overlapping samples:",length(intersected_names)))
			genotype_matcher=match(intersected_names,colnames(genotype_mat))
			expression_matcher=match(intersected_names,colnames(expression_mat_chr))
			expression_mat_sorted=expression_mat_chr[,expression_matcher]
			if(length(genotype_mats)>1){
				print("scale expression to mean 0 sd 1 per batch ...")
			}
			expression_mat_scaled=t(expression_mat_sorted)
			print("scale expression to mean 0 sd 1 ..")
			expression_mat_scaled=my_scale(expression_mat_scaled,use_minus1=FALSE)
			expression_mat_scaled=t(expression_mat_scaled)
			expression_mats_scaled[[i]]=expression_mat_scaled
			genotype_mats[[i]]=genotype_mats[[i]][,genotype_matcher]			
		}
		if(length(genotype_mats)>1){
			print("combine expression batches..")
		}
		expression_mat_scaled=do.call("cbind",expression_mats_scaled)
		outl[["expression_mat_scaled"]]=expression_mat_scaled
	}
	genotype_mat=combine_genotype_mats(genotype_mats)
	outl[["genotype_mat"]]=genotype_mat
	return(outl)
}

subset_expression_mat_chr=function(expression_mat,bagea_annotation_mat_chr){
	symbols_on_chromosome=unique(bagea_annotation_mat_chr$snp_tss_map[,symbol])
	expression_geneids=rownames(expression_mat)
	exp_on_chr_ids=which(is.element(expression_geneids,symbols_on_chromosome))
	expression_mat_chr=expression_mat[exp_on_chr_ids,]
	return(expression_mat_chr)
}

get_preprocess_expression_mats=function(set){
	if(!is.null(set$EXPRESSION_FILE_PATHS)){
		expression_mats=lapply(set$EXPRESSION_FILE_PATHS,function(filename){load_expression_mat(filename)})
		if(set$NPCS>0 && !is.null(set$GENOTYPE_PC_LIST)){
			mylen=length(expression_mats)
			expression_mats=lapply(c(1:mylen),function(i){
				out=regress_out_pcs(t(expression_mats[[i]]),set$GENOTYPE_PC_LIST[["genotype_pcs"]][[i]],set$NPCS)
				out=t(out)				
				return(out)
			})
		}
		expression_mats=intersect_expression_genes_fully(expression_mats)
	}
	return(expression_mats)
}

check_settings_list_obs=function(BAGEA_ANNOTATION_MAT,SNPS2KEEP,
	set){
	if(set$NPCS>0 & !is.null(set$GENOTYPE_PC_LIST) && !sum(class(set$GENOTYPE_PC_LIST)=="bagea_genotype_pcs_data")){
		stop("GENOTYPE_PC_LIST is set but is not a bagea_genotype_pcs_data object.")
	}
	if(!is.scalar(set$KEEP_X_y) || !is.logical(set$KEEP_X_y)){
		stop("KEEP_X_y has to be a single boolean")
	}
	if(set$NPCS>0 & !is.null(set$GENOTYPE_PC_LIST) && !sum(class(set$GENOTYPE_PC_LIST)=="bagea_genotype_pcs_data")){
		stop("GENOTYPE_PC_LIST is set but is not a bagea_genotype_pcs_data object.")
	}
	if(is.null(set$HYPER_PARAM_LIST) && is.null(set$EXPRESSION_FILE_PATHS)){
		stop("neither HYPER_PARAM_LIST for simulation supplied nor an expression file path")
	}
	if(!is.null(set$HYPER_PARAM_LIST) && !is.null(set$EXPRESSION_FILE_PATHS)){
		stop("both HYPER_PARAM_LIST for simulation  and an expression file path is set. set one to NULL")
	}
	if(!is.null(set$HYPER_PARAM_LIST) && !sum(class(set$HYPER_PARAM_LIST)=="bagea_hyperparameter_list")){
		stop("HYPER_PARAM_LIST is set but is not a bagea_hyperparameter_list object.")
	}
	if(set$ADD_SNP_FREQS2ANNOT && do_simulation){
		stop("errsmocps:cannot simulate with ADD_SNP_FREQS2ANNOT flag.")
	}
	if(sum(class(BAGEA_ANNOTATION_MAT)=="bagea_annotation_mat")==0){
		stop("BAGEA_ANNOTATION_MAT argument is not a bagea_annotation_mat object.")
	}
	### check BAGEA_ANNOTATION_MAT
	if(sum(class(BAGEA_ANNOTATION_MAT)=="bagea_annotation_mat")==0){
		stop("BAGEA_ANNOTATION_MAT argument is not a bagea_annotation_mat object.")
	}
	if(is.null(colnames(BAGEA_ANNOTATION_MAT[["annotation_mat"]]))){
		stop("BAGEA_ANNOTATION_MAT$annotation_mat has no colnames")
	}
	if(is.null(colnames(BAGEA_ANNOTATION_MAT[["annotation_shift_mat"]]))){
			stop("BAGEA_ANNOTATION_MAT$annotation_shift_mat has no colnames")
	}	
	if(!is.character(set$ANNOTATION_LIST) || !is.vector(set$ANNOTATION_LIST)){
		stop("ANNOTATION_LIST argument is not a character vector")
	}
	if(!is.character(set$ANNOTATION_SHIFT_LIST) || !is.vector(set$ANNOTATION_SHIFT_LIST)){
		stop("ANNOTATION_SHIFT_LIST argument is not a character vector")
	}
	if(!is.null(set$EXPRESSION_FILE_PATHS)){
		if(!is.vector(set$EXPRESSION_FILE_PATHS) || !is.character(set$EXPRESSION_FILE_PATHS)){
			stop("EXPRESSION_FILE_PATHS argument not a character vector")
		}
		lapply(set$EXPRESSION_FILE_PATHS, function(expresion_file_path){
			if(!file.exists(expresion_file_path)){
				stop(paste("in EXPRESSION_FILE_PATHS: ",expresion_file_path," file does not exists",sep=""))
			}
		})
	}	
	if(!is.vector(set$GENOTYPE_BASE_FILE_PATHS) || !is.character(set$GENOTYPE_BASE_FILE_PATHS)){
			stop("GENOTYPE_BASE_FILE_PATHS argument not a character vector")
	}
	### check GENOTYPE_FILE_PATHS_GZ
	if(!is.vector(set$GENOTYPE_FILE_PATHS_GZ) || !is.logical(set$GENOTYPE_FILE_PATHS_GZ)){
		stop("GENOTYPE_FILE_PATHS_GZ has to be a single boolean")
	}
	if(!is.null(set$EXPRESSION_FILE_PATHS)){
		if(length(set$EXPRESSION_FILE_PATHS)!=length(set$GENOTYPE_BASE_FILE_PATHS)){
			stop("not the same number of EXPRESSION_FILE_PATHS and GENOTYPE_BASE_FILE_PATHS supplied")
		}
	}
	if(length(set$GENOTYPE_FILE_PATHS_GZ)!=length(set$GENOTYPE_BASE_FILE_PATHS)){
		stop("not the same number of GENOTYPE_FILE_PATHS_GZ and GENOTYPE_BASE_FILE_PATHS supplied")
	}
	### check KEEP_X_y
	if(!is.scalar(set$KEEP_X_y) || !is.logical(set$KEEP_X_y)){
		stop("KEEP_X_y has to be a single boolean")
	}
	if(!is.null(set$SNPS2KEEP) && ( !is.vector(set$SNPS2KEEP) || !is.character(set$SNPS2KEEP))){
		stop("SNPS2KEEP has to be null or a character vector")
	}
	if(!is.scalar(set$ADD_SNP_FREQS2ANNOT) || !is.logical(set$ADD_SNP_FREQS2ANNOT)){
		stop("ADD_SNP_FREQS2ANNOT has to be a single boolean")
	}
	check_posinteger(set$NCORES)
	### check GENOTYPE_BASE_FILE_PATHS
	if(!is.numeric(set$CHR) || !is.vector(set$CHR)){
		stop("CHR not set to numeric.")
	}
	if(length(unique(set$CHR))!=length(set$CHR) || length(intersect(c(1:22),set$CHR))!=length(set$CHR)){
		stop("CHR has to be a unique vector of elements between 1 and 22.")
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

## aggregating settings into single object.
get_settingsout_list_obs=function(BAGEA_ANNOTATION_MAT,ANNOTATION_LIST,ANNOTATION_SHIFT_LIST,GENOTYPE_BASE_FILE_PATHS,GENOTYPE_FILE_PATHS_GZ,EXPRESSION_FILE_PATHS,VARIANCE_PERCENTAGE2KEEP,RANGE_AROUND_TSS,MAF_CUTOFF,HYPER_PARAM_LIST,CHR,KEEP_X_y,NCORES,SNPS2KEEP,ADD_SNP_FREQS2ANNOT,ADD_HERITABILTY_QUANTILE,SINGLESNP_MINPVALUE_CUTOFF,ANNOTATION_SHIFTWEIGHT_LIST,SHIFT_METAINFO_TBL,GENOTYPE_PC_LIST,NPCS){
	settings=list()
	settings[["GENOTYPE_BASE_FILE_PATHS"]]=GENOTYPE_BASE_FILE_PATHS
	settings[["GENOTYPE_FILE_PATHS_GZ"]]=GENOTYPE_FILE_PATHS_GZ
	settings[["EXPRESSION_FILE_PATHS"]]=EXPRESSION_FILE_PATHS
	settings[["BAGEA_ANNOTATION_MAT"]]=BAGEA_ANNOTATION_MAT$settings
	settings[["ANNOTATION_LIST"]]=ANNOTATION_LIST
	settings[["ANNOTATION_SHIFT_LIST"]]=ANNOTATION_SHIFT_LIST
	settings[["MAF_CUTOFF"]]=MAF_CUTOFF
	settings[["RANGE_AROUND_TSS"]]=RANGE_AROUND_TSS
	settings[["VARIANCE_PERCENTAGE2KEEP"]]=VARIANCE_PERCENTAGE2KEEP
	settings[["HYPER_PARAM_LIST"]]=HYPER_PARAM_LIST
	settings[["CHR"]]=CHR
	settings[["KEEP_X_y"]]=KEEP_X_y
	settings[["NCORES"]]=NCORES
	settings[["n_snps2keep"]]=length(SNPS2KEEP)
	settings[["ADD_SNP_FREQS2ANNOT"]]=ADD_SNP_FREQS2ANNOT
	settings[["ADD_HERITABILTY_QUANTILE"]]=ADD_HERITABILTY_QUANTILE
	settings[["SINGLESNP_MINPVALUE_CUTOFF"]]=SINGLESNP_MINPVALUE_CUTOFF
	settings[["ANNOTATION_SHIFTWEIGHT_LIST"]]=ANNOTATION_SHIFTWEIGHT_LIST
	settings[["SHIFT_METAINFO_TBL"]]=SHIFT_METAINFO_TBL
	settings[["NPCS"]]=NPCS
	return(settings)
}

regress_out_pcs=function(exp_mat,pc_mat,npcs){
	expmat_names=rownames(exp_mat)
	pc_names=rownames(pc_mat)
	if(sum(expmat_names!=pc_names)!=0){
		stop("erropmc[p:expression mat and pc mat names do not agree.")
	}
	pc_mat_to_use=pc_mat[,c(1:npcs),drop=FALSE]
	out_mat=regress_out_vects(exp_mat,pc_mat_to_use, add_intercept=FALSE)
	out_mat=my_scale(out_mat)
	return(out_mat)
}


#' set_hyperparameter_list
#' prepares the set hyperparameters that are used in bagea run.
#' @return An R-object of class bagea_hyperparameter_list; an R-list contatining all hyper-parameters.
#' @param nu_names
#' character vector. These should refer to column names in the matrix \code{BAGEA_ANNOTATION_MAT$annotation_mat} that was used. Defines which undirected annotations should be used to construct F matrix.
#' @param myc
#' either NULL, or a numerical vector, giving the prior means for the hidden variable vector nu. If NULL, 0 will be used for all elements of nu, (and  1 for the intercept term). If a vector, it has to have the same length of nu_names plus one (for the elements of nu plus the intercept term). These will give specific values to each element of nu with the first element being the intercept which is not allowed to be 0 (choose 1 for the first element for easy interpretability).
#' @param gamma1
#' scalar > 0
#' @param tau1 
#' scalar > 0
#' @param phi1 
#' scalar > 0
#' @param phi2 
#' scalar > 0
#' @param xi1 
#' scalar > 0
#' @param xi2 
#' scalar > 0
#' @param lambda1 
#' scalar > 0
#' @param rho1 
#' scalar > 0
#' @param rho2 
#' scalar > 0
#' @param chi1 
#' scalar > 0
#' @param chi2 
#' scalar > 0
#' @param chi1 
#' scalar > 0
#' @param zeta1 
#' scalar > 0
#' @param zeta2 
#' scalar > 0
#' @param p
#' scalar > 0
#' @export
set_hyperparameter_list=function(nu_names=NULL,myc=NULL,gamma1=3,tau1=2,phi1=3,phi2=3,xi1=100,xi2=0.03,lambda1=sqrt(1E5),rho1=1E5,rho2=sqrt(1E5),p=5,chi1=1,chi2=0.0003,zeta1=1,zeta2=100){
	if(is.null(nu_names)){
		stop("nu_names currently have to be set.")
	}
	if(!is.null(myc)  && length(myc)!=(length(nu_names)+1)){
		stop("myc has to be either NULL, length of nu_names plus one (for the elements of nu plut the intercept term).")
	}
	hyperparameter_list=list()
	hyperparameter_list$gamma1=gamma1
	hyperparameter_list$tau1=tau1
	hyperparameter_list$lambda1=lambda1
	hyperparameter_list$phi1=phi1
	hyperparameter_list$phi2=phi2
	hyperparameter_list$rho1=rho1
	hyperparameter_list$rho2=rho2
	hyperparameter_list$xi1=xi1
	hyperparameter_list$xi2=xi2
	hyperparameter_list$zeta1=zeta1
	hyperparameter_list$zeta2=zeta2
	if(length(p)>1 || p>0){
		if(is.null(nu_names)){
			hyperparameter_list$p=p
		}else{
			hyperparameter_list$p=c(p*50000,rep(p,length(nu_names)))
		}
		if(is.null(myc)){
			myc=rep(0,length(hyperparameter_list$p))
			myc[1]=1
			hyperparameter_list$c=myc
		}else{
			hyperparameter_list$c=myc
			if(myc[1]==0){
				stop("first element in myc cannot be 0, as it is the intercept.")
			}
		}
		updated_nu_names=c("all_id2920242",nu_names)
		names(hyperparameter_list$c)=updated_nu_names
		names(hyperparameter_list$p)=updated_nu_names
		hyperparameter_list$chi1=chi1
		hyperparameter_list$chi2=chi2
	}
	class(hyperparameter_list)=c("bagea_hyperparameter_list",class(hyperparameter_list))
	return(hyperparameter_list)
}

add_annotation_params=function(hyperparameter_list,annotation_names){
	if(sum(class(hyperparameter_list)=="bagea_hyperparameter_list")==0){
		stop("hyperparameter_list argument is not a bagea_hyperparameter_list object.")
	}
	if(!is.character(annotation_names) || !is.vector(annotation_names)){
		stop("annotation_names argument is not a character vector")
	}
	myt=length(annotation_names)
	tt=rgamma(myt,hyperparameter_list$phi1,hyperparameter_list$phi2)
	names(tt)=annotation_names
	hyperparameter_list$as=tt
	return(hyperparameter_list)
}

add_annotation_shift_params=function(hyperparameter_list,annotation_shift_names,D_mat=NULL){
	if(sum(class(hyperparameter_list)=="bagea_hyperparameter_list")==0){
		stop("hyperparameter_list argument is not a bagea_hyperparameter_list object.")
	}
	if(!is.character(annotation_shift_names) || !is.vector(annotation_shift_names)){
		stop("annotation_shift_names argument is not a character vector")
	}
	mys=length(annotation_shift_names)
	if(is.null(D_mat)){
		omega=rnorm(mys,0,sqrt(1/hyperparameter_list$delta))
	}else{
		ups_l=hyperparameter_list$upsilon_list
		my_q=dim(D_mat)[2]
		res=rep(1,mys)
		for(j in c(1:my_q)){
			res=res*ups_l[[j]][D_mat[,j]]
		}
		omega=rnorm(mys,0,1)
		omega=omega/sqrt(res)
	}
	if(!is.null(hyperparameter_list[["preset_omega"]])){
		print("presetting omega..")
		if(length(hyperparameter_list[["preset_omega"]])!=length(omega)){
			stop("Trying to preset omega. Does not have a right length.")
		}
		omega=hyperparameter_list[["preset_omega"]]
	}
	hyperparameter_list$omega=omega
	return(hyperparameter_list)
}

add_annotation_nu_params=function(hyperparameter_list,annotation_nu_names){
	if(sum(class(hyperparameter_list)=="bagea_hyperparameter_list")==0){
		stop("hyperparameter_list argument is not a bagea_hyperparameter_list object.")
	}
	if(!is.character(annotation_nu_names) || !is.vector(annotation_nu_names)){
		stop("annotation_shift_names argument is not a character vector")
	}
	annotation_nu_names=c("background",annotation_nu_names)
	mys=length(annotation_nu_names)
	p_len=length(hyperparameter_list$p)
	if(p_len!=1 && mys!=p_len){
		stop("p annot_names and  p length do not agree.")
	}
	#nu=rnorm(mys,1,sqrt(1/hyperparameter_list$p))
	nu_noshift=rnorm(mys,0,1)
	nu=(nu_noshift/sqrt(hyperparameter_list$p)+hyperparameter_list$c)
	names(nu)=annotation_nu_names
	if(!is.null(hyperparameter_list[["preset_nu"]])){
		print("presetting nu..")
		if(length(hyperparameter_list[["preset_nu"]])!=length(nu)){
			stop("Trying to preset nu. Does not have a right length.")
		}
		if(sum(names(hyperparameter_list[["preset_nu"]])!=annotation_nu_names)!=0){
			stop("Trying to preset nu. Names are incorrect.")
		}
		nu=hyperparameter_list[["preset_nu"]]		
	}
	hyperparameter_list$nu=nu
	return(hyperparameter_list)
}

get_global_param_list=function(hyperparameter_list,D_mat=NULL){
	lambda2=rgamma(1,hyperparameter_list$rho1,hyperparameter_list$rho2)
	tau2=rgamma(1,hyperparameter_list$xi1,hyperparameter_list$xi2)
	hyperparameter_list$lambda2=lambda2
	hyperparameter_list$tau2=tau2
	if(!is.null(hyperparameter_list$p)){
		if(is.null(D_mat)){
			delta=rgamma(1,hyperparameter_list$chi1,hyperparameter_list$chi2)
			hyperparameter_list$delta=delta
		}else{
			my_q=dim(D_mat)[2]
			my_s=dim(D_mat)[1]
			chi1_vec=hyperparameter_list$chi1
			if(is.null(hyperparameter_list$zeta1)){
				chi2_vec=hyperparameter_list$chi2
				if(length(chi2_vec)==1){
					chi2_vec=rep(chi2_vec,my_q)
				}
			}else{
				zeta1_vec=hyperparameter_list$zeta1
				zeta2_vec=hyperparameter_list$zeta2
				if(length(zeta1_vec)==1){
					zeta1_vec=rep(zeta1_vec,my_q)
				}
				if(length(zeta2_vec)==1){
					zeta2_vec=rep(zeta2_vec,my_q)
				}
				chi2_vec=rep(0,my_q)
				for(i in c(1:my_q)){
					chi2_vec[i]=rgamma(1,zeta1_vec[i],zeta2_vec[i])
				}
				hyperparameter_list$chi2j=chi2_vec
			}
			h_vec=apply(D_mat,2,max)
			if(length(chi1_vec)==1){
				chi1_vec=rep(chi1_vec,my_q)
			}
			if(length(chi1_vec)!=my_q){
				stop("erroaoscsdo:sl")
			}
			upsilon_list=list()
			for(j in c(1:length(h_vec))){
				upsilon_list[[j]]=rgamma(h_vec[j],chi1_vec[j],chi2_vec[j])
			}
			hyperparameter_list$upsilon_list=upsilon_list
		}
	}
	return(hyperparameter_list)
}

get_genewise_param_list=function(mapping_mat,mapping_shift_mat=NULL,hyperparameter_list,annotation_nu_names=NULL){
	kappa=rgamma(1,hyperparameter_list$tau1,hyperparameter_list$tau2)
	m=dim(mapping_mat)[1]
	aMat0=kronecker(rep(1,m),t(hyperparameter_list$as))
	gammai=exp(rowSums(log(aMat0)*mapping_mat))
	alpha=rgamma(m,hyperparameter_list$gamma1,gammai*kappa)
	if(is.null(mapping_shift_mat)){
		b=rnorm(m,0,sqrt(1/alpha))	
	}else{
		delta=hyperparameter_list$delta
		omega=hyperparameter_list$omega
		mu=matrix(mapping_shift_mat%*%omega)
		if(sum(class(hyperparameter_list)=="bagea_hyperparameter_bool_list")!=0){
			mypi=runif(m)<hyperparameter_list$p
			b=rnorm(m,mu*mypi,sqrt(1/alpha))
		}else{
			submat=cbind(TRUE,mapping_mat[,match(annotation_nu_names,colnames(mapping_mat))])
			if(m==1){
				submat=c(TRUE,mapping_mat[,match(annotation_nu_names,colnames(mapping_mat))])
				submat=as.matrix(t(submat))

			}
			Fnu=(submat%*%hyperparameter_list$nu)
			b_nonshift=rnorm(m,0,sqrt(1/alpha))
			b=b_nonshift+(mu*Fnu)[,1]
		}
	}
	lambda=rgamma(1,hyperparameter_list$lambda1,hyperparameter_list$lambda2)
	genewise_params=list()
	genewise_params$kappa=kappa
	genewise_params$gammai=gammai
	genewise_params$alpha=alpha
	genewise_params$b=b
	genewise_params$lambda=lambda
	if(sum(class(hyperparameter_list)=="bagea_hyperparameter_bool_list")!=0){
		genewise_params$mu=mu
		genewise_params$pi=mypi
	}else{
		if(!is.null(mapping_shift_mat)){
			genewise_params$mu=mu
			genewise_params$Fnu=Fnu		
		}
	}
	return(genewise_params)
}
	
simulate_y=function(X,genewise_param_list){
	n=dim(X)[2]
	lambda=genewise_param_list$lambda
	b=genewise_param_list$b	
	y=t(t(X)%*%b+(rnorm(n)/sqrt(lambda)))
	return(y)
}

add_heritability_quantile_mat=function(obs_list){
	all_he=sapply(obs_list,function(x){
			return(x$sigmasq_he)
	})
	myquantiles=quantile(all_he,c(1:4)/5)
	quantile_group=rep(0,length(all_he))
	for(i in c(1:length(myquantiles))){
		quantile_group[all_he>=myquantiles[i]]=i
	}
	quantile_inds=c(1:length(myquantiles))+1
	myquantile_names=paste("HE_quintile_",quantile_inds,sep="")
	for(i in c(1:length(quantile_group))){
		el=obs_list[[i]]
		mymapping=el$mapping_mat
		mydim1=dim(mymapping)[1]
		mydim2=length(myquantile_names)
		if(quantile_group[i]==0){
			he_quantile_mat=Matrix(0,nrow=mydim1,ncol=mydim2,sparse=TRUE)			
			}else{
			he_quantile_mat=sparseMatrix(i=c(1:mydim1),j=rep(quantile_group[i],mydim1), dims=c(mydim1,mydim2))
		}
		colnames(he_quantile_mat)=myquantile_names
		rownames(he_quantile_mat)=rownames(mymapping)
		new_mat=cbind(mymapping,he_quantile_mat)
		el$mapping_mat=new_mat==TRUE
		obs_list[[i]]=el
	}
	return(obs_list)
}

simulate_z=function(XtX_sqrt,genewise_param_list,my_n,trunc_var_scaled){
	n=my_n
	my_m=dim(XtX_sqrt)[1]
	my_mg=dim(XtX_sqrt)[2]
	lambda=genewise_param_list$lambda
	b=genewise_param_list$b
	tmp=t(XtX_sqrt)%*%b
	tmp=XtX_sqrt%*%tmp
	rhs=tmp/sqrt(my_n)
	lhs=XtX_sqrt%*%(rnorm(my_mg)/sqrt(lambda))
	#### trunc_var
	if(trunc_var_scaled!=0){
		rhs_upd=trunc_var_scaled*sqrt(my_n)*b
		rhs=rhs+rhs_upd
		lhs_upd=trunc_var_scaled*rnorm(my_m)/sqrt(lambda)
		## important: lhs_upd and lhs are using independent samples
		## therebey the covariance matrices add up.
		lhs=lhs+lhs_upd
	}
	z=rhs+lhs
	return(z)
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

load_expression_mat=function(EXPRESSION_FILE_PATH){
	###########################
	# load expression data
	###########################
	print(paste("loading expression data ",EXPRESSION_FILE_PATH,sep=""))
	print("...")
	if(!file.exists(EXPRESSION_FILE_PATH)){
		print("expression file missing:")
		print(EXPRESSION_FILE_PATH)
		stop("cant continue: abort")
	}

	is_gz=grepl(".gz$",EXPRESSION_FILE_PATH)[1]
	if(is_gz){
		expression_dat=fread(cmd=paste("gunzip -cf",EXPRESSION_FILE_PATH,sep=""),header=TRUE,sep="\t")
		}else{
		expression_dat=fread(EXPRESSION_FILE_PATH,header=TRUE,sep="\t")
	}
	expression_names_unique=length(names(expression_dat))==length(unique(names(expression_dat)))
	if(!expression_names_unique){
		stop("expression column names not unique")
		stop("abort")
	}
	if(names(expression_dat)[1]!="geneid"){
		print("first column in expression data mat not 'geneid'")
		stop("abort")
	}
	expression_geneids=expression_dat[,geneid]
	dev_nul=expression_dat[,geneid:=NULL]
	expression_mat=as.matrix(expression_dat)
	###########################
	# check gene expressions matrix
	###########################
	if(sum(is.na(expression_mat))!=0){
		stop("expression matrix was loaded with NAs. please check format")
	}
	if(!is.numeric(expression_mat)){
		stop("expression matrix couldn't be loaded as numeric matrix. please check format")
	}
	###########################
	# remove duplicate gene expressions.
	###########################
	duplicated_geneids=expression_geneids[duplicated(expression_geneids)]
	if(length(duplicated_geneids)>0){
		print(paste("found ",length(duplicated_geneids), "  duplicated gene ids in expression data",sep=""))
		print(paste("remove all duplicated genes..",sep=""))
		geneinds_tokeep=which(!is.element(expression_geneids,duplicated_geneids))
		expression_mat=expression_mat[geneinds_tokeep,]
		expression_geneids=expression_geneids[geneinds_tokeep]
		print(paste(length(expression_geneids), "  genes left in expression data",sep=""))
	}else{
		print("all genes unique")
	}
	rownames(expression_mat)=expression_geneids
	###########################
	# END: remove duplicate gene expressions.
	###########################
	###########################
	# END: expression data
	###########################
	return(expression_mat)
}


load_genotype_mat=function(chr,GENOTYPE_BASE_FILE_PATHS=NULL,GENOTYPE_FILE_PATHS_GZ=FALSE,snps2keep=NULL){
	genotype_file_path=paste(GENOTYPE_BASE_FILE_PATHS,chr,".txt",sep="")
	if(GENOTYPE_FILE_PATHS_GZ){
		genotype_file_path=paste(genotype_file_path,".gz",sep="")
	}
	if(!file.exists(genotype_file_path)){
		print("genotype file missing:")
		print(genotype_file_path)
		stop("abort")
	}
	if(GENOTYPE_FILE_PATHS_GZ){
		loadstr=paste("gunzip -cf",genotype_file_path)
		genotype_dat=fread(loadstr,header=TRUE,sep="\t")
		}else{
		genotype_dat=fread(genotype_file_path,header=TRUE,sep="\t")
	}
	genotype_names_unique=length(names(genotype_dat))==length(unique(names(genotype_dat)))
	if(!genotype_names_unique){
		stop("genotype column names not unique")
		stop("abort")
	}
	if(names(genotype_dat)[1]!="snpid"){
		print("first column in genotype data not 'snpid'")
		stop("abort")
	}
	genotype_snpids=genotype_dat[,snpid]
	dev_nul=genotype_dat[,snpid:=NULL]
	if(!is.null(snps2keep)){
		tokeep=is.element(genotype_snpids,snps2keep)
		genotype_snpids=genotype_snpids[tokeep]
		genotype_dat=genotype_dat[tokeep,]

	}
	genotype_mat=as.matrix(genotype_dat)
	if(sum(is.na(genotype_mat))!=0){
		print(paste("on chromsome ",chr,sep=""))
		stop("genotype matrix was loaded with NAs. please check format")
	}
	if(!is.numeric(genotype_mat)){
		print(paste("on chromsome ",chr,sep=""))
		stop("genotype matrix couldn't be loaded as numeric matrix. please check format")
	}
	###########################
	# remove duplicate snps from genotype.
	###########################
	duplicated_snpids=genotype_snpids[duplicated(genotype_snpids)]
	if(length(duplicated_snpids)>0){
		print(paste("found ",length(duplicated_snpids), "  duplicated snp ids in genotype data on chromosome ",chr,sep=""))
		print(paste("remove all duplicated snps..",sep=""))
		snpinds_tokeep=which(!is.element(genotype_snpids,duplicated_snpids))
		genotype_mat=genotype_mat[snpinds_tokeep,]
		genotype_snpids=genotype_snpids[snpinds_tokeep]
		print(paste(length(genotype_snpids), "  snps left in genotype data for chromosome ",chr,sep=""))
	}else{
		print(paste("all snps on chromosome ",chr," unique",sep=""))
	}
	rownames(genotype_mat)=genotype_snpids
	return(genotype_mat)
}

remove_lowmaf_snps=function(genotype_mat,MAF_CUTOFF,variance_transform=FALSE){
	###########################
	# remove low MAF snps.
	###########################
	if(variance_transform==FALSE){
		print(paste("removing  SNPs with observed maf equal or below ",MAF_CUTOFF,sep=""))
		aa=rowMeans(genotype_mat)
		snps_below_maf_cutoff=(aa<=MAF_CUTOFF*2 | (2-aa)<=MAF_CUTOFF*2)
		print(paste(sum(snps_below_maf_cutoff), " SNPs equal or below observed maf",sep=""))
		genotype_mat=genotype_mat[!(snps_below_maf_cutoff),]
	}else{
		print(paste("removing  SNPs with observed maf equal or below ",MAF_CUTOFF,sep=""))
		print(paste("assuming hardy weinberg",sep=""))
		var_cutoff=get_variance_from_maf(MAF_CUTOFF)
		aa=rowMeans(genotype_mat^2)
		if(max(rowMeans(genotype_mat))>1E-10){
			stop("demeaning didnt work")
		}
		snps_below_maf_cutoff=(aa<=var_cutoff)
		print(paste(sum(snps_below_maf_cutoff), " SNPs equal or below observed maf (variance)",sep=""))
		genotype_mat=genotype_mat[!(snps_below_maf_cutoff),]	
	}
	return(genotype_mat)
}

get_lowmaf_snps=function(genotype_mat,MAF_CUTOFF,variance_transform=FALSE){
	###########################
	# remove low MAF snps.
	###########################
	if(variance_transform==FALSE){
		print(paste("removing  SNPs with observed maf equal or below ",MAF_CUTOFF,sep=""))
		aa=rowMeans(genotype_mat)
		snps_below_maf_cutoff=(aa<=MAF_CUTOFF*2 | (2-aa)<=MAF_CUTOFF*2)
		print(paste(sum(snps_below_maf_cutoff), " SNPs equal or below observed maf",sep=""))
		removed_snps=rownames(genotype_mat)[snps_below_maf_cutoff]
	}else{
		print(paste("removing  SNPs with observed maf equal or below ",MAF_CUTOFF,sep=""))
		print(paste("assuming hardy weinberg",sep=""))
		var_cutoff=get_variance_from_maf(MAF_CUTOFF)
		aa=rowMeans(genotype_mat^2)
		if(max(rowMeans(genotype_mat))>1E-10){
			stop("demeaning didnt work")
		}
		snps_below_maf_cutoff=(aa<=var_cutoff)
		removed_snps=rownames(genotype_mat)[snps_below_maf_cutoff]
	}
	return(removed_snps)
}


get_variance_from_maf=function(MAF_CUTOFF){
	return(2*(1-MAF_CUTOFF)*MAF_CUTOFF)
}

get_snpfreq_annot=function(genotype_mat,maf_ranges_lower=c(0,0.05,0.1,0.2,0.35),maf_ranges_upper=c(0.05,0.1,0.2,0.35,0.5)){
	###########################
	# remove low MAF snps.
	###########################
	print("getting snp frequency ranges (assuming hardy weinberg)")
	aa=rowMeans(genotype_mat^2)-rowMeans(genotype_mat)^2
	lower_var=get_variance_from_maf(maf_ranges_lower)
	upper_var=get_variance_from_maf(maf_ranges_upper)
	my_snp_freqs_list=lapply(c(1:length(lower_var)),function(i){
			if(sum(aa<=upper_var[i] & aa>lower_var[i])>0){
				cbind(i,which(aa<=upper_var[i] & aa>lower_var[i]))
			}
		})
	combined=do.call("rbind",my_snp_freqs_list)
	freq_annot_mat=sparseMatrix(i=combined[,2],j=combined[,1], dims=c(length(aa),length(my_snp_freqs_list)))
	colnames(freq_annot_mat)=paste("snp_var_upto",maf_ranges_upper,sep="")
	rownames(freq_annot_mat)=rownames(genotype_mat)
	return(freq_annot_mat)

}


combine_genotype_mats=function(genotype_mats){
	###########################
	# remove low MAF snps.
	###########################
	print("center each genotype mat and combine.")
	genotype_mats=lapply(genotype_mats,function(genotype_mat){t(scale(t(genotype_mat),scale=FALSE))})	
	genotype_mat=do.call("cbind",genotype_mats)
	return(genotype_mat)
}

combine_genotype_mats_and_remove_lowmaf_snps=function(genotype_mats,MAF_CUTOFF){
	###########################
	# remove low MAF snps.
	###########################
	print("center each genotype mat and combine.")
	genotype_mats=lapply(genotype_mats,function(genotype_mat){t(scale(t(genotype_mat),scale=FALSE))})
	genotype_mat=do.call("cbind",genotype_mats)
	print(paste("removing  SNPs with observed maf equal or below ",MAF_CUTOFF,sep=""))
	print(paste("assuming hardy weinberg",sep=""))
	var_cutoff=get_variance_from_maf(MAF_CUTOFF)
	aa=rowMeans(genotype_mat^2)
	if(max(rowMeans(genotype_mat))>1E-10){
		stop("demeaning didnt work")
	}
	snps_below_maf_cutoff=(aa<=var_cutoff)
	print(paste(sum(snps_below_maf_cutoff), " SNPs equal or below observed maf (variance)",sep=""))
	genotype_mat=genotype_mat[!(snps_below_maf_cutoff),]
	return(genotype_mat)

}



intersect_expression_genes_fully=function(expression_mats){
	#############################
	# restrict to genes present in all expression matrices
	#############################
	if(length(expression_mats)>1){
			print("keep only genes that are present in all datasets ..")
			all_names_list=lapply(expression_mats,rownames)
			tabled_names=table(do.call("c",all_names_list))
			names_tokeep=names(tabled_names[tabled_names==length(expression_mats)])
			print(paste(length(names_tokeep), " genes left.",sep=""))
			if(length(names_tokeep)==0){
				stop("no gene present in all expression sets.")
			}
			expression_mats_out=lapply(c(1:length(expression_mats)),function(i){
				expression_mats[[i]][match(names_tokeep,all_names_list[[i]]),]
			})
	}else{
		expression_mats_out=expression_mats
	}
	return(expression_mats_out)
}


intersect_snps_fully=function(genotype_mats,snps2keep=NULL){
	#############################
	# restrict to genes present in all expression matrices
	#############################
	print("keep only SNPs that are present in all datasets ..")
	all_names_list=lapply(genotype_mats,rownames)
	if(!is.null(snps2keep)){
		print("overlap with external SNP list")
		all_names_list[[(length(all_names_list)+1)]]=snps2keep
	}
	tabled_names=table(do.call("c",all_names_list))
	presence_nr=(length(genotype_mats)+(!is.null(snps2keep)))
	names_tokeep=names(tabled_names[tabled_names==presence_nr])
	print(paste(length(names_tokeep), " SNPs left.",sep=""))
	if(length(names_tokeep)==0){
		stop("no SNPs present in all data sets and the external SNP list.")
	}
	genotype_mats_out=lapply(c(1:length(genotype_mats)),function(i){
		genotype_mats[[i]][match(names_tokeep,all_names_list[[i]]),]
	})
	return(genotype_mats_out)
}


prepare_data_per_chr_not1KG=function(chr,
	bagea_annotation_mat,	
	simulated_params,
	SNPS2KEEP,
	set,
	expression_mats
){
	###########################
	# restrict annotation data to single chr
	###########################
	bagea_annotation_mat=subset_annotation_mat_rows2chr(bagea_annotation_mat,chr)
	###########################
	# END: restrict annotation data to single chr
	###########################
	genotype_mats=load_all_genotype_mats(chr,bagea_annotation_mat,gwas_summary_stats_mode=FALSE,SNPS2KEEP,set=set)
	#####
	snp_tss_map_chr=bagea_annotation_mat[["snp_tss_map"]]
	annotation_mat_chr=bagea_annotation_mat[["annotation_mat"]]
	if(sum(names(bagea_annotation_mat)=="annotation_shift_mat")==1){
		shift_mat_available=TRUE
	}else{
		shift_mat_available=FALSE
	}
	if(shift_mat_available){
		annotation_shift_mat_chr=bagea_annotation_mat[["annotation_shift_mat"]]
	}
	###########################
	# restrict expression mat to single chr
	###########################
	out_mats=assemble_data_mats(expression_mats,genotype_mats,bagea_annotation_mat,do_simulation(simulated_params))
	expression_mat_scaled=out_mats$expression_mat_scaled#is empty if HYPERPARAM
	genotype_mat=out_mats$genotype_mat
	rm(out_mats)
	rm("genotype_mats")		
	###########################
	# END: order samples for genotypes.
	###########################
	my_variance_transform=TRUE
	all_snps_below_mafcutoff=get_lowmaf_snps(genotype_mat,set$MAF_CUTOFF,variance_transform=my_variance_transform)
	genotype_mat=remove_lowmaf_snps(genotype_mat,set$MAF_CUTOFF,variance_transform=my_variance_transform)
	if(set$ADD_SNP_FREQS2ANNOT){
		snpfreq_annot_mat=get_snpfreq_annot(genotype_mat,maf_ranges_lower=c(0,0.05),maf_ranges_upper=c(0.05,0.1))
	}
	###########################
	# END: remove low MAF snps.
	###########################
	###########################
	# scale genotype mat
	###########################
	print("scale genotype matrix to mean 0 sd 1 ..")
	genotype_mat_scaled=genotype_mat
	rm(genotype_mat)
	genotype_mat_scaled=t(genotype_mat_scaled)
	genotype_mat_scaled=my_scale(genotype_mat_scaled,use_minus1=FALSE)
	genotype_mat_scaled=t(genotype_mat_scaled)
	gc()
	###########################
	# END: scale genotype  mats
	###########################
	###########################
	# restrict annotation data
	###########################
	# restrict annotation data to pairs within user set RANGE
	## filter again on range as we need to filter on symbols on the same snps
	bagea_annotation_mat=subset_annotation_mat_rows2range(bagea_annotation_mat,set$RANGE_AROUND_TSS,gwas_summary_stats_mode=FALSE,symbols2filteron=rownames(expression_mat_scaled))
	bagea_annotation_mat=subset_bagea_annotation_mat_onsnpnames(bagea_annotation_mat,snpnames=rownames(genotype_mat_scaled))

	bagea_annotation_mat=add_genotypes2bagea_annotation_mat(bagea_annotation_mat,genotype_mat_scaled,set$ADD_SNP_FREQS2ANNOT)
	rm(genotype_mat_scaled)

	snp_tss_map_chr=bagea_annotation_mat[["snp_tss_map"]]
	annotation_mat_chr=bagea_annotation_mat[["annotation_mat"]]
	if(sum(names(bagea_annotation_mat)=="annotation_shift_mat")==1){
		shift_mat_available=TRUE
	}else{
		shift_mat_available=FALSE
	}
	if(shift_mat_available){
		annotation_shift_mat_chr=bagea_annotation_mat[["annotation_shift_mat"]]
	}
	genotype_mat_expanded=bagea_annotation_mat$genotype_mat_expanded
	###########################
	# END: restrict annotation data
	###########################
	### get gene list to process
	if(!do_simulation(simulated_params)){
		#genes_to_process=intersect(rownames(expression_mat_chr),snp_tss_map_chr[,symbol])
		genes_to_process=intersect(rownames(expression_mat_scaled),snp_tss_map_chr[,symbol])
	}else{
		#add global variables.
		genes_to_process=unique(snp_tss_map_chr[,symbol])
	}
	glen=length(genes_to_process)
	print("prepare observed gene list")
	print(paste(glen," to process",sep=""))
	all_Obs_list=lapply(c(1:glen),process_gene_not1KG,
		genes_to_process=genes_to_process,
		bagea_annotation_mat=bagea_annotation_mat,
		expression_mat_scaled=expression_mat_scaled,
		simulated_params=simulated_params,
		SINGLESNP_MINPVALUE_CUTOFF=set$SINGLESNP_MINPVALUE_CUTOFF,
		VARIANCE_PERCENTAGE2KEEP=set$VARIANCE_PERCENTAGE2KEEP,
		ANNOTATION_SHIFTWEIGHT_LIST=set$ANNOTATION_SHIFTWEIGHT_LIST,
		ADD_HERITABILTY_QUANTILE=set$ADD_HERITABILTY_QUANTILE,
		KEEP_X_y=set$KEEP_X_y
	)
	##remove NULL genes
	null_inds=sapply(all_Obs_list,function(x){is.null(x)})
	all_Obs_list=all_Obs_list[!null_inds]
	genes_to_process=genes_to_process[!null_inds]
	names(all_Obs_list)=genes_to_process
	attr(all_Obs_list, "all_snps_below_mafcutoff")=all_snps_below_mafcutoff
	return(all_Obs_list)	
}

process_gene_not1KG=function(i,
	genes_to_process,
	bagea_annotation_mat,
	expression_mat_scaled,
	simulated_params,
	SINGLESNP_MINPVALUE_CUTOFF,
	VARIANCE_PERCENTAGE2KEEP,
	ANNOTATION_SHIFTWEIGHT_LIST,
	ADD_HERITABILTY_QUANTILE,
	KEEP_X_y
	){
	geneid=genes_to_process[i]
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
	if(!do_simulation(simulated_params)){
		#y=expression_mat_scaled[rownames(expression_mat_chr)==geneid,,drop=F]
		y=expression_mat_scaled[rownames(expression_mat_scaled)==geneid,,drop=F]
		genewise_param_list=NULL
	}else{
		if(is.null(simulated_params[["p"]])){
			genewise_param_list=get_genewise_param_list(mapping_mat,mapping_shift_mat=NULL,hyperparameter_list=simulated_params)
		}else{
			genewise_param_list=get_genewise_param_list(mapping_mat,mapping_shift_mat,hyperparameter_list=simulated_params,annotation_nu_names=ANNOTATION_SHIFTWEIGHT_LIST)
		}
		y=simulate_y(X,genewise_param_list)
	}
	ytX=y%*%t(X)
	if(SINGLESNP_MINPVALUE_CUTOFF<1){
		pval_cutoff=pnorm(-max(abs(ytX))/sqrt(length(y[1,])))*2
		if(is.na(pval_cutoff)){stop("pval_cutoff is NA.")}
		if(pval_cutoff>SINGLESNP_MINPVALUE_CUTOFF){
			return(NULL)
		}
	}
	if(dim(X)[1]==1){
		print("currently does not support genes with only one snp. Gene is removed.")
		return(NULL)
	}
	yty=y%*%t(y)
	ns=dim(X)[2]
	eigens=rep(0,ns)
	if(ADD_HERITABILTY_QUANTILE){
		sigmasq_he=calc_he(t(y),t(X))
	}
	Xsvd=my_caught_svd(X)
	dl=length(Xsvd$d)
	eigens[c(1:dl)]=Xsvd$d^2
	eigens[eigens<0]=0
	cumeigens=cumsum(eigens)
	cumeigens_scaled=cumeigens/cumeigens[length(cumeigens)]
	selection_indices=c(1:min(which(cumeigens_scaled>=VARIANCE_PERCENTAGE2KEEP/100)))
	if(length(selection_indices)>1 && sum(selection_indices)>1){
		XtX_sqrt=Xsvd$u[,selection_indices,drop=F]%*%diag(Xsvd$d[selection_indices,drop=F])
	}else{
		XtX_sqrt=Xsvd$u[,selection_indices,drop=F]*Xsvd$d[selection_indices]
	}
	out=list()
	rownames(XtX_sqrt)=rownames(X)
	out[["chrom"]]=chr
	out[["XtX_sqrt"]]=XtX_sqrt
	if(KEEP_X_y){
		out[["X"]]=t(X)## conform to paper notation
		out[["y"]]=t(y)
	}
	out[["ytX"]]=ytX
	out["yty"]=yty
	out[["n"]]=ns
	out[["eig_values"]]=eigens
	out[["mapping_mat"]]=mapping_mat
	if(shift_mat_is_available(bagea_annotation_mat)){
		out[["mapping_shift_mat"]]=mapping_shift_mat
	}
	if(ADD_HERITABILTY_QUANTILE){
		out[["sigmasq_he"]]=sigmasq_he
	}
	out[["genewise_param_list"]]=genewise_param_list
	return(out)
}

#' prepare observed data for bagea run. 
#'
#' Takes in Genotype and expression data and combines them with a \code{bagea_annotation_mat} object to produce  an \code{bagea_observed_data} object that serves as input to \code{run_bagea}.
#' Genotype data has to be supplied indirectly by specifying the file paths of the chromosome-wise dosage files. The same goes for expression data. See Details below for a description of the procedure. 
#' @return
#' An R-object of class \code{bagea_observed_data}; an R-list contatining 4 elements. if the phenotype was simulated, the list contains 6 elements (see below).
#' \describe{\item{\code{obs_list}}{A named list of length n where n is the number of genes processed. Each element contains the data for a particular gene that is needed to run \emph{bagea}. If the phenotype was simulated, it also contains the unobserved random variables.}
#' \item{\code{settings}}{A list of parameter settings with which prepare_observed_data was run.(of bagea_annotation_mat object, only the settings-field is stored and of the SNPS2KEEP list only the length).}
#' \item{\code{simulation_unobs_list}}{Only produced if phenotype was simulated, The list contains all the unobserved random variables that are not gene specific. Additionally, it conatins all the hyper parameters.}
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
#' @param EXPRESSION_FILE_PATHS
#' Character vector or NULL. Element i gives the filepath to an expression matrix for the i-th batch in tab separated format. if not NULL, the file needs to contain a numeric tab-separated table (no NAs allowed), except the first column and row. The first row gives the column names, which are sample ids except for the first. The first column should be named 'geneid' and should contain gene hgnc symbols. All other column should contain a individuals' gene expression values. Should be NULL if the phenotype is simulated. In that case, a bagea_hyperparameter_list has to be supplied to the HYPER_PARAM_LIST argument. 
#' @param VARIANCE_PERCENTAGE2KEEP
#' Scalar between BAGEA 95 and 100. BAGEA allows to approximate the observed LD matrix with a lower rank approximation, by removing eigen vectors that do not contribute substantially to the LD matrix.  VARIANCE_PERCENTAGE2KEEP sets the amount of variance that the kept eigenvectors have to capture. has to be between 95 and 100.
#' @param RANGE_AROUND_TSS
#' Scalar between 10'000 and 500'000. 
#' Window around transcription start site for which snps are considered.
#' @param MAF_CUTOFF
#' Scalar. SNPs with observed maf below or at MAF_CUTOFF will be ignored. has to be between 0 and 0.5. The cutoff is transformed into a variance cutoff assuming hardy-weinberg equilibrium and applied to the combined dataset after mean removal of each batch.  
#' @param HYPER_PARAM_LIST
#' An object of class bagea_hyperparameter_list (produced by set_hyperparameter_list).
#' @param CHR
#' A vector of unique integer numbers between 1 and 22. Gives of chromosome to process. Default is c(1:22) 
#' @param KEEP_X_y
#' Boolean. Set TRUE if the original genotypes and phenotypes are to be kept.
#' @param NCORES 
#' Positive integer scalar. Number of cores to be used when loading the data.
#' @param SNPS2KEEP
#' Character vector or NULL. if not NULL all SNPS not in SNPS2KEEP will be removed.
#' @param ADD_SNP_FREQS2ANNOT
#' Boolean. if TRUE, SNP frequency will be added to the genomic annotations.
#' @param SINGLESNP_MINPVALUE_CUTOFF
#' A scalar set between 0 and 1. If smaller than 1, genes where no single SNP p-value reaches below this threshold is removed. 
#' @param GENOTYPE_PC_LIST
#' \code{bagea_genotype_pcs_data} object that gives genotype pc vectors used as controls (these are regressed out from the expression matrix prior to calculating the summary statistics). 
#' @param NPCS 
#' number of genotype pcs to be used in control.
#' @param SHIFT_METAINFO_TBL
#' A data frame. Defines groups shift annotations for modeling group-wise priors. The following has to hold: one column has to be named annotation_id and contain all elements in ANNOTATION_SHIFT_LIST uniquely. Every other column must contain names one particular grouping structure for directed genome annotations are used (for instance which cell type the particular annotation is mapped to).
#' @details Genotype data has to be supplied indirectly by specifying the file paths chromosome wise dosage files. The file paths are specified via 2 arguments.
#' GENOTYPE_BASE_FILE_PATHS and GENOTYPE_FILE_PATHS_GZ, where the i-th element of this vector specifies the paths for batch i.
#' Assuming that only one batch is analysed, then the paths are assembled as \preformatted{\{GENOTYPE_BASE_FILE_PATHS\}\{chr\}.txt,} where \{chr\} is a number between 1 to 22 and \code{GENOTYPE_FILE_PATHS_GZ[i]} specifies whether files in batch i are gzipped.
#' For instance, if GENOTYPE_BASE_FILE_PATHS is \preformatted{~/dosages/mygenotypes_chr,} GENOTYPE_FILE_PATHS_GZ is set to FALSE, then the function looks for 22 files of the format \preformatted{~/dosages/mygenotypes_chr1.txt} to \preformatted{~/dosages/mygenotypes_chr22.txt.} If GENOTYPE_FILE_PATHS_GZ is set to TRUE, the function will look for \preformatted{~/dosages/mygenotypes_chr22.txt.gz} etc instead.
#' All dosage files need to be numeric tab-separated tables (no NAs allowed), except the first column and row. The first row gives the column names, which are sample ids except for the first. The first column should be named 'snpid ' and should contain snp ids in rs numbers. All other column should contain a individuals' SNP dosages.
#' If multiple batches are combined in the analysis,then GENOTYPE_BASE_FILE_PATHS, GENOTYPE_FILE_PATHS_GZ,EXPRESSION_FILE_PATHS should be given as vectors, where element i in the vector is the respective parameter setting for the i-th batch and formatting follows the same rules as for one batch. 
#' To simmulate according to a model with a particular hyper-parameter setting, set HYPER_PARAM_LIST to an object of class bagea_hyperparameter_list (produced by set_hyperparameter_list), and set EXPRESSION_FILE_PATHS to NULL.
#' @export
prepare_observed_data=function(BAGEA_ANNOTATION_MAT=NULL,ANNOTATION_LIST=colnames(BAGEA_ANNOTATION_MAT$annotation_mat),ANNOTATION_SHIFT_LIST=colnames(BAGEA_ANNOTATION_MAT$annotation_shift_mat),GENOTYPE_BASE_FILE_PATHS=NULL,GENOTYPE_FILE_PATHS_GZ=rep(FALSE,length(GENOTYPE_BASE_FILE_PATHS)),EXPRESSION_FILE_PATHS=NULL,VARIANCE_PERCENTAGE2KEEP=100,RANGE_AROUND_TSS=1.5e5,MAF_CUTOFF=0.05,HYPER_PARAM_LIST=NULL,CHR=c(1:22),KEEP_X_y=FALSE,NCORES=1,SNPS2KEEP=NULL,ADD_SNP_FREQS2ANNOT=FALSE,SINGLESNP_MINPVALUE_CUTOFF=1,ANNOTATION_SHIFTWEIGHT_LIST=colnames(BAGEA_ANNOTATION_MAT$annotation_mat),GENOTYPE_PC_LIST=NULL,NPCS=10,SHIFT_METAINFO_TBL){
	###
	ADD_HERITABILTY_QUANTILE=FALSE
	settings=get_settingsout_list_obs(
		BAGEA_ANNOTATION_MAT=BAGEA_ANNOTATION_MAT,
		ANNOTATION_LIST=ANNOTATION_LIST,
		ANNOTATION_SHIFT_LIST=ANNOTATION_SHIFT_LIST,
		GENOTYPE_BASE_FILE_PATHS=GENOTYPE_BASE_FILE_PATHS,
		GENOTYPE_FILE_PATHS_GZ=GENOTYPE_FILE_PATHS_GZ,
		EXPRESSION_FILE_PATHS=EXPRESSION_FILE_PATHS,
		VARIANCE_PERCENTAGE2KEEP=VARIANCE_PERCENTAGE2KEEP,
		RANGE_AROUND_TSS=RANGE_AROUND_TSS,
		MAF_CUTOFF=MAF_CUTOFF,
		HYPER_PARAM_LIST=HYPER_PARAM_LIST,
		CHR=CHR,
		KEEP_X_y=KEEP_X_y,
		NCORES=NCORES,
		SNPS2KEEP=SNPS2KEEP,
		ADD_SNP_FREQS2ANNOT=ADD_SNP_FREQS2ANNOT,
		ADD_HERITABILTY_QUANTILE=ADD_HERITABILTY_QUANTILE,
		SINGLESNP_MINPVALUE_CUTOFF=SINGLESNP_MINPVALUE_CUTOFF,
		ANNOTATION_SHIFTWEIGHT_LIST=ANNOTATION_SHIFTWEIGHT_LIST,
		SHIFT_METAINFO_TBL=SHIFT_METAINFO_TBL,
		GENOTYPE_PC_LIST=GENOTYPE_PC_LIST,
		NPCS=NPCS
		)
	check_settings_list_obs(BAGEA_ANNOTATION_MAT,SNPS2KEEP,
	set=settings)
	if(!is.null(settings$HYPER_PARAM_LIST)){
		do_simulation=TRUE
	}else{
		do_simulation=FALSE
	}
	###########################
	# prepare annotation data
	###########################
	bagea_annotation_mat=subset_annotation_mat_cols(BAGEA_ANNOTATION_MAT,ANNOTATION_LIST,ANNOTATION_SHIFT_LIST,ANNOTATION_SHIFTWEIGHT_LIST)
	rm(BAGEA_ANNOTATION_MAT)
	###########################
	# prepare D_mat from SHIFT_METAINFO_TBL
	###########################	
	D_mat_outl=setup_D_mat_with_checks(SHIFT_METAINFO_TBL,bagea_annotation_mat)
	D_mat=D_mat_outl[["D_mat"]]
	D_names_list=D_mat_outl[["D_names_list"]]
	rm(D_mat_outl)
	###########################
	# END prepare D_mat from SHIFT_METAINFO_TBL
	###########################
	if(do_simulation(HYPER_PARAM_LIST)){
		simulated_params=simulate_params_from_hyperparams(HYPER_PARAM_LIST,D_mat,ANNOTATION_LIST,ANNOTATION_SHIFT_LIST,ANNOTATION_SHIFTWEIGHT_LIST)
		if(!is.null(HYPER_PARAM_LIST$seed_val)){
			set.seed(HYPER_PARAM_LIST$seed_val)
		}
	}else{
		simulated_params=NULL
	}
	if(!do_simulation(HYPER_PARAM_LIST)){
		expression_mats=get_preprocess_expression_mats(set=settings)
	}
	if(NCORES!=1){
		tt=mclapply(settings$CHR,prepare_data_per_chr_not1KG,
			bagea_annotation_mat=bagea_annotation_mat,	
			simulated_params=simulated_params,
			SNPS2KEEP=SNPS2KEEP,
			set=settings,
			expression_mats=expression_mats,
			mc.cores=NCORES,mc.preschedule=F)
	}else{
		tt=lapply(settings$CHR,prepare_data_per_chr_not1KG,
			bagea_annotation_mat=bagea_annotation_mat,	
			simulated_params=simulated_params,
			SNPS2KEEP=SNPS2KEEP,
			set=settings,
			expression_mats=expression_mats
			)
	}
	all_removed_snps_list=lapply(c(1:length(tt)),function(i){attr(tt[[i]], "all_snps_below_mafcutoff")})
	obs_list=unlist(tt,recursive=FALSE)
	if(ADD_HERITABILTY_QUANTILE){
		obs_list=add_heritability_quantile_mat(obs_list=obs_list)
	}
	obs_out=list()
	obs_out[["obs_list"]]=obs_list
	obs_out[["settings"]]=settings
	obs_out[["simulation_unobs_list"]]=simulated_params
	obs_out[["removed_snps"]]=all_removed_snps_list
	if(!is.null(SHIFT_METAINFO_TBL)){
		obs_out[["D_mat"]]=D_mat
		obs_out[["D_names_list"]]=D_names_list
	}
	class(obs_out)=c("bagea_observed_data",class(obs_out))
	return(obs_out)
}


get_globalE_list=function(hyperparameter_list=NULL,a_names=c("1","2"),omega_names=c("1","2"),nu_names=c("1","2"),D_mat=NULL){
	check_hyperparameter_list(hyperparameter_list)
	my_t=length(a_names)
	my_s=length(omega_names)
	my_q=length(nu_names)+1
	phi1=hyperparameter_list$phi1
	phi2=hyperparameter_list$phi2
	tau1=hyperparameter_list$tau1
	rho1=hyperparameter_list$rho1
	rho2=hyperparameter_list$rho2
	xi1=hyperparameter_list$xi1
	xi2=hyperparameter_list$xi2
	chi1=hyperparameter_list$chi1
	chi2=hyperparameter_list$chi2
	p=hyperparameter_list$p
	lambda1=hyperparameter_list$lambda1
	gamma1=hyperparameter_list$gamma1
	theta_lambda2=calc_theta_lambda2(rho1,rho2)
	theta_ak=calc_theta_ak(phi1,phi2,my_t)
	theta_tau2=calc_theta_tau2(xi1,xi2)
	Eu_tau2=calc_Eu_tau2_via_theta(theta_tau2)
	Eu_lambda2=calc_Eu_lambda2_via_theta(theta_lambda2)
	if(!is.null(hyperparameter_list$p)){
		theta_delta=calc_theta_delta(chi1,chi2)
		Eu_delta=calc_Eu_delta_via_theta(theta_delta)
		Eu_omega_list=calc_Eu_omega_via_theta(calc_theta_omegas_via_Eu(Eu_delta,my_s))
		rownames(Eu_omega_list[[1]])=omega_names
		Eu_nu_list=calc_Eu_nu_via_theta(calc_theta_nu(p,my_q))
		rownames(Eu_nu_list[[1]])=c("background",nu_names)
	}
	Eu_ak=calc_Eu_ak_via_theta(calc_theta_ak(phi1,phi2,my_t))
	rownames(Eu_ak)=a_names
	globalE_list=list()
	globalE_list$Eu_lambda2=Eu_lambda2
	globalE_list$Eu_tau2=Eu_tau2
	if(!is.null(hyperparameter_list$p)){
		if(is.null(D_mat)){
			globalE_list$Eu_delta=Eu_delta
		}else{
			my_q=dim(D_mat)[2]
			if(!is.null(hyperparameter_list$zeta1)){
				zeta1_vec=hyperparameter_list$zeta1
				zeta2_vec=hyperparameter_list$zeta2
				if(length(zeta1_vec)==1){
					zeta1_vec=rep(zeta1_vec,my_q)
					zeta2_vec=rep(zeta2_vec,my_q)
				} else if(length(zeta1_vec)!=my_q){
				stop("erromsc:chi1 has to be either length 1 or myq.")
				}
				theta_chi2j=calc_theta_chi2j(zeta1_vec,zeta2_vec)
				Eu_chi2j=calc_Eu_chi2j_via_theta(theta_chi2j)
				globalE_list$Eu_chi2j=Eu_chi2j
				chi2_vec=Eu_chi2j[,2]
				chi1_vec=chi1
			}else{
				chi2_vec=chi2
				chi1_vec=chi1
			}
			if(length(chi1_vec)==1){
				chi1_vec=rep(chi1_vec,my_q)
			}else if(length(chi1_vec)!=my_q){
				stop("erromsc:chi1 has to be either length 1 or myq.")
			}
			if(length(chi2_vec)==1){
				chi2_vec=rep(chi2_vec,my_q)
			}else if(length(chi2_vec)!=my_q){
				stop("erromsc:chi2 has to be either length 1 or myq.")
			}
			D_mat_check=apply(D_mat,2,is_ascending_integer_vec)
			h_vec=apply(D_mat,2,max)
			globalE_list$Eu_upsilon_list=calc_Eu_upsilon_via_theta_list(calc_theta_upsilon(chi1_vec,chi2_vec,h_vec=h_vec))
		}
		globalE_list$Eu_omega_list=Eu_omega_list
		globalE_list$Eu_nu_list=Eu_nu_list
	}
	globalE_list$Eu_ak=Eu_ak
	return(globalE_list)
}

get_E_list=function(hyperparameter_list=NULL,globalE_list=NULL,n=NULL,m=NULL){
	check_hyperparameter_list(hyperparameter_list)
	###
	Eu_tau2=globalE_list$Eu_tau2
	Eu_lambda2=globalE_list$Eu_lambda2
	###
	phi1=hyperparameter_list$phi1
	phi2=hyperparameter_list$phi2
	tau1=hyperparameter_list$tau1
	rho1=hyperparameter_list$rho1
	rho2=hyperparameter_list$rho2
	xi1=hyperparameter_list$xi1
	xi2=hyperparameter_list$xi2
	lambda1=hyperparameter_list$lambda1
	gamma1=hyperparameter_list$gamma1
	###
	theta_lambda2=calc_theta_lambda2(rho1,rho2)
	theta_kappa=calc_theta_kappa(tau1,Eu_tau2[2])
	theta_lambda=calc_theta_lambda(lambda1,Eu_lambda2[2])
	###
	Eu_lambda=calc_Eu_lambda_via_theta(theta_lambda)
	Eu_kappa=calc_Eu_kappa_via_theta(theta_kappa)
	Eu_alpha1=calc_Eu_gamma(gamma1,Eu_kappa[2])
	Eu_alphai=t(kronecker(t(rep(1,m)), Eu_alpha1))
	b_cur=matrix(0,m,1)
	if(exists("B_INIT_OLD_FORM") && B_INIT_OLD_FORM==TRUE){
		Eu_b_list=list()
	    Eu_b_list[["Eu_b"]]=b_cur
	    Eu_b_list[["Eu_b2"]]=b_cur%*%t(b_cur)+diag(length(b_cur))
		if(DEBUG_FORM==TRUE){
			b_cur=genelist$Unobs_list$b
			Eu_b_list[["Eu_b"]]=b_cur
			Eu_b_list[["Eu_b2"]]=b_cur%*%t(b_cur)+diag(length(b_cur))
		}
	}else{
		theta_M_inv=list()
		len=length(b_cur)
		diagDprime_alphainv=rep(2,len)
		theta_M_inv[["inversion_formula"]]="diag(diagDprime_alphainv)-XtX_sqrtAug_interimRes_sqrt%*%t(XtX_sqrtAug_interimRes_sqrt)"
		theta_M_inv[["diagDprime_alphainv"]]=diagDprime_alphainv
		Eu_b_list=list()
		Eu_b_list[["Eu_b2_formula"]]="0.5*process(theta_M_inv) + Eu_b_1 %*% t(Eu_b_1))"
		Eu_b_list[["Eu_b2_formula_expanded"]]="0.5*(diag(Eu_b_list$theta_M_inv$diagDprime_alphainv)-Eu_b_list$theta_M_inv$XtX_sqrtAug_interimRes_sqrt%*%t(Eu_b_list$theta_M_inv$XtX_sqrtAug_interimRes_sqrt)) + Eu_b_list$Eu_b_1 %*% t(Eu_b_list$Eu_b_1)"
		Eu_b_list[["Eu_b_1"]]=b_cur
		Eu_b_list[["theta_M_inv"]]=theta_M_inv
	}
	E_list=list()
	E_list$Eu_b_list=Eu_b_list
	E_list$Eu_alphai=Eu_alphai
	E_list$Eu_lambda=Eu_lambda
	E_list$Eu_kappa=Eu_kappa
	if(sum(class(hyperparameter_list)=="bagea_hyperparameter_bool_list")!=0){
		#stop("errsd,cs=mc hyperparameter_list has to be a bagea_hyperparameter_bool_list")
		if(!is.null(hyperparameter_list$p) && length(hyperparameter_list$p)==1 && hyperparameter_list$p>0){
			Eu_pi=rep(hyperparameter_list$p,m)
			E_list$Eu_pi=Eu_pi
		}
	}
	return(E_list)
}

check_paramgreaterthan0=function(param){
	if(!is.scalar(param) || !is.numeric(param)){
	stop("param should be single scalar.",sep="")
	}
	if(param<=0){
		stop("param should be larger than 0.",sep="")
	}
}



initialize_parameters=function(observed_data,hyperparameter_list,nu_names=c(""),D_mat=NULL,ak_names=NULL){
	if(sum(class(observed_data)=="bagea_observed_data")==0){
		stop("observed_data argument is not a bagea_observed_data object.")
	}
	check_hyperparameter_list(hyperparameter_list)
	obs_list=observed_data$obs_list
	my_names=names(obs_list)
	if(is.null(ak_names)){
		my_a_names=colnames(obs_list[[1]]$mapping_mat)
	}else{
		my_a_names=ak_names
	}
	mydims=dim(obs_list[[1]]$mapping_mat)
	mydims_shift=dim(obs_list[[1]]$mapping_shift_mat)
	my_t=length(my_a_names)
	my_m=mydims_shift[1]
	my_s=mydims_shift[2]
	my_omega_names=colnames(obs_list[[1]]$mapping_shift_mat)
	D_mat=observed_data$D_mat
	if(!is.null(D_mat)){
		q_names=colnames(D_mat)
		if(is.null(q_names)){
			stop("erromsdfsdf:D_mat has unnamed columns.")
		}
	}else{
		q_names=NULL
	}
	if(is.null(my_a_names)){
		stop("erromsdfsdf:mapping_mat has unnamed columns.")
	}
	if(is.null(my_omega_names)){
		stop("erra[csd,:mapping_shift_mat has unnamed columns.")
	}
	globalE_list=get_globalE_list(hyperparameter_list,a_names=my_a_names,omega_names=my_omega_names,nu_names=nu_names,D_mat=D_mat)
	gene_dat_list=lapply(c(1:length(my_names)),function(i){
		cur_obs=obs_list[[i]]
		my_n=cur_obs$n
		mydims=dim(cur_obs$mapping_mat)
		my_m=mydims[1]
		e_list=get_E_list(hyperparameter_list,globalE_list,n=my_n,m=my_m)
		out_list=list()
		out_list$Obs_list=cur_obs
		out_list$E_list=e_list
		return(out_list)
	})
	names(gene_dat_list)=my_names
	attr(gene_dat_list,"globalE_list")=globalE_list
	class(gene_dat_list)=c("bagea_initialized_gene_dat_list",class(gene_dat_list))
	return(gene_dat_list)
}
