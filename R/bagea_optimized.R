#' @import methods
#' @import data.table
#' @import parallel
#' @import pryr
#' @import Matrix
#' @importClassesFrom Matrix Matrix
#' @useDynLib bagea, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

run_gene_cycle_wrapper_lambda_b_nu=function(gene_list=NULL,Eu_omega_list=NULL,Eu_nu_list=NULL,nu_matcher=NULL,theta_lambda=NULL,Eu_lambda2=NULL,calc_L=NULL,approximation_mode=NULL,p_group_vec_gene=NULL){
	Obs_list=gene_list$Obs_list
	E_list=gene_list$E_list
	out_list=run_gene_cycle_lambda_b_nu(Eu_omega_list=Eu_omega_list,Eu_nu_list=Eu_nu_list,nu_matcher=nu_matcher,mapping_mat=Obs_list$mapping_mat,yty=Obs_list$yty,XtX_sqrt=Obs_list$XtX_sqrt,ytX=Obs_list$ytX,N=Obs_list$n,theta_lambda=theta_lambda,mapping_shift_mat=Obs_list$mapping_shift_mat,Eu_b_list=E_list$Eu_b_list,Eu_alphai=E_list$Eu_alphai,Eu_lambda=E_list$Eu_lambda,calc_L=calc_L,eig_values=Obs_list$eig_values,approximation_mode=approximation_mode)
	return(out_list)
}

process_Ls_on_node_wshift_list=function(){
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	L_Qb_part_tot=0
	Ls_list=list()
	for(i in c(1:length(gene_dat_list))){
		if(length(gene_dat_list[[i]]$E_list$Eu_b_list)==2){
			stop("errsdASfmsd:not implemented.")
		}else{
			L_b_star_part=gene_dat_list[[i]]$L_list$L_b_star_part[[1]]

		}
		L_Qb_part_tot=L_Qb_part_tot+L_b_star_part
	}
	L_obs=sum(unlist(lapply(c(1:length(gene_dat_list)),function(i){sum(unlist(gene_dat_list[[i]]$L_list$L_obs))})))
	L_lambda_star_part=sum(unlist(lapply(c(1:length(gene_dat_list)),function(i){sum(unlist(gene_dat_list[[i]]$L_list$L_lambda_star_part))})))
	L_pi_star_part=sum(unlist(lapply(c(1:length(gene_dat_list)),function(i){sum(unlist(gene_dat_list[[i]]$L_list$L_pi_star_part))})))
	Ls_dt=data.table(L_Qb_part_tot=L_Qb_part_tot,L_obs=L_obs,L_Qlambda_part=L_lambda_star_part,L_Qpi_part=L_pi_star_part)
	return(Ls_dt)
}





present_in_all=function(all_names_list,make_unique_first=TRUE){
	len=length(all_names_list)
	if(make_unique_first){
		tabled_names=table(do.call("c",lapply(all_names_list,function(x){unique(x)})))
	}else{
		tabled_names=table(do.call("c",all_names_list))
	}
	names_present_in_all=names(tabled_names[tabled_names==len])
	return(names_present_in_all)
}

get_unique_matcher_mat=function(all_names_list){
	names_present_in_all=present_in_all(all_names_list)
	get_all_matchers=lapply(all_names_list,function(mynames){
		match(names_present_in_all,mynames)
		})
	my_mat=do.call("rbind",get_all_matchers)
	rownames(my_mat)=names(all_names_list)
	colnames(my_mat)=names(names_present_in_all)
}

my_docall_rbind=function(list2bind,nsplits=1){
	if(nsplits==1){
		return(do.call("rbind",list2bind))
	}else{
		mylen=length(list2bind)
		nsplits_dist=max(floor(mylen/nsplits),1)
		myseq=seq(from=0,to=length(list2bind),by=nsplits_dist)
		if(myseq[length(myseq)]<mylen){
			myseq=c(myseq,mylen)
		}
		nested_list=list()
		for(i in c(1:(length(myseq)-1))){
			sublist=do.call("rbind",list2bind[c((myseq[i]+1):myseq[i+1])])
			nested_list[[i]]=sublist
		}

		out_list=do.call("rbind",nested_list)
		return(out_list)
	}
}

is.scalar=function(x){is.atomic(x) && length(x) == 1L}
is.pos=function(x){
	aa=is.numeric(x) && sum(x<0) == 0
	return(aa)
}

check_string=function(x,varname){
	check_scalar(x,varname)
	if(!is.character(x)){
		stop(paste(varname," in not a character scalar",sep=""))
	}
}

check_matrix=function(x,varname){
	if(!is.numeric(x) || !is.matrix(x)){
		stop(paste(varname," in not a numeric matrix",sep=""))
	}
}
check_vector=function(x,varname){
	if(!is.numeric(x) || !is.vector(x)){
		stop(paste(varname," in not a numeric vector",sep=""))
	}
}
check_logicalmatrix=function(x,varname){
	if(extends(class(x),"lsparseMatrix")){
		return()
	}
	if(class(x)=="ngCMatrix"){
		if(!is.logical(matrix(x))){
			stop(paste(varname," in not a logical matrix",sep=""))
		}
	}else{
		if(!is.logical(x) || !is.matrix(x)){
		stop(paste(varname," in not a logical matrix",sep=""))
		}
	}
}
check_path=function(x,varname){
	a=file.exists(x)
	if(!a){
		stop(paste(varname,"file  does not exist",sep=""))
	}
}
check_posinteger=function(x,varname){
	check_scalar(x,varname)
	check_numeric(x,varname)
	check_pos(x,varname)
	if((x %% 1)!=0){
		stop(paste(varname,"argument not integer",sep=""))
	}
}

check_class=function(classname,x,varname){
	if(sum(class(x)==classname)==0){
		stop(paste(varname,"argument is not a bagea_output object",sep=""))
	}
}

check_bagea_output=function(x,varname){
	if(sum(class(x)=="bagea_output")==0){
		stop(paste(varname,"argument is not a bagea_output object",sep=""))
	}
}
check_bagea_observed_data=function(x,varname){
	if(sum(class(x)=="bagea_observed_data")==0){
		stop(paste(varname,"argument is not a bagea_observed_data object",sep=""))
	}
}

check_scalar=function(x,varname){
	if(!is.scalar(x)){
		stop(paste(varname,"argument not scalar",sep=""))
	}
}
check_numeric=function(x,varname){
	if(!is.numeric(x)){
		stop(paste(varname,"argument not numeric",sep=""))
	}
}
check_pos=function(x,varname){
	if(!is.numeric(x) || !is.pos(x)){
		stop(paste(varname," argument not pos",sep=""))
	}
}
check_notnan=function(x,varname){
	if(sum(is.na(x))>0){
		stop(paste(varname," argument contains NA",sep=""))
	}
}
check_boolean=function(x,varname){
	if(!is.scalar(x) || !is.logical(x)){
		stop(paste(varname," argument non logical",sep=""))
	}
}

load_as=function(filepath){
	tt=load(filepath)
	eval(parse(text=paste("out=",tt[1],sep="")))
	return(out)
}


check_theta_gamma=function(theta_gamma,varname){
	check_numeric(theta_gamma,varname)
	check_notnan(theta_gamma,varname)
	check_pos(-theta_gamma[2],varname)
}


calc_L_b_star_part_2=function(theta_b_star_list,Eu_b_list){
	D1=theta_b_star_list$sec_diag
	U=theta_b_star_list$sec_sqrt_abs
	mm=dim(U)[2]
	D1_expanded=kronecker(rep(1,mm),t(D1))
	D2=Eu_b_list$theta_M_inv$diagDprime_alphainv
	D2_expanded=Eu_b_list$theta_M_inv$diagDprime_alphainv_expanded
	V=Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt
	gamma=Eu_b_list$Eu_b_1
	A=t(V)%*%U
	res=0.5*sum(D1*D2)+sum(gamma^2*D1)-0.5*sum(V^2*t(D1_expanded))-sum((t(gamma)%*%U)^2)-0.5*sum(D2_expanded*U^2)+0.5*sum(A^2)
	res=res*(-1)
	return(res)
}

calc_L_b_star_part_2_only_diag=function(theta_b_diag,Eu_b_list){
	D1=theta_b_diag$sec_diag
	D2=Eu_b_list$theta_M_inv$diagDprime_alphainv
	D2_expanded=Eu_b_list$theta_M_inv$diagDprime_alphainv_expanded
	V=Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt
	mm=dim(V)[2]
	D1_expanded=kronecker(rep(1,mm),t(D1))
	gamma=Eu_b_list$Eu_b_1
	res=0.5*sum(D1*D2)+sum(gamma^2*D1)-0.5*sum(V^2*t(D1_expanded))
	res=res
	return(res)
}

#tested
calculate_woodbury_inversion=function(diagDprime_alpha,eigU,eigValues){
	if(!is.null(dim(eigValues)) ||  !is.null(dim(diagDprime_alpha))){
		stop("matrix where vector demanded.")
	}
	n=dim(eigU)[1]
	k=dim(eigU)[2]

	if(length(diagDprime_alpha)!=n || length(eigValues)!=k){
		stop("dimensions do not agree.")
	}
	diagDprime_alphainv=1/diagDprime_alpha
	diagDprime_alphainv_expanded=kronecker(diagDprime_alphainv,t(rep(1,k)))
	eigUt=eigU*kronecker(rep(1,n),t(eigValues))
	eigUtAug=eigUt*diagDprime_alphainv_expanded
	eigUAug=eigU*diagDprime_alphainv_expanded
	interimRes=solve((diag(k)+t(eigUtAug)%*%eigU))
	finalRes=-eigUAug%*%interimRes%*%t(eigUtAug)
	diag(finalRes)=diag(finalRes)+diagDprime_alphainv
	return(finalRes)
}

#tested
calculate_woodbury_inversion_sym=function(diagDprime_alpha,cXtX_sqrt,only_trace=FALSE,decomposed=FALSE){
	# 
	# Note: see section 'approximating -\phi^*_b(2) for fast inversion' for context.
	# args diagDprime_alpha: follows notation of aforementioned section.
	# args cXtX_sqrt: matrix square root of cXtX
 	# (note -\phi^*_b(2) = cXtX + diagDprime_alpha)
	# value
	# finalRes: list of elements to compute processed to -\phi^*_b(2) fast, while at the same time having smaller memory footprint. (These components can be used in later steps more efficiently than the actual inversion result)
	# finalRes[["inversion_formula"]]: formula to calculate the actual inversion from the elements
	# finalRes[["cXtX_sqrtAug_interimRes_sqrt"]]: square root of \boldsymbol{{D'}_{\alpha}}^{-1}\boldsymbol{A_t}(\boldsymbol{I_t}+\boldsymbol{A^T_t}\boldsymbol{{D'}_{\alpha}}^{-1}\boldsymbol{A_t})^{-1}\boldsymbol{A^T_t}\boldsymbol{{D'}_{\alpha}}^{-1}
	# finalRes[["diagDprime_alphainv"]]: \boldsymbol{{D'}_{\alpha}}^{-1}
	# finalRes[["diagDprime_alphainv_expanded"]]: \boldsymbol{{D'}_{\alpha}}^{-1} repeated to same dimensionality as cXtX_sqrtAug_interimRes_sqrt
	n=dim(cXtX_sqrt)[1]
	k=dim(cXtX_sqrt)[2]

	diagDprime_alphainv=1/diagDprime_alpha
	diagDprime_alphainv_expanded=kronecker(diagDprime_alphainv,t(rep(1,k)))
	cXtX_sqrtAug=cXtX_sqrt*diagDprime_alphainv_expanded ##equivalent to D_alphainv%*%A_t
	interimRes=calc_fast_cov_solve(diag(k)+t(cXtX_sqrtAug)%*%cXtX_sqrt)
	if(decomposed==TRUE){
		finalRes=list()
		finalRes[["inversion_formula"]]="diag(diagDprime_alphainv)-cXtX_sqrtAug_interimRes_sqrt%*%t(cXtX_sqrtAug_interimRes_sqrt)"
		finalRes[["diagDprime_alphainv"]]=diagDprime_alphainv
		finalRes[["diagDprime_alphainv_expanded"]]=diagDprime_alphainv_expanded
		mychol=t(chol(interimRes))
		finalRes[["cXtX_sqrtAug_interimRes_sqrt"]]=cXtX_sqrtAug%*%mychol
		return(finalRes)
	}
	if(only_trace==FALSE){
		finalRes=-cXtX_sqrtAug%*%interimRes%*%t(cXtX_sqrtAug)
		diag(finalRes)=diag(finalRes)+diagDprime_alphainv
		return(finalRes)
	}else{
		RR=t(cXtX_sqrtAug)%*%cXtX_sqrtAug
		finalRes=sum(diagDprime_alphainv)-sum(RR*interimRes)
		return(finalRes)
	}
}

######## to test
get_diag_from_inv_decompose=function(my_inv){
	res=my_inv[["diagDprime_alphainv"]]-rowSums(my_inv[["cXtX_sqrtAug_interimRes_sqrt"]]^2)
	return(res)
}

calculate_determinant_sym=function(diagDprime_alpha,XtX_sqrt,mysign=1,proport=1,calc_log=TRUE){
	diagDprime_alphainv=1/diagDprime_alpha
	k=dim(XtX_sqrt)[2]
	diagDprime_alphainv_expanded=kronecker(diagDprime_alphainv,t(rep(1,k)))
	if(calc_log==FALSE){
		stop("not implemented")
	}
	mydet=determinant(diag(k)+mysign*t(XtX_sqrt*diagDprime_alphainv_expanded)%*%XtX_sqrt,logarithm=TRUE)
	mylog_det=sum(log(proport*diagDprime_alpha)) + mydet$modulus
	if(mydet$sign!=1){
			warning("degenerate")
			return(-Inf)
	}
	return(mylog_det)
}

split_along_length_lists=function(all_length_splits,mat){
	out_list=list()
	starter=1
	stopper=0
	for(j in c(1:length(all_length_splits))){
		inner_list=list()
		for(i in c(1:length(all_length_splits[[j]]))){
			stopper=all_length_splits[[j]][[i]]+stopper
			inner_list[[i]]=mat[starter:stopper,]
			starter=stopper+1
		}
		out_list[[j]]=inner_list
	}
	return(out_list)
}

split_along_lengths=function(all_length_split,vect){
	numlens=sum(unlist(all_length_split))
	if(length(vect)!=sum(numlens)){
		stop("not same length.")
	}
	out_list=list()
	starter=1
	stopper=0
	for(i in c(1:length(all_length_split))){
		print(i)
		stopper=all_length_split[i]+stopper
		out_list[[i]]=vect[starter:stopper]
		starter=stopper+1
	}
	return(out_list)
}

clusterExport_or_envExport=function(mpicl,varName_vect,envir=.GlobalEnv){
	if(class(mpicl)[1]=="SOCKcluster"){
		out=clusterExport(mpicl,varName_vect,envir=envir)
	}else{
		if(class(mpicl)[1]=="environment"){
			for(varname in varName_vect){
				local_var=get(varname,envir=envir)
				assign(varname,local_var,envir=mpicl)
			}
		}else{
			stop("neither cluster nor environment")
		}
	}
}

clusterApply_or_envApply2=function(mpicl,x,fun,...){
	if(class(mpicl)[1]=="SOCKcluster"){
		out=clusterApply(mpicl,x,fun,...)
	}else{
		if(class(mpicl)[1]=="environment"){
			arg_list=list()
			arg_list[[1]]=x[[1]]
			params=list(...)
			arg_list=append(x,params)
			fun_arg=deparse(substitute(fun))
#			out0=do.call(fun,arg_list,quote=FALSE,envir=mpicl)
			out0=do.call(fun_arg,arg_list,quote=FALSE,envir=mpicl)
			out=list()
			out[[1]]=out0
		}else{
			stop("neither cluster nor environment")
		}
	}
	return(out)
}

clusterApply_or_envApply=function(mpicl,x,fun,...){
	if(class(mpicl)[1]=="SOCKcluster"){
		out=clusterApply(mpicl,x,fun,...)
	}else{
		if(class(mpicl)[1]=="environment"){
			arg_list=list()
			arg_list[[1]]=x[[1]]
			params=list(...)
			arg_list=append(x,params)
			fun_arg=deparse(substitute(fun))
			# tt=eval("mywhere(\"IS_DATA_CONTAINER_ENVIRONMENT\")",mpicl)
			# print(tt)
			# print("ala")
			# print(fun_arg)
			# print(length(arg_list))
			# print("mm")
			# print(names(arg_list))
			# fff=paste(fun_arg,"(\"",x,"\")",sep="")
			# print(fff)
			# out0=eval(parse(text=fff),envir=mpicl)
#			out0=do.call(fun,arg_list,quote=FALSE,envir=mpicl)
			out0=do.call(fun_arg,arg_list,quote=FALSE,envir=mpicl)
			out=list()
			out[[1]]=out0
		}else{
			stop("neither cluster nor environment")
		}
	}
	return(out)
}

###
#myenvs2=clusterCall_or_envCall(mpicl, to_be_run_gamma1=gamma1,n=n,theta_lambda=theta_lambda,Eu_ak=Eu_ak,calc_L=calc_L)
clusterCall_or_envCall=function(mpicl,fun,...){
	if(class(mpicl)[1]=="SOCKcluster"){
		out=clusterCall(mpicl,fun,...)
	}else{
		if(class(mpicl)[1]=="environment"){
			arg_list=list(...)
			fun_arg=deparse(substitute(fun))
			tt=eval(parse(text="mywhere(\"IS_DATA_CONTAINER_ENVIRONMENT\")"),mpicl)
#			print(tt)
#			print("ala")
#			print(fun_arg)
#			out0=do.call(fun,arg_list,quote=FALSE,envir=mpicl)
			out0=do.call(fun_arg,arg_list,quote=FALSE,envir=mpicl)
			out=list()
			out[[1]]=out0
		}else{
			stop("neither cluster nor environment")
		}
	}
	return(out)
}

makePSOCKcluster_or_environment=function(ncores, useEnvIfOne=TRUE){
	if(ncores==1 & useEnvIfOne==TRUE){
		mpicl=new.env()
	}else{
		print("number of cluster nodes")
		print(ncores)
		mpicl=makePSOCKcluster(rep("localhost",ncores))
	}
	IS_DATA_CONTAINER_ENVIRONMENT=TRUE
	clusterExport_or_envExport(mpicl,"IS_DATA_CONTAINER_ENVIRONMENT", environment())
#		mpicl=makeSOCKcluster(rep("localhost",ncores))
	return(mpicl)
}

get_data_container_environment=function(){
	return(mywhere("IS_DATA_CONTAINER_ENVIRONMENT"))
}

clusterEvalQ_or_envEvalQ=function(mpicl,myfile="Code/random_effects/variational_aux_kappa_alphai7.R"){
	if(class(mpicl)[1]=="SOCKcluster"){
		mycmd=paste("source(\"",myfile,"\")",sep="")
		parsed=parse(text=mycmd)
		clusterCall(mpicl, eval,parsed, env = .GlobalEnv)
	}else{
		if(class(mpicl)[1]=="environment"){
			source(myfile,local=mpicl)
		}else{
			stop("neither cluster nor environment")
		}
	}
	return(mpicl)
}

#' compares two lists
#'
#' Checks whether elements in nested  named lists are the same up to a numerical margin of error. Does not check array names or attributes, etc so use with caution.
#' @return A String informing about the lists are judged the same.
#' @param listA
#' List. First list to compare.
#' @param listB
#' List. Second list to compare.
#' @param abs_err
#' Scalar. Absolute error allowed (per element.)
#' @param rel_err
#' Scalar. Error allowed relative to largest element (per element.)
#' @export
compare_lists=function(listA,listB,outStr="",abs_err=0,rel_err=1){
	namesA=names(listA)
	namesB=names(listB)
	if(!is.list(listA)){
		return("listA not a list.")
	}
	if(!is.list(listB)){
		return("listB not a list.")
	}
	if(length(namesA)!=length(namesB)){
		return(paste(outStr,"lists not equal length"))
	}
	if(sum(!is.element(namesA,namesB))!=0){
		return(paste(outStr,"lists not same names"))
	}
	if(length(namesA)==0){
		names(listA)=paste("el",c(1:length(listA)),sep="")
		names(listB)=paste("el",c(1:length(listB)),sep="")
		namesA=names(listA)
		namesB=names(listB)
	}
	for (name in namesA){
	#	print(name)
		lA=listA[[name]]
		lB=listB[[name]]
		if(length(class(lA))!=length(class(lB))){
			return(paste(outStr,":",name,"elements not same class (one has more class memberships)."))
		}
		if(sum(!is.element(class(lA),class(lB)))!=0){
			return(paste(outStr,":",name,"elements not same class"))
		}
		if(is.list(lA)){
			outval=compare_lists(lA,lB,outStr=paste(outStr,":",name),abs_err=abs_err,rel_err=rel_err)
			if(outval!="both the same"){
				return(outval)
			}else{
				next
			}
		}
		nullcounts=(is.null(dim(lA))) + (is.null(dim(lB)))
		if(nullcounts==1){
			return(paste(outStr,":",name,"one is null"))
		}	
		if(nullcounts==2){
			if(length(lA)!=length(lB)){
				return(paste(outStr,":",name,"length is not equal"))
			}
			if(abs_err==0 || !is.numeric(as.vector(lA))){
				if(sum(lA!=lB)!=0){
					return(paste(outStr,":",name,"both dim null but elements not the same"))
				}
			}else{
				# capture infs as comparison fails on inf
				same_inf_index=(lA==Inf & lB==Inf) | (lA==-Inf & lB==-Inf)
				lA[same_inf_index]=0
				lB[same_inf_index]=0
				if(sum(abs(lA-lB)>abs_err)!=0){
					return(paste(outStr,":",name,"both dim null but some numerical elements above abs_err"))
				}
				top_val=max(max(abs(lA)),max(abs(lB)))
				if(top_val>0){
					if(sum(abs(lA-lB)/top_val>rel_err)!=0){
						return(paste(outStr,":",name,"both dim null but some numerical elements above rel_err"))
					}
				}
			}
		}
		if(nullcounts==0){
			if(sum(dim(lA)!=dim(lB))!=0){
				return(paste(outStr,":",name,"dimensions not the same"))
			}
			if(abs_err==0 || !is.numeric(as.vector(lA))){
				if(sum(lA!=lB)!=0){
					return(paste(outStr,":",name,"both dims same but elements not the same"))
				}
			}else{				
				# capture infs as comparison fails on inf
				same_inf_index=(lA==Inf & lB==Inf) | (lA==-Inf & lB==-Inf)
				lA[same_inf_index]=0
				lB[same_inf_index]=0
				if(sum(abs(lA-lB)>abs_err)!=0){
					return(paste(outStr,":",name,"both dims same but some numerical elements above abs_err"))
				}
				top_val=max(max(abs(lA)),max(abs(lB)))
				if(top_val>0){
					if(sum(abs(lA-lB)/top_val>rel_err)!=0){
						return(paste(outStr,":",name,"both dim null but some numerical elements above rel_err"))
					}
				}
			}
		}
	}
	return("both the same")
}


compare_lists_test=function(){
	lA=list(A=3,B=4,C=matrix(1.5,2,2))
	lB=list(A=3,B=4,C=matrix(1.5,2,2))
	print(sum("both the same"!=compare_lists(lA,lB)))
	lB=list(A=3,B=matrix(2,2,2),C=matrix(1.5,2,2))
	print(sum(" : B elements not same class"!=compare_lists(lA,lB)))
	#lB=list(A=3,B=list(a=1,b=3),C=matrix(1.5,2,2))
	lB=list(A=3,B=4,C=matrix(1.5,4,2))
	print(sum(" : C dimensions not the same"!=compare_lists(lA,lB)))
	lB=list(A=3,B=c(4,2),C=matrix(1.5,2,2))
	print(sum(" : B length is not equal"!=compare_lists(lA,lB)))
	lA=list(A=3,B=list(a=2,b=c(2,1)),C=matrix(1.5,2,2))
	lB=list(A=3,B=list(a=2,b=c(1,1)),C=matrix(1.5,2,2))
	print(sum(" : B : b both dim null but elements not the same"!=compare_lists(lA,lB)))

	lA=list(A=3.1,B=list(a=2,b=c(2,1)),C=matrix(1.5,2,2))
	lB=list(A=3,B=list(a=2,b=c(2,1)),C=matrix(1.5,2,2))
	print(sum("both the same"!=compare_lists(lA,lB,abs_err=0.11)))
	lA=list(A=3.12,B=list(a=2,b=c(2,1)),C=matrix(1.5,2,2))
	lB=list(A=3,B=list(a=2,b=c(2,1)),C=matrix(1.5,2,2))
	print(sum(" : A both dim null but some numerical elements above abs_err"!=compare_lists(lA,lB,abs_err=0.11)))
	lA=list(A=3,B=list(a=2,b=c(2,1.1)),C=matrix(1.5,2,2))
	lB=list(A=3,B=list(a=2,b=c(2,1)),C=matrix(1.5,2,2))
	print(sum("both the same"!=compare_lists(lA,lB,abs_err=0.11)))
	lA=list(A=3,B=list(a=2,b=c(Inf,1.1)),C=matrix(1.5,2,2))
	lB=list(A=3,B=list(a=2,b=c(Inf,1)),C=matrix(1.5,2,2))
	print(sum("both the same"!=compare_lists(lA,lB,abs_err=0.11)))
	lA=list(A=3,B=list(a=2,b=c(Inf,1.2)),C=matrix(1.5,2,2))
	lB=list(A=3,B=list(a=2,b=c(Inf,1)),C=matrix(1.5,2,2))
	print(sum(" : B : b both dim null but some numerical elements above abs_err"!=compare_lists(lA,lB,abs_err=0.11)))

	lA=list(A=3,B=list(a=2,b=c(Inf,1)),C=matrix(1000.01,2,3))
	lB=list(A=3,B=list(a=2,b=c(Inf,1)),C=matrix(1000.00,2,3))
	print(sum("both the same"!=compare_lists(lA,lB,abs_err=1,rel_err=1E-4)))
	lA=list(A=3,B=list(a=2,b=c(Inf,1)),C=matrix(1000.01,2,3))
	lB=list(A=3,B=list(a=2,b=c(Inf,1)),C=matrix(1000.00,2,3))
	print(sum(" : C both dim null but some numerical elements above rel_err"!=compare_lists(lA,lB,abs_err=1,rel_err=1E-7)))
}

my_tmpdir=function(){
	return(paste(getwd(),"/interimData/tmpdir",sep=""))
}

get_tmpfile=function(pattern="",suffix=""){
	out=tempfile(pattern=pattern,tmpdir=my_tmpdir())
	if(suffix!=""){
		out=sub("$",paste(".",suffix,sep=""),out)
	}
	attributes(out)$is_deletable=TRUE
	class(out)=unique(c("filepath",class(out)))
	return(out)
}

###tested
calc_fast_cov_solve=function(A){
	#checking if all potentially all eigenvalues are non-positive (in that case, we can inverte -A instead and flip the result) 
	my_sign=1
	if(sum(diag(A))<0){
		my_sign= -1
	}
	result = tryCatch(
		{
			return(chol2inv(chol(my_sign*A))*my_sign)
		}, error = function(e) {
			N=length(A[,1])
			results = chol2inv(chol(my_sign*A+diag(N)*0.0005))*my_sign
			print("direct inversion was not possilbe tried it again with small diagonal factor")
			return(results)
		}, finally = {
		}
	)
	return(result)
}

######## tested
calc_Eu_gamma=function(shape,rate){ 
	lambda1=shape
	lambda2=rate
	Eu1_gamma=-log(lambda2)+digamma(lambda1)
	Eu2_gamma=lambda1/lambda2
	Eu_gamma=c(Eu1_gamma,Eu2_gamma)
	return(Eu_gamma)
}
######## tested indirectly
calc_Eu_gamma_func_via_theta=function(theta){
	lambda1=theta[1]+1
	lambda2=theta[2]*(-1)
	return(calc_Eu_gamma(lambda1,lambda2))
}


# tested
calc_theta_lambda=function(lambda1,lambda2,genespecific_lambda1=FALSE){
	if(genespecific_lambda1){
		theta_lambda=cbind(lambda1-1,-lambda2)
	}else{
		theta_lambda=c(lambda1-1,-lambda2)
	}
	return(theta_lambda)
}

#tested
calc_theta_kappa=function(tau1,tau2){
	theta_kappa=c(tau1-1,-tau2)
	return(theta_kappa)
}

#nottested
calc_theta_lambda2=function(rho1,rho2){
	theta_lambda2=c(rho1-1,-rho2)
	return(theta_lambda2)
}

#nottested
calc_theta_tau2=function(xi1,xi2){
	theta_tau2=c(xi1-1,-xi2)
	return(theta_tau2)
}

#tested
calc_theta_alpha=function(alpha1,alpha2){
	theta_alpha=c(alpha1-1,-alpha2)
	return(theta_alpha)
}

#tested
calc_Eu_lambda_via_theta=function(theta){
	return(calc_Eu_gamma_func_via_theta(theta))
}

calc_Eu_lambda2_via_theta=function(theta){
	return(calc_Eu_gamma_func_via_theta(theta))
}

calc_Eu_tau2_via_theta=function(theta){
	return(calc_Eu_gamma_func_via_theta(theta))
}

calc_Eu_kappa_via_theta=function(theta){
	return(calc_Eu_gamma_func_via_theta(theta))
}

######## tested
calc_Eu_alpha_via_theta=function(theta){
	return(calc_Eu_gamma_func_via_theta(theta))
}

#tested
calc_theta_ak=function(phi1,phi2,t){
	theta_phi=c(phi1-1,-phi2)
	theta_phi_mat=kronecker(rep(1,t),t(theta_phi))
	return(theta_phi_mat)
}
######## tested
calc_Eu_ak_via_theta=function(theta){
	return(t(apply(theta,1,calc_Eu_gamma_func_via_theta)))
}

######## not tested
calc_Eu_upsilon_via_theta_list=function(theta_list){
	outl=lapply(theta_list,function(theta){
		ret=t(apply(theta,1,calc_Eu_gamma_func_via_theta))
		return(ret)
	})
	return(outl)
}

calc_Eu_upsilon_via_theta=function(theta){
	ret=t(apply(theta,1,calc_Eu_gamma_func_via_theta))
	return(ret)
}


# not tested
calc_Eu_chi2j_via_theta=function(theta){
	ret=t(apply(theta,1,calc_Eu_gamma_func_via_theta))
	return(ret)
}

#not tested
calc_theta_upsilon=function(chi1_vec,chi2_vec,h_vec){
	theta_upsilon_list=list()
	for(i in c(1:length(h_vec))){
		cur_chi1=rep(chi1_vec[i],h_vec[i])
		cur_chi2=rep(chi2_vec[i],h_vec[i])
		theta_upsilon_list[[i]]=cbind(cur_chi1-1,-cur_chi2)
	}
	return(theta_upsilon_list)
}

calc_theta_upsilon_via_Eu_chi2j=function(chi1_vec,Eu_chi2j,h_vec){
	chi2_vec=Eu_chi2j[,2]
	theta_upsilon_list=list()
	if(length(chi1_vec)!=length(h_vec)){
		stop("erromsddc.")
	}
	if(length(chi1_vec)!=dim(Eu_chi2j)[1]){
		stop("errp,s]dcs")
	}
	for(i in c(1:length(h_vec))){
		cur_chi1=rep(chi1_vec[i],h_vec[i])
		cur_chi2=rep(chi2_vec[i],h_vec[i])
		theta_upsilon_list[[i]]=cbind(cur_chi1-1,-cur_chi2)
	}
	return(theta_upsilon_list)
}

#not tested
calc_theta_chi2j=function(zeta1_vec,zeta2_vec){
	calc_theta_chi2j=cbind(zeta1_vec-1,-zeta2_vec)
	return(calc_theta_chi2j)
}


######## tested
calc_Eu_alphai_via_theta=function(theta){
	return(t(apply(theta,1,calc_Eu_gamma_func_via_theta)))
}

######## to test
calc_Eu_alphai_via_theta_vectorized=function(theta){
	Eu1_gamma=-log(-theta[,2])+digamma(theta[,1]+1)
	Eu2_gamma=c(-(theta[,1]+1)/theta[,2])
	return(cbind(Eu1_gamma,Eu2_gamma))
}

######## to test
calc_Eu_alphai_via_theta_only_second=function(theta){
	Eu2_gamma=c(-(theta[,1]+1)/theta[,2])
	return(Eu2_gamma)
}

######## to test
calc_Eu_alphai_via_theta_only_second_withpreindexing=function(theta,preindex_k){
	Eu2_gamma_k=c(-(theta[preindex_k,1]+1)/theta[preindex_k,2])
	return(Eu2_gamma_k)
}

######## to test
calc_Eu_kappa_via_theta_vectorized=function(theta){
	Eu1_gamma=-log(-theta[,2])+digamma(theta[,1]+1)
	Eu2_gamma=c(-(theta[,1]+1)/theta[,2])
	return(cbind(Eu1_gamma,Eu2_gamma))
}

######## to test
calc_Eu_kappa_via_theta_only_second=function(theta){
	Eu2_gamma=c(-(theta[,1]+1)/theta[,2])
	return(Eu2_gamma)
}

######## tested
calc_Eu_gammai_via_Eu_ak=function(Eu_ak,mapping_mat){
	m=dim(mapping_mat)[1]
	aMat0=kronecker(rep(1,m),t(Eu_ak[,2]))
	gammai=exp(Matrix::rowSums(log(aMat0)*mapping_mat))
	loggammai=Matrix::rowSums(kronecker(rep(1,m),t(Eu_ak[,1]))*mapping_mat)
	Eu_gamma=cbind(loggammai,gammai)
	return(Eu_gamma)
}

######## tested
calc_Eu_gammai_via_Eu_ak_only_second=function(Eu_ak,mapping_mat){
	m=dim(mapping_mat)[1]
	Eu_gamma=matrix(0,m,2)
	aMat0=kronecker(rep(1,m),t(Eu_ak[,2]))
	gammai=exp(rowSums(log(aMat0)*mapping_mat))
	Eu_gamma[,2]=gammai
	return(Eu_gamma)
}

# ######## tested
# calc_Eu_gammai_via_Eu_ak=function(Eu_ak,mapping_mat){
# 	m=dim(mapping_mat)[1]
# 	aMat0=kronecker(rep(1,m),t(Eu_ak[,2]))
# 	gammai=exp(rowSums(log(aMat0)*mapping_mat))
# 	loggammai=rowSums(kronecker(rep(1,m),t(Eu_ak[,1]))*mapping_mat)
# 	Eu_gamma=cbind(loggammai,gammai)
# 	return(Eu_gamma)
# }


######## tested
update_Eu_gammai_second_row_via_Eu_ak=function(Eu_gammai,Eu_ak_old, Eu_ak,mapping_mat,k){
	m=dim(mapping_mat)[1]
#	Eu_gamma=matrix(0,m,2)
	indices=which(mapping_mat[,k]!=0)
	updated=mapping_mat[indices,k]*(Eu_ak[k,2]/Eu_ak_old[k,2])
	Eu_gammai[indices,2]=Eu_gammai[indices,2]*updated
	return(Eu_gammai)
}

update_Eu_gammai_second_row_via_Eu_ak_only_second_with_preindexing=function(Eu_gammai,Eu_ak_old, Eu_ak,mapping_mat,k,preindex_k){
	m=dim(mapping_mat)[1]
#	Eu_gamma=matrix(0,m,2)
#	indices=which(mapping_mat[,k]!=0)
#	updated=mapping_mat[preindex_k,k]*(Eu_ak[k,2]/Eu_ak_old[k,2])
	out=Eu_gammai[preindex_k,2]*mapping_mat[preindex_k,k]*(Eu_ak[k,2]/Eu_ak_old[k,2])
	return(out)
}

update_Eu_gammai_second_row_via_Eu_ak_only_second_with_preindexing_only_ones=function(Eu_gammai,Eu_ak_old, Eu_ak,k,preindex_k){
	out=Eu_gammai[preindex_k,2]*(Eu_ak[k,2]/Eu_ak_old[k,2])
	return(out)
}
######## tested
calc_g_alphai_via_theta=function(theta){
	g=(theta[,1]+1)*log(-theta[,2])-lgamma(theta[,1]+1)
	return(g)
}


######## tested
calc_g_alphai_via_Eu_gammai=function(gamma1,Eu_gammai){
	return((gamma1)*Eu_gammai[,1]-lgamma(gamma1))
}

######## not tested
calc_g_alphai_via_Eu_gammai_Eu_kappa=function(gamma1,Eu_gammai,Eu_kappa){
	myout=gamma1*(Eu_gammai[,1]+Eu_kappa[,1])-lgamma(gamma1)
	return(myout)
}

calc_theta_ak_via_Eu=function(Eu_ak){
	res=cbind(Eu_ak[,1]-1,-Eu_ak[,2])
	return(res)
}

calc_theta_alphai_via_Eu_gammai=function(gamma1,Eu_gammai,m){
	res=cbind(rep(gamma1-1,m),-Eu_gammai[,2])
	return(res)
}

calc_theta_kappaj_via_Eu_tau2=function(tau1,Eu_tau2,m){
	res=cbind(rep(tau1-1,m),-rep(Eu_tau2[2],m))
	return(res)
}

calc_theta_alphai_via_Eu_gammai_second=function(Eu_gammai){
	res=-Eu_gammai[,2]
	return(res)
}

calc_theta_alphai_via_Eu_gammai_Eu_kappa_second=function(Eu_gammai,Eu_kappa){
	res=-Eu_gammai[,2]*Eu_kappa[,2]
	return(res)
}

calc_theta_alphai_via_Eu_gammai_Eu_kappa_second_with_preindexing=function(Eu_gammai,Eu_kappa,preindex_k){
	res_at_preindex=-Eu_gammai[preindex_k,2]*Eu_kappa[preindex_k,2]
	return(res_at_preindex)
}


calc_theta_alphai_via_Eu_gammai_Eu_kappa=function(gamma1,Eu_gammai,Eu_kappa,m=dim(Eu_gammai)[1]){
	res=cbind(rep(gamma1-1,m),-(Eu_kappa[2]*Eu_gammai[,2]))
	return(res)
}

calc_theta_alphai_kappa_via_Eu_gammai_Eu_alphai=function(gamma1,Eu_gammai,Eu_alphai,m=length(Eu_gammai[,2])){
	res=cbind(rep(gamma1,m),-Eu_gammai[,2]*Eu_alphai[,2])
	summedres=colSums(res)
	return(summedres)
}

calc_theta_alphai_kappa_via_Eu_gammai_Eu_alphai_vect=function(gamma1,Eu_gammai,Eu_alphai,m,gene_indices){
	res_dt=data.table(a=rep(gamma1,m),b=-Eu_gammai[,2]*Eu_alphai[,2],gene_indices=gene_indices)
	summedres=res_dt[,list(asum=rep(sum(a),length(a)),bsum=rep(sum(b),length(a))),by=gene_indices]
	return(summedres)
}

calc_theta_alphai_kappa_via_Eu_gammai_Eu_alphai_from_dt=function(Eu_gammai,Eu_alphai,dt){
	dt[,b:=-Eu_gammai[,2]*Eu_alphai[,2]]
	summedres=dt[,list(asum=rep(sum(a),length(a)),bsum=rep(sum(b),length(a))),by=gene_indices.V1]
	return(summedres)
}

calc_theta_alphai_kappa_via_Eu_gammai_Eu_alphai_second_from_mat=function(Eu_gammai,Eu_alphai,index_mat){
	res=index_mat%*%(-Eu_gammai[,2]*Eu_alphai[,2])
	return(res)
}

calc_theta_alphai_kappa_via_Eu_gammai_Eu_alphai_second_from_mat_update=function(Eu_gammai,Eu_alphai,index_mat,preindex_k){
	res=index_mat[,preindex_k]%*%(-Eu_gammai[preindex_k,2]*Eu_alphai[preindex_k,2])
	return(res)
}

calc_theta_lambdai_lambda2_via_Eu_lambdai=function(lambda1,Eu_lambdai,m=dim(Eu_lambdai)[1],genespecific_lambda1=FALSE){
	if(!genespecific_lambda1){
		if(length(lambda1)!=1){
			stop("err23eodm2d")
		}
		res=cbind(rep(lambda1,m),-Eu_lambdai[,2])
	}else{
		if(length(lambda1)!=length(Eu_lambdai[,2])){
			stop("err2-0k23e")
		}
		res=cbind(lambda1,-Eu_lambdai[,2])
	}
	summedres=colSums(res)
	return(summedres)
}

calc_theta_kappaj_tau2_via_Eu_kappaj=function(tau1,Eu_kappaj,m=length(Eu_kappaj[,2])){
	res=cbind(rep(tau1,m),-Eu_kappaj[,2])
	summedres=colSums(res)
	return(summedres)
}

######## to test
calc_g_a1_via_theta=function(theta){
	g=(theta[1]+1)*log(-theta[2])-lgamma(theta[1]+1)
	return(g)
}

######## to test
calc_g_ak_via_theta=function(theta){
	g=(theta[,1]+1)*log(-theta[,2])-lgamma(theta[,1]+1)
	return(g)
}


######## tested
calc_theta_alphai_ak_via_Eu=function(gamma1,Eu_alphai,Eu_ak,mapping_mat,m,t,i,k){
	if(mapping_mat[i,k]==0){
		return(c(0,0))
	}
	out=c(gamma1*mapping_mat[i,k],-exp(sum(mapping_mat[i,-k]*log(Eu_ak[-k,2])))*Eu_alphai[i,2])
	return(out)
}

calc_theta_alphai_ak_kappa_via_Eu=function(gamma1,Eu_alphai,Eu_ak,Eu_kappa,mapping_mat,m,t,i,k){
	if(mapping_mat[i,k]==0){
		return(c(0,0))
	}
	out=c(gamma1*mapping_mat[i,k],-exp(sum(mapping_mat[i,-k]*log(Eu_ak[-k,2])))*Eu_alphai[i,2]*Eu_kappa[i,2])
	return(out)
}


########  tested
calc_theta_alpha_ak_via_Eugammai=function(Eu_gammai,gamma1,Eu_alphai,Eu_kappa,Eu_ak,mapping_mat,k){
	divider=mapping_mat[,k]/Eu_ak[k,2]
	sec_col=-sum(Eu_gammai[,2]*divider*Eu_alphai[,2]*Eu_kappa[2])
	first_col=gamma1*sum(mapping_mat[,k])
	return(c(first_col,sec_col))
}

########  tested
calc_theta_alpha_ak_via_Eugammai_vect=function(Eu_gammai,gamma1,Eu_alphai,Eu_kappa,Eu_ak,mapping_mat,k){
	divider=mapping_mat[,k]/Eu_ak[k,2]
	sec_col=-sum(Eu_gammai[,2]*divider*Eu_alphai[,2]*Eu_kappa[,2])
	first_col=gamma1*sum(mapping_mat[,k])
	return(c(first_col,sec_col))
}

calc_theta_alpha_ak_via_Eugammai_vect_only_second=function(Eu_gammai,gamma1,Eu_alphai,Eu_kappa,Eu_ak,mapping_mat,k){
#	divider=mapping_mat[,k]/Eu_ak[k,2]
#	sec_col=-sum(Eu_gammai[,2]*divider*Eu_alphai[,2]*Eu_kappa[,2])
#	divider=mapping_mat[,k]/Eu_ak[k,2]
	sec_col=-sum(Eu_gammai[,2]*mapping_mat[,k]*Eu_alphai[,2]*Eu_kappa[,2]/Eu_ak[k,2])
	return(sec_col)
}

calc_theta_alpha_ak_via_Eugammai_vect_only_second_with_preindexing=function(Eu_gammai,gamma1,Eu_alphai,Eu_kappa,Eu_ak,mapping_mat,k,non_null_mapping_mat_inds){
	sec_col=-sum(Eu_gammai[non_null_mapping_mat_inds,2]*mapping_mat[non_null_mapping_mat_inds,k]*Eu_alphai[non_null_mapping_mat_inds,2]*Eu_kappa[non_null_mapping_mat_inds,2])/Eu_ak[k,2]
	return(sec_col)
}

calc_theta_alpha_ak_via_Eugammai_vect_only_second_with_preindexing_only_ones=function(Eu_gammai,gamma1,Eu_alphai,Eu_kappa,Eu_ak,k,non_null_mapping_mat_inds){

	sec_col=-sum(Eu_gammai[non_null_mapping_mat_inds,2]*Eu_alphai[non_null_mapping_mat_inds,2]*Eu_kappa[non_null_mapping_mat_inds,2])/Eu_ak[k,2]
	return(sec_col)
}



########  tested
calc_theta_alpha_aks_via_Eu=function(gamma1,Eu_alphai,Eu_kappa,Eu_ak,mapping_mat,m,t){
	out=matrix(0,m,t)
##	tt=kronecker(rep(1,m),t(Eu_ak[,2]))
	Eu_ak_mat=kronecker(rep(1,m),t(Eu_ak[,2]))
	print("a3")
	Eu_ak_prod=exp(rowSums(log(Eu_ak_mat)*mapping_mat))
	Eu_ak_prod_mat=kronecker(Eu_ak_prod,t(rep(1,t)))
	Eu_ak_prod_mat_minus_k=(Eu_ak_prod_mat/Eu_ak_mat)
	Eu_alphai_mat=kronecker(Eu_alphai[,2],t(rep(1,t)))
	print("a4")
	Eu_alphai_mat=Eu_alphai_mat*mapping_mat
	upper=-Eu_ak_prod_mat_minus_k*Eu_alphai_mat
	theta_alphai_ak=cbind(gamma1*c(colSums(mapping_mat)),colSums(upper)*Eu_kappa[2])
	return(theta_alphai_ak)
}


######## to test
calc_theta_alpha_ak_via_Eu=function(gamma1,Eu_alphai,Eu_ak,mapping_mat,m,t){
	out=matrix(0,m,t)
##	tt=kronecker(rep(1,m),t(Eu_ak[,2]))
	Eu_ak_mat=kronecker(rep(1,m),t(Eu_ak[,2]))
	Eu_ak_prod=exp(rowSums(log(Eu_ak_mat)*mapping_mat))
	print("a1")
	Eu_ak_prod_mat=kronecker(Eu_ak_prod,t(rep(1,t)))
	Eu_ak_prod_mat_minus_k=(Eu_ak_prod_mat/Eu_ak_mat)
	Eu_alphai_mat=kronecker(Eu_alphai[,2],t(rep(1,t)))
	print("a2")
	Eu_alphai_mat=Eu_alphai_mat*mapping_mat
	upper=-Eu_ak_prod_mat_minus_k*Eu_alphai_mat
	theta_alphai_ak=cbind(gamma1*c(colSums(mapping_mat)),colSums(upper))
	return(theta_alphai_ak)
}

######## to test
calc_theta_alpha_ak_kappa_via_Eu=function(gamma1,Eu_alphai,Eu_ak,Eu_kappa,mapping_mat,m,t){
	out=matrix(0,m,t)
##	tt=kronecker(rep(1,m),t(Eu_ak[,2]))
	Eu_ak_mat=kronecker(rep(1,m),t(Eu_ak[,2]))
	Eu_ak_prod=exp(rowSums(log(Eu_ak_mat)*mapping_mat))
	Eu_ak_prod_mat=kronecker(Eu_ak_prod,t(rep(1,t)))
	Eu_ak_prod_mat_minus_k=(Eu_ak_prod_mat/Eu_ak_mat)
	Eu_alphai_mat=kronecker(Eu_alphai[,2],t(rep(1,t)))
	Eu_alphai_mat=Eu_alphai_mat*mapping_mat
	stop("not ready1")
	upper=-Eu_ak_prod_mat_minus_k*Eu_alphai_mat
	theta_alphai_ak=cbind(gamma1*c(colSums(mapping_mat)),colSums(upper))
	return(theta_alphai_ak)
}

######## tested
calc_Eu_alphai_via_theta=function(theta){
	m=length(theta[,1])
	out=theta-theta
	for(i in c(1:m)){
		out[i,]=calc_Eu_gamma_func_via_theta(theta[i,])
	}
	return(out)
}

######## tested
calc_Eu_b_via_theta=function(theta_list){
	Eu_b_list=list()
	theta_M=-2*(theta_list[[2]])
	theta_M_inv=calc_fast_cov_solve(theta_M)
	Eu_b_1=theta_M_inv%*%theta_list[[1]]
	Eu_b_2=theta_M_inv + Eu_b_1%*%t(Eu_b_1)
	Eu_b_list[["Eu_b"]]=Eu_b_1
	Eu_b_list[["Eu_b2"]]=Eu_b_2
	return(Eu_b_list)
}

######## tested
calc_Eu_b_via_theta4eigen=function(theta_list,decomposed=FALSE,missing_variance=0,myc=-1){
	## (explains decomposed=TRUE)
	## 
	## theta_list: decomposed version to calculate  Eu_b2_formula:
	## TODO: fill in 
	## arguments theta_list
	##
	##
	## value
	## Eu_b_list: decomposed version to calculate  Eu_b2_formula as well as Eu_b1 (allows to compute Eu_b2_ fast, while at the same time having smaller memory footprint. (These components can be used in later steps more efficiently than the actual Eu_b2 result)):
	## Eu_b_list[["Eu_b1"]]: Eu_b1
	## Eu_b_list[["Eu_b2_formula_expanded"]]: formula to calculate  Eu_b2 from the list elements.
	## Eu_b_list[["Eu_b2_formula"]]: collapsed version of Eu_b2_formula_expanded
	## Eu_b_list[["theta_M_inv"]]: output from calculate_woodbury_inversion_sym(), i. e.: elements to calculate --\phi^*_b(2)
	Eu_b_list=list()
	diag2inverse=(-1)*theta_list[["sec_diag"]]
	cXtXsqrt2inverse=theta_list[["sec_sqrt_abs"]]
	if(missing_variance!=0){
		mym=length(diag2inverse)
		missing_frac=missing_variance/mym
		diag2inverse=diag2inverse+missing_frac*(-myc)
	}
	if(decomposed==FALSE){
		my_inv=(0.5)*calculate_woodbury_inversion_sym(diagDprime_alpha=diag2inverse,cXtX_sqrt=cXtXsqrt2inverse,only_trace=FALSE,decomposed=FALSE)
		theta_M_inv=my_inv
		Eu_b_1=theta_M_inv%*%theta_list[["first"]]
		Eu_b_2=theta_M_inv + Eu_b_1%*%t(Eu_b_1)
		Eu_b_list[["Eu_b"]]=Eu_b_1
		Eu_b_list[["Eu_b2"]]=Eu_b_2
	}else{
		Eu_b_list[["Eu_b2_formula"]]="0.5*process(theta_M_inv) + Eu_b_1 %*% t(Eu_b_1))"
		Eu_b_list[["Eu_b2_formula_expanded"]]="0.5*(diag(Eu_b_list$theta_M_inv$diagDprime_alphainv)-Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt%*%t(Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt)) + Eu_b_list$Eu_b_1 %*% t(Eu_b_list$Eu_b_1)"
		my_inv_decomposed=calculate_woodbury_inversion_sym(diagDprime_alpha=diag2inverse,cXtX_sqrt=cXtXsqrt2inverse,only_trace=FALSE,decomposed=TRUE)
		LL=t(my_inv_decomposed[["cXtX_sqrtAug_interimRes_sqrt"]])%*%theta_list[["first"]]
		LL=my_inv_decomposed[["cXtX_sqrtAug_interimRes_sqrt"]]%*%LL
		Eu_b_1=0.5*(my_inv_decomposed[["diagDprime_alphainv"]]*theta_list[["first"]]-LL)
		Eu_b_list[["Eu_b_1"]]=Eu_b_1
		Eu_b_list[["theta_M_inv"]]=my_inv_decomposed
	}
	return(Eu_b_list)
}

######## tested
calc_theta_b=function(alpha,m){
	theta_b_list=list()
	theta_b_1=matrix(0,m,1)
	theta_b_2= -alpha*diag(m)/2
	theta_b_list[[1]]=theta_b_1
	theta_b_list[[2]]=theta_b_2
	return(theta_b_list)
}

######## tested
calc_theta_b_from_natural=function(mu,Sigma){
	SigmaInv=calc_fast_cov_solve(Sigma)
	theta_b_list=list()
	theta_b_1=SigmaInv%*%mu
	theta_b_2= -SigmaInv/2
	theta_b_list[[1]]=theta_b_1
	theta_b_list[[2]]=theta_b_2
	return(theta_b_list)
}

######## tested
calc_natural_from_theta_b=function(theta_b_list){
	SigmaInv=-theta_b_list[[2]]*2
	Sigma=calc_fast_cov_solve(SigmaInv)
	mu=Sigma%*%theta_b_list[[1]]
	natural_list=list()
	natural_list[[1]]=mu
	natural_list[[2]]=Sigma
	return(natural_list)
}

######## tested
calc_g_b_from_natural=function(Sigma,mu){
	m=length(Sigma[,1])
	if(m==1){
		M=-0.5*log(Sigma)
	}else{
		mydet=determinant(Sigma,logarithm=TRUE)
		if(mydet$sign!=1){stop("degenerate")}
		M=-0.5*mydet$modulus
	}
	N=-m/2*log(2*pi)
	F=-0.5*t(mu)%*%calc_fast_cov_solve(Sigma)%*%mu
	res=M+N+F
	return(res)
}

######## tested
calc_f_b_from_natural=function(){
	return(0)
}

######## tested
calc_g_b=function(alpha,m){
	return(m*log(alpha)/2-m/2*log(2*pi))
}

######## tested
calc_g_lambda=function(lambda1,lambda2){
	return(lambda1*log(lambda2)-lgamma(lambda1))
}

######## not tested
calc_g_lambda_via_Eu_lambda2=function(lambda1,Eu_lambda2){
	return(lambda1*Eu_lambda2[1]-lgamma(lambda1))
}

calc_g_kappa=function(tau1,tau2){
	return(tau1*log(tau2)-log(gamma(tau1)))
}

######## tested
calc_g_alpha=function(alpha1,alpha2){
	return(alpha1*log(alpha2)-lgamma(alpha1))
}

######## tested
calc_theta_y_lambda=function(n,b,ytX, XtX,yty){
	theta_y_lambda_1=n/2
	theta_y_lambda_2=-0.5*(yty-2*ytX%*%b+t(b)%*%XtX%*%b)
	theta_y_lambda=c(theta_y_lambda_1,theta_y_lambda_2)
	return(theta_y_lambda)
}

######## not tested
calc_theta_z_lambda=function(myn,b,z,Sigma){
	mym=dim(Sigma)[1]
	theta_z_lambda_1=mym/2
	theta_z_lambda_2=-0.5*(t(z)%*%Sigma%*%z-2*z%*%b*sqrt(myn)+myn*t(b)%*%Sigma%*%b)
	theta_z_lambda=c(theta_z_lambda_1,theta_z_lambda_2)
	return(theta_z_lambda)
}

######## tested
calc_g_y_lambda=function(n){
	return(-n/2*log(2*pi))
}

######## tested
calc_theta_y_b=function(lambda,ytX, XtX){
	theta_y_b_list=list()
	theta_y_b_1=lambda*t(ytX)
	theta_y_b_2=-lambda*XtX/2
	theta_y_b_list[[1]]=theta_y_b_1
	theta_y_b_list[[2]]=theta_y_b_2
	return(theta_y_b_list)
}

######## not tested
calc_theta_z_b=function(lambda,z,myn,Sigma){
	theta_z_b_list=list()
	theta_z_b_1=lambda*t(z)*sqrt(myn)
	theta_z_b_2=-lambda*Sigma*myn/2
	theta_z_b_list[[1]]=theta_z_b_1
	theta_z_b_list[[2]]=theta_z_b_2
	return(theta_z_b_list)
}

######## tested
calc_g_y_b=function(lambda,yty,n){
	a=-n/2*log(2*pi)+n/2*log(lambda)-lambda/2*yty
	return(a)
}

######## tested
calc_theta_b_alpha=function(m,b){
	theta_b_alpha_1=m/2
	theta_b_alpha_2=-sum(b^2)/2
	theta_b_alpha=c(theta_b_alpha_1,theta_b_alpha_2)
	return(theta_b_alpha)
}

######## tested
calc_theta_b_via_Eu=function(Eu_alpha,m){
	theta_b_list=list()
	theta_b_1=matrix(0,m,1)
	theta_b_2= -Eu_alpha[2]*(diag(m))/2
	theta_b_list[[1]]=theta_b_1
	theta_b_list[[2]]=theta_b_2
	return(theta_b_list)
}

######## tested
calc_theta_bi_via_Eu=function(Eu_alphai){
	m=length(Eu_alphai[,1])
	theta_b_list=list()
	theta_b_1=matrix(0,m,1)
	if(m>1){
	theta_b_2=-diag(Eu_alphai[,2])/2
	}else{
		theta_b_2=as.matrix(-Eu_alphai[,2]/2)
	}
	theta_b_list[[1]]=theta_b_1
	theta_b_list[[2]]=theta_b_2
	return(theta_b_list)
}

######## to test
calc_theta_bi_via_Eu_diag=function(Eu_alphai){
	m=length(Eu_alphai[,1])
	theta_b_list=list()
	theta_b_1=matrix(0,m,1)
	theta_b_2=-Eu_alphai[,2]/2
	theta_b_list[["first"]]=theta_b_1
	theta_b_list[["sec_diag"]]=theta_b_2
	return(theta_b_list)
}


######## tested
calc_theta_lambda_via_Eu=function(lambda1,lambda2){
	theta_lambda=c(lambda1-1,-lambda2)
	return(theta_lambda)
}


######## tested
calc_theta_alpha_via_Eu=function(alpha1,alpha2){
	theta_alpha=c(alpha1-1,-alpha2)
	return(theta_alpha)
}

######## tested
calc_theta_y_lambda_via_Eu=function(myn,Eu_b_list,ytX, XtX,yty){
	theta_y_lambda_1=myn/2
	theta_y_lambda_2=-0.5*(yty-2*ytX%*%Eu_b_list[[1]]+sum(XtX*Eu_b_list[[2]]))
	theta_y_lambda=c(theta_y_lambda_1,theta_y_lambda_2)
	return(theta_y_lambda)
}

calc_pseudo_inverse=function(XtX,ndim){
	myeig=eigen(XtX)
	val=myeig$values
	val[c(1:length(val))>ndim]=0
	val[c(1:length(val))<=ndim]=1/val[c(1:length(val))<=ndim]
	pseudoinv=myeig$vectors%*%diag(val)%*%t(myeig$vectors)
	return(pseudoinv)
}

######## not tested
calc_theta_z_lambda_via_Eu=function(myn=NULL,Eu_b_list,ytX, XtX,yty,adjust_yty=FALSE){
	#z=ytX/sqrt(n)
	#Sigma=XtX/(n)
	mym=dim(XtX)[1]
	theta_z_lambda_1=mym/2
	###calculate pseudo-inverse
	ndim=rankMatrix(XtX)
	pseudoinv=calc_pseudo_inverse(XtX,ndim)
	XtX_inv=pseudoinv
	###END: calculate pseudo-inverse###ytX%*%XtX_inv%*%t(ytX)-yty
	if(adjust_yty){
		theta_z_lambda_2=-0.5*(yty-2*ytX%*%Eu_b_list[[1]]+sum(XtX*Eu_b_list[[2]]))
	}else{
		theta_z_lambda_2=-0.5*(ytX%*%XtX_inv%*%t(ytX)-2*ytX%*%Eu_b_list[[1]]+sum(XtX*Eu_b_list[[2]]))
	}
	theta_z_lambda=c(theta_z_lambda_1,theta_z_lambda_2[1])
	return(theta_z_lambda)
}

calc_Eb2_XtX_trace=function(XtX_sqrt,Eu_b_list,missing_variance=0){
#	tr(Eb2%*%XtX)
#	= tr(XtX%*%(solve(-phi_st_b2)/2+Eb1%*%t(Eb1))
#	=tr(XtX%*%Eb1%*%t(Eb1))+tr(XtX%*%(solve(-phi_st_b2)/2)
#
# 	tr(XtX%*%Eb1%*%t(Eb1))=tr(t(Eb1)%*%XtX%*%Eb1)
#	= tr(t(Eb1)%*%XtX_sqrt%*%t(XtX_sqrt)%*%Eb1)
# 	= sum((t(XtX_sqrt)%*%Eb1)^2)
#
#  tr(XtX%*%(solve(-phi_st_b2)/2)
# since -phi_st_b2=diagDprime-cXtX_sqrtAug_interimRes_sqrt%*%t(cXtX_sqrtAug_interimRes_sqrt)
# we have
#  tr(XtX%*%(solve(-phi_st_b2)/2)=tr(diag(XtX)%*%diagDprime/2)
# -sum((XtX_sqrt%*%cXtX_sqrtAug_interimRes_sqrt)^2)/2
#
#  Note: if the missing_variance != 0, then XtX_sqrt is the variance truncated matrix square root. We therefore need to add back in the diagonal. (thats what the if block is for. )


	first=sum((t(XtX_sqrt)%*%Eu_b_list$Eu_b_1)^2)
	second=sum(XtX_sqrt^2*Eu_b_list$theta_M_inv$diagDprime_alphainv_expanded)
	third=sum((t(XtX_sqrt)%*%Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt)^2)
	res=first+0.5*second-0.5*third
	# print("first")
	# print(first)
	# print("second")
	# print(second)
	# print("third")
	# print(third)
	if(missing_variance!=0){
		mym = length(XtX_sqrt[,1])
		missing_frac = missing_variance/mym
		first_missingvar = sum(Eu_b_list$Eu_b_1^2)*missing_frac
		second_missingvar = missing_frac*sum(Eu_b_list$theta_M_inv$diagDprime_alphainv)
		third_missingvar= sum(Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt^2)*missing_frac
		res_missingvar=first_missingvar+0.5*second_missingvar-0.5*third_missingvar
		calc_direct=FALSE
		if(calc_direct==TRUE){
			TT = 0.5*(diag(Eu_b_list$theta_M_inv$diagDprime_alphainv)-Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt%*%t(Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt)) + Eu_b_list$Eu_b_1 %*% t(Eu_b_list$Eu_b_1)
			res_missing_var=missing_frac*sum(diag(TT))
		}
		res=res+res_missingvar
	}
	return(res)
}

######## SPEEDUP
calc_theta_y_lambda_via_Eu_svd=function(myn,Eu_b_list,ytX,XtX_sqrt,yty,decomposed=FALSE,missing_variance=0){
	if(decomposed==FALSE){
		mym=length(XtX_sqrt[,1])
		XtX=(XtX_sqrt)%*%t(XtX_sqrt)
		if(missing_variance!=0){
			diag(XtX)=diag(XtX)+missing_variance/mym
		}
		theta_y_lambda_1=myn/2
		theta_y_lambda_2=-0.5*(yty-2*ytX%*%Eu_b_list[[1]]+sum(XtX*Eu_b_list[[2]]))
	}else{
		theta_y_lambda_1=myn/2
		temp_res=calc_Eb2_XtX_trace(XtX_sqrt,Eu_b_list,missing_variance=missing_variance)
		theta_y_lambda_2=-0.5*(yty-2*ytX%*%Eu_b_list$Eu_b_1+temp_res)
	}
	theta_y_lambda=c(theta_y_lambda_1,theta_y_lambda_2)	
	return(theta_y_lambda)
}

######## SPEEDUP
calc_theta_z_lambda_via_Eu_svd=function(myn,Eu_b_list,ytX,XtX_sqrt,ZtXtX_invZ,decomposed=FALSE,missing_variance=0){
	mym=length(XtX_sqrt[,1])
	if(decomposed==FALSE){
		XtX=(XtX_sqrt)%*%t(XtX_sqrt)
		stop("errndfkgsp:not implemented.")
		if(missing_variance!=0){
			diag(XtX)=diag(XtX)+missing_variance/mym
		}
		theta_y_lambda_1=mym/2
		theta_y_lambda_2=-0.5*(yty-2*ytX%*%Eu_b_list[[1]]+sum(XtX*Eu_b_list[[2]]))
	}else{
		theta_z_lambda_1=mym/2
		temp_res=calc_Eb2_XtX_trace(XtX_sqrt,Eu_b_list,missing_variance=missing_variance)
		theta_z_lambda_2=-0.5*(ZtXtX_invZ-2*ytX%*%Eu_b_list$Eu_b_1+temp_res)
	}
	theta_z_lambda=c(theta_z_lambda_1,theta_z_lambda_2)	
	return(theta_z_lambda)
}

######## tested
calc_theta_y_b_via_Eu=function(Eu_lambda,ytX, XtX){
	theta_y_b_list=list()
	theta_y_b_1=Eu_lambda[[2]]*t(ytX)
	theta_y_b_2=-Eu_lambda[[2]]*XtX/2
	theta_y_b_list[[1]]=theta_y_b_1
	theta_y_b_list[[2]]=theta_y_b_2
	return(theta_y_b_list)
}

######## not tested
calc_theta_z_b_via_Eu=function(Eu_lambda,ytX, XtX){
	theta_z_b_list=list()
	theta_z_b_1=Eu_lambda[[2]]*t(ytX)
	theta_z_b_2=-Eu_lambda[[2]]*XtX/2
	theta_z_b_list[[1]]=theta_z_b_1
	theta_z_b_list[[2]]=theta_z_b_2
	return(theta_z_b_list)
}

######## SPEEDUP
calc_theta_y_b_via_Eu_sqrt_abs=function(Eu_lambda,ytX,XtX_sqrt,missing_variance=0){
	theta_y_b_list=list()
	theta_y_b_1=Eu_lambda[[2]]*t(ytX)
	mym=length(XtX_sqrt[,1])
	theta_y_b_2_sqrt_abs=sqrt(Eu_lambda[[2]]/2)*XtX_sqrt
	theta_y_b_list[["first"]]=theta_y_b_1
	theta_y_b_list[["sec_sqrt_abs"]]=theta_y_b_2_sqrt_abs
	return(theta_y_b_list)
}

######## tested
calc_theta_b_alpha_via_Eu=function(m,Eu_b_list){
	theta_b_alpha_1=m/2
	theta_b_alpha_2=-sum(diag(Eu_b_list[[2]]))/2
	if(dim(Eu_b_list[[2]])[1]==1){
		theta_b_alpha_2=-sum(Eu_b_list[[2]])/2
	}
	theta_b_alpha=c(theta_b_alpha_1,theta_b_alpha_2)
	return(theta_b_alpha)
}



######## tested
calc_theta_bi_alphai_via_Eu=function(m,Eu_b_list,decomposed=FALSE){
	theta_b_alpha_1=rep(1/2,m)
	if(decomposed==FALSE){
		theta_b_alpha_2=-(diag(Eu_b_list[[2]]))/2
		if(dim(Eu_b_list[[2]])[1]==1){
			theta_b_alpha_2=-(Eu_b_list[[2]])/2
		}
	}else{
		theta_b_alpha_2=-(0.5*get_diag_from_inv_decompose(Eu_b_list$theta_M_inv)+(Eu_b_list$Eu_b_1)^2)/2
	}
	theta_b_alpha=cbind(theta_b_alpha_1,theta_b_alpha_2)
	return(theta_b_alpha)
}


######## tested
calc_g_lambda_via_theta=function(theta){
	return((theta[1]+1)*log(-theta[2])-lgamma(theta[1]+1))
}

######## to test
calc_g_lambda2_via_theta=function(theta){
	return((theta[1]+1)*log(-theta[2])-lgamma(theta[1]+1))
}

calc_g_lambda2=function(rho1,rho2){
	myg=rho1*log(rho2)-lgamma(rho1)
	return(myg)
}

######## to test
calc_g_tau2_via_theta=function(theta){
	return((theta[1]+1)*log(-theta[2])-lgamma(theta[1]+1))
}

calc_g_tau2=function(tau1,tau2){
	return(tau1*log(tau2)-lgamma(tau1))
}


######## to test
calc_g_kappa_via_theta=function(theta){
	return((theta[1]+1)*log(-theta[2])-lgamma(theta[1]+1))
}

######## to test
calc_g_kappa_via_Eu_tau2=function(tau1,Eu_tau2){
	return(tau1*Eu_tau2[1]-lgamma(tau1))
}
######## tested
calc_g_alpha_via_theta=function(theta){
	return((theta[1]+1)*log(-theta[2])-lgamma(theta[1]+1))
}

######## tested
calc_g_alphai_via_theta=function(theta){
	return((theta[,1]+1)*log(-theta[,2])-lgamma(theta[,1]+1))
}

######## to test
calc_g_b_via_Eu_alphai=function(Eu_alphai){
	m=dim(Eu_alphai)[1]
	return(-0.5*log(2*pi)*m + sum(0.5*Eu_alphai[,1]))
}


######## tested
calc_g_b_via_thetaOld=function(theta_list){
	m=length(theta_list[[1]])
	theta_M=-2*theta_list[[2]]
	theta_M_inv=calc_fast_cov_solve(theta_M)
	mu=theta_M_inv%*%theta_list[[1]]
	N=-m/2*log(2*pi)
	if(dim(theta_M)[1]==1){
		M=0.5*log(theta_M)
	}else{
		mydet=determinant(theta_M,logarithm=TRUE)
		if(mydet$sign!=1){stop("degenerate")}
		M=0.5*mydet$modulus
	}
	F=-0.5*t(mu)%*%theta_list[[1]]
	res=M+N+F
	return(res)
}

######## tested
calc_g_b_via_theta=function(theta_list){
	m=length(theta_list[[1]])
	N=-m/2*log(2*pi)
	if(dim(theta_list[[2]])[1]==1){
		M=0.5*log(-2*theta_list[[2]])
	}else{
		mydet=determinant((-2*theta_list[[2]]),logarithm=TRUE)
		if(mydet$sign!=1){
			warning("degenerate")
			return(-Inf)
		}
		M=0.5*mydet$modulus
	}
	F=0.25*t(theta_list[[1]])%*%solve(theta_list[[2]],theta_list[[1]])
	res=M+N+F
	return(res)
}

######## to test
calc_g_b_via_theta4eigen=function(theta_list){
	if(!is.null(theta_list[["first"]])){
		m=length(theta_list[["first"]])
	}else{
		m=length(theta_list[[1]])
	}
	theta2_inv=-calculate_woodbury_inversion_sym(diagDprime_alpha=-theta_list[["sec_diag"]],cXtX_sqrt=theta_list[["sec_sqrt_abs"]],only_trace=FALSE)
	N=-m/2*log(2*pi)
	if(dim(theta2_inv)[1]==1){
		sec=(theta_list[["sec_diag"]])-theta_list[["sec_sqrt_abs"]]%*%t(theta_list[["sec_sqrt_abs"]])
		M=0.5*log(-2*sec)
	}else{
		mydet=tryCatch(
			{
				result=(calculate_determinant_sym(diagDprime_alpha=theta_list[["sec_diag"]],XtX_sqrt=theta_list[["sec_sqrt_abs"]],mysign=-1,proport=-2,calc_log=TRUE))
			},error=function(e){
				warning("degenerate")
				return(-Inf)
			},finally = {
			}
		)
		M=0.5*mydet
	}
	F=0.25*sum(theta2_inv*(theta_list[["first"]]%*%t(theta_list[["first"]])))
	res=M+N+F
	return(res)
}


######## to test
calc_g_b_via_theta4eigen_only_diag=function(theta_list){
	if(!is.null(theta_list[["first"]])){
		m=length(theta_list[["first"]])
	}else{
		m=length(theta_list[[1]])
	}
	theta2_inv=1/theta_list[["sec_diag"]]

	N=-m/2*log(2*pi)
	if(length(theta2_inv)==1){
		sec=theta_list[["sec_diag"]]
		M=0.5*log(-2*sec)
	}else{
		mydet=sum(log(-2*theta_list[["sec_diag"]]))
		M=0.5*mydet
	}
	F=0.25*sum(theta2_inv*(theta_list[["first"]]%*%t(theta_list[["first"]])))
	res=M+N+F
	res=res
	return(res)
}

######## tested (indirectly)
calc_y_canonical=function(y){
	return(cbind(y,y^2))
}

######## tested (indirectly)
calc_theta_y_canonical=function(Eu_b_list,Eu_lambda,X){
	theta1=Eu_lambda[[2]]*X%*%Eu_b_list[[1]]
	theta2=-Eu_lambda[[2]]/2
	theta_list=list()
	theta_list[[1]]=theta1
	theta_list[[2]]=theta2
	return(theta_list)
}

######## tested (indirectly)
calc_g_y_canonical=function(Eu_b_list,Eu_lambda,X){
	n=length(X[,1])
	g_ys=rep(0,n)
	for(i in c(1:n)){
		g_ys[i]=0.5*(Eu_lambda[[1]]-sum(Eu_lambda[[2]]*(t(X[i,,drop=F])%*%X[i,,drop=F])*Eu_b_list[[2]])-log(2*pi))
	}
	return(g_ys)
}

######## tested
calc_L_y_slow=function(y,Eu_b_list,Eu_lambda,X){
	ss=calc_theta_y_canonical(Eu_b_list,Eu_lambda,X)
	aa=calc_g_y_canonical(Eu_b_list,Eu_lambda,X)
	bb=rep(0,length(aa))
	for(i in c(1:length(y))){
		bb[i]=sum(calc_y_canonical(y[i])*c(ss[[1]][i],ss[[2]]))
	}
	cc=aa+bb
}


######## to test
calc_L_y_observed=function(yty,ytX,XtX,Eu_b_list,Eu_lambda,n){
	a=-n*log(2*pi)/2+n*Eu_lambda[[1]]/2
	b=-Eu_lambda[[2]]*(yty-2*ytX%*%Eu_b_list[[1]]+sum(Eu_b_list[[2]]*XtX))/2
	ret=a+b
	return(ret)
}

######## to test
calc_L_z_observed=function(ytX,XtX,Eu_b_list,Eu_lambda,myn,yty,adjust_yty=FALSE){
	mym=dim(XtX)[1]
	ndim=rankMatrix(XtX)
	a2=sum(log(eigen(XtX/myn)$values[1:ndim]))
	pseudoinv=calc_pseudo_inverse(XtX,ndim)
	XtX_inv=pseudoinv
	a=-mym*log(2*pi)/2+mym*Eu_lambda[[1]]/2-a2/2
	#b=-Eu_lambda[[2]]*(ytX%*%XtX_inv%*%t(ytX)-2*ytX%*%Eu_b_list[[1]]+sum(Eu_b_list[[2]]*XtX))/2
	if(adjust_yty){
		b=-Eu_lambda[[2]]*(yty-2*ytX%*%Eu_b_list[[1]]+sum(Eu_b_list[[2]]*XtX))/2
	}else{
		b=-Eu_lambda[[2]]*(ytX%*%XtX_inv%*%t(ytX)-2*ytX%*%Eu_b_list[[1]]+sum(Eu_b_list[[2]]*XtX))/2
	}
	ret=a+b
	return(ret)
}

######## to test
calc_L_y_observed4eigen=function(yty,ytX,XtX_sqrt,Eu_b_list,Eu_lambda,n,decomposed=FALSE,missing_variance=0){
	a=-n*log(2*pi)/2+n*Eu_lambda[[1]]/2
	if(decomposed==FALSE){
		XtX=XtX_sqrt%*%t(XtX_sqrt)
		if(missing_variance!=0){
			mym=length(XtX_sqrt[,1])
			diag(XtX)=diag(XtX)+missing_variance/mym
		}
		b=-Eu_lambda[[2]]*(yty-2*ytX%*%Eu_b_list[[1]]+sum(Eu_b_list[[2]]*XtX))/2
	}else{
		mytr=calc_Eb2_XtX_trace(XtX_sqrt,Eu_b_list,missing_variance=missing_variance)
		b=-Eu_lambda[[2]]*(yty-2*ytX%*%Eu_b_list[["Eu_b_1"]]+mytr)/2
	}
	ret=a+b
	return(ret)
}

######## to test
calc_L_z_observed4eigen=function(ZtSigma_invZ,ytX,XtX_sqrt,Eu_b_list,Eu_lambda,n,missing_variance=0){
	mym=dim(XtX_sqrt)[1]
	myd=svd(XtX_sqrt/sqrt(n))$d
	myeigen=myd^2
	missing_var_element_contrib=missing_variance/(n*mym)
	all_eigen=rep(missing_var_element_contrib,mym)
	all_eigen[c(1:length(myd))]=all_eigen[c(1:length(myd))]+myeigen
	a2=sum(log(all_eigen))
	a=-mym*log(2*pi)/2+mym*Eu_lambda[[1]]/2-a2/2
	mytr=calc_Eb2_XtX_trace(XtX_sqrt,Eu_b_list,missing_variance=missing_variance)
	b=-Eu_lambda[[2]]*(ZtSigma_invZ-2*ytX%*%Eu_b_list[["Eu_b_1"]]+mytr)/2
	ret=a+b
	return(ret)
}

######## to test
calc_L_b_observed=function(Eu_b_list,Eu_alphai){
	if(dim(Eu_b_list[[2]])[1]==1){
		ret=-0.5*log(2*pi)+0.5*Eu_alphai[,1]-0.5*Eu_alphai[,2]*Eu_b_list[[2]][1,1]
	}else{
		ret=-0.5*log(2*pi)+0.5*Eu_alphai[,1]-0.5*Eu_alphai[,2]*diag(Eu_b_list[[2]])
	}
	return(ret)
}

######## to test
calc_L_alpha_observed=function(gamma1,Eu_alphai,Eu_ak,mapping_mat){
	m=dim(mapping_mat)[1]
	Eu_ak1_mat=t(kronecker(t(rep(1,m)),Eu_ak[,1]))
	Eu_ak2_mat=t(kronecker(t(rep(1,m)),Eu_ak[,2]))
	a=gamma1*rowSums(mapping_mat*Eu_ak1_mat)
	b=-lgamma(gamma1)+(gamma1-1)*Eu_alphai[,1]
	d=-Eu_alphai[,2]*exp(rowSums(mapping_mat*log(Eu_ak2_mat)))
	ret=a+b+d
	return(ret)
}

calculate_missing_variance=function(eig_values,XtX_sqrt){
	res=sum(eig_values)-sum(XtX_sqrt^2)
	res=max(0,res)
	return(res)
}

expand_Eu_b_list=function(Eu_b_list,k){
	if(is.null(Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt)){
		len=length(Eu_b_list$theta_M_inv$diagDprime_alphainv)
		Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt=matrix(0,len,k)
	}
	if(k!=dim(Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt)[2]){
		stop("errospamc=:dimension error")
	}
	Eu_b_list$theta_M_inv$diagDprime_alphainv_expanded=kronecker(Eu_b_list$theta_M_inv$diagDprime_alphainv,t(rep(1,k)))
	return(Eu_b_list)
}

run_gene_cycle_only_b=function(yty,XtX_sqrt,ytX,N,theta_lambda,mapping_mat,Eu_b_list,Eu_alphai,Eu_lambda,calc_L=TRUE,eig_values,approximation_mode,Eu_lambda2){
	if(is.null(N)){
		stop("N is missing")
	}
	if(length(Eu_b_list)==2){
		M=length(Eu_b_list[[1]])
	}else{
		M=length(Eu_b_list[[3]])
	}
	nodenr=get("nodenr",get_data_container_environment())
	missing_var=calculate_missing_variance(eig_values,XtX_sqrt)
	if(length(Eu_b_list)==2){
		if(approximation_mode){
			stop("erraox,ps:approx not implemented.")
		}
		theta_y_lambda=calc_theta_y_lambda_via_Eu_svd(myn=N,Eu_b_list=Eu_b_list,ytX=ytX, XtX_sqrt=XtX_sqrt,yty=yty,decomposed=FALSE,missing_variance=missing_var)
	}else{
		Eu_b_list=expand_Eu_b_list(Eu_b_list,k=dim(XtX_sqrt)[2])
		if(approximation_mode){
			theta_y_lambda=calc_theta_z_lambda_via_Eu_svd(myn=N,Eu_b_list=Eu_b_list,ytX=ytX,XtX_sqrt,ZtXtX_invZ=yty,decomposed=TRUE,missing_variance=missing_var)
			}else{
		theta_y_lambda=calc_theta_y_lambda_via_Eu_svd(myn=N,Eu_b_list=Eu_b_list,ytX=ytX, XtX_sqrt=XtX_sqrt,yty=yty,decomposed=TRUE,missing_variance=missing_var)
		}
	}
	current_round=get("current_round",get_data_container_environment())
	theta_lambda_star=theta_lambda+theta_y_lambda
	Eu_lambda=calc_Eu_lambda_via_theta(theta_lambda_star)
	if(sum(is.na(Eu_lambda))){
		print("Eu_lambda")
		print(Eu_lambda)
		browser()
		stop("NA in Eu_lambda")
	}
	theta_b_diag=calc_theta_bi_via_Eu_diag(Eu_alphai=Eu_alphai)
##	theta_b_diag[["sec_diag"]]=-D_alpha=-0.5
	theta_y_b_sqrt_abs=calc_theta_y_b_via_Eu_sqrt_abs(Eu_lambda,ytX,XtX_sqrt)
	theta_b_star_list=list()
	theta_b_star_list[["formula"]]="diag(sec_diag)-sec_sqrt_abs%*%t(sec_sqrt_abs)"
	theta_b_star_list[["formula_expanded"]]="diag(theta_b_star_list$sec_diag)-theta_b_star_list$sec_sqrt_abs%*%t(theta_b_star_list$sec_sqrt_abs)"
	theta_b_star_list[["first"]]=theta_b_diag[["first"]]+theta_y_b_sqrt_abs[["first"]]
	theta_b_star_list[["sec_diag"]]=theta_b_diag[["sec_diag"]]
	theta_b_star_list[["sec_sqrt_abs"]]=theta_y_b_sqrt_abs[["sec_sqrt_abs"]]
	##
	Eu_b_list=calc_Eu_b_via_theta4eigen(theta_b_star_list,decomposed=TRUE,missing_variance=missing_var,myc=-Eu_lambda[2]/2)
	theta_b_alphai=calc_theta_bi_alphai_via_Eu(M,Eu_b_list,decomposed=TRUE)
	if(calc_L){
			L_lambda_star_part=-sum((theta_lambda_star)*Eu_lambda)-calc_g_lambda_via_theta(theta_lambda_star)
			if(length(Eu_b_list)==2){
				stop("not aligned anymore with length(Eu_b_list)==4")
				theta_b_star_list[["sec"]]=diag(theta_b_star_list[["sec_diag"]])-theta_b_star_list[["sec_sqrt_abs"]]%*%t(theta_b_star_list[["sec_sqrt_abs"]])
				L_b_star_part=sum(-(theta_b_star_list[["first"]])*(Eu_b_list[[1]]))+sum(-unlist(theta_b_star_list[["sec"]])*unlist(Eu_b_list[[2]]))-calc_g_b_via_theta4eigen(theta_b_star_list)				
				L_obs=calc_L_y_observed4eigen(yty,ytX,XtX_sqrt=XtX_sqrt,Eu_b_list,Eu_lambda,N,decomposed=FALSE,missing_variance=missing_var)
			}else{

				fake_B_list=list()
				fake_B_list[[1]]=Eu_b_list[["Eu_b_1"]]
				fake_B_list[[2]]=0.5*(diag(Eu_b_list$theta_M_inv$diagDprime_alphainv)-Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt%*%t(Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt)) + Eu_b_list$Eu_b_1 %*% t(Eu_b_list$Eu_b_1)
				fake_theta_list=list()
				fake_theta_list[[1]]=theta_b_star_list[["first"]]
				fake_theta_list[[2]]=diag(theta_b_star_list$sec_diag)-theta_b_star_list$sec_sqrt_abs%*%t(theta_b_star_list$sec_sqrt_abs)
				if(missing_var!=0){
					mym=length(fake_theta_list[[1]])
					myadder=(missing_var/mym)*Eu_lambda[2]/2
					diag(fake_theta_list[[2]])=diag(fake_theta_list[[2]])-myadder
					theta_b_star_list[["sec_diag"]]=theta_b_star_list[["sec_diag"]]-myadder
				}
				L_b_star_part=sum(-unlist(fake_theta_list)*unlist(fake_B_list))-calc_g_b_via_theta4eigen(theta_b_star_list)
				mym=dim(fake_B_list[[2]])[1]
				if(mym==1){
					stop("nsnps has to be larger than 1")
				}
				L_b_part=sum(unlist(theta_b_diag[["sec_diag"]])*diag(fake_B_list[[2]]))+ calc_g_b_via_Eu_alphai(Eu_alphai)
				L_b=L_b_part+L_b_star_part
				if(approximation_mode){
					L_obs=calc_L_z_observed4eigen(ZtSigma_invZ=yty,ytX,XtX_sqrt=XtX_sqrt,Eu_b_list,Eu_lambda,N,missing_variance=missing_var)
				}else{
					L_obs=calc_L_y_observed4eigen(yty,ytX,XtX_sqrt=XtX_sqrt,Eu_b_list,Eu_lambda,N,decomposed=TRUE,missing_variance=missing_var)
				}
			}
	}
	Eu_b_list$theta_M_inv$diagDprime_alphainv_expanded=NULL
############
	out_list=list(
		E_list=list(
			Eu_b_list=list(
				Eu_b="",
				Eu_b2=""
			),
			Eu_lambda=""
		),
		L_list=list(
			L_lambda_star_part="",
			L_b_star_part="",
			L_b="",
			L_y_obs=""
		),
		theta_b_alphai="",
		theta_b_star_list=""
	)
	out_list$E_list$Eu_b_list=Eu_b_list
	out_list$E_list$Eu_lambda=Eu_lambda
	out_list$theta_b_alphai=theta_b_alphai
	if(calc_L){	
		out_list$L_list$L_lambda_star_part=L_lambda_star_part
		out_list$L_list$L_b_star_part=L_b_star_part
		out_list$L_list$L_b=L_b
		out_list$L_list$L_y_obs=L_obs
	}
	return(out_list)
}


### runs on node.
calculate_ZtSigma_invZ_list=function(){
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	for(i in c(1:length(gene_dat_list))){
		XtX_sqrt=gene_dat_list[[i]]$Obs_list$XtX_sqrt
		eig_values=gene_dat_list[[i]]$Obs_list$eig_values
		ytX=gene_dat_list[[i]]$Obs_list$ytX
		mym=dim(XtX_sqrt)[1]
		myn=dim(XtX_sqrt)[1]
		missing_var=calculate_missing_variance(eig_values,XtX_sqrt)
		mydiag=rep(missing_var/mym,mym)
		if(sum(mydiag)<1E-10){
			print(paste("for gene",names(gene_dat_list)[i]))
			print("no variance removed. weak manual regularization necessary.")
			#mydiag=rep(missing_var/mym,mym)
			tot2add=1E-9*mym
			mydiag=mydiag+rep(tot2add,mym)/mym
			if(length(eig_values)>length(mydiag)){
				eig_values[1:length(mydiag)]=eig_values[1:length(mydiag)]+rep(tot2add,mym)/mym
			}else{
				eig_values=eig_values+tot2add/myn
			}
			gene_dat_list[[i]]$Obs_list$eig_values=eig_values
		}
		myinv=calculate_woodbury_inversion_sym(diagDprime_alpha=mydiag,cXtX_sqrt=XtX_sqrt,only_trace=FALSE,decomposed=TRUE)
		Xw=myinv$cXtX_sqrtAug_interimRes_sqrt
		tmp=ytX%*%Xw
		a=tmp%*%t(tmp)
		b=sum(myinv$diagDprime_alphainv*(ytX^2))
		if(mym > myn*1.2){
			browser
		}
		myret=b-a
		gene_dat_list[[i]]$Obs_list$yty=myret
	}
	assign("gene_dat_list",gene_dat_list,envir=get_data_container_environment())
}

check_yty_on_node_list=function(){
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	for(i in c(1:length(gene_dat_list))){
		XtX_sqrt=gene_dat_list[[i]]$Obs_list$XtX_sqrt
		eig_values=gene_dat_list[[i]]$Obs_list$eig_values
		ytX=gene_dat_list[[i]]$Obs_list$ytX
		yty=gene_dat_list[[i]]$Obs_list$yty
		mym=dim(XtX_sqrt)[1]
		missing_var=calculate_missing_variance(eig_values,XtX_sqrt)
		mydiag=rep(missing_var/mym,mym)
		if(missing_var>1E-9){
			myinv=calculate_woodbury_inversion_sym(diagDprime_alpha=mydiag,cXtX_sqrt=XtX_sqrt,only_trace=FALSE,decomposed=TRUE)
			Xw=myinv$cXtX_sqrtAug_interimRes_sqrt
			tmp=ytX%*%Xw
			a=tmp%*%t(tmp)
			b=sum(myinv$diagDprime_alphainv*(ytX^2))
			myret=b-a
		}else{
			mysvd=svd(XtX_sqrt)
			tmp=ytX%*%mysvd$u
			tmp=tmp*(1/mysvd$d)
			myret=sum(tmp^2)
		}		
		gene_dat_list[[i]][["Obs_list"]][["ytyold"]]=yty
		gene_dat_list[[i]][["Obs_list"]][["ZtSigma_invZ"]]=myret
		gene_dat_list[[i]][["Obs_list"]][["yty"]]=max(myret,yty)
	}
	assign("gene_dat_list",gene_dat_list,envir=get_data_container_environment())
}

#runs on main with 
calc_ZtSigma_invZ=function(gene_dat_list){
	all_ZtSigma_invZ=rep(0,length(gene_dat_list))
	for(i in c(1:length(gene_dat_list))){
		XtX_sqrt=gene_dat_list[[i]]$Obs_list$XtX_sqrt
		eig_values=gene_dat_list[[i]]$Obs_list$eig_values
		ytX=gene_dat_list[[i]]$Obs_list$ytX
		mym=dim(XtX_sqrt)[1]
		missing_var=calculate_missing_variance(eig_values,XtX_sqrt)
		mydiag=rep(missing_var/mym,mym)
		if(missing_var>1E-9){
			myinv=calculate_woodbury_inversion_sym(diagDprime_alpha=mydiag,cXtX_sqrt=XtX_sqrt,only_trace=FALSE,decomposed=TRUE)
			Xw=myinv$cXtX_sqrtAug_interimRes_sqrt
			tmp=ytX%*%Xw
			a=tmp%*%t(tmp)
			b=sum(myinv$diagDprime_alphainv*(ytX^2))
			myret=b-a
		}else{
			mysvd=svd(XtX_sqrt)
			tmp=ytX%*%mysvd$u
			tmp=tmp*(1/mysvd$d)
			myret=sum(tmp^2)
		}
		all_ZtSigma_invZ[i]=myret
	}
	return(all_ZtSigma_invZ)
}

process_Ls_on_node_list=function(Eu_alphai){
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	L_btot=0
	L_alphaitot=0
	aa=20
#	print("tmptmofw")
	Ls_list=list()
	for(i in c(1:length(gene_dat_list))){
		Eu_b_list=gene_dat_list[[i]]$E_list$Eu_b_list
		if(length(Eu_b_list)==2){
			theta_b=calc_theta_bi_via_Eu(Eu_alphai=Eu_alphai[[i]])
			L_b=sum(unlist(theta_b)*unlist(Eu_b_list))+sum(0.5*Eu_alphai[[i]][,1]-0.5*log(2*pi))
		}else{
			L_b=gene_dat_list[[i]]$L_list$L_b

		}
		L_btot=L_btot+L_b
	}
	L_bs=L_btot
	L_y_obs=sum(unlist(lapply(c(1:length(gene_dat_list)),function(i){sum(unlist(gene_dat_list[[i]]$L_list$L_y_obs))})))
	L_lambda_star_part=sum(unlist(lapply(c(1:length(gene_dat_list)),function(i){sum(unlist(gene_dat_list[[i]]$L_list$L_lambda_star_part))})))
	Ls_list=list(L_bs=L_bs,L_y_obs=L_y_obs,L_lambda_star_part=L_lambda_star_part)
	return(Ls_list)
}


get_L_y_obs_from_node=function(){
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	L_y_obs_tot=sum(unlist(lapply(c(1:length(gene_dat_list)),function(i){sum(unlist(gene_dat_list[[i]]$L_list$L_y_obs))})))
	return(L_y_obs_tot)
}

get_L_lambda_from_node=function(){
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	L_lambda_tot=sum(unlist(lapply(c(1:length(gene_dat_list)),function(i){sum(unlist(gene_dat_list[[i]]$L_list$L_lambda))})))
	return(L_lambda_tot)
}

run_gene_cycle_only_alpha=function(gene_list,yty,XtX,ytX,N,theta_lambda,theta_kappa,mapping_mat,Eu_b_list,Eu_ak,Eu_alphai,Eu_lambda,Eu_kappa,gamma1,theta_b_alphai,theta_alphai_ak,calc_L=TRUE,calc_B=TRUE,k_minus1,Eu_ak_old,Eu_gammai_old){
	K=dim(Eu_ak)[1]
	M=length(Eu_b_list[[1]])
	#Eu_gammai=calc_Eu_gammai_via_Eu_ak_only_second(Eu_ak,mapping_mat)# 0.7
	if(k_minus1>0){
		Eu_gammai=update_Eu_gammai_second_row_via_Eu_ak(Eu_gammai_old,Eu_ak_old, Eu_ak,mapping_mat,k_minus1)
		}else{
			Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat)
		}
	theta_alphai=calc_theta_alphai_via_Eu_gammai_Eu_kappa(gamma1,Eu_gammai,Eu_kappa,M)#0.08
#	print("tt")
	if(k_minus1==2){
	}
	theta_alphai_star=theta_alphai+theta_b_alphai		
	Eu_alphai=calc_Eu_alphai_via_theta_vectorized(theta_alphai_star)#0.1
	theta_alphai_kappa=calc_theta_alphai_kappa_via_Eu_gammai_Eu_alphai(gamma1,Eu_gammai,Eu_alphai,M)# 0.14
	theta_kappa_star=theta_kappa+theta_alphai_kappa
	theta_alphai_ak[k_minus1+1,]=calc_theta_alpha_ak_via_Eugammai(Eu_gammai,gamma1,Eu_alphai,Eu_kappa,Eu_ak,mapping_mat,k=(k_minus1+1))
	out_list=list(
			E_list=list(
				Eu_ak="",
				Eu_alphai="",
				Eu_gammai=""
				),
			theta_alphai_ak="",
			theta_alphai_star="",
			theta_b_alphai="")
	out_list$E_list$Eu_alphai=Eu_alphai
	out_list$E_list$Eu_gammai=Eu_gammai
	out_list$theta_alphai_ak=theta_alphai_ak
	out_list$theta_alphai_star=theta_alphai_star
	out_list$theta_b_alphai=theta_b_alphai	
	return(out_list)
}


run_gene_cycle_wrapper_only_b=function(gene_list=NULL,theta_lambda=NULL,Eu_lambda2=NULL,calc_L=NULL,approximation_mode=NULL){
	Obs_list=gene_list$Obs_list
	E_list=gene_list$E_list
	out_list=run_gene_cycle_only_b(yty=Obs_list$yty,XtX_sqrt=Obs_list$XtX_sqrt,ytX=Obs_list$ytX,N=Obs_list$n,theta_lambda=theta_lambda,mapping_mat=Obs_list$mapping_mat,Eu_b_list=E_list$Eu_b_list,Eu_alphai=E_list$Eu_alphai,Eu_lambda=E_list$Eu_lambda,calc_L=calc_L,eig_values=Obs_list$eig_values,approximation_mode=approximation_mode,Eu_lambda2=Eu_lambda2)
	return(out_list)
}

run_gene_cycle_wrapper=function(gene_list=NULL,gamma1=NULL,n=NULL,theta_lambda=NULL,theta_kappa=NULL,Eu_ak=NULL,calc_L=TRUE,calc_B=TRUE,k_minus1=NULL,Eu_ak_old=NULL,calc_Kappa=TRUE){
	Obs_list=gene_list$Obs_list
	E_list=gene_list$E_list
	if(calc_B==TRUE || calc_Kappa==TRUE){
		out_list=run_gene_cycle(yty=Obs_list$yty,XtX=Obs_list$XtX,ytX=Obs_list$ytX,N=n,theta_lambda=theta_lambda,theta_kappa=theta_kappa,mapping_mat=Obs_list$mapping_mat,Eu_b_list=E_list$Eu_b_list,Eu_ak=Eu_ak,Eu_alphai=E_list$Eu_alphai,Eu_lambda=E_list$Eu_lambda,Eu_kappa=E_list$Eu_kappa,gamma1=gamma1,calc_L=calc_L,calc_B=calc_B,calc_Kappa=calc_Kappa)
	}else{
		out_list=run_gene_cycle_only_alpha(gene_list=gene_list,yty=Obs_list$yty,XtX=Obs_list$XtX,ytX=Obs_list$ytX,N=n,theta_lambda=theta_lambda,theta_kappa=theta_kappa,mapping_mat=Obs_list$mapping_mat,Eu_b_list=E_list$Eu_b_list,Eu_ak=Eu_ak,Eu_alphai=E_list$Eu_alphai,Eu_lambda=E_list$Eu_lambda,Eu_kappa=E_list$Eu_kappa,gamma1=gamma1,theta_b_alphai=gene_list$theta_b_alphai,theta_alphai_ak=gene_list$theta_alphai_ak,calc_L=calc_L,calc_B=calc_B,k_minus1=k_minus1,Eu_ak_old=Eu_ak_old,Eu_gammai=E_list$Eu_gammai)
	}
	return(out_list)
}



### this function is used for memory purposes:: if an element already exists in a list, then copying happens in place.
### important the output still has two refs so outside of the function we need to call 
### out=make_empty_list(c("aa","bb"))
### print(refs(out))
### out=lapply(out,function(x){x})
### print(refs(out))

make_empty_list=function(names){
	empty_list=list()
	for(name in names){
		 empty_list[[name]]=""
	}
	return(empty_list)
}

#### asigns gene_dat_list into the environment into which variational_aux was sourced
load_data_to_cluster_node=function(gene_dat_list_split){
#	assign("gene_dat_list",gene_dat_list_split,envir=get_data_container_environment())
	print("kkom")
	assign("gene_dat_list",gene_dat_list_split,envir=get_data_container_environment())
#	myenv=parent.
	return(TRUE)
}

load_data_from_cluster_node=function(){
	mygenedat=get("gene_dat_list",envir=get_data_container_environment())
	return(mygenedat)
}


to_be_run_on_node_only_b_part=function(){
	### Assumes that EU_alphai was updated
	current_round=get("current_round",get_data_container_environment())
	calc_L=get("calc_L",get_data_container_environment())
	calc_L_rounds=get("calc_L_rounds",get_data_container_environment())
	calc_B_rounds=get("calc_B_rounds",get_data_container_environment())
	gamma1=get("gamma1",get_data_container_environment())
	nodenr=get("nodenr",get_data_container_environment())
	theta_lambda=get("theta_lambda",get_data_container_environment())
	Eu_lambda2=get("Eu_lambda2",get_data_container_environment())
	approximation_mode=get("approximation_mode",get_data_container_environment())
	############
	############
	if(current_round>1){
		current_round=current_round+calc_B_rounds
	}
	if(current_round==1 && calc_B_rounds==1){
		current_round=2
	}
	if(current_round==1 && calc_B_rounds>1){
		current_round=calc_B_rounds
	}
	if(current_round==0){
		current_round=1
	}
		
	############
	############
	if((current_round %% calc_L_rounds)==0 || current_round==1){
		local_calc_L=TRUE
	}else{
		local_calc_L=FALSE
	}
	if(calc_L==FALSE){
		local_calc_L=FALSE
	}
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	if(approximation_mode){
		theta_index=match(names(gene_dat_list),rownames(theta_lambda))
		out_list=lapply(c(1:length(gene_dat_list)),function(i){			
			run_gene_cycle_wrapper_only_b(gene_dat_list[[i]],theta_lambda=theta_lambda[theta_index[i],],Eu_lambda2=Eu_lambda2,calc_L=local_calc_L,approximation_mode=approximation_mode)})
		names(out_list)=names(gene_dat_list)
	}else{
		out_list=lapply(gene_dat_list,run_gene_cycle_wrapper_only_b,theta_lambda=theta_lambda,Eu_lambda2=Eu_lambda2,calc_L=local_calc_L,approximation_mode=approximation_mode)	
	}
	for(myname in names(out_list)){
		gene_dat_list[[myname]]$E_list$Eu_b_list=out_list[[myname]]$E_list$Eu_b_list
		gene_dat_list[[myname]]$E_list$Eu_lambda=out_list[[myname]]$E_list$Eu_lambda
		gene_dat_list[[myname]]$theta_b_alphai=out_list[[myname]]$theta_b_alphai
		if(local_calc_L){	
			gene_dat_list[[myname]]$L_list$L_lambda_star_part=out_list[[myname]]$L_list$L_lambda_star_part
			gene_dat_list[[myname]]$L_list$L_b_star_part=out_list[[myname]]$L_list$L_b_star_part
			gene_dat_list[[myname]]$L_list$L_b=out_list[[myname]]$L_list$L_b
			gene_dat_list[[myname]]$L_list$L_y_obs=out_list[[myname]]$L_list$L_y_obs
		}		
	}
	all_theta_b_alphai_second=lapply(gene_dat_list,function(x){return(x$theta_b_alphai[,2,drop=F])})
	all_theta_b_alphai_vect=get_all_theta_b_alphai_vect(all_theta_b_alphai_second)
	assign("all_theta_b_alphai_vect",all_theta_b_alphai_vect,envir=get_data_container_environment())
	assign("current_round",current_round,envir=get_data_container_environment())
	assign("gene_dat_list",gene_dat_list,envir=get_data_container_environment())

	########################################
	## pull Eu_lambda_vect as an attribute
	Eu_lambda_list=lapply(names(out_list),function(myname){
		gene_dat_list[[myname]]$E_list$Eu_lambda}
		)
	Eu_lambda_vect=my_docall_rbind(Eu_lambda_list,nsplits=10)
	attr(all_theta_b_alphai_vect,"Eu_lambda_vect")=Eu_lambda_vect
	########################################	
	return(all_theta_b_alphai_vect)
}



#### gets all gene_dat_list$theta_alphai_ak from the environment into which variational_aux was sourced
get_all_theta_alphai_ak=function(){
	out=lapply(get("gene_dat_list",envir=get_data_container_environment()),function(x){return(x$theta_alphai_ak)})
	return(out)
}

### gets all gene_dat_list$theta_alphai_ak from the environment into which variational_aux was sourced
get_all_theta_alphai_ak_k=function(k){
	out=lapply(get("gene_dat_list",envir=get_data_container_environment()),function(x){return(x$theta_alphai_ak[k,,drop=F])})
	return(out)
}

### gets all gene_dat_list$theta_b from the environment into which variational_aux was sourced
get_all_theta_b_alphai_second=function(){
	out=lapply(get("gene_dat_list",envir=get_data_container_environment()),function(x){return(x$theta_b_alphai[,2,drop=F])})
	return(out)
}

set_all_theta_alphai_star_second=function(){
	mylen=length(gene_dat_list)
	get("gene_dat_list",envir=get_data_container_environment())
	if(mylen!=length(all_theta_b_alphai_second)){
		stop("lengths do not correspond.")
	}
	for(i in c(1:mylen)){
		gene_dat_list[[i]]$theta_alphai_star=all_theta_b_alphai_second[[i]]
	}
	assign("gene_dat_list",envir=get_data_container_environment())
}

get_all_theta_b_alphai_gene_indices=function(all_theta_b_alphai_second){
	mylist=all_theta_b_alphai_second
	a=unlist(lapply(mylist,length))
	has_list_el_of_length=sum(a==0)>0
	if(has_list_el_of_length){
		stop("there seems to be a gene without SNPs. remove gene prior to running algorithm.")
	}
	out=unlist(lapply(c(1:length(mylist)),function(x){rep(x,length(mylist[[x]]))}))
	return(out)
}

get_all_theta_b_alphai_vect=function(all_theta_b_alphai_second){
	out=unlist(all_theta_b_alphai_second)
	return(out)
}
	
##TODO: rename
get_gene_indices_vec_boundaries=function(gene_indices){
	uniq_gene_inds=unique(gene_indices)
	res_min=unlist(lapply(uniq_gene_inds,function(x){min(which(gene_indices==x))}))
	res_max=unlist(lapply(uniq_gene_inds,function(x){max(which(gene_indices==x))}))
	gene_indices_vec_boundaries=cbind(res_min,res_max)
	return(gene_indices_vec_boundaries)
}

get_all_theta_alphai_star_from_vect_form=function(vect,gene_indices){
	uniq_gene_inds=unique(gene_indices)
	res=lapply(uniq_gene_inds,function(x){min(which(gene_indices==x))})
	return(res)
}

load_alphai_star_second_to_cluster_node=function(alphai_star_second_vect){
	assign("alphai_star_second_vect",alphai_star_second_vect,envir=get_data_container_environment())		
	myenv=get_data_container_environment()
	return(myenv)
}

all_theta_b_alphai_gene_indices_to_cluster_node=function(all_theta_b_alphai_gene_indices){
	assign("all_theta_b_alphai_gene_indices",all_theta_b_alphai_gene_indices,envir=get_data_container_environment())	
	return(TRUE)
}

load_Eu_alphai_second_list_to_cluster_node=function(Eu_alphai_second_vect){
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	gene_indices_vec_boundaries=get("all_theta_b_alphai_gene_indices",get_data_container_environment())
	mylen=length(gene_indices_vec_boundaries[,1])
	for(i in c(1:mylen)){
		down=gene_indices_vec_boundaries[i,1]
		upper=gene_indices_vec_boundaries[i,2]
		gene_dat_list[[i]]$E_list$Eu_alphai[,2]=Eu_alphai_second_vect[c(down:upper)]
	}
	assign("gene_dat_list",gene_dat_list,get_data_container_environment())
	return(NULL)
}

load_Eu_alphai_second_to_cluster_node=function(Eu_alphai_second_vect){
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	gene_indices_vec_boundaries=get("all_theta_b_alphai_gene_indices",get_data_container_environment())
	mylen=length(gene_indices_vec_boundaries[,1])
	for(i in c(1:mylen)){
		down=gene_indices_vec_boundaries[i,1]
		upper=gene_indices_vec_boundaries[i,2]
		gene_dat_list[[i]]$E_list$Eu_alphai[,2]=Eu_alphai_second_vect[c(down:upper)]
	}
	assign("gene_dat_list",gene_dat_list,get_data_container_environment())
	return(NULL)
}

load_theta_b_alphai_vect_from_cluster_node=function(){
	theta_b_alphai_vect=get("theta_b_alphai_vect",envir=get_data_container_environment())
	return(theta_b_alphai_vect)
}

calculate_L_ak=function(Eu_ak,theta_ak,theta_ak_star){
	L_ak=sum((theta_ak-theta_ak_star)*Eu_ak)+sum(calc_g_ak_via_theta(theta_ak))-sum(calc_g_ak_via_theta(theta_ak_star))
	return(L_ak)
}

get_all_Ls_from_node=function(){
	Lgs=sum(unlist(lapply(get("gene_dat_list",envir=get_data_container_environment()),function(x){sum(unlist(x$L_list))})))
	return(Lgs)
}


calculate_L_minus_ak_alpha_split=function(gene_dat_list,gamma1,Eu_ak,theta_alphai_star,mapping_mat,Eu_alphai,Eu_kappa,theta_alphai_first_row){
	L_alphai=calc_L_alpha(gamma1=gamma1,Eu_ak=Eu_ak,theta_alphai_star=theta_alphai_star,mapping_mat=mapping_mat,Eu_alphai=Eu_alphai,Eu_kappa=Eu_kappa,theta_alphai_first_row=theta_alphai_first_row)
	Lgs=sum(unlist(lapply(gene_dat_list,function(x){sum(unlist(x$L_list))})))
	Ltot_minus_ak=Lgs+sum(L_alphai)
	return(Ltot_minus_ak)
}

calc_L_alpha=function(gamma1,Eu_ak,theta_alphai_star,mapping_mat,Eu_alphai,Eu_kappa,theta_alphai_first_row){
	Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat)
	E_g_alphai_par=calc_g_alphai_via_Eu_gammai_Eu_kappa(gamma1,Eu_gammai,Eu_kappa)
	theta_alphai_vect_second=calc_theta_alphai_via_Eu_gammai_Eu_kappa_second(Eu_gammai,Eu_kappa)
	theta_alphai=cbind(theta_alphai_first_row,theta_alphai_vect_second)
	L_alphai=sum(rowSums((theta_alphai-theta_alphai_star)*Eu_alphai)+E_g_alphai_par-calc_g_alphai_via_theta(theta_alphai_star))
	return(L_alphai)
}
	
calc_L_lambda=function(lambda1,Eu_lambda2,approximation_mode,Eu_lambda_vect,Ls_from_node_dt){
	theta_lambda=calc_theta_lambda(lambda1,Eu_lambda2[2],genespecific_lambda1=approximation_mode)
	#Eu_lambda_vect
	g_lambda_Epar=calc_g_lambda_via_Eu_lambda2(lambda1,Eu_lambda2)
	theta_lambda=calc_theta_lambda(lambda1,Eu_lambda2[2],genespecific_lambda1=approximation_mode)
	if(approximation_mode){
		g_lambda_Epar=sum(g_lambda_Epar)
	}else{
		g_lambda_Epar=g_lambda_Epar*dim(Eu_lambda_vect)[1]
		theta_lambda=rep(1,dim(Eu_lambda_vect)[1])%*%t(theta_lambda)
	}
	Ltot_lambda_Epar=sum(theta_lambda*Eu_lambda_vect)+g_lambda_Epar
	Ltot_lambda=Ltot_lambda_Epar+Ls_from_node_dt[,L_Qlambda_part]
	return(Ltot_lambda)
}

calc_L_b=function(Eu_alphai_vect,Eu_omega_list,mapping_shift_mat_expanded,Eu_nu_list,mapping_mat_expanded_nu,Ls_from_node_dt,all_Eu_b1,all_Eu_b2_diag){
	theta_bi_list_diag=calc_theta_bi_via_Eu_wshift_nu_diag(Eu_alphai=Eu_alphai_vect,Eu_omega_list=Eu_omega_list,mapping_shift_mat=mapping_shift_mat_expanded,Eu_nu_list=Eu_nu_list,mapping_mat_expanded_nu,add_intercept=TRUE,mapping_mat_nu_sorted=TRUE)
	g_b_Epar=calc_g_b_via_Eu_wshift_nu_fast(Eu_alphai=Eu_alphai_vect,Eu_omega_list=Eu_omega_list,V=mapping_shift_mat_expanded,F=mapping_mat_expanded_nu,Eu_nu_list=Eu_nu_list,add_intercept=TRUE,mapping_mat_nu_sorted=TRUE)
	Ltot_b_Epar=sum(theta_bi_list_diag[["first"]]*all_Eu_b1)+sum(theta_bi_list_diag[["sec_diag"]]*all_Eu_b2_diag)+g_b_Epar
		Ltot_b=Ltot_b_Epar+Ls_from_node_dt[,L_Qb_part_tot]
	return(Ltot_b)
}

run_cycle_across_allgenes_mpi_only_b_part=function(mpicl=NULL){
	if(is.null(mpicl)){
		stop("no mpi cluster provided.")
	}
	all_theta_b_alphai_vect_list=clusterCall_or_envCall(mpicl, to_be_run_on_node_only_b_part)
	return(all_theta_b_alphai_vect_list)
}

has_unique_names=function(mylist){
	out=TRUE
	mynames=names(mylist)
	if(is.null(mynames)){
		out=FALSE
	}else{
		if(length(unique(mynames))!=length(mynames)){
			out=FALSE
		}
	}
	return(out)
}

has_unique_names=function(mylist){
	out=TRUE
	mynames=names(mylist)
	if(is.null(mynames)){
		out=FALSE
	}else{
		if(length(unique(mynames))!=length(mynames)){
			out=FALSE
		}
	}
	return(out)
}

get_nsnps_per_gene=function(gene_list){
	len=length(gene_list$Obs_list$ytX)
	return(len)
}

### splits the data in gene_dat_list such that the gene size distribution is balanced between different splits.
split_gene_dat_list=function(mpicl,gene_dat_list,balanceManually=FALSE,verbose=TRUE){
	nsnps_vec=unlist(lapply(gene_dat_list,get_nsnps_per_gene))
	ngenes=length(nsnps_vec)
	if(sum(is.na(nsnps_vec))>0 || sum(nsnps_vec==0)>0){
		stop("some genes don't have any snps attached, abort")
	}
	if(balanceManually==TRUE){
		if(!has_unique_names(gene_dat_list)){
			stop("gene_dat_list needs unique identifying names to do manual split")
		}
		##prepare_empty gene_dat_list
		gene_dat_list_splits=list()
		if(class(mpicl)[1]=="SOCKcluster"){
			cl_len=length(mpicl)
		}else{
			if(class(mpicl)[1]=="environment"){
				cl_len=1
			}else{
				stop("neither cluster nor envir.")
			}
		}
		for(i in c(1:cl_len)){
			gene_dat_list_splits[[i]]=list()			
		}
		genenames=names(gene_dat_list)
		geneindices_descending_by_genesize=order(nsnps_vec,decreasing=TRUE)
		counter=0
		while(counter<ngenes){
			print("genes processed")
			print(counter)
			limit=min(counter+cl_len,length(nsnps_vec))
			indices_to_assign=c((counter+1):limit)
			if(length(indices_to_assign)>1){
				indices_to_assign_shuffled=sample(indices_to_assign,length(indices_to_assign))# this is done so that the first node in the list always gets the biggest genes
			}else{
				indices_to_assign_shuffled=indices_to_assign
			}
			for(i in c(1:length(indices_to_assign_shuffled))){
				current_geneind=geneindices_descending_by_genesize[indices_to_assign_shuffled[i]]
				gene_dat_list_splits[[i]][[genenames[current_geneind]]]=genenames[current_geneind]
			}
			counter=limit
		}
		## currently only gene names are saved in gene_dat_list_splits (to avoid R copynightmare)
		copy_over_for_split=function(split){
			lapply(split, function(myname){split[[myname]]=gene_dat_list[[myname]]})
		}
		gene_dat_list_splits=lapply(gene_dat_list_splits,copy_over_for_split)
		return(gene_dat_list_splits)
	}else{
		gene_dat_list_splits=clusterSplit(mpicl,gene_dat_list)
		return(gene_dat_list_splits)
	}				
}

get_alphasplits=function(mat_to_split,index_range_mappings){
	split_mat=list()
	for(i in c(1:length(index_range_mappings))){
		if(is.null(index_range_mappings[[i]])){
			next
		}
		cur_len=dim(index_range_mappings[[i]])[1]
		if(cur_len==0){
			next
		}
		cur_split=list()
		for(j in c(1:cur_len)){
			cur_index_range=index_range_mappings[[i]][j,]
			lower=cur_index_range[1]
			upper=cur_index_range[2]
			cur_split[[j]]=mat_to_split[c(lower:upper),,drop=FALSE]			
		}
		split_mat[[i]]=cur_split
	}
	return(split_mat)
}

pull_out_dat_list_for_alphasplit=function(mpicl,Eu_ak,Eu_lambda2,Eu_tau2,Ltot,all_Ls,Eu_ak_trace,Eu_kappa_vect,Eu_alphai_vect,index_range_mappings,gene_names,Eu_omega_list=NULL,Eu_delta=NULL,Eu_upsilon_list=NULL,Eu_omega1_trace=NULL,Eu_nu_list=NULL,Eu_nu1_trace=NULL,Eu_chi2j=NULL){
	gene_dat_list_splits=clusterCall_or_envCall(mpicl,load_data_from_cluster_node)
	for(i in c(1:length(gene_dat_list_splits))){
		if(is.null(gene_dat_list_splits[[i]])){
			next
		}
		cur_len=length(gene_dat_list_splits[[i]])
		if(cur_len==0){
			next
		}
		for(j in c(1:cur_len)){
			cur_index_range=index_range_mappings[[i]][j,]
			lower=cur_index_range[1]
			upper=cur_index_range[2]
			gene_dat_list_splits[[i]][[j]]$E_list$Eu_kappa=Eu_kappa_vect[c(lower:upper),,drop=FALSE][1,]
			gene_dat_list_splits[[i]][[j]]$E_list$Eu_alphai=Eu_alphai_vect[c(lower:upper),,drop=FALSE]
		}
	}
	gene_dat_list=unlist(gene_dat_list_splits,recursive=FALSE)
	unsorted_names=names(gene_dat_list)
	matcher=match(gene_names,unsorted_names)
	gene_dat_list=gene_dat_list[matcher]
	out_dat_list=list(Eu_ak="",Ltot="",gene_dat_list="",Eu_ak_trace=Eu_ak_trace)
	out_dat_list$Eu_ak=Eu_ak
	out_dat_list$Eu_lambda2=Eu_lambda2
	out_dat_list$Eu_tau2=Eu_tau2
	out_dat_list$Eu_ak_trace=Eu_ak_trace
	out_dat_list$Ltot=Ltot
	out_dat_list$all_Ls=all_Ls
	out_dat_list$gene_dat_list=gene_dat_list
	out_dat_list$Eu_omega_list=Eu_omega_list
	out_dat_list$Eu_delta=Eu_delta
	out_dat_list$Eu_upsilon_list=Eu_upsilon_list
	out_dat_list$Eu_omega1_trace=Eu_omega1_trace
	out_dat_list$Eu_nu_list=Eu_nu_list
	out_dat_list$Eu_nu1_trace=Eu_nu1_trace
	out_dat_list$Eu_chi2j=Eu_chi2j
	return(out_dat_list)
}

load_node_nr_to_cluster=function(i){
#	assign("nodenr",i,get_data_container_environment())
#	assign("nodenr",i,get_data_container_environment())
	assign("nodenr",i,get_data_container_environment())
	return()
}


check_core_nrs=function(ncores){
	if(ncores<1){
		stop("no reasonable core nr specified.")
	}
}

check_round_indices=function(nrounds,calc_L_rounds,calc_B_rounds){
	if(calc_L_rounds %% calc_B_rounds!=0){
		stop("calc_L_rounds not a multiple of calc_B_rounds. This is not allowed")
	}
	if(nrounds %% calc_L_rounds!=0){
		stop("nrounds not a multiple of calc_L_rounds. This is not allowed")
	}
}

rm_cluster=function(mpicl){
	if(class(mpicl)[1]=="SOCKcluster"){
			stopCluster(mpicl)
	}else{
		if(class(mpicl)[1]=="environment"){
			rm(mpicl)
		}else{
			stop("neither cluster nor environment")
		}
	}
}

source_on_nodes=function(mpicl,to_be_sourced_on_nodes=NULL){
	if(is.null(to_be_sourced_on_nodes)){
		stop("no sourcing script to source on nodes")
	}
	print("loading functions on each node")
	clusterEvalQ_or_envEvalQ(mpicl,myfile=to_be_sourced_on_nodes)
}

get_cluster_size=function(mpicl){
	if(class(mpicl)[1]=="SOCKcluster"){
		cl_len=length(mpicl)
	}else{
		if(class(mpicl)[1]=="environment"){
			cl_len=1
		}else{
			stop("neither cluster nor envir.")
		}
	}
	return(cl_len)
}

paste_theta_b_alphai_vect_second=function(all_theta_b_alphai_vect_list,all_theta_b_alphai_gene_indices,Mtot){
	all_theta_b_alphai_vect_second=rep(NA,Mtot)
	tt=unlist(lapply(all_theta_b_alphai_gene_indices,function(x){length(x[,1])}))
	for(j in c(1:length(tt))){
		cur_starter=1
		cur_stopper=0
		for(jj in c(1:tt[j])){
			cur=all_theta_b_alphai_gene_indices[[j]][jj,]
			cur_len=cur[2]-cur[1]+1					
			cur_stopper=cur_stopper+cur_len
			all_theta_b_alphai_vect_second[cur[1]:cur[2]]=all_theta_b_alphai_vect_list[[j]][cur_starter:cur_stopper]
			cur_starter=cur_stopper+1
		}
	}
	return(all_theta_b_alphai_vect_second)
}

collapse_Eu_kappa_vect=function(index_mat_t,Eu_kappa_pergene_vect){
	if(dim(Eu_kappa_pergene_vect)[1]==1){
		Eu_kappa_vect=as.matrix(t(index_mat_t))%*%Eu_kappa_pergene_vect
	}else{
		Eu_kappa_vect=index_mat_t%*%Eu_kappa_pergene_vect
	}
	return(Eu_kappa_vect)
}


check_obs_el=function(Obs_list){
	objmissing=setdiff(c("yty","ytX","XtX_sqrt","n","mapping_mat","eig_values"), names(Obs_list))
	if(length(objmissing)!=0){
		stop(paste("Obs_list missing, first obj missing: ",objmissing[1]),sep="")
	}
	check_scalar(Obs_list$yty,"Obs_list$yty")
	check_pos(Obs_list$yty,"Obs_list$yty")
	check_scalar(Obs_list$n,"Obs_list$n")
	check_posinteger(Obs_list$n,"Obs_list$n")
	if(!is.matrix(Obs_list$XtX_sqrt)){
		stop("Obs_list$XtX_sqrt is not a matrix")
	}
	check_numeric(Obs_list$XtX_sqrt)
	check_notnan(Obs_list$XtX_sqrt)
	check_matrix(Obs_list$XtX_sqrt,"XtX_sqrt")
	check_notnan(Obs_list$XtX_sqrt,"XtX_sqrt")
	if(dim(Obs_list$ytX)[1]!=1){
		stop("Obs_list$ytX is not one dim")
	}
	check_matrix(Obs_list$ytX,"ytX")
	check_notnan(Obs_list$ytX,"ytX")
	check_logicalmatrix(Obs_list$mapping_mat,"mapping_mat")
	check_notnan(Obs_list$mapping_mat,"mapping_mat")
	check_vector(Obs_list$eig_values,"eig_values")
	check_notnan(Obs_list$eig_values,"eig_values")
	if(dim(Obs_list$mapping_mat)[1]!=dim(Obs_list$XtX_sqrt)[1]){
		stop("Obs_list$mapping_mat dim not equal to Obs_list$XtX_sqrt")
	}
	if(dim(Obs_list$mapping_mat)[1]!=dim(Obs_list$XtX_sqrt)[1]){
		stop("Obs_list$mapping_mat dim not equal to Obs_list$XtX_sqrt")
	}
	if(dim(Obs_list$mapping_mat)[1]!=dim(Obs_list$ytX)[2]){
		stop("Obs_list$mapping_mat dim not equal to Obs_list$ytX")
	}
}

check_E_el=function(E_list){
	objmissing=setdiff(c("Eu_b_list","Eu_alphai","Eu_lambda","Eu_kappa"), names(E_list))
	if(length(objmissing)!=0){
		stop(paste("E_list missing, first obj missing: ",objmissing[1]),sep="")
	}
	check_vector(E_list$Eu_kappa,"Eu_kappa")
	check_notnan(E_list$Eu_kappa,"Eu_kappa")
	check_pos(E_list$Eu_kappa[2],"Eu_kappa[2]")
	check_vector(E_list$Eu_lambda,"Eu_lambda")
	check_notnan(E_list$Eu_lambda,"Eu_lambda")
	check_pos(E_list$Eu_lambda[2],"Eu_lambda[,2]")
	if(length(E_list$Eu_b_list)!=2 && length(E_list$Eu_b_list)!=4){
		stop("Eu_b_list length has to be 2 or 4.")
	}
	if(length(E_list$Eu_b_list)==2){
		check_vector(E_list$Eu_b_list[[1]],"Eu_b_list[[1]]")
		check_notnan(E_list$Eu_b_list[[1]],"Eu_b_list[[1]]")
		check_matrix(E_list$Eu_b_list[[2]],"Eu_b_list[[2]]")
		check_notnan(E_list$Eu_b_list[[2]],"Eu_b_list[[2]]")
	}
	check_matrix(E_list$Eu_alphai,"E_list$Eu_alphai")
	check_notnan(E_list$Eu_alphai,"E_list$Eu_alphai")
	check_pos(E_list$Eu_alphai[,2],"E_list$Eu_alphai[,2]")
}

check_gene_dat_list=function(gene_dat_list=NULL){
	tester=unlist(lapply(c(1:length(gene_dat_list)),function(i){dim(gene_dat_list[[i]]$Obs_list$XtX_sqrt)[2]}))
	if(sum(tester==0)>0){
		print("some gene data seem to have empty XtX_sqrt matrices. This is not allowed")
		stop("abort")
	}
	tester=unlist(lapply(c(1:length(gene_dat_list)),function(i){dim(gene_dat_list[[i]]$Obs_list$XtX_sqrt)[1]}))
	if(sum(tester==0)>0){
		print("some gene data seem to have empty XtX_sqrt matrices. This is not allowed")
		stop("abort")
	}
	lapply(gene_dat_list,function(gene_dat_el){
		if(sum(is.element(c("Obs_list","E_list"),names(gene_dat_el)))!=2){
			stop("not all elements in gene_dat_el")
		}
	})
	lapply(gene_dat_list,function(gene_dat_el){check_obs_el(gene_dat_el$Obs_list)})
	lapply(gene_dat_list,function(gene_dat_el){check_E_el(gene_dat_el$E_list)})
	lapply(gene_dat_list,function(gene_dat_el){
		aa=length(gene_dat_el$Obs_list$ytX)
		bb=dim(gene_dat_el$E_list$Eu_alphai)[1]
		if(aa!=bb){
			stop("Obs_list and E_list lengths are inconsistent.")
		}
	})
}

check_hyperparameter_list=function(hyperparameter_list){
	if(sum(class(hyperparameter_list)=="bagea_hyperparameter_list")==0){
		stop("hyperparameter_list argument is not a bagea_hyperparameter_list object.")
	}
	allparams=c("gamma1","tau1","lambda1","phi1","phi2","rho1","rho2","xi1","xi2")
	for(i in c(1:length(allparams))){
		check_pos(hyperparameter_list[[allparams[i]]],allparams[i])
	}
	if(!is.null(hyperparameter_list$p)){
		shift_params=c("p","chi1","chi2")
		for(i in c(1:length(shift_params))){
			check_pos(hyperparameter_list[[shift_params[i]]],shift_params[i])
		}		
		check_pos(hyperparameter_list[["p"]],"p")
	}
}

check_input=function(hyperparameter_list=NULL,calc_L=NULL,calc_B=NULL,ncores=NULL,nrounds=NULL,calc_L_rounds=NULL,calc_B_rounds=NULL,write_logfile=NULL,approximation_mode=NULL,nu_names=NULL,ak_names=NULL){
	check_hyperparameter_list(hyperparameter_list)
	check_boolean(calc_L,"calc_L")
	check_boolean(calc_B,"calc_B")
	check_posinteger(ncores,"ncores")
	check_posinteger(nrounds,"nrounds")
	check_posinteger(calc_L_rounds,"calc_L_rounds")
	check_posinteger(calc_B_rounds,"calc_B_rounds")
	check_boolean(write_logfile,"write_logfile")
	check_boolean(approximation_mode,"approximation_mode")
	sapply(nu_names,function(x){check_string(x,"nu_names")})
	sapply(ak_names,function(x){check_string(x,"ak_names")})
	check_core_nrs(ncores)
	check_round_indices(nrounds,calc_L_rounds,calc_B_rounds)
}


mywhere2=function(name, env=parent.frame(1)){
    stopifnot(is.character(name), length(name) == 1)
    browser()
    while(TRUE){
    	if (identical(env, emptyenv())){
        	stop("Can't find ", name, call. = FALSE)
        }
    	if (exists(name, env, inherits = FALSE)){
        	return(env)
    	}
    	env=parent.env(env)
    }
}


mywhere=function(name, mylevel=1){
    stopifnot(is.character(name), length(name) == 1)
    while(TRUE){
    	env=parent.frame(mylevel)
    	if (identical(env, emptyenv())){
        	stop("Can't find ", name, call. = FALSE)
        }
    	if (exists(name, env, inherits = FALSE)){
        	return(env)
    	}
    	mylevel=mylevel+1
    }
}



