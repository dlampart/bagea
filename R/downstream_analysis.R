
#' show E_ak_convergence
#' tries to show how much Eu_ak values have converged. For each Eu_ak parameter, the trace of values
#' across the entire run is split into segments and for each segment the variances are calculated.
#
#' param  bagea_result
#' A bagea output object
#' param  partition_nr
#' number of segments into which the Eu_ak trace should be split.
#' return a matrix of size n,m where n equals the number of ak estimates and m equals partition_nr
#' export giving the variance of  the Eu_ak trace in the respective trace partition.
get_bagea_E_ak_convergence_matrix=function(bagea_result,partition_nr){
	check_bagea_output(bagea_result,"bagea_result")
	check_posinteger(partition_nr,"partition_nr")
	m=partition_nr
	n=dim(bagea_result$Eu_ak_trace[[1]])[1]
	mylen=length(bagea_result$Eu_ak_trace)
	myby=floor(mylen/partition_nr)
	myseq=ceiling(seq(1,mylen,by=myby)-1)
	mydiff=myseq[2]-myseq[1]
	a=matrix(NA,n,m)
	mycolnames=rep("",m)
	for(k in c(1:n)){

		tt=unlist(lapply(c(1:mylen),function(i){bagea_result$Eu_ak_trace[[i]][k,2]}))
		a[k,]=unlist(lapply(c(1:m),function(i){var(tt[c(1:mydiff)+myseq[i]])}))
	}
	mycolnames=unlist(lapply(c(1:m),function(i){paste("var",paste(c(1,mydiff)+myseq[i],collapse="-"),sep="")}))
	colnames(a)=mycolnames
	return(a)
}


#' Remove genewise matrices from bagea_output object to reduce size.
#'
#' @return
#' An R-object of class \code{bagea_reduced_output}; This is \code{bagea_output} object with genewise matrices removed to save space.
#' @param bagea_result
#' An R-object of class \code{bagea}; This is \code{bagea_output} object with genewise matrices removed to save space.
#' keep_only_essential_bagea_results
#' removes all observed data objects from the bagea_result object
#' @export
reduce_bagea_results=function(bagea_result){
	check_bagea_output(bagea_result,"bagea_output")
	mylen=length(bagea_result$gene_dat_list)
	for(i in c(1:mylen)){
		bagea_result[["gene_dat_list"]][[i]][["E_list"]][["Eu_b_list"]][["theta_M_inv"]]=NULL
		bagea_result[["gene_dat_list"]][[i]][["E_list"]][["Eu_b_list"]][["Eu_b2_formula"]]=NULL
		bagea_result[["gene_dat_list"]][[i]][["E_list"]][["Eu_b_list"]][["Eu_b2_formula_expanded"]]=NULL
		bagea_result[["gene_dat_list"]][[i]][["Obs_list"]]=NULL
	}
	class(bagea_result)=c("bagea_reduced_output",class(bagea_result))
	return(bagea_result)
}


#' removes all observed data objects from the bagea_result object
get_Eu_b_table=function(bagea_result){
	check_bagea_output(bagea_result,"bagea_result")
	mylen=length(bagea_result$gene_dat_list)
	mynames=names(bagea_result$gene_dat_list)
	all_lists=lapply(c(1:mylen),function(i){
		eu_b1=bagea_result[["gene_dat_list"]][[i]][["E_list"]][["Eu_b_list"]][["Eu_b_1"]]
		rs_names=rownames(eu_b1)
		dt=data.table(gname=mynames[i],rs=rs_names,b_hat=eu_b1[,1])
		return(dt)
	})
	dt_out=my_docall_rbind(all_lists, 10)
	return(dt_out)
}

#' extracts snp names
#'
#' Returns all SNP names of SNPs that are present in the bagea_oberved data object.
#' @return
#' A vector snp names.
#' @param bagea_obs_list
#' An R-object of class bagea_observed_data producted by \code{\link{prepare_observed_data}} or \code{\link{prepare_observed_1KG_data}}.
#' @export
extract_snp_names=function(bagea_obs_list){
	check_bagea_observed_data(bagea_obs_list,"bagea_obs_list")
	snp_lists=lapply(bagea_obs_list$obs_list,function(x){return(colnames(x$ytX))})
	all_snps=unlist(snp_lists)
	snps_uniq=unique(all_snps)
	return(snps_uniq)
}


#' rescaled effect size for omegas based on the variance and bias of directed annotation. i.e.
# \omega_i_sc=\omega_i*sqrt(\sum_j((V^i_j)^2))/sqrt(n_genes); if use_Fnu_scaling=TRUE, we replace each V^i_j with (V^i_j)*(Fnu)_j. 
##
calc_rescaled_omega_estimates=function(bagea_result,use_Fnu_scaling=TRUE,supress_check=FALSE){
	if(!supress_check){
		check_bagea_output(bagea_result,"bagea_result")
	}
	omega1=bagea_result$Eu_omega_list[[1]]
	i=1
	if(use_Fnu_scaling){
		nu1=bagea_result$Eu_nu_list[[1]]
		nu1=nu1[2:dim(nu1)[1],,drop=FALSE]

	}
	ngenes=length(bagea_result$gene_dat_list)
	all_vals=lapply(c(1:ngenes),function(i){
		cur_obs=bagea_result$gene_dat_list[[i]][["Obs_list"]]
		if(use_Fnu_scaling){
			mmap=cur_obs$mapping_mat
			mymatch=match(rownames(nu1),colnames(mmap))
			resc=mmap[,mymatch,drop=FALSE]%*%nu1
			resc=resc+1
			mmap2=cur_obs$mapping_shift_mat
			rescmat=resc%*%t(rep(1,dim(mmap2)[2]))
			ret=colSums((mmap2*rescmat)^2)
		}else{
			ret=colSums(cur_obs$mapping_shift_mat^2)
		}
		return(ret)
	})
	all_vals=colSums(do.call("rbind",all_vals))
	all_n=sapply(c(1:ngenes),function(i){
		return(dim(bagea_result$gene_dat_list[[i]][["Obs_list"]]$mapping_shift_mat)[1])
	})
	rescaling_fact_sq=all_vals/ngenes
	rescaling_fact=sqrt(rescaling_fact_sq)
	out=omega1*rescaling_fact_sq
	return(out)
}




get_XtX=function(x){
	XtX_sqrt=x[["XtX_sqrt"]]
	XtX=XtX_sqrt%*%t(XtX_sqrt)
	a=sum(XtX_sqrt^2)
	b=sum(x$eig_values)
	b-a
	mym=dim(XtX)[1]
	if(!(mym>1)){
		stop("errmv3-c")
	}
	XtX=XtX+diag(mym)*(b-a)/mym
	return(XtX)
}


#' prepares data table to calculate average MSE_dir .
#'
#' Using omega estimates extracted from \code{bagea_result}, calculates director predictors mu for genes in \code{observed_data} and derives their \code{MSE_dir}.
#' @return data.table. Contains MSE_dir and S for all genes in \code{observed_data}.
#' @param observed_data
#' An R-object of class bagea_observed_data, producted for instance by \code{\link{prepare_observed_data}}
#' @param bagea_result
#' An R-object of class \code{bagea_output}, producted by \code{\link{run_bagea}}
#' @export
predict_directed_fast=function(observed_data,bagea_result){
	calc_eta=function(F,nu,V,omega){
		eta = (V%*%omega)*(F%*%nu)
	}
	calc_mu=function(X,eta){
		mu = X%*%eta
	}
	calc_mse_dir=function(mu,y){
		n=length(y)
		sum((y-mu)^2)/n
	}
	calc_mse_dir_approx=function(eta,z,Sigma,n){
		val1=t(eta)%*%Sigma%*%eta
		val2=-2*t(eta)%*%z/sqrt(n)
		val3=1
		return(sum(val1+val2+val3))
	}
	get_results_params=function(bagea_result){
		omega1=bagea_result$Eu_omega_list[[1]]
		nu1=bagea_result$Eu_nu_list[[1]]
		params=list(omega1=omega1,nu1=nu1)
	}
	get_obs_params=function(obs_list){
		XtX=get_XtX(obs_list)
		V=obs_list$mapping_shift_mat
		F=cbind(TRUE,obs_list$mapping_mat)
		ytX=obs_list$ytX
		n_indiv=obs_list$n
		z=t(ytX)/sqrt(n_indiv)
		Sigma=XtX/n_indiv
		list(V=V,F=F,n_indiv=n_indiv,z=z,Sigma=Sigma)
	}
	get_fullobs_params=function(obs_list){
		myl=get_obs_params(obs_list)
		myl[["X"]]=obs_list$X
		myl[["y"]]=obs_list$y
		return(myl)
	}
	calc_S=function(mu){
		t(mu)%*%mu/length(mu)
	}
	calc_yty_sc=function(y,n){
		t(y)%*%y/n
	}
	calc_S_approx=function(eta,Sigma,n){
		(t(eta)%*%Sigma%*%eta)
	}
	fitted_params=get_results_params(bagea_result)
	myKEEP_X_y=observed_data$settings$KEEP_X_y
	only_approx=FALSE
	if(is.null(myKEEP_X_y) || myKEEP_X_y!=TRUE){
		print("KEEP_X_y not TRUE. only approximation values computed.")
		only_approx=TRUE
	}
	gene_names=names(observed_data$obs_list)
	all_vals=lapply(c(1:length(gene_names)),function(i){
		if(i %% 100==0){
			print(i)
		}
		cur_obs=observed_data$obs_list[[i]]
		if(only_approx){
			obs_param=get_obs_params(cur_obs)
		}else{
			obs_param=get_fullobs_params(cur_obs)
		}
		eta=calc_eta(obs_param$F,fitted_params$nu,obs_param$V,fitted_params$omega)
		mse_dir_approx=calc_mse_dir_approx(eta=eta,z=obs_param$z,Sigma=obs_param$Sigma,obs_param$n)
		S_approx=calc_S_approx(eta=eta,Sigma=obs_param$Sigma,obs_param$n)[1]
		maxpval=pnorm(-max(abs(obs_param$z)))*2
		out=data.table(gname=gene_names[i],S_approx=S_approx,mse_dir_approx=mse_dir_approx,maxpval=maxpval)
		if(!only_approx){
			mu=calc_mu(obs_param$X,eta)
			mse_dir=calc_mse_dir(mu,obs_param$y)
			myS=calc_S(mu)[1]
			yty_sc=calc_yty_sc(obs_param$y,obs_param$n)[1]
			out[,mse_dir:=mse_dir]
			out[,myS:=myS]
			out[,yty_sc:=yty_sc]			
		}
		return(out)
	})
	return(do.call("rbind",all_vals))
}





