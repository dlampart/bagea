get_Eu_kappa_vect_from_update_outl=function(update_out_list){
	Eu_kappa_vect=do.call("rbind",lapply(update_out_list$Eu_alphaiTimesEu_kappa,function(x){
					attr(x,"Eu_kappa_vect")
		}))
	return(Eu_kappa_vect)
}

pull_out_dat=function(mpicl,update_out_list,traces_l,index_range_mappings,gene_names,settings,D_mat,D_names_list
	){
	Eu_nu1_trace=traces_l$Eu_nu1_trace
	Eu_omega1_trace=traces_l$Eu_omega1_trace
	Eu_ak_trace=traces_l$Eu_ak_trace
	Ltot=traces_l$Ltot
	all_Ls=traces_l$all_Ls
	####
	global_E_list=update_out_list$global_E_list
	Eu_alphai_vect=update_out_list$Eu_alphai_vect
	Eu_kappa_vect=get_Eu_kappa_vect_from_update_outl(update_out_list)
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
	out_dat_list$Eu_ak=global_E_list$Eu_ak
	out_dat_list$Eu_lambda2=global_E_list$Eu_lambda2
	out_dat_list$Eu_tau2=global_E_list$Eu_tau2
	out_dat_list$Eu_ak_trace=Eu_ak_trace
	out_dat_list$Ltot=Ltot
	out_dat_list$all_Ls=all_Ls
	out_dat_list$gene_dat_list=gene_dat_list
	out_dat_list$Eu_omega_list=global_E_list$Eu_omega_list
	out_dat_list$Eu_delta=global_E_list$Eu_delta
	out_dat_list$Eu_upsilon_list=global_E_list$Eu_upsilon_list
	out_dat_list$Eu_omega1_trace=Eu_omega1_trace
	out_dat_list$Eu_nu_list=global_E_list$Eu_nu_list
	out_dat_list$Eu_nu1_trace=Eu_nu1_trace
	out_dat_list$Eu_chi2j=global_E_list$Eu_chi2j
	out_dat_list$settings=settings
	out_dat_list$D_mat=D_mat
	out_dat_list$D_names_list=D_names_list
	return(out_dat_list)
}

update_params=function(i,mpicl,Eu_ak_old,global_E_list,mysettings,theta_root_list,global_theta_star_list,mapping_mat_ak_nonnull_element_list,hyperparameter_list,mapping_mat_expanded_ak,theta_alphai_vect_first_row,D_mat,theta_alphai_ak_first_col,Eu_gammai_vect,nu_names){
	gamma1=hyperparameter_list$gamma1
	chi1_vec=hyperparameter_list$chi1_vec
	if(is.null(chi1_vec)){
		stop("err2emcom: hyperparameter_list not updated yet.") 
	}
	lambda1=hyperparameter_list$lambda1
	#####
	local_calc_L= mysettings$calc_L && (((i %% mysettings$calc_L_rounds)==0) || i==1)
	local_calc_B= mysettings$calc_B && (((i %% mysettings$calc_B_rounds)==0) || i==1)
	print(Sys.time())
	print("round")
	print(i)
	if(local_calc_B){
		##############################
		## update Eu_lambdaj;Eu_bj
		##############################
		Eu_alphaiTimesEu_kappa=run_cycle_across_allgenes_mpi_loc(mpicl=mpicl)
		Eu_alphai_vectTimesEu_kappa_vect=my_docall_rbind(Eu_alphaiTimesEu_kappa)
		theta_lambdai_lambda2_mat=do.call("rbind",lapply(Eu_alphaiTimesEu_kappa,function(x){
				attr(x,"theta_lambdai_lambda2")
			}))
		all_Eu_b1=do.call("rbind",lapply(Eu_alphaiTimesEu_kappa,function(x){
				attr(x,"all_Eu_b1")
		}))
		Eu_alphai_vect=do.call("rbind",lapply(Eu_alphaiTimesEu_kappa,function(x){
				attr(x,"Eu_alphai_vect")
			}))
			theta_lambdai_lambda2=colSums(theta_lambdai_lambda2_mat)
			theta_b_nu=list()
			theta_b_nu[[1]]=Reduce("+",lapply(Eu_alphaiTimesEu_kappa,function(x){
					attr(x,"theta_b_nu")[[1]]}))
			theta_b_nu[[2]]=Reduce("+",lapply(Eu_alphaiTimesEu_kappa,function(x){
					attr(x,"theta_b_nu")[[2]]}))
			theta_kappaj_tau2_mat=do.call("rbind",lapply(Eu_alphaiTimesEu_kappa,function(x){
					attr(x,"theta_kappaj_tau2")
				}))
			theta_kappaj_tau2=colSums(theta_kappaj_tau2_mat)
			##############################
			## END: update Eu_lambdaj;Eu_bj
			##############################
		}
		outl=calc_Eu_ak_optimized(i,global_E_list$Eu_ak,Eu_ak_old,mapping_mat_ak_nonnull_element_list=mapping_mat_ak_nonnull_element_list,mapping_mat_expanded=mapping_mat_expanded_ak,local_calc_L,theta_ak=theta_root_list$theta_ak,Eu_alphai_vectTimesEu_kappa_vect,gamma1,theta_alphai_ak_first_col=theta_alphai_ak_first_col,
			Eu_gammai_vect=Eu_gammai_vect
			)
		global_E_list$Eu_ak=outl[["Eu_ak"]]
		Eu_ak_old=outl[["Eu_ak_old"]]
		global_theta_star_list$theta_ak_star=outl[["theta_ak_star"]]
		Eu_gammai_vect=outl$Eu_gammai_vect
		##############################
		## update Eu_lambda2
		##############################		
		global_theta_star_list$theta_lambda2_star=theta_root_list$theta_lambda2+theta_lambdai_lambda2
		global_E_list$Eu_lambda2=calc_Eu_lambda2_via_theta(global_theta_star_list$theta_lambda2_star)
		##############################
		## END: update Eu_lambda2
		##############################
		theta_lambda=calc_theta_lambda(lambda1,global_E_list$Eu_lambda2[2],genespecific_lambda1=mysettings$approximation_mode)
		clusterExport_or_envExport(mpicl,c("theta_lambda"),envir=environment())
		########################
		global_theta_star_list$theta_tau2_star=theta_root_list$theta_tau2+theta_kappaj_tau2
		###########################
		## update Eu_tau2
		###########################
		global_E_list$Eu_tau2=calc_Eu_tau2_via_theta(global_theta_star_list$theta_tau2_star)
		###########################
		## END: update Eu_tau2
		###########################
		###########################
		## update Eu_nu_list
		###########################
		theta_nu_star_list=theta_root_list$theta_nu_list
		theta_nu_star_list[[1]]=theta_nu_star_list[[1]]+theta_b_nu[[1]]
		theta_nu_star_list[[2]]=theta_nu_star_list[[2]]+theta_b_nu[[2]]
		global_theta_star_list$theta_nu_star_list=theta_nu_star_list
		global_E_list$Eu_nu_list=calc_Eu_nu_via_theta(global_theta_star_list$theta_nu_star_list)
		print("Eu_nu1")
		print(global_E_list$Eu_nu_list[[1]])
		upload_global_Elist(mpicl,global_E_list)
		###########################
		## END: update Eu_nu_list
		###########################
		###########################
		## update Eu_omega_list
		###########################
		theta_b_omega_lists=clusterCall_or_envCall(mpicl, run_calc_theta_b_omega_via_Eu_nu_on_node)
		theta_b_omega=list()
		theta_b_omega[[1]]=Reduce("+",lapply(theta_b_omega_lists,function(x){x[[1]]}))
		theta_b_omega[[2]]=Reduce("+",lapply(theta_b_omega_lists,function(x){x[[2]]}))
		theta_omega_list=calc_theta_omegas_via_Eu_upsilon(global_E_list$Eu_upsilon_list,D_mat)
		theta_omega_star_list=theta_omega_list
		theta_omega_star_list[[1]]=theta_omega_star_list[[1]]+theta_b_omega[[1]]
		theta_omega_star_list[[2]]=theta_omega_star_list[[2]]+theta_b_omega[[2]]
		global_theta_star_list$theta_omega_star_list=theta_omega_star_list
		global_E_list$Eu_omega_list=calc_Eu_omega_via_theta(global_theta_star_list$theta_omega_star_list)
		outl=process_upsilon(chi1_vec,global_E_list$Eu_chi2j,global_E_list$Eu_omega_list,global_E_list$Eu_upsilon_list,D_mat)
		global_E_list$Eu_upsilon_list=outl[["Eu_upsilon_list"]]
		global_theta_star_list$theta_upsilon_star_list=outl[["theta_upsilon_star_list"]]
		theta_upsilon_chi2j=calc_theta_upsilon_chi2j_via_Eu(global_E_list$Eu_upsilon_list,chi1_vec)
		global_theta_star_list$theta_chi2j_star=theta_root_list$theta_chi2j+theta_upsilon_chi2j
		global_E_list$Eu_chi2j=calc_Eu_chi2j_via_theta(global_theta_star_list$theta_chi2j_star)
		print("Eu_omega1")
		print(global_E_list$Eu_omega_list[[1]])
		upload_global_Elist(mpicl,global_E_list)
		rownames(global_E_list$Eu_nu_list[[1]])=c("background",nu_names)
		outl=list()
		outl$Eu_alphaiTimesEu_kappa=Eu_alphaiTimesEu_kappa
		outl$Eu_alphai_vect=Eu_alphai_vect
		outl$global_E_list=global_E_list
		outl$global_theta_star_list=global_theta_star_list
		outl$all_Eu_b1=all_Eu_b1
		outl$Eu_ak_old=Eu_ak_old
		outl$Eu_gammai_vect=Eu_gammai_vect
		return(outl)
}

calc_L_global=function(mpicl,mapping_mats_l,update_out_list,g_root_list,theta_root_list,hyperparameter_list,mapping_mat_expanded_ak,theta_alphai_vect_first_row,D_mat,mapping_mat_expanded_nu,approximation_mode,mapping_shift_mat_expanded){
	gamma1=hyperparameter_list$gamma1
	chi1_vec=hyperparameter_list$chi1_vec
	if(is.null(chi1_vec)){
		stop("err2emcom: hyperparameter_list not updated yet.") 
	}
	lambda1=hyperparameter_list$lambda1
	######
	mapping_mat_expanded_nu=mapping_mats_l$mapping_mat_expanded_nu
	mapping_shift_mat_expanded=mapping_mats_l$mapping_shift_mat_expanded
	mapping_mat_expanded_ak=mapping_mats_l$mapping_mat_expanded_ak

	Eu_alphaiTimesEu_kappa=update_out_list[["Eu_alphaiTimesEu_kappa"]]
	global_E_list=update_out_list[["global_E_list"]]
	global_theta_star_list=update_out_list[["global_theta_star_list"]]
	Eu_alphai_vect=update_out_list[["Eu_alphai_vect"]]
	all_Eu_b1=update_out_list[["all_Eu_b1"]]

	Eu_omega_list=global_E_list$Eu_omega_list
	Eu_nu_list=global_E_list$Eu_nu_list
	Eu_chi2j=global_E_list$Eu_chi2j
	Eu_ak=global_E_list$Eu_ak
	Eu_lambda2=global_E_list$Eu_lambda2
	Eu_tau2=global_E_list$Eu_tau2
	Eu_upsilon_list=global_E_list$Eu_upsilon_list
	####
	theta_lambda2_star=global_theta_star_list$theta_lambda2_star
	theta_tau2_star=global_theta_star_list$theta_tau2_star
	theta_omega_star_list=global_theta_star_list$theta_omega_star_list
	theta_nu_star_list=global_theta_star_list$theta_nu_star_list
	theta_upsilon_star_list=global_theta_star_list$theta_upsilon_star_list
	theta_ak_star=global_theta_star_list$theta_ak_star
	theta_chi2j_star=global_theta_star_list$theta_chi2j_star
	####
	all_Eu_b2_diag=do.call("rbind",lapply(Eu_alphaiTimesEu_kappa,function(x){
			attr(x,"all_Eu_b2_diag")
	}))
	Eu_kappa_pergene_vect=do.call("rbind",lapply(Eu_alphaiTimesEu_kappa,function(x){
			attr(x,"Eu_kappa_pergene_vect")
	}))
	Eu_lambda_vect=do.call("rbind",lapply(Eu_alphaiTimesEu_kappa,function(x){
			attr(x,"Eu_lambda_vect")
	}))
	Eu_kappa_vect=do.call("rbind",lapply(Eu_alphaiTimesEu_kappa,function(x){
			attr(x,"Eu_kappa_vect")
	}))
	L_kappa_star_part=Reduce("+",lapply(Eu_alphaiTimesEu_kappa,function(x){
			attr(x,"L_kappa_star_part")
	}))
	L_alphai_star_part=Reduce("+",lapply(Eu_alphaiTimesEu_kappa,function(x){
			attr(x,"L_alphai_star_part")
	}))
	n_genes=dim(Eu_kappa_pergene_vect)[1]
	L_lambda2=sum((theta_root_list$theta_lambda2-theta_lambda2_star)*Eu_lambda2)+g_root_list$g_lambda2-calc_g_lambda2_via_theta(theta_lambda2_star)
	L_tau2=sum((theta_root_list$theta_tau2-theta_tau2_star)*Eu_tau2)+g_root_list$g_tau2-calc_g_tau2_via_theta(theta_tau2_star)
	g_kappa=calc_g_kappa_via_Eu_tau2(hyperparameter_list[["tau1"]],Eu_tau2)
	theta_kappa_pergene_vect=calc_theta_kappaj_via_Eu_tau2(hyperparameter_list[["tau1"]],Eu_tau2,n_genes)
	L_kappa=sum((theta_kappa_pergene_vect)*Eu_kappa_pergene_vect)+g_kappa*dim(theta_kappa_pergene_vect)[1]-L_kappa_star_part
	##########
	#L_alphai
	##########
	L_alphai_minus_star=sum(calc_L_alpha_minus_star(gamma1,Eu_ak,mapping_mat=mapping_mat_expanded_ak,Eu_alphai=Eu_alphai_vect,Eu_kappa=Eu_kappa_vect,theta_alphai_first_row=theta_alphai_vect_first_row))
	L_alphai=L_alphai_minus_star-L_alphai_star_part
	##########
	#L_ak
	##########
	L_ak=calculate_L_ak(Eu_ak,theta_root_list$theta_ak,theta_ak_star)
	##########
	#L_omega
	##########
	g_omega_pa=sum(calc_g_omega_via_Eu_upsilon(Eu_upsilon_list,D_mat))
	theta_omega_list=calc_theta_omegas_via_Eu_upsilon(Eu_upsilon_list,D_mat)
	L_omega=sum((theta_omega_list[[1]]-theta_omega_star_list[[1]])*Eu_omega_list[[1]])+sum((theta_omega_list[[2]]-theta_omega_star_list[[2]])*Eu_omega_list[[2]])+g_omega_pa-calc_g_omega_via_theta(theta_omega_star_list)
	##########
	#L_upsilon
	##########
	L_upsilon=calc_L_upsilon(chi1_vec,D_mat,Eu_chi2j,theta_upsilon_star_list,Eu_upsilon_list)
	##########
	#L_chi2j
	##########
	g_chi2j_pa=g_root_list[["g_chi2j"]]
	g_chi2j_theta_star=calc_g_chi2j_via_theta(theta_chi2j_star)
	L_chi2j=sum((theta_root_list[["theta_chi2j"]]-theta_chi2j_star)*Eu_chi2j)+sum(g_chi2j_pa)-sum(g_chi2j_theta_star)
	##########
	#L_nu
	##########
	g_nu_pa=calc_g_nu_via_theta(theta_root_list$theta_nu_list)[1]
	L_nu=(sum((theta_root_list$theta_nu_list[[1]]-theta_nu_star_list[[1]])*Eu_nu_list[[1]])+sum((theta_root_list$theta_nu_list[[2]]-theta_nu_star_list[[2]])*Eu_nu_list[[2]])+g_root_list$g_nu_list-calc_g_nu_via_theta(theta_nu_star_list))[1]
	##########
	#process Ls on node
	##########
	Ls_from_node_dt=data.table(t(colSums(do.call("rbind",clusterCall_or_envCall(mpicl,process_Ls_on_node_wshift_list)))))
	Ltot_b=calc_L_b(Eu_alphai_vect=Eu_alphai_vect,Eu_omega_list=Eu_omega_list,mapping_shift_mat_expanded=mapping_shift_mat_expanded,Eu_nu_list=Eu_nu_list,mapping_mat_expanded_nu=mapping_mat_expanded_nu,Ls_from_node_dt=Ls_from_node_dt,all_Eu_b1=all_Eu_b1,all_Eu_b2_diag=all_Eu_b2_diag)
	Ltot_lambda=calc_L_lambda(lambda1=lambda1,Eu_lambda2=Eu_lambda2,approximation_mode=approximation_mode,Eu_lambda_vect=Eu_lambda_vect,Ls_from_node_dt=Ls_from_node_dt)
	Ltot_obs=Ls_from_node_dt[,L_obs]
	all_Ls=c(
		Ltot_obs,
		Ltot_b,
		Ltot_lambda,			
		L_lambda2,
		L_alphai,
		L_kappa,
		L_tau2,
		L_ak,
		L_omega,
		L_nu,
		L_upsilon,
		L_chi2j
	)
	return(all_Ls)
}

process_upsilon=function(chi1_vec,Eu_chi2j,Eu_omega_list,Eu_upsilon_list,D_mat){
	myq=dim(D_mat)[2]
	h_vec=apply(D_mat,2,max)
	theta_upsilon_list=calc_theta_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec=h_vec)
	theta_upsilon_star_list=list()
	for(myj in c(1:myq)){
		theta_omega_upsilon_j=calc_theta_omegas_upsilon_via_Eu(Eu_omega_list,Eu_upsilon_list,D_mat,myj)
		theta_upsilon_j=theta_upsilon_list[[myj]]
		theta_upsilon_star_j=theta_omega_upsilon_j+theta_upsilon_j
		theta_upsilon_star_list[[myj]]=theta_upsilon_star_j
		Eu_upsilon_j=calc_Eu_upsilon_via_theta(theta_upsilon_star_j)
		Eu_upsilon_list[[myj]]=Eu_upsilon_j
	}
	outl=list()
	outl[["Eu_upsilon_list"]]=Eu_upsilon_list
	outl[["theta_upsilon_star_list"]]=theta_upsilon_star_list
	return(outl)
}

calc_L_upsilon=function(chi1_vec,D_mat,Eu_chi2j,theta_upsilon_star_list,Eu_upsilon_list){
	myq=dim(D_mat)[2]
	h_vec=apply(D_mat,2,max)
	theta_upsilon_list=calc_theta_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec=h_vec)
	g_upsilon_pa=calc_g_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec)
	g_upsilon_theta_star=calc_g_upsilon_via_theta_list(theta_upsilon_star_list)
	L_upsilons=rep(0,myq)
	for(j in c(1:myq)){
		L_ups_j=(theta_upsilon_list[[j]]-theta_upsilon_star_list[[j]])*Eu_upsilon_list[[j]]
		L_ups_j=sum(L_ups_j)
		L_ups_j=L_ups_j+g_upsilon_pa[j]-g_upsilon_theta_star[j]
		L_upsilons[j]=L_ups_j
	}
	L_upsilon=sum(L_upsilons)
}

calc_Eu_ak_optimized=function(i,Eu_ak,Eu_ak_old,mapping_mat_ak_nonnull_element_list,mapping_mat_expanded,local_calc_L,theta_ak,Eu_alphai_vectTimesEu_kappa_vect,gamma1,theta_alphai_ak_first_col,Eu_gammai_vect){
	myK=dim(Eu_ak)[1]
	theta_ak_star_update=Eu_ak-Eu_ak
	if(i>1){
		Eu_ak_old_last=Eu_ak_old
	}
	Eu_ak_old=Eu_ak
	for(k in c(1:myK)){
		preindex_k=mapping_mat_ak_nonnull_element_list[[k]]
		if(k==1){
			if(local_calc_L){
				Eu_gammai_vect=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat_expanded)
			}else{
				if(i==1){
					Eu_gammai_vect[,2]=calc_Eu_gammai_via_Eu_ak_only_second(Eu_ak,mapping_mat_expanded)[,2]
				}else{
					preindex_k_minus_1=mapping_mat_ak_nonnull_element_list[[myK]]
					Eu_gammai_vect[preindex_k_minus_1,2]=update_Eu_gammai_second_row_via_Eu_ak_only_second_with_preindexing_only_ones(Eu_gammai_vect,Eu_ak_old_last, Eu_ak,myK,preindex_k_minus_1)	
				}
			}
		}else{
			preindex_k_minus_1=mapping_mat_ak_nonnull_element_list[[k-1]]
			Eu_gammai_vect[preindex_k_minus_1,2]=update_Eu_gammai_second_row_via_Eu_ak_only_second_with_preindexing_only_ones(Eu_gammai_vect,Eu_ak_old, Eu_ak,k-1,preindex_k_minus_1)
		}		
 		theta_alphai_ak=theta_ak-theta_ak

		theta_alphai_ak_k_sec=calc_theta_alpha_ak_via_EugammaiTimesEu_alphai_vect_only_second_with_preindexing_only_ones(Eu_gammai_vect,gamma1,Eu_alphai_vectTimesEu_kappa_vect,Eu_ak,k=k, mapping_mat_ak_nonnull_element_list[[k]])
		theta_alphai_ak[k,]=c(theta_alphai_ak_first_col[k],theta_alphai_ak_k_sec)			
		theta_ak_star=theta_ak+theta_alphai_ak
		##############################
		## update Eu_ak
		##############################
		Eu_ak_k=calc_Eu_ak_via_theta(theta_ak_star[k,,drop=F])
		theta_ak_star_update[k,]=theta_ak_star[k,]
		Eu_ak[k,]=Eu_ak_k
		##############################
		## END: update Eu_ak
		##############################
	}
	outl=list()
	outl$Eu_ak=Eu_ak
	outl$Eu_ak_old=Eu_ak_old
	outl$theta_ak_star=theta_ak_star_update
	outl$Eu_gammai_vect=Eu_gammai_vect
	return(outl)
}

initialize_lambda1nrho1_for_approx_mode=function(gene_dat_list=gene_dat_list,lambda1=lambda1,rho1=rho1,constant_noise_prior_pergene=constant_noise_prior_pergene,approximation_mode=approximation_mode){
		if(!approximation_mode){
			outl=list()
			outl[["lambda1"]]=lambda1
			outl[["rho1"]]=rho1
			return(outl)
		}
		my_ms=sapply(gene_dat_list,function(x){dim(x$Obs_list$XtX_sqrt)[1]})
		my_ns=sapply(gene_dat_list,function(x){x$Obs_list$n})
		if(constant_noise_prior_pergene){
			print("shifting parameters in approximation mode..")
			print("constant noise prior per gene.")
			lambda1old=lambda1
			lambda1=lambda1old+my_ns/2-my_ms/2
			rho1old=rho1
			rhodiff=lambda1old*length(gene_dat_list)-sum(lambda1)
			rho1=rho1old+rhodiff
			print(paste("min lambda1 after shift:",min(lambda1)))
			if(min(lambda1)<0){
				stop("min lambda1 below 0. increase base lambda1.")
			}
			print(paste("rho1 after shift:",rho1))
			if(rho1<0){
				stop("rho1 below 0. increase base rho1.")
			}
		}else{
			print("constant noise prior per snp.")
			lambda1old=lambda1
			lambda1=rep(lambda1old,length(my_ns))
			names(lambda1)=names(my_ns)
		}
		outl=list()
		outl[["lambda1"]]=lambda1
		outl[["rho1"]]=rho1
		return(outl)
}

get_matcher_vectors=function(gene_dat_list,ak_names,nu_names){
	if(is.null(ak_names)){
		ak_names=get_mapping_mat_colnames(gene_dat_list)
	}
	if(is.null(nu_names)){
		nu_names=get_mapping_mat_colnames(gene_dat_list)
	}
	ak_matcher=match(ak_names,get_mapping_mat_colnames(gene_dat_list))
	nu_matcher=match(nu_names,get_mapping_mat_colnames(gene_dat_list))
	if(sum(is.na(ak_matcher))>0){
		stop("ak_names not all available.")
	}
	if(sum(is.na(nu_matcher))>0){
		stop("nu_names not all available.")
	}
	outl=list()
	outl[["ak_matcher"]]=ak_matcher
	outl[["nu_matcher"]]=nu_matcher
	outl[["ak_names"]]=ak_names
	outl[["nu_names"]]=nu_names
	return(outl)
}


upload_index_mat=function(mpicl,index_mat,num_genes_per_split,gene_split_index,cluster_indices){
	index_mat_t=t(index_mat)
	##############################
	index_mat_list=list()
	index_mat_t_list=list()
	index_mat_list=lapply(c(1:length(num_genes_per_split)),function(ii){
		return(index_mat[gene_split_index==ii,cluster_indices==ii])
	})
	index_mat_t_list=lapply(c(1:length(num_genes_per_split)),function(ii){
		return(t(index_mat_list[[ii]]))
	})
	aa=clusterApply_or_envApply(mpicl,index_mat_list,load_index_mat_list_to_cluster_node)
	aa=clusterApply_or_envApply(mpicl,index_mat_t_list,load_index_mat_t_list_to_cluster_node)
}

split_and_upload=function(mpicl,gene_dat_list){
	print("split data into chunks..")
	gene_dat_list_splits=split_gene_dat_list(mpicl,gene_dat_list,balanceManually=TRUE)
	print("upload chunks to nodes..")
	myenvs1=clusterApply_or_envApply(mpicl,gene_dat_list_splits,load_data_to_cluster_node)
	print("finished uploading chunks to nodes..")
	gene_names_splits=lapply(gene_dat_list_splits,function(x){names(x)})
	all_lengths=lapply(gene_dat_list_splits,function(x){lapply(x,function(xx){length(xx$Obs_list$mapping_mat[,1])})})
	cluster_indices=do.call("c",lapply(c(1:length(all_lengths)),function(i){rep(i,sum(unlist(all_lengths[[i]])))}))
	num_genes_per_split=lapply(all_lengths,function(x){length(unlist(x))})
	ncores=length(num_genes_per_split)
	gene_split_index=unlist(lapply(c(1:ncores),function(i){rep(i,unlist(num_genes_per_split)[i])}))
	outl=list()
	outl[["gene_names_splits"]]=gene_names_splits
	outl[["cluster_indices"]]=cluster_indices
	return(outl)
}


process_mapping_mats=function(mpicl,gene_dat_list,snp_index_tbl,num_genes_per_split,ak_matcher,nu_matcher,use_compiled){
	cluster_indices=snp_index_tbl[,cluster_index]
	dd=c(1:length(gene_dat_list))
	aa=split(dd, ceiling(seq_along(dd)/130))
	#####
	mapping_mat_expanded_list=lapply(c(1:length(aa)),function(j){do.call("rbind",lapply(aa[[j]],function(i){gene_dat_list[[i]]$Obs_list$mapping_mat}))})
	mapping_mat_expanded=do.call("rbind",mapping_mat_expanded_list)
	mapping_mat_expanded_list=lapply(c(1:length(num_genes_per_split)),function(ii){
		return(mapping_mat_expanded[cluster_indices==ii,])
	})
	a=clusterApply_or_envApply(mpicl,mapping_mat_expanded_list,load_mapping_data_to_cluster_node)
	rm(mapping_mat_expanded_list)
	###### if ak_matcher has different order or has columns missing rerun mapping_expanded creation with correct order (we need to upload all vectors to slave nodes because ak_names may not include all nu_names).rerunning because subsetting sparse large Matrices can take a lot of resources 
	mapping_mat_expanded_full=mapping_mat_expanded 
	if(length(ak_matcher)!=dim(mapping_mat_expanded)[2] || sum(ak_matcher!=colnames(mapping_mat_expanded))!=0){
		print("prepare for main node ..")
		mapping_mat_expanded=mapping_mat_expanded[,ak_matcher]
	}
	inds=which(mapping_mat_expanded!=0,arr.ind=TRUE)
#	inds2=which(mapping_mat_expanded2!=0,arr.ind=TRUE)
	mapping_mat_expanded=sparseMatrix(i=inds[,1],j=inds[,2],dims=dim(mapping_mat_expanded))
	inds=which(mapping_mat_expanded_full!=0,arr.ind=TRUE)
#	inds2=which(mapping_mat_expanded2!=0,arr.ind=TRUE)
	mapping_mat_expanded_full=sparseMatrix(i=inds[,1],j=inds[,2],dims=dim(mapping_mat_expanded_full))
	mapping_mat_expanded_nu=mapping_mat_expanded_full[,nu_matcher,drop=FALSE]
	mapping_mat_expanded_ak=mapping_mat_expanded_full[,ak_matcher,drop=FALSE]
	rm(mapping_mat_expanded_full)
	rm(inds)
	########################
	## create mapping_shift_mat_expanded
	########################
	mapping_shift_mat_expanded_list=lapply(c(1:length(aa)),function(j){do.call("rbind",lapply(aa[[j]],function(i){gene_dat_list[[i]]$Obs_list$mapping_shift_mat}))})
	mapping_shift_mat_expanded=do.call("rbind",mapping_shift_mat_expanded_list)
	mapping_shift_mat_expanded_list=lapply(c(1:length(num_genes_per_split)),function(ii){
		return(mapping_shift_mat_expanded[cluster_indices==ii,])
	})
	a=clusterApply_or_envApply(mpicl,mapping_shift_mat_expanded_list,load_mapping_shift_data_to_cluster_node)
	rm(mapping_shift_mat_expanded_list)
	outl=list()
	###########
	###########
	if(use_compiled){
		clusterCall_or_envCall(mpicl,get_shift_mat_melted_on_node)
		}else{
		clusterCall_or_envCall(mpicl,get_shift_mat_melted_on_node_sham)
	}
	###########
	###########
	mapping_mat_ak_nonnull_element_list=list()
	myK=dim(mapping_mat_expanded_ak)[2]
	for(k in c(1:myK)){
	mapping_mat_ak_nonnull_element_list[[k]]=which(mapping_mat_expanded_ak[,k]!=0)
	}
	outl[["mapping_shift_mat_expanded"]]=mapping_shift_mat_expanded
	outl[["mapping_mat_expanded_ak"]]=mapping_mat_expanded_ak
	outl[["mapping_mat_expanded_nu"]]=mapping_mat_expanded_nu
	outl[["mapping_mat_ak_nonnull_element_list"]]=mapping_mat_ak_nonnull_element_list
	return(outl)
}

upload_global_Elist=function(mpicl,global_E_list){
	Eu_tau2=global_E_list$Eu_tau2
	Eu_nu_list=global_E_list$Eu_nu_list
	Eu_lambda2=global_E_list$Eu_lambda2
	Eu_omega_list=global_E_list$Eu_omega_list
	Eu_ak=global_E_list$Eu_ak
	clusterExport_or_envExport(mpicl,c("Eu_tau2","Eu_nu_list","Eu_omega_list","Eu_lambda2","Eu_ak"),envir=environment())
}


initialize_alphai_kappa=function(mpicl,gamma1,gene_split_index,concat_local_E_l,gene_indices_vec,ncores){
	Eu_gammai_vect=concat_local_E_l$Eu_gammai_vect
	Eu_alphai_vect=concat_local_E_l$Eu_alphai_vect
	theta_alphai_kappa_vect_dt=data.table(a=gamma1,b=-Eu_gammai_vect[,2]*Eu_alphai_vect[,2],gene_indices=gene_indices_vec)
	theta_alphai_kappa_vect_first=theta_alphai_kappa_vect_dt[,list(asum=sum(a)),by=gene_indices][,asum]
	theta_alphai_kappa_vect_first_list=list()
	for(ii in c(1:ncores)){
		theta_alphai_kappa_vect_first_list[[ii]]=theta_alphai_kappa_vect_first[gene_split_index==ii]
	}
	aa=clusterApply_or_envApply(mpicl,theta_alphai_kappa_vect_first_list,load_takvf_list_to_cluster_node)
	return(TRUE)
}

get_gene_lengths=function(gene_dat_list){
	mylen=length(gene_dat_list)
	myMs=unlist(lapply(c(1:mylen),function(i){length(gene_dat_list[[i]]$Obs_list$mapping_mat[,1])}))
}

get_gene_indices_vec=function(gene_dat_list){
	myMs=get_gene_lengths(gene_dat_list)
	Mtot=sum(myMs)
	gene_indices_vec=matrix(0,Mtot,1)
	starter=1
	stopper=0
	mylen=length(gene_dat_list)
	stoppers=cumsum(myMs)
	starters=c(0,stoppers[1:(length(stoppers)-1)])+1	
	for(i in c(1:mylen)){
		print(i)
		starter=starters[i]
		stopper=stoppers[i]
		gene_indices_vec[starter:stopper]=i
	}
	return(gene_indices_vec)
}



initialize_concat_local_E_lists=function(global_E_list,mapping_mat_expanded_ak,gene_dat_list){
	myMs=get_gene_lengths(gene_dat_list)
	Mtot=sum(myMs)
	Eu_alphai_vect=matrix(0,Mtot,2)
	Eu_kappa_vect=matrix(0,Mtot,2)	
	starter=1
	stopper=0
	mylen=length(gene_dat_list)
	stoppers=cumsum(myMs)
	starters=c(0,stoppers[1:(length(stoppers)-1)])+1
	Eu_kappa_pergene_vect=matrix(0,mylen,2)
	for(i in c(1:mylen)){
		print(i)
		starter=starters[i]
		stopper=stoppers[i]
		Eu_alphai_vect[starter:stopper,]=gene_dat_list[[i]]$E_list$Eu_alphai
		Eu_kappa_vect[starter:stopper,1]=gene_dat_list[[i]]$E_list$Eu_kappa[1]
		Eu_kappa_vect[starter:stopper,2]=gene_dat_list[[i]]$E_list$Eu_kappa[2]
		Eu_kappa_pergene_vect[i,]=gene_dat_list[[i]]$E_list$Eu_kappa
	}
	Eu_gammai_vect=calc_Eu_gammai_via_Eu_ak(global_E_list$Eu_ak,mapping_mat_expanded_ak)
	outl=list()
	outl[["Eu_alphai_vect"]]=Eu_alphai_vect
	outl[["Eu_kappa_vect"]]=Eu_kappa_vect
	outl[["Eu_kappa_pergene_vect"]]=Eu_kappa_pergene_vect
	outl[["Eu_gammai_vect"]]=Eu_gammai_vect
	return(outl)
}

initialize_index_range_mappings=function(mpicl,gene_indices_vec,gene_names_splits){
	gene_names=unlist(gene_names_splits)
	gene_indices_vec_boundaries=get_gene_indices_vec_boundaries(gene_indices_vec)
	all_theta_b_alphai_gene_indices=list()
	print(Sys.time())
	print("prepare boundaries ..")
	for(i in c(1:length(gene_names_splits))){
		cur=gene_indices_vec_boundaries[match(gene_names_splits[[i]],gene_names),,drop=F]
		rownames(cur)=gene_names_splits[[i]]
		all_theta_b_alphai_gene_indices[[i]]=cur
	}
	mytrues=clusterApply_or_envApply(mpicl,all_theta_b_alphai_gene_indices,all_theta_b_alphai_gene_indices_to_cluster_node)
	index_range_mappings=all_theta_b_alphai_gene_indices
	return(index_range_mappings)
}

initialize_index_mats=function(mpicl,snp_index_tbl,index_range_mappings){
	num_genes_per_split=get_numgenes_persplit_fromrange_obj(index_range_mappings)
	gene_split_index=get_gene_cluster_index_fromrange_obj(index_range_mappings)
	print("prepare index_mat")
	gene_indices_vec=snp_index_tbl[,gene_index]
	cluster_indices=snp_index_tbl[,cluster_index]
	gene_indices_vec_boundaries=get_gene_indices_vec_boundaries(gene_indices_vec)
	num_genes=dim(gene_indices_vec_boundaries)[1]
	Mtot=length(gene_indices_vec)[1]
	inds_list=lapply(c(1:num_genes),function(i){cbind(c(gene_indices_vec_boundaries[i,1]:gene_indices_vec_boundaries[i,2]),i)})
	inds_concat=do.call("rbind",inds_list)
	index_mat = sparseMatrix(i = inds_concat[,2], j = inds_concat[,1], x = 1,dims=c(num_genes,Mtot))

	print(Sys.time())
	print("transpose index_mat.")
	print("upload index mats.")
	upload_index_mat(mpicl=mpicl,index_mat=index_mat,num_genes_per_split=num_genes_per_split,gene_split_index=gene_split_index,cluster_indices=cluster_indices)
	return(TRUE)
}

initialize_traces=function(nrounds){
	Eu_ak_trace=make_empty_list(paste("round",c(1:nrounds),sep=""))
	Eu_omega1_trace=make_empty_list(paste("round",c(1:nrounds),sep=""))
	Eu_nu1_trace=make_empty_list(paste("round",c(1:nrounds),sep=""))
	Ltot=rep(-Inf,nrounds)
	L_names=c("L_obs","L_b","L_lambda","L_lambda2","L_alphai","L_kappa","L_tau2","L_ak","L_omega","L_upsilon","L_nu","L_chi2j")
	all_Ls=matrix(-Inf,nrounds,length(L_names))
	colnames(all_Ls)=L_names
	traces_l=list()
	traces_l[["Eu_ak_trace"]]=Eu_ak_trace
	traces_l[["Eu_omega1_trace"]]=Eu_omega1_trace
	traces_l[["Ltot"]]=Ltot
	traces_l[["all_Ls"]]=all_Ls
	return(traces_l)
}

update_traces=function(traces_l,i,all_Ls_cur,global_E_list){
	if(!is.null(all_Ls_cur)){
		traces_l[["all_Ls"]][i,]=all_Ls_cur
		traces_l[["Ltot"]][i]=sum(all_Ls_cur)
	}
	traces_l[["Eu_ak_trace"]][[paste("round",i,sep="")]]=global_E_list$Eu_ak
	traces_l[["Eu_omega1_trace"]][[paste("round",i,sep="")]]=global_E_list$Eu_omega_list[[1]]
	traces_l[["Eu_nu1_trace"]][[paste("round",i,sep="")]]=as.matrix(global_E_list$Eu_nu_list[[1]])
	print(traces_l[["Ltot"]][i])
	return(traces_l)
}

inititalize_constant_theta_first_cols=function(gamma1,mapping_mat_expanded_ak){
	Mtot=dim(mapping_mat_expanded_ak)[1]
	theta_alphai_ak_first_col=gamma1*colSums(mapping_mat_expanded_ak)
	theta_alphai_vect_first_col=rep(gamma1-1,Mtot)
	constant_theta_first_cols_l=list()
	constant_theta_first_cols_l[["theta_alphai_ak_first_col"]]=theta_alphai_ak_first_col
	constant_theta_first_cols_l[["theta_alphai_vect_first_col"]]=theta_alphai_vect_first_col
	return(constant_theta_first_cols_l)
}

initialize_theta_alphai_star_vect_l=function(mpicl,gamma1,concat_local_E_l,num_genes_per_split,all_M_split_indices){
	Mtot=dim(concat_local_E_l$Eu_gammai_vect)[1]

	theta_alphai_vect=calc_theta_alphai_via_Eu_gammai(gamma1,concat_local_E_l$Eu_gammai_vect,Mtot)
	theta_alphai_star_vect=theta_alphai_vect-theta_alphai_vect
	theta_alphai_star_vect[,1]=(gamma1-1)+0.5
	theta_alphai_star_vect_list=lapply(c(1:length(num_genes_per_split)),function(ii){
		return(theta_alphai_star_vect[all_M_split_indices==ii,])
	})
	aa=clusterApply_or_envApply(mpicl,theta_alphai_star_vect_list,load_theta_alphai_star_vect_to_cluster_node)
	return(TRUE)
}

upload_lambda1=function(hyperparameter_list,mpicl,index_range_mappings,approximation_mode,gene_names_orig,gene_names){
	lambda1=hyperparameter_list$lambda1
	if(approximation_mode){
		gene_split_index=get_gene_cluster_index_fromrange_obj(index_range_mappings)
		ncores=length(index_range_mappings)
		### this  reordering only matters constant_noise_prior_pergene=TRUE, which is only used for testing purposes
		lambda1=lambda1[match(gene_names,gene_names_orig)]
		lambda1_list=list()
		for(ii in c(1:ncores)){
			lambda1_list[[ii]]=lambda1[gene_split_index==ii]
		}
		aa=clusterApply_or_envApply(mpicl,lambda1_list,load_lambda1_list_to_cluster_node)
	}else{
		clusterExport_or_envExport(mpicl,"lambda1",envir=environment())
	}
	hyperparameter_list$lambda1=lambda1
	return(hyperparameter_list)
}

initialize_contiguous_vectors=function(mpicl,global_E_list,mapping_mats_l,gene_dat_list,hyperparameter_list,snp_index_tbl,ncores,index_range_mappings){
	gene_indices_vec=snp_index_tbl[,gene_index]
	cluster_indices=snp_index_tbl[,cluster_index]

	gene_split_index=get_gene_cluster_index_fromrange_obj(index_range_mappings)
	num_genes_per_split=get_numgenes_persplit_fromrange_obj(index_range_mappings)

	concat_local_E_l=initialize_concat_local_E_lists(global_E_list=global_E_list,mapping_mat_expanded_ak=mapping_mats_l$mapping_mat_expanded_ak,gene_dat_list=gene_dat_list)
	Eu_alphai_vect=concat_local_E_l[["Eu_alphai_vect"]]
	Eu_kappa_vect=concat_local_E_l[["Eu_kappa_vect"]]
	Eu_kappa_pergene_vect=concat_local_E_l[["Eu_kappa_pergene_vect"]]
	###################
	initialize_theta_alphai_star_vect_l(mpicl,hyperparameter_list[["gamma1"]],concat_local_E_l,num_genes_per_split,cluster_indices)
	##############################
	initialize_alphai_kappa(mpicl=mpicl,gamma1=hyperparameter_list[["gamma1"]],gene_split_index,concat_local_E_l=concat_local_E_l,gene_indices_vec=gene_indices_vec,ncores=ncores)
	constant_theta_first_cols_l=inititalize_constant_theta_first_cols(hyperparameter_list[["gamma1"]],mapping_mat_expanded_ak=mapping_mats_l$mapping_mat_expanded_ak)
	Eu_alphai_list=get_alphasplits(mat_to_split=Eu_alphai_vect,index_range_mappings=index_range_mappings)
	a=clusterApply_or_envApply(mpicl,Eu_alphai_list,load_Eu_alphai_to_cluster_node)
	Eu_gammai_vect=concat_local_E_l[["Eu_gammai_vect"]]
	outl=list()
	outl[["Eu_gammai_vect"]]=Eu_gammai_vect
	outl[["constant_theta_first_cols_l"]]=constant_theta_first_cols_l
	return(outl)
}

upload_aux_data=function(mpicl,theta_root_list,hyperparameter_list,settings,global_E_list,ak_matcher,nu_matcher,calc_B,write_logfile){
	lambda1=hyperparameter_list$lambda1
	calc_B=calc_B
	write_logfile=write_logfile
	approximation_mode=settings$approximation_mode
	calc_L=settings$calc_L
	calc_B_rounds=settings$calc_B_rounds
	calc_L_rounds=settings$calc_L_rounds
	current_round=0
	theta_ak=theta_root_list$theta_ak
	theta_lambda=calc_theta_lambda(lambda1,global_E_list$Eu_lambda2[2],genespecific_lambda1=approximation_mode)
	theta_kappa=calc_theta_lambda(hyperparameter_list[["tau1"]],global_E_list$Eu_tau2[2])
	gamma1=hyperparameter_list[["gamma1"]]
	tau1=hyperparameter_list[["tau1"]]
	lambda1=hyperparameter_list[["lambda1"]]
	upload_global_Elist(mpicl,global_E_list)
	clusterExport_or_envExport(mpicl,c("tau1","nu_matcher","ak_matcher","theta_lambda","theta_kappa","theta_ak","gamma1","calc_L","calc_B","calc_L_rounds","calc_B_rounds","current_round","write_logfile","approximation_mode"),envir=environment())
	cl_len=get_cluster_size(mpicl)
	aa=clusterApply_or_envApply(mpicl,c(1:cl_len),load_node_nr_to_cluster)
}



get_numgenes_persplit=function(gene_names_splits){
	lapply(gene_names_splits,length)
}

get_numgenes_persplit_fromrange_obj=function(index_range_mappings){
	get_numgenes_persplit(get_gene_names_splits_fromrange_obj(index_range_mappings))
}

get_gene_cluster_index_fromrange_obj=function(index_range_mappings){
	get_gene_cluster_index(get_gene_names_splits_fromrange_obj(index_range_mappings))
}

get_gene_names_splits_fromrange_obj=function(index_range_mappings){
	lapply(index_range_mappings,rownames)
}

get_gene_cluster_index=function(gene_names_splits){
	num_genes_per_split=get_numgenes_persplit(gene_names_splits)
	gene_split_index=unlist(lapply(c(1:length(gene_names_splits)),function(i){rep(i,unlist(num_genes_per_split)[i])}))
	return(gene_split_index)
}

update_hyperparameter_list=function(gene_dat_list,hyperparameter_list,settings,D_mat){
	outl=initialize_lambda1nrho1_for_approx_mode(gene_dat_list=gene_dat_list,lambda1=hyperparameter_list[["lambda1"]],rho1=hyperparameter_list[["rho1"]],constant_noise_prior_pergene=settings$constant_noise_prior_pergene,approximation_mode=settings$approximation_mode)
	lambda1=outl[["lambda1"]]
	rho1=outl[["rho1"]]
	hyperparameter_list$lambda1Old=hyperparameter_list$lambda1
	hyperparameter_list$rho11Old=hyperparameter_list$rho1
	hyperparameter_list$lambda1=lambda1
	hyperparameter_list$rho1=rho1
	myq=dim(D_mat)[2]
	chi1_vec=hyperparameter_list$chi1
	if(length(chi1_vec)==1){
		chi1_vec=rep(chi1_vec,myq)
	}
	hyperparameter_list$chi1_vec=chi1_vec
	return(hyperparameter_list)
}

get_nu_names=function(hyperparameter_list){
	nu_names_extended=names(hyperparameter_list$c)
	if(is.null(nu_names_extended)){
		stop("hyperparameter_list$c is unnamed.")
	}
	if(nu_names_extended[1]!="all_id2920242"){
		stop("hyperparameter_list doesn't look like it was made by function set_hyperparameter_list.")
	}
	nu_names=nu_names_extended[2:length(nu_names_extended)]
	return(nu_names)
}

#' Runs variational updates.
#'
#' Runs variational updates according to the bagea model.
#' Data is provided by the object \code{gene_dat_list}, hyperparameter settings of the model is by provided by \code{hyperparameter_list}.
#' @param gene_dat_list
#' An R-object of class \code{bagea_observed_data} producted by \code{\link{prepare_observed_data}}
#' @param hyperparameter_list
#' An R-object of class bagea_hyperparameter_list produced by \code{\link{set_hyperparameter_list}}
#' @param calc_L
#' Boolean. Set true if the trace of L (bound on the full likelihood) should be calculated. 
#' @param ncores
#' Positive integer scalar. Number of parallel cores should be used.
#' @param nrounds
#' Positive integer scalar. Number of rounds the update should be run.
#' @param calc_L_rounds
#' Positive integer scalar. number of rounds before a new value of L should be calculated. nrounds has to be an integer multiple of  calc_L_rounds.
#' @param  approximation_mode
#' Boolean. Should approximation mode be used? Should be switched on if summary statistics are used.
#' @param ak_names
#' Character vector. Which undirected annotations should be used to construct C matrix? if  NULL, all undirected annotations available in gene_dat_list will be used.
#' @param  use_compiled
#' Boolean. Should compiled code be used? Set this to FALSE if compilation and linking failed.
#' @return
#' An R-object of class \code{bagea_output}; an \code{R}-list containing multiple objects.
#' \describe{\item{global variable estimates}{Expectations of the sufficient statistics of global variable estimates are provided. For gamma distributed variables they are given as follows: c(E[log(X)],E[X]). Gamma distributed variable outputs are: \code{Eu_ak}, \code{Eu_lambda2}, \code{Eu_tau2},  \code{Eu_upsilon_list}, \code{Eu_chi2j}.For multivariate normal the output format is list(E[X],E[XX^T]),Multivariate normally distributed variables  are: \code{Eu_omega_list}, \code{Eu_nu_list}.}
#' \item{\code{gene_dat_list}}{A named list containing all the relevant estimates per gene.}
#' \item{\code{Ltot}}{A trace of the data likelihood bound along the iterations. If not computed for the round in question, this is set to -Inf.}
#' \item{\code{D_mat}}{A matrix that maps the index of the meta-annotation group indices to the indices of the directed annotations.}
#' \item{\code{D_names_list}}{A list of names that contains to the meta-annotation names.}
#' \item{\code{D_names_list}}{A list of names that contains to the meta-annotation names.}
#' \item{selected global variable traces}{For particular global variabels, the whole trace of estimates across iterations are given. These are given as \code{Eu_nu1_trace}, \code{Eu_omega1_trace},\code{Eu_ak_trace}}
#' \item{\code{settings}}{A list of parameter settings with which prepare_annotation_mapping_mat was run.}}
#' @export
run_bagea=function(gene_dat_list=NULL,hyperparameter_list=NULL,calc_L=TRUE,ncores=NULL,nrounds=1000,calc_L_rounds=50,approximation_mode=TRUE,ak_names=NULL,use_compiled=TRUE){
	####################
	## some imput checks
	####################
	Eu_ak=NULL## init Eu_ak value
	calc_B_rounds=1 ## number of rounds before b estimates should be updated.
	calc_B=TRUE
	write_logfile=FALSE
	constant_noise_prior_pergene=FALSE
	gene_dat_list_settings=gene_dat_list$settings
	####################
	outl=extract_D_mat_info(gene_dat_list)
	D_mat=outl[["D_mat"]]
	D_names_list=outl[["D_names_list"]]
	nu_names=get_nu_names(hyperparameter_list)
	check_input(hyperparameter_list=hyperparameter_list,calc_L=calc_L,calc_B=calc_B,ncores=ncores,nrounds=nrounds,calc_L_rounds=calc_L_rounds,calc_B_rounds=calc_B_rounds,write_logfile=write_logfile,approximation_mode=approximation_mode,nu_names=nu_names,ak_names=ak_names)


	outl=get_matcher_vectors(gene_dat_list=gene_dat_list,ak_names=ak_names,nu_names=nu_names)
	ak_matcher=outl[["ak_matcher"]]
	nu_matcher=outl[["nu_matcher"]]
	nu_names=outl[["nu_names"]]
	ak_names=outl[["ak_names"]]
	rm(outl)
	if((length(hyperparameter_list$p)-1)!=length(nu_names)){
		print(length(hyperparameter_list$p))
		print(nu_names)
		print(length(nu_names))
		stop("length of hyperparameter vector p has to be one larger than nu_names.")
	}
	if((length(hyperparameter_list$c)-1)!=length(nu_names)){
		print(length(hyperparameter_list$c))
		print(nu_names)
		print(length(nu_names))
		stop("length of hyperparameter vector c has to be one larger than nu_names.")
	}
	outl=prepare_inputdata_lists(gene_dat_list,hyperparameter_list,nu_names,ak_names)
	gene_dat_list=outl[["gene_dat_list"]]
	global_E_list=outl[["global_E_list"]]
	input_settings=outl[["inputsettings"]]
	rm(outl)
	check_gene_dat_list(gene_dat_list=gene_dat_list)
	check_global_E_list(global_E_list)
	global_theta_star_list=list()
	print("checks done")	
	####################
	## END: some imput checks
	####################
	mysettings=get_settingsout_list_bagea(gene_dat_list_settings=input_settings,hyperparameter_list=hyperparameter_list,calc_L=calc_L,ncores=ncores,nrounds=nrounds,calc_L_rounds=calc_L_rounds,Eu_ak=global_E_list$Eu_ak,calc_B_rounds=calc_B_rounds,approximation_mode=approximation_mode,constant_noise_prior_pergene=constant_noise_prior_pergene,nu_names=nu_names,ak_names=ak_names)
	hyperparameter_list_updated=update_hyperparameter_list(gene_dat_list=gene_dat_list,hyperparameter_list=hyperparameter_list,settings=mysettings,D_mat=D_mat)
	theta_root_list=prepare_theta_root_list(my_t=dim(global_E_list$Eu_ak)[1],hyperparameter_list=hyperparameter_list_updated,D_mat=D_mat,nu_names=nu_names)
	g_root_list=prepare_g_root_list(theta_root_list)
	print("making cluster ..")
	mpicl=makePSOCKcluster_or_environment(ncores,useEnvIfOne=TRUE)
	###################
	print("uploading gene wise data ..")
	gene_names_orig=names(gene_dat_list)
	outl=split_and_upload(mpicl,gene_dat_list)
	gene_names_splits=outl[["gene_names_splits"]]
	cluster_indices=outl[["cluster_indices"]]
	num_genes_per_split=get_numgenes_persplit(gene_names_splits)
	gene_split_index=get_gene_cluster_index(gene_names_splits)
	rm(outl)
	gene_names=unlist(gene_names_splits)
	gene_dat_list=gene_dat_list[match(gene_names,gene_names_orig)]
	gene_indices_vec=get_gene_indices_vec(gene_dat_list)
	snp_index_tbl=data.table(gene_index=gene_indices_vec[,1],cluster_index=cluster_indices)
	index_range_mappings=initialize_index_range_mappings(mpicl=mpicl,gene_indices_vec=snp_index_tbl$gene_index,gene_names_splits=gene_names_splits)
	initialize_index_mats(mpicl=mpicl,snp_index_tbl=snp_index_tbl,index_range_mappings=index_range_mappings)
	print(Sys.time())
	hyperparameter_list_updated=upload_lambda1(hyperparameter_list=hyperparameter_list_updated,mpicl=mpicl,index_range_mappings=index_range_mappings,approximation_mode=approximation_mode,gene_names_orig=gene_names_orig,gene_names=gene_names)
	####################
	####################
	print("upload auxiliary data ..")
	upload_aux_data(mpicl=mpicl,theta_root_list=theta_root_list,hyperparameter_list=hyperparameter_list_updated,settings=mysettings,global_E_list=global_E_list,ak_matcher=ak_matcher,nu_matcher=nu_matcher,calc_B=calc_B,write_logfile=write_logfile)
	mapping_mats_l=process_mapping_mats(mpicl=mpicl,
		gene_dat_list=gene_dat_list,
		snp_index_tbl=snp_index_tbl,
		num_genes_per_split=get_numgenes_persplit_fromrange_obj(index_range_mappings),
		ak_matcher=ak_matcher,
		nu_matcher=nu_matcher,
		use_compiled=use_compiled
		)
	###################
	traces_l=initialize_traces(nrounds)
	##### split this out next
	print(Sys.time())
		print("preparing contiguous vectors..")
	outl=initialize_contiguous_vectors(mpicl=mpicl,
		global_E_list=global_E_list,
		mapping_mats_l=mapping_mats_l,
		gene_dat_list=gene_dat_list,
		hyperparameter_list=hyperparameter_list_updated,
		snp_index_tbl=snp_index_tbl,
		ncores=ncores,
		index_range_mappings=index_range_mappings)
	Eu_gammai_vect=outl[["Eu_gammai_vect"]]
	constant_theta_first_cols_l=outl[["constant_theta_first_cols_l"]]
	##############################
	## compute t(z)%*%Sigma^-1%*%z
	##############################	
	if(approximation_mode){
		print("in approximation mode: overwrite yty with t(z)%*%Sigma^-1z")
		clusterCall_or_envCall(mpicl,calculate_ZtSigma_invZ_list)
	}
	##############################
	## END: setup 
	##############################
	for(i in c(1:nrounds)){
		outl=update_params(
			i=i,
			mpicl=mpicl,
			Eu_ak_old=Eu_ak_old,
			global_E_list=global_E_list,
			mysettings=mysettings,
			theta_root_list=theta_root_list,
			global_theta_star_list=global_theta_star_list,
			mapping_mat_ak_nonnull_element_list=mapping_mats_l$mapping_mat_ak_nonnull_element_list,
			hyperparameter_list=hyperparameter_list_updated,
			mapping_mat_expanded_ak=mapping_mats_l$mapping_mat_expanded_ak,
			theta_alphai_vect_first_row=constant_theta_first_cols_l$theta_alphai_vect_first_col,
			D_mat=D_mat,
			theta_alphai_ak_first_col=constant_theta_first_cols_l$theta_alphai_ak_first_col,
			Eu_gammai_vect=Eu_gammai_vect,
			nu_names=nu_names)
		Eu_ak_old=outl[["Eu_ak_old"]]
		Eu_gammai_vect=outl[["Eu_gammai_vect"]]
		global_E_list=outl[["global_E_list"]]		
		local_calc_L= mysettings$calc_L && (((i %% mysettings$calc_L_rounds)==0) || i==1)
		if(local_calc_L){
			all_Ls_cur=calc_L_global(
				mpicl=mpicl,
				update_out_list=outl,
				mapping_mats_l=mapping_mats_l,
				g_root_list=g_root_list,
				theta_root_list=theta_root_list,
				hyperparameter_list=hyperparameter_list_updated,	
				theta_alphai_vect_first_row=constant_theta_first_cols_l$theta_alphai_vect_first_col,
				D_mat=D_mat,
				approximation_mode=approximation_mode
			)			
		}else{
			all_Ls_cur=NULL
		}
		traces_l=update_traces(
			traces_l=traces_l,
			i=i,
			all_Ls_cur=all_Ls_cur,
			global_E_list=global_E_list)
	}
	#Rprof(NULL)
	out_dat_list=pull_out_dat(
		mpicl=mpicl,
		update_out_list=outl,
		traces_l=traces_l,
		index_range_mappings=index_range_mappings,
		gene_names=gene_names_orig,
		settings=mysettings,
		D_mat=D_mat,
		D_names_list=D_names_list
	)
	rm_cluster(mpicl)
	out_dat_list$omegahat_resc=calc_rescaled_omega_estimates(out_dat_list,use_Fnu_scaling=TRUE,supress_check=TRUE)
	out_dat_list$all_Ls=NULL # remove for release
	out_dat_list$omegahat_resc=NULL # remove for release
	class(out_dat_list)=c("bagea_output",class(out_dat_list))
	return(out_dat_list)
}