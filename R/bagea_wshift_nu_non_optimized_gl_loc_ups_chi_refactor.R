simple_subsub_chi=function(upsilon_list,chi1_vec,zeta1_vec,zeta2_vec){
	myq=length(chi1_vec)
	Eu_upsilon_list=list()
	for(j in c(1:myq)){
		Eu_upsilon_list[[j]]=cbind(log(upsilon_list[[j]]),upsilon_list[[j]])
	}
	theta_chi2j=calc_theta_chi2j(zeta1_vec,zeta2_vec)
	Eu_chi2j=calc_Eu_chi2j_via_theta(theta_chi2j)
	for(round in c(1:200)){
		theta_upsilon_star_list=list()
		theta_chi2j=calc_theta_chi2j(zeta1_vec,zeta2_vec)
		theta_upsilon_chi2j=calc_theta_upsilon_chi2j_via_Eu(Eu_upsilon_list,chi1_vec)
		theta_chi2j_star=theta_chi2j+theta_upsilon_chi2j
		Eu_chi2j=calc_Eu_chi2j_via_theta(theta_chi2j_star)
	}
	outl=list()
	outl[["Eu_upsilon_list"]]=Eu_upsilon_list
	outl[["Eu_chi2j"]]=Eu_chi2j
	return(outl)
}

simple_sub_chi=function(omega,D_mat,chi1_vec,zeta1_vec,zeta2_vec,nrounds=200){
	Eu_omega_list=list()
	Eu_omega_list[[1]]=omega
	Eu_omega_list[[2]]=omega%*%t(omega)
	D_mat_check=apply(D_mat,2,is_ascending_integer_vec)
	myq=dim(D_mat)[2]
	h_vec=apply(D_mat,2,max)
	theta_chi2j=calc_theta_chi2j(zeta1_vec,zeta2_vec)
	Eu_chi2j=calc_Eu_chi2j_via_theta(theta_chi2j)
	theta_upsilon_list=calc_theta_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec=h_vec)
	Eu_upsilon_list=calc_Eu_upsilon_via_theta_list(theta_upsilon_list)	
	for(round in c(1:nrounds)){
		theta_upsilon_star_list=list()
		for(myj in c(1:myq)){
			theta_omega_upsilon_j=calc_theta_omegas_upsilon_via_Eu(Eu_omega_list,Eu_upsilon_list,D_mat,myj)
			theta_upsilon_j=theta_upsilon_list[[myj]]
			theta_upsilon_star_j=theta_omega_upsilon_j+theta_upsilon_j
			theta_upsilon_star_list[[myj]]=theta_upsilon_star_j
			Eu_upsilon_j=calc_Eu_upsilon_via_theta(theta_upsilon_star_j)
			Eu_upsilon_list[[myj]]=Eu_upsilon_j
		}
		theta_chi2j=calc_theta_chi2j(zeta1_vec,zeta2_vec)
		theta_upsilon_chi2j=calc_theta_upsilon_chi2j_via_Eu(Eu_upsilon_list,chi1_vec)
		theta_chi2j_star=theta_chi2j+theta_upsilon_chi2j		
		Eu_chi2j=calc_Eu_chi2j_via_theta(theta_chi2j_star)
		theta_upsilon_list=calc_theta_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec=h_vec)
	}
	outl=list()
	outl[["Eu_upsilon_list"]]=Eu_upsilon_list
	outl[["Eu_chi2j"]]=Eu_chi2j
	return(outl)
}

get_mapping_mat_colnames=function(gene_dat_list){
	input_is_bagea_output=sum(class(gene_dat_list)=="bagea_output")!=0
	input_is_bagea_initalized_input=sum(class(gene_dat_list)=="bagea_initialized_gene_dat_list")!=0
	if(!input_is_bagea_output && !input_is_bagea_initalized_input){
		return(colnames(gene_dat_list$obs_list[[1]]$mapping_mat))
	}
	if(input_is_bagea_output){
		return(colnames(gene_dat_list$gene_dat_list[[1]]$Obs_list$mapping_mat))
	}
}

prepare_nonoptim_outl=function(gene_dat_list,mysettings,global_E_list,all_Ls,ssMat,D_mat,D_names_list){
	out_list=list()
	out_list[["gene_dat_list"]]=gene_dat_list
	out_list[["settings"]]=mysettings
	out_list[["Eu_ak"]]=global_E_list$Eu_ak
	out_list[["Eu_upsilon_list"]]=global_E_list$Eu_upsilon_list
	out_list[["Eu_omega_list"]]=global_E_list$Eu_omega_list
	out_list[["Eu_nu_list"]]=global_E_list$Eu_nu_list
	out_list[["Eu_tau2"]]=global_E_list$Eu_tau2
	out_list[["Eu_lambda2"]]=global_E_list$Eu_lambda2
	out_list[["Eu_chi2j"]]=global_E_list$Eu_chi2j
	out_list[["all_Ls"]]=all_Ls
	ssMat[1,]=-Inf
	out_list[["ssMat"]]=ssMat
	out_list[["D_mat"]]=D_mat
	out_list[["D_names_list"]]=D_names_list
	class(out_list)=c("bagea_output",class(out_list))
	return(out_list)
}

check_E_gamma_el=function(E_el,strname){
	if(is.null(dim(E_el))){
		if(length(E_el)!=2){
			stop(paste(strname,": not correct length."))
		}
		E_el=t(matrix(E_el))
	}
	check_matrix(E_el,strname)
	check_pos(E_el[,2],paste(strname,"[,2]",sep=""))
	violates_jensen=sum(E_el[,1]>log(E_el[,2]))
	if(violates_jensen){
		stop(paste(strname," : jensen ineq E[log(x)]<=log(E[x]) violated."))
	}
}

check_E_gaussian_el=function(E_el,strname){
	mu=E_el[[1]]
	e2=E_el[[2]]
	myv=e2-mu%*%t(mu)
	check_matrix(mu,paste(strname,": mu not formatted correctly."))
	check_matrix(e2,paste(strname,": e2 not formatted correctly."))
	smallest_eigenval=tail(eigen(myv)$values,1)	
	if(smallest_eigenval< -1E-9){
		stop(paste(strname,": variance not semidefinite."))
	}	
}

check_global_E_list=function(global_E_list){
	expected_E_names=c("Eu_lambda2","Eu_tau2","Eu_chi2j","Eu_upsilon_list","Eu_omega_list","Eu_nu_list","Eu_ak")
	is_gaussian=rep(FALSE,length(expected_E_names))
	is_gaussian[c(5,6)]=TRUE
	correct_nrs=length(union(names(global_E_list),expected_E_names))==length(intersect(names(global_E_list),expected_E_names))
	if(!correct_nrs){
		stop("global E_list not right number of elements")
	}
	for(i in c(1:length(expected_E_names))){
		cur_name=expected_E_names[i]
		if(is_gaussian[i]){			
			check_E_gaussian_el(global_E_list[[cur_name]],cur_name)
		}else{
			if(cur_name=="Eu_upsilon_list"){
				lapply(global_E_list[[cur_name]],check_E_gamma_el,cur_name)
			}else{
				check_E_gamma_el(global_E_list[[cur_name]],cur_name)
			}
		}
	}
}

extract_D_mat_info=function(gene_dat_list){
	D_mat=gene_dat_list$D_mat
	if(is.null(D_mat)){
		stop("D_mat not allowed to be NULL.")
	}
	####################
	## check D_mat
	####################
	D_mat_check=apply(D_mat,2,is_ascending_integer_vec)
	####################
	## END: check D_mat
	####################
	D_names_list=gene_dat_list$D_names_list
	outl=list()
	outl[["D_mat"]]=D_mat
	outl[["D_names_list"]]=D_names_list
	return(outl)
}

prepare_inputdata_lists=function(gene_dat_list,hyperparameter_list,nu_names,ak_names){
	inputsettings=gene_dat_list$settings
	input_is_bagea_output=sum(class(gene_dat_list)=="bagea_output")!=0
	input_is_bagea_initalized_input=sum(class(gene_dat_list)=="bagea_initialized_gene_dat_list")!=0
	if(!input_is_bagea_output && !input_is_bagea_initalized_input){
		gene_dat_list=initialize_parameters(gene_dat_list,hyperparameter_list,nu_names=nu_names,ak_names=ak_names)
		global_E_list=attributes(gene_dat_list)$globalE_list
		outl=list()
		outl[["gene_dat_list"]]=gene_dat_list
		outl[["global_E_list"]]=global_E_list
	}
	if(input_is_bagea_output){
		outl=reformat_bagea_output2input(gene_dat_list,ak_names)
	}
	outl[["inputsettings"]]=inputsettings
	return(outl)
}

prepare_g_root_list=function(theta_root_list){
	g_root_list=list()
	g_lambda2=calc_g_lambda2_via_theta(theta_root_list[["theta_lambda2"]])
	g_tau2=calc_g_tau2_via_theta(theta_root_list[["theta_tau2"]])
	g_chi2j=calc_g_chi2j_via_theta(theta_root_list[["theta_chi2j"]])
	g_nu=calc_g_nu_via_theta(theta_root_list[["theta_nu_list"]])
	g_ak=calc_g_ak_via_theta(theta_root_list[["theta_ak"]])
	g_root_list[["g_lambda2"]]=g_lambda2
	g_root_list[["g_tau2"]]=g_tau2
	g_root_list[["g_ak"]]=g_ak
	g_root_list[["g_chi2j"]]=g_chi2j
	g_root_list[["g_nu_list"]]=g_nu
	return(g_root_list)

}

prepare_theta_root_list=function(my_t,hyperparameter_list,D_mat,nu_names){
	myq=dim(D_mat)[2]
	h_vec=apply(D_mat,2,max)
	chi1_vec=hyperparameter_list$chi1
	if(length(chi1_vec)==1){
		chi1_vec=rep(chi1_vec,myq)
	}
	zeta1_vec=hyperparameter_list$zeta1
	zeta2_vec=hyperparameter_list$zeta2
	if(length(zeta1_vec)==1){
		zeta1_vec=rep(zeta1_vec,myq)
		zeta2_vec=rep(zeta2_vec,myq)
	}
	theta_lambda2=calc_theta_lambda2(hyperparameter_list[["rho1"]],hyperparameter_list[["rho2"]])
	theta_tau2=calc_theta_tau2(hyperparameter_list[["xi1"]],hyperparameter_list[["xi2"]])
	theta_ak=kronecker(t(calc_theta_lambda(hyperparameter_list[["phi1"]],hyperparameter_list[["phi2"]])),rep(1,my_t))
	theta_nu_list=calc_theta_nu(p=hyperparameter_list$p,q=(1+length(nu_names)),myc=hyperparameter_list$c)
	theta_chi2j=calc_theta_chi2j(zeta1_vec,zeta2_vec)
	theta_root_list=list()
	theta_root_list[["theta_lambda2"]]=theta_lambda2
	theta_root_list[["theta_tau2"]]=theta_tau2
	theta_root_list[["theta_ak"]]=theta_ak
	theta_root_list[["theta_chi2j"]]=theta_chi2j
	theta_root_list[["theta_nu_list"]]=theta_nu_list
	return(theta_root_list)
}

check_Eu_omega_list_fromoutput=function(gene_dat_list,Eu_omega_list){
	my_s=dim(gene_dat_list[[1]][["Obs_list"]][["mapping_shift_mat"]])[2]
	if(input_is_bagea_output){
		if(sum(dim(Eu_omega_list[[1]])!=c(my_s,1))!=0){
			stop("Eu_omega not the right dimension")
		}
		if(sum(dim(Eu_omega_list[[2]])!=c(my_s,my_s))!=0){
			stop("Eu_omega not the right dimension")
		}
		check_matrix(Eu_omega_list[[2]],"Eu_omega_list[[2]]")
	}
}

reformat_bagea_output2input=function(bagea_output,ak_names){
	global_E_list=list()
	global_E_list[["Eu_lambda2"]]=bagea_output$Eu_lambda2
	global_E_list[["Eu_tau2"]]=bagea_output$Eu_tau2
	if(dim(bagea_output$Eu_ak)[1]!=length(ak_names) || sum(rownames(bagea_output$Eu_ak)!=ak_names)!=0){
		mymatcher=match(ak_names,rownames(bagea_output$Eu_ak))
		if(sum(is.na(mymatcher))!=0){
			stop("erraomscomsdc")
		}
		global_E_list[["Eu_ak"]]=bagea_output$Eu_ak[mymatcher,]
		}else{
		global_E_list[["Eu_ak"]]=bagea_output$Eu_ak
	}
	global_E_list[["Eu_upsilon_list"]]=bagea_output$Eu_upsilon_list
	global_E_list[["Eu_chi2j"]]=bagea_output$Eu_chi2j
	if(is.null(global_E_list[["Eu_upsilon_list"]])){
		stop("errmwd-w")
	}
	global_E_list[["Eu_omega_list"]]=bagea_output$Eu_omega_list
	global_E_list[["Eu_nu_list"]]=bagea_output$Eu_nu_list
	gene_dat_list=bagea_output$gene_dat_list
	for(i in c(1:length(gene_dat_list))){
		gene_dat_list[[i]]$L_list=list()
	}
	check_Eu_omega_list_fromoutput(gene_dat_list,global_E_list$Eu_omega_list)
	outl=list()
	outl[["gene_dat_list"]]=gene_dat_list
	outl[["global_E_list"]]=global_E_list
	return(outl)
}


get_settingsout_list_bagea=function(gene_dat_list_settings,hyperparameter_list,calc_L,ncores,nrounds,calc_L_rounds,Eu_ak,calc_B_rounds,approximation_mode,constant_noise_prior_pergene,nu_names,ak_names){
	mysettings=list()
	mysettings[["hyperparameter_list"]]=hyperparameter_list
	mysettings[["observed_data_settings"]]=gene_dat_list_settings
	mysettings[["calc_L"]]=calc_L
	mysettings[["ncores"]]=ncores
	mysettings[["nrounds"]]=nrounds
	mysettings[["calc_L_rounds"]]=calc_L_rounds
	mysettings[["Eu_ak_init"]]=Eu_ak
	mysettings[["calc_B_rounds"]]=calc_B_rounds
	mysettings[["approximation_mode"]]=approximation_mode
	mysettings[["constant_noise_prior_pergene"]]=constant_noise_prior_pergene
	mysettings[["nu_names"]]=nu_names
	mysettings[["ak_names"]]=ak_names
	return(mysettings)
}

update_Eu_alphai=function(gene_dat_list,Eu_ak,ak_matcher,hyperparameter_list,Eu_omega_list,Eu_nu_list){
	###########################
	## update Eu_alphai
	###########################	
	for(j in c(1:length(gene_dat_list))){
		cur_e=gene_dat_list[[j]]$E_list
		cur_obs=gene_dat_list[[j]]$Obs_list
		Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat=cur_obs$mapping_mat[,ak_matcher])
		theta_alphai=calc_theta_alphai_via_Eu_gammai_Eu_kappa(hyperparameter_list[["gamma1"]],Eu_gammai,cur_e$Eu_kappa)
		cur_mapping_mat=cbind(TRUE,cur_obs$mapping_mat[,match(nu_names,colnames(cur_obs$mapping_mat))])
		theta_bi_alphai=calc_theta_bi_alphai_via_Eu_wshift_nu(m=dim(theta_alphai)[1],Eu_b_list=cur_e$Eu_b_list,mapping_shift_mat=cur_obs$mapping_shift_mat,Eu_omega_list=Eu_omega_list,mapping_mat=cur_mapping_mat,Eu_nu_list=Eu_nu_list)
		theta_alphai_star=theta_alphai+theta_bi_alphai
		Eu_alphai=calc_Eu_alphai_via_theta(theta_alphai_star)
		gene_dat_list[[j]]$E_list$Eu_alphai=Eu_alphai
		###########################
		gene_dat_list[[j]]$theta_list$theta_alphai_star=theta_alphai_star
		gene_dat_list[[j]]$theta_list$theta_alphai=theta_alphai
	}
	###########################
	## END: update Eu_alphai
	###########################
	###########################
	## update relevant L_s: L_alphai;L_b
	###########################
	for(j in c(1:length(gene_dat_list))){
		cur_theta=gene_dat_list[[j]]$theta_list
		cur_e=gene_dat_list[[j]]$E_list
		cur_obs=gene_dat_list[[j]]$Obs_list
		Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat=cur_obs$mapping_mat[,ak_matcher])
		theta_alphai=calc_theta_alphai_via_Eu_gammai_Eu_kappa(hyperparameter_list[["gamma1"]],Eu_gammai,cur_e$Eu_kappa)
		Eu_kappa_expanded=rep(1,dim(Eu_gammai)[1])%*%t(cur_e$Eu_kappa)
		g_alphai_pa=calc_g_alphai_via_Eu_gammai_Eu_kappa(hyperparameter_list[["gamma1"]],Eu_gammai,Eu_kappa_expanded)
		L_alphai=sum((theta_alphai-cur_theta$theta_alphai_star)*cur_e$Eu_alphai)+sum(g_alphai_pa)-sum(calc_g_alphai_via_theta(cur_theta$theta_alphai_star))
		gene_dat_list[[j]][["L_list"]][["L_alphai"]]=L_alphai
		cur_mapping_mat=cbind(TRUE,cur_obs$mapping_mat[,match(nu_names,colnames(cur_obs$mapping_mat))])
		theta_b=calc_theta_bi_via_Eu_wshift_nu(cur_e$Eu_alphai,Eu_omega_list,cur_obs$mapping_shift_mat,Eu_nu_list,cur_mapping_mat)
		g_b_pa=calc_g_b_via_Eu_wshift_nu(Eu_alphai=cur_e$Eu_alphai,Eu_omega_list=Eu_omega_list,V=cur_obs$mapping_shift_mat,F=cur_mapping_mat,Eu_nu_list=Eu_nu_list)
		L_b=sum((theta_b[[1]]-cur_theta$theta_b_star_list[[1]])*cur_e$Eu_b_list[[1]])+sum((theta_b[[2]]-cur_theta$theta_b_star_list[[2]])*cur_e$Eu_b_list[[2]])+g_b_pa-calc_g_b_via_theta(cur_theta$theta_b_star_list)
		gene_dat_list[[j]][["L_list"]][["L_b"]]=L_b[1]
	}
	###########################
	## END: update relevant L_s: L_alphai;L_b(offspring)
	###########################
	return(gene_dat_list)
}

update_Eu_kappai=function(i,gene_dat_list,Eu_tau2,hyperparameter_list,ak_matcher,Eu_ak){
		###########################
		## update Eu_kappai
		###########################		
		theta_kappa=calc_theta_kappaj_via_Eu_tau2(hyperparameter_list[["tau1"]],Eu_tau2)
		for(j in c(1:length(gene_dat_list))){
			cur_obs=gene_dat_list[[j]]$Obs_list
			cur_e=gene_dat_list[[j]]$E_list
			Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat=cur_obs$mapping_mat[,ak_matcher])
			theta_alphai_kappa=calc_theta_alphai_kappa_via_Eu_gammai_Eu_alphai(gamma1=hyperparameter_list[["gamma1"]],Eu_gammai=Eu_gammai,Eu_alphai=cur_e$Eu_alphai)
			theta_kappa_star=theta_alphai_kappa+theta_kappa
			Eu_kappa=calc_Eu_kappa_via_theta(theta_kappa_star)
			gene_dat_list[[j]]$E_list$Eu_kappa=Eu_kappa
			###########################
			gene_dat_list[[j]]$theta_list$theta_kappa=theta_kappa
			gene_dat_list[[j]]$theta_list$theta_kappa_star=theta_kappa_star
		}		
		###########################
		## END: update Eu_kappai
		###########################
		###########################
		## update relevant L_s: L_kappai;L_alphai(offspring)
		###########################
		theta_kappa=calc_theta_kappaj_via_Eu_tau2(hyperparameter_list[["tau1"]],Eu_tau2)
		for(j in c(1:length(gene_dat_list))){
			print(j)
			cur_theta=gene_dat_list[[j]]$theta_list
			cur_e=gene_dat_list[[j]]$E_list
			cur_obs=gene_dat_list[[j]]$Obs_list
			Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat=cur_obs$mapping_mat[,ak_matcher])
			theta_alphai=calc_theta_alphai_via_Eu_gammai_Eu_kappa(hyperparameter_list[["gamma1"]],Eu_gammai,cur_e$Eu_kappa)
			Eu_kappa_expanded=rep(1,dim(Eu_gammai)[1])%*%t(cur_e$Eu_kappa)
			g_alphai_pa=calc_g_alphai_via_Eu_gammai_Eu_kappa(hyperparameter_list[["gamma1"]],Eu_gammai,Eu_kappa_expanded)
			if(i==1){
				cur_theta$theta_alphai_star=theta_alphai
			}
			L_alphai=sum((theta_alphai-cur_theta$theta_alphai_star)*cur_e$Eu_alphai)+sum(g_alphai_pa)-sum(calc_g_alphai_via_theta(cur_theta$theta_alphai_star))
			theta_kappa=calc_theta_kappaj_via_Eu_tau2(hyperparameter_list[["tau1"]],Eu_tau2)
			g_kappa_pa=calc_g_kappa_via_Eu_tau2(hyperparameter_list[["tau1"]],Eu_tau2)
			L_kappa=sum((theta_kappa-cur_theta$theta_kappa_star)*cur_e$Eu_kappa)+g_kappa_pa-calc_g_kappa_via_theta(cur_theta$theta_kappa_star)
			gene_dat_list[[j]][["L_list"]][["L_kappa"]]=L_kappa
			gene_dat_list[[j]][["L_list"]][["L_alphai"]]=L_alphai
		}	
		###########################
		## END: update relevant L_s: L_kappai;L_alphai(offspring)
		###########################
		return(gene_dat_list)
}



update_Eu_bj_Eu_lambdaj=function(gene_dat_list,Eu_lambda2,lambda1,approximation_mode,Eu_omega_list,Eu_nu_list){
		###########################
		## update Eu_bj,Eu_lambdaj
		###########################
		for(j in c(1:length(gene_dat_list))){
			cur_obs=gene_dat_list[[j]]$Obs_list
			cur_e=gene_dat_list[[j]]$E_list
			if(approximation_mode){	
				lambda1j=lambda1[j]
				theta_lambda=calc_theta_lambda(lambda1j,Eu_lambda2[2])
				theta_z_lambda=calc_theta_z_lambda_via_Eu(myn=cur_obs$n,Eu_b_list=cur_e$Eu_b_list,ytX=cur_obs$ytX,XtX=cur_obs$XtX,yty=cur_obs$yty)
				theta_lambda_star=theta_lambda+theta_z_lambda
			}else{
				(print("ss"))
				print(j)
				theta_lambda=calc_theta_lambda(lambda1,Eu_lambda2[2])
				theta_y_lambda=calc_theta_y_lambda_via_Eu(myn=cur_obs$n,Eu_b_list=cur_e$Eu_b_list,ytX=cur_obs$ytX,XtX=cur_obs$XtX,yty=cur_obs$yty)
				theta_lambda_star=theta_lambda+theta_y_lambda
			}
			Eu_lambda=calc_Eu_lambda_via_theta(theta_lambda_star)
			cur_mapping_mat=cbind(TRUE,cur_obs$mapping_mat[,match(nu_names,colnames(cur_obs$mapping_mat))])
			theta_b=calc_theta_bi_via_Eu_wshift_nu(cur_e$Eu_alphai,Eu_omega_list,cur_obs$mapping_shift_mat,Eu_nu_list,cur_mapping_mat)
			theta_b_star_list=list()
			if(approximation_mode){
				theta_z_b=calc_theta_z_b_via_Eu(Eu_lambda=Eu_lambda,ytX=cur_obs$ytX,XtX=cur_obs$XtX)
				theta_b_star_list[[1]]=theta_b[[1]]+theta_z_b[[1]]
				theta_b_star_list[[2]]=theta_b[[2]]+theta_z_b[[2]]
			}else{
				theta_y_b=calc_theta_y_b_via_Eu(Eu_lambda=Eu_lambda,ytX=cur_obs$ytX,XtX=cur_obs$XtX)
				theta_b_star_list[[1]]=theta_b[[1]]+theta_y_b[[1]]
				theta_b_star_list[[2]]=theta_b[[2]]+theta_y_b[[2]]
			}
			Eu_b_list=calc_Eu_b_via_theta(theta_b_star_list)	
			######################
			gene_dat_list[[j]][["E_list"]][["Eu_b_list"]]=Eu_b_list
			gene_dat_list[[j]][["E_list"]][["Eu_lambda"]]=Eu_lambda
			#######################
			gene_dat_list[[j]][["theta_list"]][["theta_lambda"]]=theta_lambda
			gene_dat_list[[j]][["theta_list"]][["theta_lambda_star"]]=theta_lambda_star
			gene_dat_list[[j]][["theta_list"]][["theta_b_star_list"]]=theta_b_star_list
			if(approximation_mode){			
				gene_dat_list[[j]][["theta_list"]][["theta_b_z_or_y_list"]]=theta_z_b
			}else{
				gene_dat_list[[j]][["theta_list"]][["theta_b_z_or_y_list"]]=theta_y_b
			}
			gene_dat_list[[j]][["theta_list"]][["theta_b"]]=theta_b
		}
		###########################
		## END: update Eu_bj,Eu_lambdaj
		###########################
		###########################
		## update relevant Ls:L_obs;L_b;L_lambda;
		###########################
		## to update all Ls that are relevant from the updated variables:
		## to be updated:
		##		L_obs
		##		L_b
		##		L_lambda
		for(j in c(1:length(gene_dat_list))){
			cur_theta=gene_dat_list[[j]]$theta_list
			cur_e=gene_dat_list[[j]]$E_list
			cur_obs=gene_dat_list[[j]]$Obs_list
			###### calculate L_b ######
			theta_b=cur_theta$theta_b
			cur_mapping_mat=cbind(TRUE,cur_obs$mapping_mat[,match(nu_names,colnames(cur_obs$mapping_mat))])
			g_b_pa=calc_g_b_via_Eu_wshift_nu(Eu_alphai=cur_e$Eu_alphai,Eu_omega_list=Eu_omega_list,V=cur_obs$mapping_shift_mat,F=cur_mapping_mat,Eu_nu_list=Eu_nu_list)
			L_b=sum((theta_b[[1]]-cur_theta$theta_b_star_list[[1]])*cur_e$Eu_b_list[[1]])+sum((theta_b[[2]]-cur_theta$theta_b_star_list[[2]])*cur_e$Eu_b_list[[2]])+g_b_pa-calc_g_b_via_theta(cur_theta$theta_b_star_list)
			###### calculate L_b ######	
			if(approximation_mode){
				L_z_obs=calc_L_z_observed(cur_obs$ytX,cur_obs$XtX,cur_e$Eu_b_list,cur_e$Eu_lambda,cur_obs$n)
			}else{
				L_y_obs=calc_L_y_observed(cur_obs$yty,cur_obs$ytX,cur_obs$XtX,cur_e$Eu_b_list,cur_e$Eu_lambda,cur_obs$n)
			}
			if(approximation_mode){	
				lambda1j=lambda1[j]
				theta_lambda=calc_theta_lambda(lambda1j,Eu_lambda2[2])
				g_lambda_pa=calc_g_lambda_via_Eu_lambda2(lambda1j,Eu_lambda2)
			}else{
				theta_lambda=calc_theta_lambda(lambda1,Eu_lambda2[2])
				g_lambda_pa=calc_g_lambda_via_Eu_lambda2(lambda1,Eu_lambda2)
			}
			L_lambda=sum((theta_lambda-cur_theta$theta_lambda_star)*cur_e$Eu_lambda)+g_lambda_pa-calc_g_lambda_via_theta(cur_theta$theta_lambda_star)
			#######################
			if(approximation_mode){
				gene_dat_list[[j]][["L_list"]][["L_z_obs"]]=L_z_obs[1]
			}else{
				gene_dat_list[[j]][["L_list"]][["L_y_obs"]]=L_y_obs[1]
			}		
			gene_dat_list[[j]][["L_list"]][["L_lambda"]]=L_lambda
			gene_dat_list[[j]][["L_list"]][["L_b"]]=L_b[1]
		}
		###########################
		## END: update relevant Ls:L_obs;L_b;L_lambda;
		###########################
		return(gene_dat_list)
}

reassemble_Eu_b_list=function(gene_dat_list){
	for (i in c(1:length(gene_dat_list))){
		Eu_b_list=gene_dat_list[[i]][["E_list"]][["Eu_b_list"]]
		new_list=list()
		new_list[["Eu_b"]]=Eu_b_list[["Eu_b_1"]]
		if(is.null(Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt)){
			## use formula : 'Eu_b2_formula_expanded' with XtX_sqrtAug_interimRes_sqrt set to 0.
			new_list[["Eu_b2"]]=0.5*(diag(Eu_b_list$theta_M_inv$diagDprime_alphainv)+ Eu_b_list$Eu_b_1 %*% t(Eu_b_list$Eu_b_1))
		}else{	
			new_list[["Eu_b2"]]=0.5*(diag(Eu_b_list$theta_M_inv$diagDprime_alphainv)-Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt%*%t(Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt)) + Eu_b_list$Eu_b_1 %*% t(Eu_b_list$Eu_b_1)
		}
		gene_dat_list[[i]][["E_list"]][["Eu_b_list"]]=new_list
	}
	return(gene_dat_list)
}

reassemble_XtX=function(gene_dat_list){
	for (i in c(1:length(gene_dat_list))){
		XtX_sqrt=gene_dat_list[[i]][["Obs_list"]][["XtX_sqrt"]]
		XtX=XtX_sqrt%*%t(XtX_sqrt)
		a=sum(svd(gene_dat_list[[i]][["Obs_list"]]$XtX_sqrt)$d^2)
		b=sum(gene_dat_list[[i]][["Obs_list"]]$eig_values)
		b-a
		mym=dim(XtX)[1]
		if(!(mym>1)){
			stop("errmv3-c")
		}
		XtX=XtX+diag(mym)*(b-a)/mym
		gene_dat_list[[i]][["Obs_list"]][["XtX"]]=XtX
	}
	return(gene_dat_list)
}

run_globalE_list=function(i,L_global_list,global_E_list,gene_dat_list,ak_matcher,hyperparameter_list,approximation_mode,lambda1,D_mat,theta_root_list){
	theta_ak=theta_root_list$theta_ak
	theta_lambda2=theta_root_list$theta_lambda2
	theta_tau2=theta_root_list$theta_tau2
	theta_chi2j=theta_root_list$theta_chi2j
	theta_nu_list=theta_root_list$theta_nu_list
	######
	Eu_ak=global_E_list$Eu_ak
	Eu_upsilon_list=global_E_list$Eu_upsilon_list
	Eu_lambda2=global_E_list$Eu_lambda2
	Eu_tau2=global_E_list$Eu_tau2
	Eu_chi2j=global_E_list$Eu_chi2j
	Eu_omega_list=global_E_list$Eu_omega_list
	Eu_nu_list=global_E_list$Eu_nu_list
	######
	myq=dim(D_mat)[2]
	h_vec=apply(D_mat,2,max)
	chi1_vec=hyperparameter_list$chi1
	if(length(chi1_vec)==1){
		chi1_vec=rep(chi1_vec,myq)
	}
	ssvec=rep(0,6)
	outl=update_Eu_ak(theta_ak=theta_ak,gene_dat_list=gene_dat_list,Eu_ak=Eu_ak,hyperparameter_list=hyperparameter_list,L_global_list=L_global_list,ak_matcher=ak_matcher)
	L_global_list=outl[["L_global_list"]]
	Eu_ak=outl[["Eu_ak"]]
	gene_dat_list=outl[["gene_dat_list"]]
	rm(outl)
	aa=get_global_L_sum(L_global_list)+get_local_L_sum(gene_dat_list)
	ssvec[1]=aa
	print(aa)
	print("update Eu_lambda2")
	###########################
	## update Eu_lambda2
	###########################
	outl=update_Eu_lambda2(theta_lambda2=theta_lambda2,gene_dat_list=gene_dat_list,lambda1=lambda1,approximation_mode=approximation_mode,L_global_list=L_global_list)
	L_global_list=outl[["L_global_list"]]
	Eu_lambda2=outl[["Eu_lambda2"]]
	gene_dat_list=outl[["gene_dat_list"]]
	###########################
	## END: update relevant L_s: L_lambda2;L_lambda
	###########################	
	aa=get_global_L_sum(L_global_list)+get_local_L_sum(gene_dat_list)
	ssvec[2]=aa
	print(aa)
	print("update Eu_tau2")
	outl=update_Eu_tau2(i=i,gene_dat_list=gene_dat_list,L_global_list=L_global_list,hyperparameter_list=hyperparameter_list,theta_tau2=theta_tau2)
	L_global_list=outl[["L_global_list"]]
	Eu_tau2=outl[["Eu_tau2"]]
	gene_dat_list=outl[["gene_dat_list"]]
	###########################
	## update Eu_tau2
	###########################
	aa=get_global_L_sum(L_global_list)+get_local_L_sum(gene_dat_list)
	ssvec[3]=aa
	print(aa)
	print("update Eu_nu_list")
	outl=update_Eu_nu_list(gene_dat_list=gene_dat_list,L_global_list=L_global_list,Eu_omega_list=Eu_omega_list,hyperparameter_list=hyperparameter_list)
	L_global_list=outl[["L_global_list"]]
	Eu_nu_list=outl[["Eu_nu_list"]]
	gene_dat_list=outl[["gene_dat_list"]]
	###########################
	## update Eu_omega_list
	###########################
	outl=update_Eu_omega_list(gene_dat_list=gene_dat_list,L_global_list=L_global_list,hyperparameter_list=hyperparameter_list,Eu_nu_list=Eu_nu_list,Eu_upsilon_list=Eu_upsilon_list,D_mat=D_mat)
	L_global_list=outl[["L_global_list"]]
	Eu_omega_list=outl[["Eu_omega_list"]]
	gene_dat_list=outl[["gene_dat_list"]]	
	theta_omega_star_list=outl[["theta_omega_star_list"]]
	###########################
	## END: update relevant L_s: L_omega;L_b
	###########################
	aa=get_global_L_sum(L_global_list)+get_local_L_sum(gene_dat_list)
	ssvec[4]=aa
	print(aa)
	print("update Eu_upsilon_list")
	outl=update_Eu_upsilon_list(theta_omega_star_list=theta_omega_star_list,gene_dat_list=gene_dat_list,Eu_omega_list=Eu_omega_list,Eu_upsilon_list=Eu_upsilon_list,L_global_list=L_global_list,h_vec=h_vec,myq=myq,chi1_vec=chi1_vec,Eu_chi2j=Eu_chi2j,D_mat=D_mat)
	L_global_list=outl[["L_global_list"]]
	Eu_upsilon_list=outl[["Eu_upsilon_list"]]
	gene_dat_list=outl[["gene_dat_list"]]
	theta_upsilon_star_list=outl[["theta_upsilon_star_list"]]
	aa=get_global_L_sum(L_global_list)+get_local_L_sum(gene_dat_list)
	ssvec[5]=aa
	print(aa)
	print("update Eu_chi2j")
	outl=update_Eu_chi2j(L_global_list=L_global_list,theta_upsilon_star_list=theta_upsilon_star_list,theta_chi2j=theta_chi2j,Eu_upsilon_list=Eu_upsilon_list,chi1_vec=chi1_vec,h_vec=h_vec,myq=myq)
	L_global_list=outl[["L_global_list"]]
	Eu_chi2j=outl[["Eu_chi2j"]]
	rm(outl)
	aa=get_global_L_sum(L_global_list)+get_local_L_sum(gene_dat_list)
	ssvec[6]=aa
	print(aa)
	global_E_list$Eu_ak=Eu_ak
	global_E_list$Eu_upsilon_list=Eu_upsilon_list
	global_E_list$Eu_lambda2=Eu_lambda2
	global_E_list$Eu_tau2=Eu_tau2
	global_E_list$Eu_chi2j=Eu_chi2j
	
	global_E_list$Eu_omega_list=Eu_omega_list
	global_E_list$Eu_nu_list=Eu_nu_list

	outl=list()
	outl[["L_global_list"]]=L_global_list
	outl[["global_E_list"]]=global_E_list
	outl[["gene_dat_list"]]=gene_dat_list
	outl[["ssvec"]]=ssvec
	return(outl)
}

update_Eu_ak=function(theta_ak,gene_dat_list,Eu_ak,hyperparameter_list,L_global_list,ak_matcher){
	###########################
	## update E_ak
	###########################	
	theta_ak_star=theta_ak-theta_ak
	my_t=dim(theta_ak_star)[1]
	for(k in c(1:my_t)){
		theta_alpha_ak_l=list()
		for(j in c(1:length(gene_dat_list))){
			cur_e=gene_dat_list[[j]]$E_list
			cur_obs=gene_dat_list[[j]]$Obs_list
			Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat=cur_obs$mapping_mat[,ak_matcher])
			theta_alpha_ak=calc_theta_alpha_ak_via_Eugammai(Eu_gammai=Eu_gammai,gamma1=hyperparameter_list[["gamma1"]],Eu_alphai=cur_e$Eu_alphai,Eu_kappa=cur_e$Eu_kappa,Eu_ak=Eu_ak,mapping_mat=cur_obs$mapping_mat[,ak_matcher],k=k)
			theta_alpha_ak_l[[j]]=theta_alpha_ak
		}
		theta_alpha_ak_k=do.call("rbind",theta_alpha_ak_l)
		theta_ak_k_star=colSums(theta_alpha_ak_k)+theta_ak[k,,drop=FALSE]
		Eu_ak[k,]=calc_Eu_ak_via_theta(theta_ak_k_star)
		theta_ak_star[k,]=theta_ak_k_star
	}
	###########################
	## END: update E_ak
	###########################
	###########################
	## update relevant L_s: L_ak;L_alphai
	###########################
	for(j in c(1:length(gene_dat_list))){
		cur_theta=gene_dat_list[[j]]$theta_list
		cur_e=gene_dat_list[[j]]$E_list
		cur_obs=gene_dat_list[[j]]$Obs_list
		Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat=cur_obs$mapping_mat[,ak_matcher])
		theta_alphai=calc_theta_alphai_via_Eu_gammai_Eu_kappa(hyperparameter_list[["gamma1"]],Eu_gammai,cur_e$Eu_kappa)
		Eu_kappa_expanded=rep(1,dim(Eu_gammai)[1])%*%t(cur_e$Eu_kappa)
		g_alphai_pa=calc_g_alphai_via_Eu_gammai_Eu_kappa(hyperparameter_list[["gamma1"]],Eu_gammai,Eu_kappa_expanded)
		L_alphai=sum((theta_alphai-cur_theta$theta_alphai_star)*cur_e$Eu_alphai)+sum(g_alphai_pa)-sum(calc_g_alphai_via_theta(cur_theta$theta_alphai_star))
		gene_dat_list[[j]][["L_list"]][["L_alphai"]]=L_alphai
	}
	L_ak=calculate_L_ak(Eu_ak,theta_ak,theta_ak_star)
	L_global_list[["L_ak"]]=L_ak
	###########################
	## END: update relevant L_s: L_ak;L_alphai(offspring)
	###########################
	outl=list()
	outl[["L_global_list"]]=L_global_list
	outl[["Eu_ak"]]=Eu_ak
	outl[["gene_dat_list"]]=gene_dat_list
	return(outl)
}


update_Eu_nu_list=function(gene_dat_list,L_global_list=L_global_list,Eu_omega_list,hyperparameter_list){
	theta_b_nu_l=list()
	for(j in c(1:length(gene_dat_list))){
		cur_e=gene_dat_list[[j]]$E_list
		cur_obs=gene_dat_list[[j]]$Obs_list
		cur_mapping_mat=cbind(TRUE,cur_obs$mapping_mat[,match(nu_names,colnames(cur_obs$mapping_mat))])
		theta_b_nu=calc_theta_bi_nu_via_Eu_wshift(Eu_b_list=cur_e$Eu_b_list,mapping_shift_mat=cur_obs$mapping_shift_mat,Eu_alphai=cur_e$Eu_alphai,Eu_omega_list=Eu_omega_list,mapping_mat=cur_mapping_mat)
		theta_b_nu_l[[j]]=theta_b_nu
	}
	theta_nu_list=calc_theta_nu(p=hyperparameter_list$p,q=dim(cur_mapping_mat)[2],myc=hyperparameter_list$c)
	theta_nu_star_list=theta_nu_list
	for(j in c(1:length(gene_dat_list))){
		theta_nu_star_list[[1]]=theta_nu_star_list[[1]]+theta_b_nu_l[[j]][[1]]
		theta_nu_star_list[[2]]=theta_nu_star_list[[2]]+theta_b_nu_l[[j]][[2]]
	}
	Eu_nu_list=calc_Eu_nu_via_theta(theta_nu_star_list)
	###########################
	## END: update Eu_nu_list;
	###########################
	###########################
	## update relevant L_s: L_nu; L_bs
	###########################
	for(j in c(1:length(gene_dat_list))){
		cur_e=gene_dat_list[[j]]$E_list
		cur_obs=gene_dat_list[[j]]$Obs_list
		cur_theta=gene_dat_list[[j]]$theta_list	
		cur_mapping_mat=cbind(TRUE,cur_obs$mapping_mat[,match(nu_names,colnames(cur_obs$mapping_mat))])
		theta_b=calc_theta_bi_via_Eu_wshift_nu(cur_e$Eu_alphai,Eu_omega_list,cur_obs$mapping_shift_mat,Eu_nu_list,cur_mapping_mat)
		g_b_pa=calc_g_b_via_Eu_wshift_nu(Eu_alphai=cur_e$Eu_alphai,Eu_omega_list=Eu_omega_list,V=cur_obs$mapping_shift_mat,F=cur_mapping_mat,Eu_nu_list=Eu_nu_list)
		L_b=sum((theta_b[[1]]-cur_theta$theta_b_star_list[[1]])*cur_e$Eu_b_list[[1]])+sum((theta_b[[2]]-cur_theta$theta_b_star_list[[2]])*cur_e$Eu_b_list[[2]])+g_b_pa-calc_g_b_via_theta(cur_theta$theta_b_star_list)
		gene_dat_list[[j]][["L_list"]][["L_b"]]=L_b[1]
	}

	g_nu_pa=calc_g_nu_via_theta(theta_nu_list)
	L_nu=sum((theta_nu_list[[1]]-theta_nu_star_list[[1]])*Eu_nu_list[[1]])+sum((theta_nu_list[[2]]-theta_nu_star_list[[2]])*Eu_nu_list[[2]])+g_nu_pa-calc_g_nu_via_theta(theta_nu_star_list)
	L_global_list[["L_nu"]]=L_nu[1]
	###########################	
	outl=list()
	outl[["L_global_list"]]=L_global_list
	outl[["Eu_nu_list"]]=Eu_nu_list
	outl[["gene_dat_list"]]=gene_dat_list
	return(outl)
}


update_Eu_omega_list=function(gene_dat_list,L_global_list=L_global_list,Eu_nu_list,hyperparameter_list,Eu_upsilon_list,D_mat){
	theta_b_omega_l=list()
	for(j in c(1:length(gene_dat_list))){
		cur_e=gene_dat_list[[j]]$E_list
		cur_obs=gene_dat_list[[j]]$Obs_list
		cur_mapping_mat=cbind(TRUE,cur_obs$mapping_mat[,match(nu_names,colnames(cur_obs$mapping_mat))])
		theta_b_omega=calc_theta_b_omega_via_Eu_nu(Eu_b_list=cur_e$Eu_b_list,Eu_alphai=cur_e$Eu_alphai,mapping_shift_mat=cur_obs$mapping_shift_mat,mapping_mat=cur_mapping_mat,Eu_nu_list=Eu_nu_list)	
		theta_b_omega_l[[j]]=theta_b_omega
	}
	theta_omega_list=calc_theta_omegas_via_Eu_upsilon(Eu_upsilon_list,D_mat)
	theta_omega_star_list=theta_omega_list
	for(j in c(1:length(gene_dat_list))){
		theta_omega_star_list[[1]]=theta_omega_star_list[[1]]+theta_b_omega_l[[j]][[1]]
		theta_omega_star_list[[2]]=theta_omega_star_list[[2]]+theta_b_omega_l[[j]][[2]]
	}
	Eu_omega_list=calc_Eu_omega_via_theta(theta_omega_star_list)
	###########################
	## END: update Eu_omega_list
	###########################
	###########################
	## update relevant L_s: L_omega;L_b
	###########################
	for(j in c(1:length(gene_dat_list))){
		cur_e=gene_dat_list[[j]]$E_list
		cur_obs=gene_dat_list[[j]]$Obs_list
		cur_theta=gene_dat_list[[j]]$theta_list	
		cur_mapping_mat=cbind(TRUE,cur_obs$mapping_mat[,match(nu_names,colnames(cur_obs$mapping_mat))])
		theta_b=calc_theta_bi_via_Eu_wshift_nu(cur_e$Eu_alphai,Eu_omega_list,cur_obs$mapping_shift_mat,Eu_nu_list,cur_mapping_mat)
		g_b_pa=calc_g_b_via_Eu_wshift_nu(Eu_alphai=cur_e$Eu_alphai,Eu_omega_list=Eu_omega_list,V=cur_obs$mapping_shift_mat,F=cur_mapping_mat,Eu_nu_list=Eu_nu_list)
		L_b=sum((theta_b[[1]]-cur_theta$theta_b_star_list[[1]])*cur_e$Eu_b_list[[1]])+sum((theta_b[[2]]-cur_theta$theta_b_star_list[[2]])*cur_e$Eu_b_list[[2]])+g_b_pa-calc_g_b_via_theta(cur_theta$theta_b_star_list)
		gene_dat_list[[j]][["L_list"]][["L_b"]]=L_b[1]
	}
	####
	g_omega_pa=sum(calc_g_omega_via_Eu_upsilon(Eu_upsilon_list,D_mat))
	theta_omega_list=calc_theta_omegas_via_Eu_upsilon(Eu_upsilon_list,D_mat)
	L_omega=sum((theta_omega_list[[1]]-theta_omega_star_list[[1]])*Eu_omega_list[[1]])+sum((theta_omega_list[[2]]-theta_omega_star_list[[2]])*Eu_omega_list[[2]])+g_omega_pa-calc_g_omega_via_theta(theta_omega_star_list)
	L_global_list[["L_omega"]]=L_omega[1]
	###########################	
	outl=list()
	outl[["L_global_list"]]=L_global_list
	outl[["Eu_omega_list"]]=Eu_omega_list
	outl[["gene_dat_list"]]=gene_dat_list
	outl[["theta_omega_star_list"]]=theta_omega_star_list
	return(outl)
}

update_Eu_upsilon_list=function(theta_omega_star_list,gene_dat_list,Eu_omega_list,Eu_upsilon_list,L_global_list,h_vec,myq,chi1_vec,Eu_chi2j,D_mat){
	theta_upsilon_star_list=list()
	theta_upsilon_list=calc_theta_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec=h_vec)
	for(myj in c(1:myq)){
		theta_omega_upsilon_j=calc_theta_omegas_upsilon_via_Eu(Eu_omega_list,Eu_upsilon_list,D_mat,myj)
		theta_upsilon_j=theta_upsilon_list[[myj]]
		theta_upsilon_star_j=theta_omega_upsilon_j+theta_upsilon_j
		theta_upsilon_star_list[[myj]]=theta_upsilon_star_j
		Eu_upsilon_j=calc_Eu_upsilon_via_theta(theta_upsilon_star_j)
		Eu_upsilon_list[[myj]]=Eu_upsilon_j
	}
	###########################
	## END: Eu_upsilon_list
	###########################	
	###########################
	## update relevant L_s: L_upsilon_list;L_omega
	###########################
	g_upsilon_pa=calc_g_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec)
	g_upsilon_theta_star=calc_g_upsilon_via_theta_list(theta_upsilon_star_list)
	L_upsilons=rep(0,myq)
	for(j in c(1:myq)){
		L_ups_j=(theta_upsilon_list[[j]]-theta_upsilon_star_list[[j]])*Eu_upsilon_list[[j]]
		L_ups_j=sum(L_ups_j)
		L_ups_j=L_ups_j+g_upsilon_pa[j]-g_upsilon_theta_star[j]
		L_upsilons[j]=L_ups_j
	}
	g_omega_pa=sum(calc_g_omega_via_Eu_upsilon(Eu_upsilon_list,D_mat))
	theta_omega_list=calc_theta_omegas_via_Eu_upsilon(Eu_upsilon_list,D_mat)
	L_omega=sum((theta_omega_list[[1]]-theta_omega_star_list[[1]])*Eu_omega_list[[1]])+sum((theta_omega_list[[2]]-theta_omega_star_list[[2]])*Eu_omega_list[[2]])+g_omega_pa-calc_g_omega_via_theta(theta_omega_star_list)
	L_global_list[["L_upsilon"]]=sum(L_upsilons)
	L_global_list[["L_omega"]]=L_omega[1]
	###########################	
	outl=list()
	outl[["L_global_list"]]=L_global_list
	outl[["Eu_upsilon_list"]]=Eu_upsilon_list
	outl[["gene_dat_list"]]=gene_dat_list
	outl[["theta_upsilon_star_list"]]=theta_upsilon_star_list
	return(outl)
}


update_Eu_lambda2=function(theta_lambda2,gene_dat_list,lambda1,approximation_mode,L_global_list){
	all_Eu_lambda_l=lapply(c(1:length(gene_dat_list)),function(j){
			return(gene_dat_list[[j]][["E_list"]][["Eu_lambda"]])
	})
	all_Eu_lambda=do.call("rbind",all_Eu_lambda_l)
	theta_lambdai_lambda2=calc_theta_lambdai_lambda2_via_Eu_lambdai(lambda1=lambda1,Eu_lambdai=all_Eu_lambda,genespecific_lambda1=approximation_mode)
	theta_lambda2_star=theta_lambdai_lambda2+theta_lambda2
	Eu_lambda2=calc_Eu_lambda2_via_theta(theta_lambda2_star)
	###########################
	## END: update Eu_lambda2
	###########################
	###########################
	## update relevant L_s: L_lambda2;L_lambda
	###########################
	L_lambda2=sum((theta_lambda2-theta_lambda2_star)*Eu_lambda2)+calc_g_lambda_via_theta(theta_lambda2)-calc_g_lambda_via_theta(theta_lambda2_star)
	L_global_list[["L_lambda2"]]=L_lambda2
	for(j in c(1:length(gene_dat_list))){
		cur_theta=gene_dat_list[[j]]$theta_list
		cur_e=gene_dat_list[[j]]$E_list
		cur_obs=gene_dat_list[[j]]$Obs_list
		if(approximation_mode){	
			lambda1j=lambda1[j]
			theta_lambda=calc_theta_lambda(lambda1j,Eu_lambda2[2])
			g_lambda_pa=calc_g_lambda_via_Eu_lambda2(lambda1j,Eu_lambda2)
		}else{
			theta_lambda=calc_theta_lambda(lambda1,Eu_lambda2[2])
			g_lambda_pa=calc_g_lambda_via_Eu_lambda2(lambda1,Eu_lambda2)
		}
		L_lambda=sum((theta_lambda-cur_theta$theta_lambda_star)*cur_e$Eu_lambda)+g_lambda_pa-calc_g_lambda_via_theta(cur_theta$theta_lambda_star)
		gene_dat_list[[j]][["L_list"]][["L_lambda"]]=L_lambda
	}
	outl=list()
	outl[["L_global_list"]]=L_global_list
	outl[["Eu_lambda2"]]=Eu_lambda2
	outl[["gene_dat_list"]]=gene_dat_list
	return(outl)
}

update_Eu_tau2=function(i,gene_dat_list,L_global_list,hyperparameter_list,theta_tau2){
	all_Eu_kappa_l=lapply(c(1:length(gene_dat_list)),function(j){
			return(gene_dat_list[[j]][["E_list"]][["Eu_kappa"]])
	})
	all_Eu_kappa=do.call("rbind",all_Eu_kappa_l)
	theta_kappaj_tau2=calc_theta_kappaj_tau2_via_Eu_kappaj(hyperparameter_list[["tau1"]],Eu_kappaj=all_Eu_kappa)
	theta_tau2_star=theta_tau2+theta_kappaj_tau2
	Eu_tau2=calc_Eu_tau2_via_theta(theta_tau2_star)
	###########################
	## END: update Eu_tau2
	###########################
	###########################
	## update relevant L_s: L_tau2;L_kappa
	###########################
	L_tau2=sum((theta_tau2-theta_tau2_star)*Eu_tau2)+calc_g_tau2_via_theta(theta_tau2)-calc_g_tau2_via_theta(theta_tau2_star)
	L_global_list[["L_tau2"]]=L_tau2
	for(j in c(1:length(gene_dat_list))){
		cur_theta=gene_dat_list[[j]]$theta_list
		cur_e=gene_dat_list[[j]]$E_list
		cur_obs=gene_dat_list[[j]]$Obs_list
		theta_kappa=calc_theta_kappaj_via_Eu_tau2(hyperparameter_list[["tau1"]],Eu_tau2)
		g_kappa_pa=calc_g_kappa_via_Eu_tau2(hyperparameter_list[["tau1"]],Eu_tau2)
		theta_kappa_star=cur_theta$theta_kappa_star
		L_kappa=sum((theta_kappa-theta_kappa_star)*cur_e$Eu_kappa)+g_kappa_pa-calc_g_kappa_via_theta(theta_kappa_star)
		gene_dat_list[[j]][["L_list"]][["L_kappa"]]=L_kappa
	}
	outl=list()
	outl[["L_global_list"]]=L_global_list
	outl[["Eu_tau2"]]=Eu_tau2
	outl[["gene_dat_list"]]=gene_dat_list
	return(outl)
}


update_Eu_chi2j=function(L_global_list,theta_upsilon_star_list,theta_chi2j,Eu_upsilon_list,chi1_vec,h_vec,myq){
	###########################
	## Eu_chi2j
	###########################
	theta_upsilon_chi2j=calc_theta_upsilon_chi2j_via_Eu(Eu_upsilon_list,chi1_vec)
	theta_chi2j_star=theta_chi2j+theta_upsilon_chi2j
	Eu_chi2j=calc_Eu_chi2j_via_theta(theta_chi2j_star)
	###########################
	## END: Eu_chi2j
	###########################	
	###########################
	## update relevant L_s: L_upsilon_list;L_chi2j
	###########################
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
	g_chi2j_pa=calc_g_chi2j_via_theta(theta_chi2j)
	g_chi2j_theta_star=calc_g_chi2j_via_theta(theta_chi2j_star)
	L_chi2j=sum((theta_chi2j-theta_chi2j_star)*Eu_chi2j)+sum(g_chi2j_pa)-sum(g_chi2j_theta_star)
	L_global_list[["L_upsilon"]]=sum(L_upsilons)
	L_global_list[["L_chi2j"]]=L_chi2j
	###########################
	## END update relevant L_s: L_upsilon_list;L_chi2j
	###########################
	outl=list()
	outl[["L_global_list"]]=L_global_list
	outl[["Eu_chi2j"]]=Eu_chi2j
	return(outl)
}

### update ordering simplification. run_bagea_wshift_non_optim
### uses the following update strategy.
### if we want to update X, we check whether any  coparents of X co(X) (say Z)were updated. If so, the child nodes that intersect  are updated first, i.e. if  Y=ch(X) \intersect ch(co(X)) then Y has to be updated first.
### However, this is not really necessary. Rather,  m_yx is a function of Z as well as Y, i.e.  m_yx(Y,Z). We therefore send m_yx(Y,Z*) where Z* is the updated sufficient statistic of Z. When calculating L we have to add <P(Y|Z,X)> plugging in Z* and X*. However, Q(Y) is not updated yet. 


#' run nonoptimized variational updates.
#'
#' This function runs the same variational updates as \code{run_bagea} but updating steps are non optimized. It therefore can be used for integration testing. Note that ,while the various estimates should be in very good agreement, They will most likely not be identical due to numerical precision of the various algorithms employed.
#' @export
run_bagea_nonoptimized=function(gene_dat_list=NULL,hyperparameter_list=NULL,calc_L=TRUE,ncores=NULL,nrounds=1000,calc_L_rounds=50,approximation_mode=TRUE,ak_names=NULL){
	####################
	## some imput checks
	####################
	Eu_ak=NULL## init Eu_ak value
	calc_B_rounds=1 ## number of rounds before b estimates should be updated.
	calc_B=TRUE
	constant_noise_prior_pergene=FALSE
	if(ncores!=1){
		stop("ncores has to be 1 for run_bagea_non_optim")
	}
	if(calc_B_rounds!=1){
		stop("calc_B_rounds has to be 1 for run_bagea_non_optim")
	}
	nu_names=get_nu_names(hyperparameter_list)
	outl=extract_D_mat_info(gene_dat_list)
	D_mat=outl[["D_mat"]]
	D_names_list=outl[["D_names_list"]]
	if(is.null(ak_names)){
		ak_names=get_mapping_mat_colnames(gene_dat_list)
	}
	ak_matcher=match(ak_names,get_mapping_mat_colnames(gene_dat_list))
	check_input(hyperparameter_list=hyperparameter_list,calc_L=calc_L,calc_B=TRUE,ncores=ncores,nrounds=nrounds,calc_L_rounds=calc_L_rounds,calc_B_rounds=calc_B_rounds,write_logfile=FALSE,approximation_mode=approximation_mode,nu_names=nu_names,ak_names)
	################
	outl=prepare_inputdata_lists(gene_dat_list,hyperparameter_list,nu_names,ak_names)
	gene_dat_list=outl[["gene_dat_list"]]
	global_E_list=outl[["global_E_list"]]
	input_settings=outl[["inputsettings"]]
	rm(outl)
	check_gene_dat_list(gene_dat_list=gene_dat_list)
	check_global_E_list(global_E_list)
	print("checks done")
	####################
	## END: some imput checks
	####################
	mysettings=get_settingsout_list_bagea(gene_dat_list_settings=input_settings,hyperparameter_list=hyperparameter_list,calc_L=calc_L,ncores=ncores,nrounds=nrounds,calc_L_rounds=calc_L_rounds,Eu_ak=global_E_list$Eu_ak,calc_B_rounds=calc_B_rounds,approximation_mode=approximation_mode,constant_noise_prior_pergene=constant_noise_prior_pergene,nu_names=nu_names,ak_names=ak_names)
	gene_dat_list=reassemble_Eu_b_list(gene_dat_list)
	gene_dat_list=reassemble_XtX(gene_dat_list)
	outl=initialize_lambda1nrho1_for_approx_mode(gene_dat_list,hyperparameter_list[["lambda1"]],hyperparameter_list[["rho1"]],constant_noise_prior_pergene,approximation_mode)
	lambda1=outl[["lambda1"]]
	rho1=outl[["rho1"]]
	####################
	## prepare (unchanging) theta for root nodes
	####################
	theta_root_list=prepare_theta_root_list(my_t=dim(global_E_list$Eu_ak)[1],hyperparameter_list,D_mat,nu_names)
	all_Ls_list=list()
	L_global_list=list()
	ssMat=matrix(0,nrounds,9)
	for(i in c(1:nrounds)){	
		####################
		####################
		####################
		## update order:
		## 	local:		
		## 		genewise(loop):		
		##		Eu_lambdaj
		##		Eu_bj
		##	local:
		##		Eu_kappa(all):
		##		Eu_alphai(all):
		##	global:
		##		Eu_ak: 
		##	global:
		##		Eu_lambda2
		##		Eu_tau2	
		##	global:
		##		Eu_nu_list	
		##	global:
		##		Eu_omega_list
		## global:
		##		Eu_upsilon_list
		####################
		local_calc_L=calc_L && (((i %% calc_L_rounds)==0) || i==1)
		print(Sys.time())
		print("round")
		print(i)
		gene_dat_list=update_Eu_bj_Eu_lambdaj(gene_dat_list,global_E_list$Eu_lambda2,lambda1,approximation_mode,Eu_omega_list=global_E_list$Eu_omega_list,Eu_nu_list=global_E_list$Eu_nu_list)
		print("ss L")
		aa=get_global_L_sum(L_global_list)+get_local_L_sum(gene_dat_list)
		ssMat[i,1]=aa
		print(aa)
		gene_dat_list=update_Eu_kappai(i=i,gene_dat_list=gene_dat_list,Eu_tau2=global_E_list$Eu_tau2,hyperparameter_list=hyperparameter_list,ak_matcher=ak_matcher,Eu_ak=global_E_list$Eu_ak)
		aa=get_global_L_sum(L_global_list)+get_local_L_sum(gene_dat_list)
		ssMat[i,2]=aa
		print(aa)
		gene_dat_list=update_Eu_alphai(
			gene_dat_list=gene_dat_list,
			Eu_ak=global_E_list$Eu_ak,
			ak_matcher=ak_matcher,
			hyperparameter_list=hyperparameter_list,
			Eu_omega_list=global_E_list$Eu_omega_list,
			Eu_nu_list=global_E_list$Eu_nu_list)
		print("ss3 L")
		aa=get_global_L_sum(L_global_list)+get_local_L_sum(gene_dat_list)
		ssMat[i,3]=aa
		print(aa)
		###########################
		global_ret=run_globalE_list(i=i,L_global_list=L_global_list,global_E_list,gene_dat_list,ak_matcher=ak_matcher,hyperparameter_list=hyperparameter_list,approximation_mode=approximation_mode,lambda1=lambda1,D_mat=D_mat,theta_root_list=theta_root_list)
		L_global_list=global_ret$L_global_list
		global_E_list=global_ret$global_E_list
		gene_dat_list=global_ret$gene_dat_list
		ssvec=global_ret$ssvec
		ssMat[i,(3+(1:length(ssvec)))]=ssvec
		###########################
		gene_L_list=lapply(gene_dat_list,function(x){
			return(x$L_list)
		})
	}
	all_Ls=list()
	all_Ls[["L_global_list"]]=L_global_list
	all_Ls[["L_local_list"]]=get_local_L_sum(gene_dat_list)
	out_list=prepare_nonoptim_outl(gene_dat_list=gene_dat_list,mysettings=mysettings,global_E_list=global_E_list,all_Ls=all_Ls,ssMat=ssMat,D_mat=D_mat,D_names_list=D_names_list)	
	return(out_list)	
}
