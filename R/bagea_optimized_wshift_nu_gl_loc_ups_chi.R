get_objsize_on_node=function(){
	### Assumes that EU_alphai was updated
	current_round=get("current_round",get_data_container_environment())#update local
	calc_L=get("calc_L",get_data_container_environment())#const
	calc_L_rounds=get("calc_L_rounds",get_data_container_environment())#const
	calc_B_rounds=get("calc_B_rounds",get_data_container_environment())#const
	gamma1=get("gamma1",get_data_container_environment())#const
	nodenr=get("nodenr",get_data_container_environment())#const
	theta_lambda=get("theta_lambda",get_data_container_environment())#update global
	Eu_lambda2=get("Eu_lambda2",get_data_container_environment())#update global
	approximation_mode=get("approximation_mode",get_data_container_environment())#const
	Eu_omega_list=get("Eu_omega_list",get_data_container_environment())
	Eu_nu_list=get("Eu_nu_list",get_data_container_environment())
	nu_matcher=get("nu_matcher",get_data_container_environment())
	ak_matcher=get("ak_matcher",get_data_container_environment())
	#update_global
#	Eu_alphai_vect=get("Eu_alphai",get_data_container_environment())
	tau1=get("tau1",get_data_container_environment())#const
	lambda1=get("lambda1",get_data_container_environment())#const
	mapping_mat_expanded=get("mapping_mat_expanded",get_data_container_environment())#const
	index_mat=get("index_mat",get_data_container_environment())#const
	index_mat_t=get("index_mat_t",get_data_container_environment())#const
	mapping_shift_mat_expanded=get("mapping_shift_mat_expanded",get_data_container_environment())#const
	theta_alphai_kappa_vect_first=get("theta_alphai_kappa_vect_first",get_data_container_environment())#const
	theta_alphai_star_vect=get("theta_alphai_star_vect",get_data_container_environment())#
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
	ret=get_objsizes()
	return(ret)
}

#get_shift_mat_index_list
get_shift_mat_melted_on_node=function(){
	mapping_shift_mat_expanded=get("mapping_shift_mat_expanded",get_data_container_environment())
	mapping_shift_mat_expanded=t(mapping_shift_mat_expanded)
	tt=which(mapping_shift_mat_expanded!=0)
	vals=mapping_shift_mat_expanded[tt]
	tt=which(mapping_shift_mat_expanded!=0,2)
	shift_mat_melted=cbind(tt,vals)
	assign("shift_mat_melted",shift_mat_melted,envir=get_data_container_environment())
	return(NULL)
}

#set shift_mat_melted to NULL this will prevent running the compiled code
get_shift_mat_melted_on_node_sham=function(){
	shift_mat_melted=NULL
	assign("shift_mat_melted",shift_mat_melted,envir=get_data_container_environment())
	return(NULL)
}


get_objsize_on_node2=function(){
	### Assumes that EU_alphai was updated
	all_obj=ls(get_data_container_environment())
	all_sizes=sapply(all_obj,function(x){
		obj=get(x,get_data_container_environment())
		size=pryr::object_size(obj)
	})
	out=data.table(name=all_obj,size=all_sizes)
	setkey(out,size)
	return(out)
}



# not tested
calc_theta_delta=function(chi1,chi2){
	theta_delta=c(chi1-1,-chi2)
	return(theta_delta)
}

# not tested
calc_Eu_delta_via_theta=function(theta){
	return(calc_Eu_gamma_func_via_theta(theta))
}


# not tested
calc_Eu_omega_via_theta=function(theta_omega_list){
	theta_omega_list2_inv=solve(theta_omega_list[[2]])
	Eu_omega1=-theta_omega_list2_inv%*%theta_omega_list[[1]]/2
	Eu_omega2=-theta_omega_list2_inv/2+Eu_omega1%*%t(Eu_omega1)
	Eu_omega_list=list()
	Eu_omega_list[[1]]=Eu_omega1
	Eu_omega_list[[2]]=Eu_omega2
	return(Eu_omega_list)
}

calc_Eu_nu_via_theta=function(theta_nu_list){
	return(calc_Eu_omega_via_theta(theta_nu_list))
}

calc_theta_omegas_via_Eu=function(Eu_delta,s){
	theta_list=list()
	theta_list[[1]]=rep(0,s)
	if(s<1){
		traceback()
		stop("errm2-ce:s has to be larger than 1.")
	}
	theta_list[[2]]=-Eu_delta[2]*diag(rep(1,s))/2
	return(theta_list)
}

calc_theta_omegas_via_Eu_upsilon=function(Eu_upsilon_list,D_mat){
	theta_list=list()
	s=dim(D_mat)[1]
	myq=length(Eu_upsilon_list)
	theta_list[[1]]=rep(0,s)
	if(s<1){
		traceback()
		stop("errm2-ce:s has to be larger than 1.")
	}
	res=rep(1,s)
	for(j in c(1:myq)){
		res=res*Eu_upsilon_list[[j]][D_mat[,j],2]
	}
	theta_list[[2]]=-diag(res)/2
	return(theta_list)
}


calc_theta_nu=function(p,q=length(p),myc=NULL){
	theta_list=list()
	if(length(p)==1 && is.null(myc)){
		theta_list[[1]]=p*rep(1,q)
		if(q<1){
			traceback()
			stop("erre-sms-ce:s has to be larger than 1.")
		}
		theta_list[[2]]=-diag(rep(p,q))/2
	}else{
		if(!is.null(myc)){
			if(length(p)!=length(myc)){
				stop("errpmasdco:p and myc have to have same length.")
			}
			theta_list[[1]]=myc*p
		}else{
			theta_list[[1]]=p			
		}
		theta_list[[2]]=-diag(p)/2
	}
	return(theta_list)
}

calc_theta_omegas_deltas_via_Eu=function(Eu_omega_list){
	s=length(Eu_omega_list[[1]])
	theta_delta1=s/2
	theta_delta2=-sum(diag(Eu_omega_list[[2]]))/2
	theta_delta=c(theta_delta1,theta_delta2)
	return(theta_delta)
}

calc_theta_omegas_upsilon_via_Eu=function(Eu_omega_list,Eu_upsilon_list,D_mat,myj=myj){
	theta_list=list()
	myq=dim(D_mat)[2]
	h_vec=apply(D_mat,2,max)
	h_vec_j=h_vec[myj]
	theta_first=rep(0,h_vec_j)
	for(k in c(1:h_vec_j)){
		theta_first[k]=0.5*sum(D_mat[,myj]==k)
	}
	j_primes=c(1:myq)[-myj]
	res=rep(1,dim(D_mat)[1])
	if(length(j_primes)>0){
		for(j_prime in j_primes){
			res=res*Eu_upsilon_list[[j_prime]][D_mat[,j_prime],2]
		}
	}
	res=-res*diag(Eu_omega_list[[2]])/2
	theta_sec=rep(0,h_vec_j)
	for(k in (1:h_vec_j)){
		theta_sec[k]=sum(res[D_mat[,myj]==k])
	}
	theta_out=cbind(theta_first,theta_sec)
	return(theta_out)
}


calc_theta_upsilon_chi2j_via_Eu=function(Eu_upsilon_list,chi1_vec){
	theta_list=list()
	myq=length(chi1_vec)
	theta_out=matrix(0,myq,2)
	for(j in c(1:myq)){
		h_vec_j=dim(Eu_upsilon_list[[j]])[1]
		theta1=chi1_vec[j]*h_vec_j
		theta2=-sum(Eu_upsilon_list[[j]][,2])
		theta_out[j,1]=theta1
		theta_out[j,2]=theta2
	}
	return(theta_out)
}


calc_theta_b_omega_via_Eu_nu=function(Eu_b_list,Eu_alphai,mapping_shift_mat,mapping_mat,Eu_nu_list){
	s=length(mapping_shift_mat[1,])
	first=mapping_mat%*%Eu_nu_list[[1]]
	first=(Eu_alphai[,2]*Eu_b_list[[1]]*first)
	first_mat=first%*%t(rep(1,s))
	theta_b_omega_1=Matrix::colSums(first_mat*mapping_shift_mat)
	second=mapping_mat%*%Eu_nu_list[[2]]
	second=Matrix::rowSums(second*mapping_mat)
	second_mat=-((Eu_alphai[,2]*second)/2)%*%t(rep(1,s))
	theta_b_omega_2=t(second_mat*mapping_shift_mat)%*%mapping_shift_mat
	theta_b_omega_list=list()
	theta_b_omega_list[[1]]=as.matrix(theta_b_omega_1)
	theta_b_omega_list[[2]]=as.matrix(theta_b_omega_2)
	return(theta_b_omega_list)
}



calc_theta_b_omega_via_Eu_nu_vectorized_fast=function(all_Eu_b1,all_Eu_alphai,mapping_shift_mat_expanded,mapping_mat_expanded,Eu_nu_list,nu_matcher,add_intercept=FALSE){
	non_zeros2=which(mapping_shift_mat_expanded!=0,2)
	if(add_intercept){
		mapping_mat_expanded=cbind(TRUE,mapping_mat_expanded[,nu_matcher])
	}
	s=length(mapping_shift_mat_expanded[1,])
	first=mapping_mat_expanded%*%Eu_nu_list[[1]]
	first=(all_Eu_alphai[,2]*all_Eu_b1*first)
#	first_mat=first%*%t(rep(1,s))
	first_mat=sparseMatrix(i=non_zeros2[,1],j=non_zeros2[,2],x=first[non_zeros2[,1]],dims=dim(mapping_shift_mat_expanded))

	theta_b_omega_1=colSums(first_mat*mapping_shift_mat_expanded)
	second=mapping_mat_expanded%*%Eu_nu_list[[2]]
	second=rowSums(second*mapping_mat_expanded)
	#second_mat
	sec=sqrt((all_Eu_alphai[,2]*second)/2)
	preadder=sparseMatrix(i=non_zeros2[,1],j=non_zeros2[,2],x=sec[non_zeros2[,1]],dims=dim(mapping_shift_mat_expanded))
	Vprime=preadder*mapping_shift_mat_expanded
	theta_b_omega_2=-t(Vprime)%*%Vprime
	theta_b_omega_list=list()
	theta_b_omega_list[[1]]=as.matrix(theta_b_omega_1)
	theta_b_omega_list[[2]]=as.matrix(theta_b_omega_2)
	return(theta_b_omega_list)
}


######## not tested
calc_theta_bi_alphai_via_Eu_wshift_nu=function(m,Eu_b_list,mapping_shift_mat,Eu_omega_list,decomposed=FALSE,mapping_mat,Eu_nu_list,add_intercept=FALSE){
	if(add_intercept){
		mapping_mat=cbind(TRUE,mapping_mat)
	}
	theta_bi_alpha_1=rep(1/2,m)
	if(decomposed){
		first=(0.5*get_diag_from_inv_decompose(Eu_b_list$theta_M_inv)+(Eu_b_list$Eu_b_1)^2)
	}else{
		first=diag(Eu_b_list[[2]])
	}
 	sec1=mapping_shift_mat%*%Eu_omega_list[[2]]
 	sec1=rowSums(sec1*mapping_shift_mat)
 	sec2=mapping_mat%*%Eu_nu_list[[2]]
 	sec2=Matrix::rowSums(sec2*mapping_mat)
 	sec=sec1*sec2
 	if(decomposed){
 		third=(Eu_b_list[[3]])
 	}else{
 		third=(Eu_b_list[[1]])
 	}
 	third=-2*third*(mapping_shift_mat%*%Eu_omega_list[[1]])
 	third=third*(mapping_mat%*%Eu_nu_list[[1]])
 	theta_bi_alphai_2=-(first+sec+third)/2
 	theta_bi_alphai=as.matrix(cbind(theta_bi_alpha_1,theta_bi_alphai_2))
	return(theta_bi_alphai)
}


######## not tested
calc_theta_bi_nu_via_Eu_wshift=function(Eu_b_list,mapping_shift_mat,Eu_alphai,Eu_omega_list,decomposed=FALSE,mapping_mat){
	if(decomposed){
 		first=(Eu_b_list[[3]])
 	}else{
 		first=(Eu_b_list[[1]])
 	}
 	first=first*Eu_alphai[,2]*(mapping_shift_mat%*%Eu_omega_list[[1]])
 	first_mat=first%*%t(rep(1,dim(mapping_mat)[2]))
	theta_b_nu_1=Matrix::colSums(first_mat*mapping_mat)
	###
	sec=mapping_shift_mat%*%Eu_omega_list[[2]]
	sec=Matrix::rowSums(mapping_shift_mat*sec)
	sec=sqrt(sec*Eu_alphai[,2]/2)
	###
	non_zeros2=Matrix::which(mapping_mat!=0,2)
	preadder=sparseMatrix(i=non_zeros2[,1],j=non_zeros2[,2],x=sec[non_zeros2[,1]],dims=dim(mapping_mat))
	Fprime=preadder*mapping_mat
	FptFp=-Matrix::t(Fprime)%*%Fprime
	out_l=list()
	out_l[[1]]=theta_b_nu_1
	out_l[[2]]=FptFp
	return(out_l)
}

calc_diag_mtCm_flipped=function(C_mat,mapping_shift_mat){
	mapping_shift_mat_t=t(mapping_shift_mat)
	myn=dim(mapping_shift_mat)[1]
	mym=dim(mapping_shift_mat)[2]
	myeig=eigen(C_mat)
	vec_t=t(myeig$vectors)
	out_m=rep(0,myn)
	for(i in c(1:mym)){
		if(i%%200==0){
			print(i)
		}
		rr=vec_t[i,]%*%mapping_shift_mat_t
		rr=rr^2*myeig$values[i]
		out_m=rr+out_m
	}
	out_m=t(out_m)
	return(out_m)
}

calc_diag_mtCm_flipped2=function(C_mat,mapping_shift_mat){
	mapping_shift_mat_t=t(mapping_shift_mat)
	myn=dim(mapping_shift_mat)[1]
	mym=dim(mapping_shift_mat)[2]
	myeig=eigen(C_mat)
	vec_t=t(myeig$vectors)
	out_m=rep(0,myn)
	svdvalmat=sqrt(myeig$values)%*%t(rep(1,length(myeig$values)))
	vec_t=vec_t*svdvalmat
	for(i in c(1:mym)){
		if(i%%200==0){
			print(i)
		}
		rr=vec_t[i,]%*%mapping_shift_mat_t
		rr=rr^2
		out_m=rr+out_m
	}
	out_m=t(out_m)
	return(out_m)
}

calc_diag_mtCm_flipped_preindexed=function(C_mat,mapping_shift_mat,shift_mat_index_list){
	mapping_shift_mat_t=t(mapping_shift_mat)
	myn=dim(mapping_shift_mat)[1]
	mym=dim(mapping_shift_mat)[2]
	out_m=rep(0,myn)
	mylist=list()
	for(i in c(1:mym)){
		if(i%%200==0){
			print(i)
		}
		if(length(shift_mat_index_list[[i]])>0){
			mylist[[i]]=as.vector(C_mat[i,]%*%mapping_shift_mat_t[,shift_mat_index_list[[i]]])
		}
	}
	out_m=tryCatch(
			{
				out_m=do.call("c",mylist)				
			},error=function(){
				mylen=length(mylist)
				dd=ceiling(seq(from=0,to=mylen,length.out=20))
				myupper=dd[2:20]
				mylower=dd[1:19]+1
				out_m=lapply(c(1:19),function(i){
					cur_o=mylist[mylower[i]:myupper[i]]
					out_el=do.call("c",cur_o)
					return(out_el)
				})
				out_m=do.call("c",out_m)
				return(out_m)
			},finally={
			}
		)
	non_zeros2=which(mapping_shift_mat!=0,2)
	preadder=sparseMatrix(i=non_zeros2[,1],j=non_zeros2[,2],x=out_m,dims=dim(mapping_shift_mat))
	ret=rowSums(preadder*mapping_shift_mat)
	return(ret)
}

calc_diag_mtCm_flipped_preindexed2=function(C_mat,mapping_shift_mat,shift_mat_index_list){
	mapping_shift_mat_t=t(mapping_shift_mat)
	myn=dim(mapping_shift_mat)[1]
	mym=dim(mapping_shift_mat)[2]
	out_m=rep(0,myn)
	for(i in c(1:mym)){
		if(i%%200==0){
			print(i)
		}
		if(length(shift_mat_index_list[[i]])>0){
			gg=C_mat[i,]%*%mapping_shift_mat_t[,shift_mat_index_list[[i]]]
			out_m[shift_mat_index_list[[i]]]=out_m[shift_mat_index_list[[i]]]+(gg*mapping_shift_mat_t[i,shift_mat_index_list[[i]]])
		}
	}
	out_m=t(out_m)
	return(out_m)
}

calc_diag_mtCm=function(C_mat,mapping_shift_mat){
	myn=dim(mapping_shift_mat)[1]
	mym=dim(mapping_shift_mat)[2]
	myeig=eigen(C_mat)
	out_m=rep(0,myn)
	for(i in c(1:mym)){
		if(i%%200==0){
			print(i)
		}
		rr=mapping_shift_mat%*%myeig$vectors[,i]
		rr=rr^2*myeig$values[i]
		out_m=rr+out_m
	}
	return(out_m)
}


calc_theta_b_nu_via_Eu_vectorized_fast=function(all_Eu_b1,all_Eu_alphai,Eu_omega_list,mapping_shift_mat_expanded,mapping_mat_expanded,nu_matcher,add_intercept=FALSE,low_memory=TRUE,shift_mat_melted=NULL){	
	mapping_mat_expanded=mapping_mat_expanded[,nu_matcher]
	if(add_intercept){
		mapping_mat_expanded=cbind(TRUE,mapping_mat_expanded)
	}
 	first=all_Eu_b1*all_Eu_alphai[,2]*(mapping_shift_mat_expanded%*%Eu_omega_list[[1]])
 	first_mat=first%*%t(rep(1,dim(mapping_mat_expanded)[2]))
	theta_b_nu_1=Matrix::colSums(first_mat*mapping_mat_expanded)
	###
	print(Sys.time())
	if(low_memory){
		if(is.null(shift_mat_melted)){
			sec=calc_diag_mtCm_flipped(C_mat=Eu_omega_list[[2]],mapping_shift_mat=mapping_shift_mat_expanded)
		}else{
			lenM=dim(mapping_shift_mat_expanded)[1]
			sec=getDiagMtCM_optim(Cmat=as.matrix(Eu_omega_list[[2]]),shift_mat_melted,lenM)
		}
	}else{		
		sec=mapping_shift_mat_expanded%*%Eu_omega_list[[2]]
		sec=Matrix::rowSums(mapping_shift_mat_expanded*sec)
	}
	print(Sys.time())
	###
	sec=sqrt(sec*all_Eu_alphai[,2]/2)
	###
	non_zeros2=which(mapping_mat_expanded!=0,2)
	preadder=sparseMatrix(i=non_zeros2[,1],j=non_zeros2[,2],x=sec[non_zeros2[,1]],dims=dim(mapping_mat_expanded))
	Fprime=preadder*mapping_mat_expanded
	FptFp=-t(Fprime)%*%Fprime
	out_l=list()
	out_l[[1]]=theta_b_nu_1
	out_l[[2]]=FptFp
	return(out_l)
}

calc_theta_b_nu_via_Eu_vectorized_slow=function(all_Eu_b1,all_Eu_alphai,Eu_omega_list,mapping_shift_mat_expanded,mapping_mat_expanded,nu_matcher,add_intercept=FALSE){
	mapping_mat_expanded=mapping_mat_expanded[,nu_matcher]
	if(add_intercept){
		mapping_mat_expanded=cbind(TRUE,mapping_mat_expanded)
	}
 	first=all_Eu_b1*all_Eu_alphai[,2]*(mapping_shift_mat_expanded%*%Eu_omega_list[[1]])
 	first_mat=first%*%t(rep(1,dim(mapping_mat_expanded)[2]))
	theta_b_nu_1=Matrix::colSums(first_mat*mapping_mat_expanded)
	###
	sec=mapping_shift_mat_expanded%*%Eu_omega_list[[2]]
	sec=Matrix::rowSums(mapping_shift_mat_expanded*sec)
	sec=sqrt(sec*all_Eu_alphai[,2]/2)
	###
	non_zeros2=which(mapping_mat_expanded!=0,2)
	preadder=sparseMatrix(i=non_zeros2[,1],j=non_zeros2[,2],x=sec[non_zeros2[,1]],dims=dim(mapping_mat_expanded))
	Fprime=preadder*mapping_mat_expanded
	FptFp=-t(Fprime)%*%Fprime
	out_l=list()
	out_l[[1]]=theta_b_nu_1
	out_l[[2]]=FptFp
	out_l[[3]]=sec
	out_l[[4]]=Fprime
	return(out_l)
}

####### not tested
calc_theta_bi_via_Eu_wshift_nu=function(Eu_alphai,Eu_omega_list,mapping_shift_mat,Eu_nu_list,mapping_mat_sub){	
	tmp=mapping_shift_mat%*%Eu_omega_list[[1]]
	if(is.null(dim(mapping_mat_sub)[2]) ||dim(mapping_mat_sub)[2]!=length(Eu_nu_list[[1]])){
		stop("errwocew:lengths do not agree.")
	}
	tmp2=mapping_mat_sub%*%Eu_nu_list[[1]]
	theta_first=matrix(Eu_alphai[,2]*tmp2*tmp)
	if(length(Eu_alphai[,2])==1){
		stop("err20kdssf")
	}
	theta_sec=diag(-Eu_alphai[,2]/2)
	theta_out_list=list()
	theta_out_list[[1]]=theta_first
	theta_out_list[[2]]=theta_sec
	return(theta_out_list)
}



calc_g_delta_via_theta=function(theta){
	g=(theta[1]+1)*log(-theta[2])-lgamma(theta[1]+1)
	return(g)
}

calc_g_omega_via_theta=function(theta_list){
	g=calc_g_b_via_theta(theta_list)
	return(g)
}

calc_g_nu_via_theta=function(theta_list){
	g=calc_g_b_via_theta(theta_list)
	return(g)
}


calc_g_b_via_Eu_wshift_nu=function(Eu_alphai,Eu_omega_list,V,F,Eu_nu_list){
	F_first=F%*%Eu_nu_list[[2]]
	F_first=Matrix::rowSums(F_first*F)
	firster=sqrt((Eu_alphai[,2]*F_first)/2)
	preadder=firster%*%t(rep(1,dim(V)[2]))
	Vprime=preadder*V
	VptVp=t(Vprime)%*%Vprime
	my_trace=sum(Eu_omega_list[[2]]*as.matrix(VptVp))
	myn=dim(Eu_alphai)[1]
	my_g=-0.5*log(2*pi)*myn+sum(0.5*Eu_alphai[,1])-my_trace
	return(my_g)
}

calc_g_b_via_Eu_wshift_nu_fast=function(Eu_alphai,Eu_omega_list,V,F,Eu_nu_list,nu_matcher=NULL,add_intercept,mapping_mat_nu_sorted=FALSE){
	if(mapping_mat_nu_sorted==TRUE && !is.null(nu_matcher)){
		stop("mapping_mat_nu_sorted and nu_matcher should not be set together.")
	}
	if(mapping_mat_nu_sorted==FALSE && is.null(nu_matcher)){
		stop("either mapping_mat_nu_sorted or nu_matcher should be set.")
	}
	if(!mapping_mat_nu_sorted){
		F=F[,nu_matcher,drop=FALSE]
	}
	if(add_intercept){
		F=cbind(TRUE,F)
	}
	F_first=F%*%Eu_nu_list[[2]]
	F_first=rowSums(F_first*F)
	firster=sqrt((Eu_alphai[,2]*F_first)/2)	
	non_zeros2=which(V!=0,2)
	preadder=sparseMatrix(i=non_zeros2[,1],j=non_zeros2[,2],x=firster[non_zeros2[,1]],dims=dim(V))
	Vprime=preadder*V
	VptVp=t(Vprime)%*%Vprime
	my_trace=sum(Eu_omega_list[[2]]*as.matrix(VptVp))
	myn=dim(Eu_alphai)[1]
	my_g=-0.5*log(2*pi)*myn+sum(0.5*Eu_alphai[,1])-my_trace
	return(my_g)
}

calc_g_nu_via_Eu=function(p,q=length(p),myc=NULL){
	if(length(p)==1 && is.null(myc)){
		myg= -0.5*log(2*pi)*q + (0.5*q)*log(p)-0.5*p*q
	}else{
		if(!is.null(myc)){
			myg= -0.5*log(2*pi)*q+0.5*sum(log(p))-0.5*sum(p*myc^2)
		}else{
			myg= -0.5*log(2*pi)*q+0.5*sum(log(p))-0.5*sum(p)
		}
	}
	return(myg)
}

calc_g_omega_via_Eu=function(mys,Eu_delta){
	myg= -0.5*log(2*pi)*mys + (0.5*Eu_delta[1])*mys
	return(myg)
}

calc_g_omega_via_Eu_upsilon=function(Eu_upsilon_list,D_mat){
	mys=dim(D_mat)[1]
	myq=dim(D_mat)[2]
	myg= -0.5*log(2*pi)
	res=rep(0,mys)
	for(j in c(1:myq)){
		res=res+Eu_upsilon_list[[j]][D_mat[,j],1]
	}
	res=0.5*res+myg
	return(res)
}

calc_g_delta_via_Eu=function(chi1,chi2){
	myg=chi1*log(chi2)-lgamma(chi1)
	return(myg)
}

calc_g_upsilon_via_Eu=function(chi1_vec,chi2_vec,h_vec){
	all_g=rep(0,length(chi2_vec))
	for(j in c(1:length(chi2_vec))){
		myg=chi1_vec[j]*log(chi2_vec[j])-lgamma(chi1_vec[j])
		all_g[j]=myg*h_vec[j]
	}
	return(all_g)
}

calc_g_chi2j_via_Eu=function(zeta1_vec,zeta2_vec){
	all_g=rep(0,length(zeta2_vec))
	for(j in c(1:length(zeta2_vec))){
		myg=zeta1_vec[j]*log(zeta2_vec[j])-lgamma(zeta1_vec[j])
		all_g[j]=myg
	}
	return(all_g)
}

calc_g_upsilon_via_Eu_chi2j=function(chi1_vec,Eu_chi2j,h_vec){
	logchi2_vec=Eu_chi2j[,1]
	all_g=rep(0,length(logchi2_vec))
	for(j in c(1:length(logchi2_vec))){
		myg=chi1_vec[j]*logchi2_vec[j]-lgamma(chi1_vec[j])
		all_g[j]=myg*h_vec[j]
	}
	return(all_g)
}

calc_g_upsilon_via_theta_list=function(theta_list){
	all_g=rep(0,length(theta_list))
	for(j in c(1:length(theta_list))){
		cur_theta=theta_list[[j]]
		all_gs=(cur_theta[,1]+1)*log(-cur_theta[,2])-lgamma(cur_theta[,1]+1)
		gout_j=sum(all_gs)
		all_g[j]=gout_j
	}
	return(all_g)
}


calc_g_chi2j_via_theta=function(theta){
	all_gs=(theta[,1]+1)*log(-theta[,2])-lgamma(theta[,1]+1)
	gout=sum(all_gs)
	return(gout)
}

get_local_L_sum=function(gene_dat_list){
	if(length(gene_dat_list)==0){
		return(0)
	}
	gene_L_list=lapply(gene_dat_list,function(x){
		return(x$L_list)
	})
	return(sum(unlist(gene_L_list)))
}


get_global_L_sum=function(L_global_list){
	if(length(L_global_list)==0){
		return(0)
	}
	return(sum(unlist(L_global_list)))
}


calc_L_alpha_minus_star=function(gamma1,Eu_ak,mapping_mat,Eu_alphai,Eu_kappa,theta_alphai_first_row){
	Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat)
	E_g_alphai_par=calc_g_alphai_via_Eu_gammai_Eu_kappa(gamma1,Eu_gammai,Eu_kappa)
	theta_alphai_vect_second=calc_theta_alphai_via_Eu_gammai_Eu_kappa_second(Eu_gammai,Eu_kappa)
	theta_alphai=cbind(theta_alphai_first_row,theta_alphai_vect_second)
	L_alphai_minus_star=sum(rowSums(theta_alphai*Eu_alphai)+E_g_alphai_par)
	return(L_alphai_minus_star)
}


calc_theta_alpha_ak_via_EugammaiTimesEu_alphai_vect_only_second_with_preindexing_only_ones=function(Eu_gammai,gamma1,Eu_alphaiTimesEu_kappa,Eu_ak,k,non_null_mapping_mat_inds){

	sec_col=-sum(Eu_gammai[non_null_mapping_mat_inds,2]*Eu_alphaiTimesEu_kappa[non_null_mapping_mat_inds,2])/Eu_ak[k,2]
	return(sec_col)
}

run_cycle_across_allgenes_mpi_loc=function(mpicl=NULL){
	if(is.null(mpicl)){
		stop("no mpi cluster provided.")
	}
	Eu_alphaiTimesEu_kappa=clusterCall_or_envCall(mpicl,to_be_run_on_node_loc)
	return(Eu_alphaiTimesEu_kappa)
}


load_mapping_shift_data_to_cluster_node=function(mapping_shift_mat_expanded){
	assign("mapping_shift_mat_expanded",mapping_shift_mat_expanded,envir=get_data_container_environment())
	return(TRUE)
}

load_mapping_data_to_cluster_node=function(mapping_mat_expanded){
	assign("mapping_mat_expanded",mapping_mat_expanded,envir=get_data_container_environment())
	return(TRUE)
}

load_Eu_alphai_to_cluster_node=function(Eu_alphai){
	assign("Eu_alphai",Eu_alphai,envir=get_data_container_environment())
	return(TRUE)
}
load_lambda1_list_to_cluster_node=function(lambda1){
	assign("lambda1",lambda1,envir=get_data_container_environment())
	return(TRUE)
}
load_index_mat_list_to_cluster_node=function(index_mat){
	assign("index_mat",index_mat,envir=get_data_container_environment())
	return(TRUE)
}
load_index_mat_t_list_to_cluster_node=function(index_mat_t){
	assign("index_mat_t",index_mat_t,envir=get_data_container_environment())
	return(TRUE)
}

load_takvf_list_to_cluster_node=function(theta_alphai_kappa_vect_first){
	assign("theta_alphai_kappa_vect_first",theta_alphai_kappa_vect_first,envir=get_data_container_environment())
	return(TRUE)
}

load_theta_alphai_star_vect_to_cluster_node=function(theta_alphai_star_vect){
	assign("theta_alphai_star_vect",theta_alphai_star_vect,envir=get_data_container_environment())
}

to_be_run_on_node_loc=function(){
	### Assumes that EU_alphai was updated
	current_round=get("current_round",get_data_container_environment())#update local
	calc_B_rounds=get("calc_B_rounds",get_data_container_environment())#const
	nodenr=get("nodenr",get_data_container_environment())#const
	######################
	######################
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
	######################
	######################
	calc_L=get("calc_L",get_data_container_environment())#const
	calc_L_rounds=get("calc_L_rounds",get_data_container_environment())#const
	gamma1=get("gamma1",get_data_container_environment())#const
	theta_lambda=get("theta_lambda",get_data_container_environment())#update global
	Eu_lambda2=get("Eu_lambda2",get_data_container_environment())#update global
	approximation_mode=get("approximation_mode",get_data_container_environment())#const
	Eu_omega_list=get("Eu_omega_list",get_data_container_environment())
	Eu_nu_list=get("Eu_nu_list",get_data_container_environment())
	nu_matcher=get("nu_matcher",get_data_container_environment())
	ak_matcher=get("ak_matcher",get_data_container_environment())
	#update_global
#	Eu_alphai_vect=get("Eu_alphai",get_data_container_environment())
	Eu_tau2=get("Eu_tau2",get_data_container_environment())#update global
	tau1=get("tau1",get_data_container_environment())#const
	lambda1=get("lambda1",get_data_container_environment())#const
	mapping_mat_expanded=get("mapping_mat_expanded",get_data_container_environment())#const
	index_mat=get("index_mat",get_data_container_environment())#const
	index_mat_t=get("index_mat_t",get_data_container_environment())#const
	mapping_shift_mat_expanded=get("mapping_shift_mat_expanded",get_data_container_environment())#const
	Eu_ak=get("Eu_ak",get_data_container_environment())#update global
	theta_alphai_kappa_vect_first=get("theta_alphai_kappa_vect_first",get_data_container_environment())#const
	theta_alphai_star_vect=get("theta_alphai_star_vect",get_data_container_environment())#first column const second gets overwritten everytime
	shift_mat_melted=get("shift_mat_melted",get_data_container_environment())
	######################
	######################
	if((current_round %% calc_L_rounds)==0 || current_round==1){
		local_calc_L=TRUE
	}else{
		local_calc_L=FALSE
	}
	if(calc_L==FALSE){
		local_calc_L=FALSE
	}
	gene_dat_list=get("gene_dat_list",get_data_container_environment())
#	rm("gene_da_list",envir = get_data_container_environment())
#	mygc_res=gc()
	##############
	##############
	##############
	##############
	Eu_alphai_vect=my_docall_rbind(lapply(gene_dat_list,function(x){x$E_list$Eu_alphai}))
	if(approximation_mode){
		theta_index=match(names(gene_dat_list),rownames(theta_lambda))
		out_list=lapply(c(1:length(gene_dat_list)),function(i){
			run_gene_cycle_wrapper_lambda_b_nu(gene_dat_list[[i]],Eu_omega_list,Eu_nu_list=Eu_nu_list,nu_matcher=nu_matcher,theta_lambda=theta_lambda[theta_index[i],],Eu_lambda2=Eu_lambda2,calc_L=local_calc_L,approximation_mode=approximation_mode)
		})
		names(out_list)=names(gene_dat_list)
	}else{		
		#stop("erropxwoexsd:not yet implemented.")
		out_list=lapply(gene_dat_list,run_gene_cycle_wrapper_lambda_b_nu,Eu_omega_list,Eu_nu_list=Eu_nu_list,nu_matcher=nu_matcher,theta_lambda=theta_lambda,Eu_lambda2=Eu_lambda2,calc_L=local_calc_L,approximation_mode=approximation_mode)	
	}
	##############
	##############
	##############
	names_outl=names(out_list)
	for(myname in names(out_list)){
		gene_dat_list[[myname]]$E_list$Eu_b_list=out_list[[myname]]$E_list$Eu_b_list
		gene_dat_list[[myname]]$E_list$Eu_lambda=out_list[[myname]]$E_list$Eu_lambda
		gene_dat_list[[myname]]$theta_b_alphai=out_list[[myname]]$theta_b_alphai
		gene_dat_list[[myname]]$theta_b_nu=out_list[[myname]]$theta_b_nu
		if(local_calc_L){
			gene_dat_list[[myname]]$L_list$L_lambda_star_part=out_list[[myname]]$L_list$L_lambda_star_part
			gene_dat_list[[myname]]$L_list$L_b_star_part=out_list[[myname]]$L_list$L_b_star_part
			gene_dat_list[[myname]]$L_list$L_obs=out_list[[myname]]$L_list$L_y_obs
			gene_dat_list[[myname]]$L_aux$processed_Eu_b2_diag=out_list[[myname]]$L_aux$processed_Eu_b2_diag
		}
	}
	rm(list="out_list")
	mygc_res=gc()
	########################################
	## pull Eu_lambda_vect as an attribute
	Eu_lambda_list=lapply(names_outl,function(myname){
		gene_dat_list[[myname]]$E_list$Eu_lambda}
		)
	Eu_lambda_vect=my_docall_rbind(Eu_lambda_list,nsplits=10)
	########################################
	all_theta_b_alphai_second=lapply(gene_dat_list,function(x){return(x$theta_b_alphai[,2,drop=F])})
	all_theta_b_alphai_vect=get_all_theta_b_alphai_vect(all_theta_b_alphai_second)
	if(length(gene_dat_list[[1]]$E_list$Eu_b_list)==2){
		all_Eu_b1=lapply(gene_dat_list,function(x){return(x$E_list$Eu_b_list[[1]])})
	}else{
		all_Eu_b1=lapply(gene_dat_list,function(x){return(x$E_list$Eu_b_list[[3]])})
	}
	all_Eu_b1=my_docall_rbind(all_Eu_b1,nsplits=10)
	if(local_calc_L){
		all_Eu_b2_diag=lapply(gene_dat_list,function(x){return(as.matrix(x$L_aux$processed_Eu_b2_diag))})
		all_Eu_b2_diag=my_docall_rbind(all_Eu_b2_diag,nsplits=10)
	}
	######################
	######################
	##todo Mtot?
	myMs=unlist(lapply(c(1:length(gene_dat_list)),function(i){length(gene_dat_list[[i]]$Obs_list$mapping_mat[,1])}))
	Mtot=sum(myMs)
	all_theta_b_alphai_gene_indices=get("all_theta_b_alphai_gene_indices",get_data_container_environment())
	all_theta_b_alphai_second_list=list()
	all_theta_b_alphai_second_list[[1]]=all_theta_b_alphai_vect
	all_theta_b_alphai_gene_indices_list=list()
	first_ind=all_theta_b_alphai_gene_indices[1,1]
	all_theta_b_alphai_gene_indices_list[[1]]=all_theta_b_alphai_gene_indices-first_ind+1

	all_theta_b_alphai_vect_second=paste_theta_b_alphai_vect_second(all_theta_b_alphai_second_list,all_theta_b_alphai_gene_indices_list,Mtot)
	###################
	###################
	Eu_gammai_vect=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat_expanded[,ak_matcher])
	n_genes=length(gene_dat_list)
	###################
	theta_kappa_pergene_vect=calc_theta_kappaj_via_Eu_tau2(tau1,Eu_tau2,n_genes)
	###################
	theta_alphai_kappa_vect_second=calc_theta_alphai_kappa_via_Eu_gammai_Eu_alphai_second_from_mat(Eu_gammai_vect,Eu_alphai_vect,index_mat)	
	###################
	theta_alphai_kappa_vect=cbind(theta_alphai_kappa_vect_first,theta_alphai_kappa_vect_second)
	theta_kappa_star_vect=theta_kappa_pergene_vect+theta_alphai_kappa_vect
	##############################
	## update Eu_kappa
	##############################
	Eu_kappa_pergene_vect=calc_Eu_kappa_via_theta_vectorized(theta_kappa_star_vect)
	Eu_kappa_vect=collapse_Eu_kappa_vect(index_mat_t,Eu_kappa_pergene_vect)
	##############################
	## END: update Eu_kappa
	##############################
	theta_alphai_vect_second=calc_theta_alphai_via_Eu_gammai_Eu_kappa_second(Eu_gammai_vect,Eu_kappa_vect)
	theta_alphai_star_vect[,2]=theta_alphai_vect_second+all_theta_b_alphai_vect_second
	Eu_alphai_vect=calc_Eu_alphai_via_theta_vectorized(theta_alphai_star_vect)
	######################
	######################
	#######	#######	#######	#######	#######
	first_ind=all_theta_b_alphai_gene_indices[1,1]
	gene_indices_vec_boundaries=all_theta_b_alphai_gene_indices-first_ind+1
	#return(gene_indices_vec_boundaries)
	mylen=length(gene_indices_vec_boundaries[,1])
	for(i in c(1:mylen)){
		down=gene_indices_vec_boundaries[i,1]
		upper=gene_indices_vec_boundaries[i,2]
		gene_dat_list[[i]]$E_list$Eu_alphai[,2]=Eu_alphai_vect[c(down:upper),2]
	}
	#######	#######	#######	#######	#######
	mylen=dim(Eu_kappa_pergene_vect)[1]
	for(i in c(1:mylen)){
		gene_dat_list[[i]]$E_list$Eu_kappa=Eu_kappa_pergene_vect[i,]
	}
	######################
	#######	#######	#######	#######	#######
	theta_lambdai_lambda2=calc_theta_lambdai_lambda2_via_Eu_lambdai(lambda1=lambda1,Eu_lambdai=Eu_lambda_vect,m=dim(Eu_lambda_vect)[1],genespecific_lambda1=approximation_mode)
	n_genes=dim(Eu_kappa_pergene_vect)[1]
	theta_kappaj_tau2=calc_theta_kappaj_tau2_via_Eu_kappaj(tau1,Eu_kappaj=Eu_kappa_pergene_vect,m=n_genes)
	##############
	##############
	theta_b_nu=calc_theta_b_nu_via_Eu_vectorized_fast(all_Eu_b1=all_Eu_b1,all_Eu_alphai=Eu_alphai_vect,Eu_omega_list=Eu_omega_list,mapping_shift_mat_expanded=mapping_shift_mat_expanded,mapping_mat_expanded=mapping_mat_expanded,nu_matcher=nu_matcher,add_intercept=TRUE,shift_mat_melted=shift_mat_melted)
	######################
	######################
	Eu_alphai_vectTimesEu_kappa_vect=Eu_kappa_vect*Eu_alphai_vect
	attr(Eu_alphai_vectTimesEu_kappa_vect,"theta_b_nu")=theta_b_nu
	attr(Eu_alphai_vectTimesEu_kappa_vect,"theta_kappaj_tau2")=theta_kappaj_tau2
	attr(Eu_alphai_vectTimesEu_kappa_vect,"theta_lambdai_lambda2")=theta_lambdai_lambda2
	attr(Eu_alphai_vectTimesEu_kappa_vect,"all_Eu_b1")=all_Eu_b1
	attr(Eu_alphai_vectTimesEu_kappa_vect,"Eu_alphai_vect")=Eu_alphai_vect

	if(local_calc_L){
		L_kappa_star_part=sum(theta_kappa_star_vect*Eu_kappa_pergene_vect)+sum(apply(theta_kappa_star_vect,1,calc_g_kappa_via_theta))
		L_alphai_star_part=sum(theta_alphai_star_vect*Eu_alphai_vect)+sum(calc_g_alphai_via_theta(theta_alphai_star_vect))
		attr(Eu_alphai_vectTimesEu_kappa_vect,"L_kappa_star_part")=L_kappa_star_part
		attr(Eu_alphai_vectTimesEu_kappa_vect,"L_alphai_star_part")=L_alphai_star_part
		attr(Eu_alphai_vectTimesEu_kappa_vect,"all_Eu_b2_diag")=all_Eu_b2_diag
		attr(Eu_alphai_vectTimesEu_kappa_vect,"Eu_kappa_pergene_vect")=Eu_kappa_pergene_vect
		attr(Eu_alphai_vectTimesEu_kappa_vect,"Eu_kappa_vect")=Eu_kappa_vect
		attr(Eu_alphai_vectTimesEu_kappa_vect,"Eu_lambda_vect")=Eu_lambda_vect
	}
	assign("current_round",current_round,envir=get_data_container_environment())
	assign("gene_dat_list",gene_dat_list,envir=get_data_container_environment())
	assign("all_Eu_b1",all_Eu_b1,envir=get_data_container_environment())
	assign("Eu_alphai_vect",Eu_alphai_vect,envir=get_data_container_environment())
	######################
	######################
	return(Eu_alphai_vectTimesEu_kappa_vect)
}


run_calc_theta_b_omega_via_Eu_nu_on_node=function(){
	all_Eu_b1=get("all_Eu_b1",get_data_container_environment())
	Eu_alphai_vect=get("Eu_alphai_vect",get_data_container_environment())
	mapping_shift_mat_expanded=get("mapping_shift_mat_expanded",get_data_container_environment())
	mapping_mat_expanded=get("mapping_mat_expanded",get_data_container_environment())
	mapping_mat_expanded=get("mapping_mat_expanded",get_data_container_environment())
	Eu_nu_list=get("Eu_nu_list",get_data_container_environment())
	nu_matcher=get("nu_matcher",get_data_container_environment())
	theta_b_omega=calc_theta_b_omega_via_Eu_nu_vectorized_fast(all_Eu_b1=all_Eu_b1,all_Eu_alphai=Eu_alphai_vect,mapping_shift_mat_expanded=mapping_shift_mat_expanded,mapping_mat_expanded=mapping_mat_expanded,Eu_nu_list=Eu_nu_list,nu_matcher,add_intercept=TRUE)
	return(theta_b_omega)
}

####### not tested
calc_theta_bi_via_Eu_wshift_nu_diag=function(Eu_alphai,Eu_omega_list,mapping_shift_mat,Eu_nu_list,mapping_mat,nu_matcher=NULL,add_intercept=FALSE,mapping_mat_nu_sorted=FALSE){
	if(mapping_mat_nu_sorted==TRUE && !is.null(nu_matcher)){
		stop("mapping_mat_nu_sorted and nu_matcher should not be set together.")
	}
	if(mapping_mat_nu_sorted==FALSE && is.null(nu_matcher)){
		stop("either mapping_mat_nu_sorted or nu_matcher should be set.")
	}
	if(mapping_mat_nu_sorted){
		nu_matcher=c(1:dim(mapping_mat)[2])
	}
	tmp=mapping_shift_mat%*%Eu_omega_list[[1]]
	mylen=length(Eu_nu_list[[1]])
	if(add_intercept){
		tmp2=(as.matrix(mapping_mat[,nu_matcher])%*%Eu_nu_list[[1]][2:mylen])+Eu_nu_list[[1]][1]
	}else{
		tmp2=as.matrix(mapping_mat[,nu_matcher]%*%Eu_nu_list[[1]])
	}
	theta_first=as.matrix(Eu_alphai[,2]*tmp2*tmp)
	if(length(Eu_alphai[,2])==1){
		stop("err20kdssf")
	}
	theta_sec=-Eu_alphai[,2]/2
	theta_b_list=list()
	theta_b_list[["first"]]=theta_first
	theta_b_list[["sec_diag"]]=theta_sec
	return(theta_b_list)
}


run_gene_cycle_lambda_b_nu=function(Eu_omega_list,Eu_nu_list,nu_matcher,mapping_mat,yty,XtX_sqrt,ytX,N,theta_lambda,mapping_shift_mat,Eu_b_list,Eu_alphai,Eu_lambda,calc_L=TRUE,eig_values,approximation_mode,p_group_vec_gene=NULL){
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
		Eu_b_list=expand_Eu_b_list(Eu_b_list,k=dim(XtX_sqrt)[2])##prepares empty Eu_b matrix if doesnt exist yet. otherwise makes dimensionaliy check 
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
	#
	theta_b_diag=calc_theta_bi_via_Eu_wshift_nu_diag(Eu_alphai=Eu_alphai,Eu_omega_list,mapping_shift_mat=mapping_shift_mat,Eu_nu_list=Eu_nu_list,mapping_mat=mapping_mat,nu_matcher=nu_matcher,add_intercept=TRUE)
	theta_y_b_sqrt_abs=calc_theta_y_b_via_Eu_sqrt_abs(Eu_lambda,ytX,XtX_sqrt)
	theta_b_star_list=list()
	theta_b_star_list[["formula"]]="diag(sec_diag)-sec_sqrt_abs%*%t(sec_sqrt_abs)"
	theta_b_star_list[["formula_expanded"]]="diag(theta_b_star_list$sec_diag)-theta_b_star_list$sec_sqrt_abs%*%t(theta_b_star_list$sec_sqrt_abs)"
	theta_b_star_list[["first"]]=theta_b_diag[["first"]]+theta_y_b_sqrt_abs[["first"]]
	theta_b_star_list[["sec_diag"]]=theta_b_diag[["sec_diag"]]
	theta_b_star_list[["sec_sqrt_abs"]]=theta_y_b_sqrt_abs[["sec_sqrt_abs"]]
	##
	Eu_b_list=calc_Eu_b_via_theta4eigen(theta_b_star_list,decomposed=TRUE,missing_variance=missing_var,myc=-Eu_lambda[2]/2)
	##########################
	##########################	
	theta_b_alphai=calc_theta_bi_alphai_via_Eu_wshift_nu(m=length(Eu_b_list[[1]]),Eu_b_list=Eu_b_list,mapping_shift_mat=mapping_shift_mat,Eu_omega_list=Eu_omega_list,decomposed=TRUE,mapping_mat=mapping_mat[,nu_matcher],Eu_nu_list=Eu_nu_list,add_intercept=TRUE)

	if(calc_L){
			L_lambda_star_part=-sum((theta_lambda_star)*Eu_lambda)-calc_g_lambda_via_theta(theta_lambda_star)
			if(length(Eu_b_list)==2){
				stop("not aligned anymore with length(Eu_b_list)==4")
				theta_b_star_list[["sec"]]=diag(theta_b_star_list[["sec_diag"]])-theta_b_star_list[["sec_sqrt_abs"]]%*%t(theta_b_star_list[["sec_sqrt_abs"]])
				L_b_star_part=sum(-(theta_b_star_list[["first"]])*(Eu_b_list[[1]]))+sum(-unlist(theta_b_star_list[["sec"]])*unlist(Eu_b_list[[2]]))-calc_g_b_via_theta4eigen(theta_b_star_list)	
				L_obs=calc_L_y_observed4eigen(yty,ytX,XtX_sqrt=XtX_sqrt,Eu_b_list,Eu_lambda,N,decomposed=FALSE,missing_variance=missing_var)
			}else{
				processed_Eu_b_list=list()
				processed_Eu_b_list[[1]]=as.matrix(Eu_b_list[["Eu_b_1"]])
				processed_Eu_b_list[[2]]=as.matrix(0.5*(diag(Eu_b_list$theta_M_inv$diagDprime_alphainv)-Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt%*%t(Eu_b_list$theta_M_inv$cXtX_sqrtAug_interimRes_sqrt)) + Eu_b_list$Eu_b_1 %*% t(Eu_b_list$Eu_b_1))
				fake_theta_list=list()
				fake_theta_list[[1]]=as.matrix(theta_b_star_list[["first"]])
				fake_theta_list[[2]]=diag(theta_b_star_list$sec_diag)-theta_b_star_list$sec_sqrt_abs%*%t(theta_b_star_list$sec_sqrt_abs)
				if(missing_var!=0){
					mym=length(fake_theta_list[[1]])
					myadder=(missing_var/mym)*Eu_lambda[2]/2
					diag(fake_theta_list[[2]])=diag(fake_theta_list[[2]])-myadder
					theta_b_star_list[["sec_diag"]]=theta_b_star_list[["sec_diag"]]-myadder
				}
				#####
				L_b_star_part=sum(-unlist(fake_theta_list)*unlist(processed_Eu_b_list))-calc_g_b_via_theta4eigen(theta_b_star_list)
				#L_Q=L_b_star_part
				mym=dim(processed_Eu_b_list[[2]])[1]
				if(mym==1){
					stop("nsnps has to be larger than 1")
				}
				if(approximation_mode){
					L_obs=calc_L_z_observed4eigen(ZtSigma_invZ=yty,ytX,XtX_sqrt=XtX_sqrt,Eu_b_list,Eu_lambda,N,missing_variance=missing_var)
				}else{
					L_obs=calc_L_y_observed4eigen(yty,ytX,XtX_sqrt=XtX_sqrt,Eu_b_list,Eu_lambda,N,decomposed=TRUE,missing_variance=missing_var)
				}
			}
	}
############
	##clean up Eu_b_list=
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
			L_lambda="",
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
		out_list$L_list$L_y_obs=L_obs
		out_list$L_aux$processed_Eu_b2_diag=diag(processed_Eu_b_list[[2]])
	}
	return(out_list)
}
