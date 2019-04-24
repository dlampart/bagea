context("gamma nodes.")
#########################################
## gamma test helper functions
#########################################
## out:(E[log(x)],E[x])
get_gammaE_from_pars=function(shape,rate){
	elogx=digamma(shape)-log(rate)
	ex=shape/rate
	cbind(elogx,ex)
}
## out:(E_x[log(P(x|shape,rate)])
get_gamma_Entropy_from_pars=function(shape,rate){
	entropy=shape-log(rate)+lgamma(shape)+(1-shape)*digamma(shape)
	return(entropy)
}

## set y = rate for var x then
## out E_y[log(P(x|shape,y))]
get_E_yLogPx=function(x,shape,E_y){
	out=shape*E_y[1]-lgamma(shape)+(shape-1)*log(x)-E_y[2]*x
	return(out)
}
#########################################
## END: gamma test helper functions
#########################################

#########################################
## gamma root test functions
#########################################
# testing Eu_x;theta_x
test_root_gamma=function(calc_Eu_via_theta_x,calc_theta_x){
	ret1=calc_Eu_via_theta_x(calc_theta_x(0.3,2.4))
	ret2=get_gammaE_from_pars(shape=0.3,rate=2.4)
	test_res=sum(abs(ret1-ret2))
	return(test_res)
}

# testing g_x g_theta
test_root_gamma2=function(calc_g_x_via_Eu,calc_g_via_theta_x,calc_theta_x){	
	a1=calc_g_via_theta_x(calc_theta_x(0.3,2.6))
	a2=calc_g_x_via_Eu(0.3,2.6)
	ret=sum(abs(a1-a2))
	return(ret)
}

# testing logP
test_root_gamma3=function(calc_theta_x,calc_g_x_via_Eu){
	myshape=0.3
	myrate=2.6
	myx=6.6
	myvec=c(log(myx),myx)
	tt=sum(myvec*calc_theta_x(myshape,myrate))+calc_g_x_via_Eu(myshape,myrate)
	tt2=dgamma(x=myx,shape=myshape,rate=myrate,log=TRUE)
	ret=tt-tt2
	return(ret)
}
#########################################
## END: gamma root test functions
#########################################


#########################################
##  gamma child test functions
#########################################

### testing calc_theta_xy
test_child_gamma=function(calc_theta_xy){
	shape_x=0.5
	rate_x=2	
	Eu_x=get_gammaE_from_pars(shape_x,rate_x)
	negEntropy_x=-get_gamma_Entropy_from_pars(shape_x,rate_x)
	theta_xy=calc_theta_xy(shape_x=shape_x,Eu_x=Eu_x)
	theta_xy=as.vector(theta_xy)
	negEntropy_indir=theta_xy[1]*log(rate_x)-lgamma(shape_x)+(shape_x-1)*Eu_x[1,1]+rate_x*theta_xy[2]
	ret=negEntropy_indir-negEntropy_x
	return(ret)
}

# testing calc_Eu_x_via_theta_x;calc_theta_x_via_Eu_x
test_child_gamma2=function(calc_Eu_via_theta_x,calc_theta_x){
	shape_x=0.4
	Eu_y=cbind(0.5,2)
	ret1=calc_Eu_via_theta_x(calc_theta_x(shape_x=shape_x,Eu_y=Eu_y))
	ret2=get_gammaE_from_pars(shape=shape_x,rate=Eu_y[2])
	ret=sum(abs(ret1-ret2))
	return(ret)
}

# testing calc_g_theta via logP_theta
test_child_gamma3=function(calc_theta_x,calc_g_via_theta_x){
	shape_x=0.4
	Eu_y=cbind(0.5,2)
	myx=0.6
	tot_res=get_E_yLogPx(myx,shape_x,Eu_y)
	theta_x=calc_theta_x(shape_x,Eu_y)
	calc_g_x_theta=calc_g_via_theta_x(theta_x)
	myx=0.2
	myvec=c(log(myx),myx)
	tt=sum(myvec*theta_x)+calc_g_x_theta
	tt2=dgamma(x=myx,shape=shape_x,rate=Eu_y[,2],log=TRUE)
	ret=tt-tt2
	return(ret)
}

# testing g_Ex via ElogP
test_child_gamma4=function(calc_g_x_via_Eu){
	shape_x=0.4
	Eu_y=cbind(0.5,2)
	myx=0.2
	ret1=get_E_yLogPx(x=myx,shape=shape_x,E_y=Eu_y)
	ret2_tmp=calc_g_x_via_Eu(shape_x,Eu_y)
	ret2=ret2_tmp+(shape_x-1)*log(myx)-Eu_y[2]*myx
	out=ret1-ret2
	return(out)
}
#########################################
## END: gamma child test functions
#########################################
#########################################
## combined test functions.
#########################################
#### test_gamma_root_node
#### required function arguments:
#### calc_theta_x: shape_x(1),Eu_y(1*2)
#### calc_Eu_via_theta_x: theta(2)
#### calc_g_via_theta_x: theta(2)
#### calc_g_x_via_Eu: shape_x(1),Eu_y(1*2)
test_gamma_root_node=function(calc_theta_x,calc_Eu_via_theta_x,calc_g_via_theta_x,calc_g_x_via_Eu){
	ret1=test_root_gamma(calc_Eu_via_theta_x,calc_theta_x)
	ret2=test_root_gamma2(calc_g_x_via_Eu,calc_g_via_theta_x,calc_theta_x)
	ret3=test_root_gamma3(calc_theta_x,calc_g_x_via_Eu)
	out=c(ret1,ret2,ret3)
	names(out)=paste("test",c(1:length(out)))
	return(out)
}

#### test_gamma_child_node
#### required function arguments:
#### calc_theta_xy: shape_x(1),Eu_x(1*2)
#### calc_theta_x: shape_x(1),Eu_y(1*2)
#### calc_Eu_via_theta_x: theta(2)
#### calc_g_via_theta_x: theta(2)
#### calc_g_x_via_Eu: shape_x(1),Eu_y(1*2)
test_gamma_child_node=function(calc_theta_xy,calc_theta_x,calc_Eu_via_theta_x,calc_g_via_theta_x,calc_g_x_via_Eu){
	ret1=test_child_gamma(calc_theta_xy)
	ret2=test_child_gamma2(calc_Eu_via_theta_x,calc_theta_x)
	ret3=test_child_gamma3(calc_theta_x,calc_g_via_theta_x)
	ret4=test_child_gamma4(calc_g_x_via_Eu)
	out=c(ret1,ret2,ret3,ret4)
	names(out)=paste("test",c(1:length(out)))
	return(out)
}
#########################################
###################################
###### test node alphai (full)
###################################
test_alphai_node=function(){
	set.seed(12)
	nalphais=10
	naks=5
	gamma1=0.5
	mapping_mat=Matrix(matrix(rnorm(nalphais*naks)>0.5,nalphais,naks),sparse=TRUE)
	Eu_aks_1=rgamma(naks,shape=1,rate=1)
	Eu_ak=cbind(log(Eu_aks_1)-0.1,Eu_aks_1)
	Eu_kappa_1=rgamma(1,shape=1,rate=1)
	Eu_kappa=cbind(log(Eu_kappa_1)-0.1,Eu_kappa_1)
	Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat=mapping_mat)
	theta_alphai=calc_theta_alphai_via_Eu_gammai_Eu_kappa(gamma1,Eu_gammai,Eu_kappa)
	Eu_alphai=calc_Eu_alphai_via_theta(theta_alphai)
	Eu_kappa_expanded=rep(1,nalphais)%*%Eu_kappa
	g_alphai_pa=calc_g_alphai_via_Eu_gammai_Eu_kappa(gamma1,Eu_gammai,Eu_kappa_expanded)
	E_full_lPb_1=sum(theta_alphai*Eu_alphai)+sum(g_alphai_pa)
	theta_alphai_kappa=calc_theta_alphai_kappa_via_Eu_gammai_Eu_alphai(gamma1=gamma1,Eu_gammai=Eu_gammai,Eu_alphai=Eu_alphai)
	calc_g_alpha_kappa=function(gamma1,Eu_gammai,Eu_alphai){
		ret=gamma1*Eu_gammai[,1]-lgamma(gamma1)+(gamma1-1)*Eu_alphai[,1]
		return(ret)
	}
	g_alpha_kappa=calc_g_alpha_kappa(gamma1,Eu_gammai,Eu_alphai)
	E_full_lPb_2=sum(theta_alphai_kappa*Eu_kappa)+sum(g_alpha_kappa)
	out1=E_full_lPb_1-E_full_lPb_2
	myk=1
	out2=0
	for(myk in c(1:naks)){
		theta_alpha_ak=calc_theta_alpha_ak_via_Eugammai(Eu_gammai=Eu_gammai,gamma1=gamma1,Eu_alphai=Eu_alphai,Eu_kappa=Eu_kappa,Eu_ak=Eu_ak,mapping_mat=mapping_mat,k=myk)
		calc_g_alpha_ak=function(gamma1,mapping_mat,Eu_ak,Eu_alphai,Eu_kappa,myk){
			Eu_gamma_minusi=calc_Eu_gammai_via_Eu_ak(Eu_ak[-myk,],mapping_mat[,-myk])
			ret=gamma1*(Eu_kappa[1]+Eu_gamma_minusi[,1])-lgamma(gamma1)+(gamma1-1)*Eu_alphai[,1]
			ret[mapping_mat[,myk]==0]=0
			return(ret)
		}
		g_alpha_ak=calc_g_alpha_ak(gamma1,mapping_mat,Eu_ak,Eu_alphai,Eu_kappa,myk)
		##### !!! note: only active variables in myk contribute
		acvtive_vars=mapping_mat[,myk]
		E_full_lPb_1_k=sum(theta_alphai[acvtive_vars,]*Eu_alphai[acvtive_vars,])+sum(g_alphai_pa[acvtive_vars])
		E_full_lPb_1_k
		E_full_lPb_3=sum(theta_alpha_ak*Eu_ak[myk,])+sum(g_alpha_ak)
		out2=out2+sum(abs(E_full_lPb_3-E_full_lPb_1_k))
	}
	#### 2. use deterministic parents to compare Expectation with negative entropy.	#### 
	#### step1 make parent nodes deterministic
	Eu_kappa[1]=log(Eu_kappa[2])
	Eu_ak[,1]=log(Eu_ak[,2])
	Eu_gammai=calc_Eu_gammai_via_Eu_ak(Eu_ak,mapping_mat=mapping_mat)
	#if Eu_ak is deterministic(i.e. known) then Eu_gammai is deterministic. We check this
	out3=sum(abs(Eu_gammai[,1]-log(Eu_gammai[,2])))
	logX=Eu_gammai[,1]+Eu_kappa[1]
	X=Eu_gammai[,2]*Eu_kappa[2]
	negEntropy_x=-get_gamma_Entropy_from_pars(shape=gamma1,rate=X)
	#### step3 compare to results from functions.
	theta_alphai=calc_theta_alphai_via_Eu_gammai_Eu_kappa(gamma1,Eu_gammai,Eu_kappa)
	Eu_alphai=calc_Eu_alphai_via_theta(theta_alphai)
	Eu_kappa_expanded=rep(1,nalphais)%*%Eu_kappa
	g_alphai_pa=calc_g_alphai_via_Eu_gammai_Eu_kappa(gamma1,Eu_gammai,Eu_kappa_expanded)
	negEntropy_x_indir=rowSums(theta_alphai*Eu_alphai)+g_alphai_pa
	out4=sum(abs(negEntropy_x_indir-negEntropy_x))
	myout=c(out1,out2,out3,out4)
	return(myout)
}
###################################
###### END: test node alphai (full)
###################################
###################################
###### test node upsilon (full)
###################################
test_upsilon_node=function(){
	set.seed(11)
	h_vec=c(10,6)
	myq=length(h_vec)
	chi1_vec=c(0.3,0.4)
	Eu_chi2j_2=rgamma(myq,1,1)
	Eu_chi2j=cbind(log(Eu_chi2j_2)-0.2,Eu_chi2j_2)
	theta_upsilon_list=calc_theta_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec=h_vec)
	Eu_upsilon_list=lapply(theta_upsilon_list,function(x)calc_Eu_upsilon_via_theta(x))
	g_upsilon_pa=calc_g_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec)
	theta_upsilon_chi2j=calc_theta_upsilon_chi2j_via_Eu(Eu_upsilon_list,chi1_vec)
	calc_g_upsilon_chi2j=function(chi1_vec,Eu_upsilon_list){
		myq=length(chi1_vec)
		out=lapply(c(1:myq),function(i){
			return(-lgamma(chi1_vec[i])+(chi1_vec[i]-1)*Eu_upsilon_list[[i]][,1])
			})
		return(out)
	}
	g_upsilon_chi2j=calc_g_upsilon_chi2j(chi1_vec,Eu_upsilon_list)
	g_upsilon_chi2j=sapply(g_upsilon_chi2j,sum)
	E_full_lPb_1=sapply(c(1:myq),function(myk){
	sum(theta_upsilon_list[[myk]]*Eu_upsilon_list[[myk]])+g_upsilon_pa[myk]})
	E_full_lPb_2=rowSums(theta_upsilon_chi2j*Eu_chi2j)+g_upsilon_chi2j
	out1=sum(abs(E_full_lPb_1-E_full_lPb_2))
	#### 2. use deterministic parents to compare Expectation with negative entropy.	#### 
	#### step1 make parent nodes deterministic
	Eu_chi2j[,1]=log(Eu_chi2j[,2])
	theta_upsilon_list=calc_theta_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec=h_vec)
	Eu_upsilon_list=lapply(theta_upsilon_list,function(x)calc_Eu_upsilon_via_theta(x))
	g_upsilon_pa=calc_g_upsilon_via_Eu_chi2j(chi1_vec,Eu_chi2j,h_vec)
	E_full_lPb_1=sapply(c(1:myq),function(myk){
	sum(theta_upsilon_list[[myk]]*Eu_upsilon_list[[myk]])+g_upsilon_pa[myk]})
	#print(E_full_lPb_1)
	X=Eu_chi2j[,2]
	negEntropy_x=-get_gamma_Entropy_from_pars(shape=chi1_vec,rate=X)
	negent=negEntropy_x*h_vec
	out2=sum(abs(negent-E_full_lPb_1))
	out=c(out1,out2)
	return(out)
}
###################################
###### END: test node upsilon (full)
###################################
###################################
###### test node chi2j 
###################################
test_that("chi2j node functions are consistent.",{
	calc_theta_x=calc_theta_chi2j
	calc_Eu_via_theta_x=calc_Eu_chi2j_via_theta
	calc_g_via_theta_x=calc_g_chi2j_via_theta
	calc_g_x_via_Eu=calc_g_lambda2
	ret=test_gamma_root_node(calc_theta_x,calc_Eu_via_theta_x,calc_g_via_theta_x,calc_g_x_via_Eu)
	names(ret)=NULL
	expect_equal(ret[1],0,tolerance=1e-14)
	expect_equal(ret[2],0,tolerance=1e-14)
	expect_equal(ret[3],0,tolerance=1e-14)
})
###################################
###### END: test node chi2j 
###################################
###################################
###### test node kappaj
###################################
test_that("kappa node functions are consistent.",{
	calc_theta_xy=function(shape_x,Eu_x){
		return(calc_theta_kappaj_tau2_via_Eu_kappaj(tau1=shape_x,Eu_kappaj=Eu_x))
	}
	calc_theta_x=function(shape_x,Eu_y){
		calc_theta_kappaj_via_Eu_tau2(shape_x,Eu_y,1)
	}
	calc_Eu_via_theta_x=calc_Eu_kappa_via_theta
	calc_g_via_theta_x=calc_g_kappa_via_theta
	calc_g_x_via_Eu=calc_g_kappa_via_Eu_tau2
	ret=test_gamma_child_node(calc_theta_xy,calc_theta_x,calc_Eu_via_theta_x,calc_g_via_theta_x,calc_g_x_via_Eu)
	names(ret)=NULL
	expect_equal(ret[1],0,tolerance=1e-14)
	expect_equal(ret[2],0,tolerance=1e-14)
	expect_equal(ret[3],0,tolerance=1e-14)
	expect_equal(ret[4],0,tolerance=1e-14)
})
###################################
###### END: test node kappaj
###################################
###################################
###### test node lambda
###################################
test_that("lambda node functions are consistent.",{
	calc_theta_xy=function(shape_x,Eu_x){
		return(calc_theta_lambdai_lambda2_via_Eu_lambdai(lambda1=shape_x,Eu_lambdai=Eu_x))
	}
	calc_theta_x=function(shape_x,Eu_y){
		calc_theta_lambda(lambda1=shape_x,lambda2=Eu_y[2])
	}
	calc_Eu_via_theta_x=calc_Eu_lambda_via_theta
	calc_g_via_theta_x=calc_g_lambda_via_theta
	calc_g_x_via_Eu=calc_g_lambda_via_Eu_lambda2
	ret=test_gamma_child_node(calc_theta_xy,calc_theta_x,calc_Eu_via_theta_x,calc_g_via_theta_x,calc_g_x_via_Eu)
	names(ret)=NULL
	expect_equal(ret[1],0,tolerance=1e-14)
	expect_equal(ret[2],0,tolerance=1e-14)
	expect_equal(ret[3],0,tolerance=1e-14)
	expect_equal(ret[4],0,tolerance=1e-14)
})
###################################
###### END: test node lambda
###################################
###################################
###### test node lambda2
###################################
test_that("lambda2 node functions are consistent.",{
	calc_theta_x=calc_theta_lambda2
	calc_Eu_via_theta_x=calc_Eu_lambda2_via_theta
	calc_g_via_theta_x=calc_g_lambda2_via_theta
	calc_g_x_via_Eu=calc_g_lambda2
	ret=test_gamma_root_node(calc_theta_x,calc_Eu_via_theta_x,calc_g_via_theta_x,calc_g_x_via_Eu)
	names(ret)=NULL
	expect_equal(ret[1],0,tolerance=1e-12)
	expect_equal(ret[2],0,tolerance=1e-12)
	expect_equal(ret[3],0,tolerance=1e-12)
})
###################################
###### END: test node lambda2
###################################
###################################
###### test node upsilon (edge case)
###################################
test_that("upsilon node functions are consistent. (edge case)",{
	calc_theta_xy=function(shape_x,Eu_x){
		Eu_upsilon_list=list()
		Eu_upsilon_list[[1]]=Eu_x
		ret=calc_theta_upsilon_chi2j_via_Eu(Eu_upsilon_list,chi1_vec=shape_x)
		return(ret)
	}
	calc_theta_x=function(shape_x,Eu_y){
		out=calc_theta_upsilon_via_Eu_chi2j(chi1_vec=shape_x,Eu_chi2j=Eu_y,h_vec=1)
		return(out[[1]])
	}
	calc_Eu_via_theta_x=calc_Eu_upsilon_via_theta
	calc_g_via_theta_x=function(theta){
		theta_list=list()
		theta_list[[1]]=theta
		calc_g_upsilon_via_theta_list(theta_list)
	}
	calc_g_x_via_Eu=function(shape_x,Eu_y){
		return(calc_g_upsilon_via_Eu_chi2j(chi1_vec=shape_x,Eu_chi2j=Eu_y,h_vec=1))
	}
	ret=test_gamma_child_node(calc_theta_xy,calc_theta_x,calc_Eu_via_theta_x,calc_g_via_theta_x,calc_g_x_via_Eu)
	names(ret)=NULL
	expect_equal(ret[1],0,tolerance=1e-14)
	expect_equal(ret[2],0,tolerance=1e-14)
	expect_equal(ret[3],0,tolerance=1e-14)
	expect_equal(ret[4],0,tolerance=1e-14)
})
###################################
###### END: test node upsilon (edge case)
###################################
###################################
###### test node tau2
###################################
test_that("tau2 node functions are consistent.",{
	calc_theta_x=calc_theta_tau2
	calc_Eu_via_theta_x=calc_Eu_tau2_via_theta
	calc_g_via_theta_x=calc_g_tau2_via_theta
	calc_g_x_via_Eu=calc_g_tau2
	ret=test_gamma_root_node(calc_theta_x,calc_Eu_via_theta_x,calc_g_via_theta_x,calc_g_x_via_Eu)
	names(ret)=NULL
	expect_equal(ret[1],0,tolerance=1e-12)
	expect_equal(ret[2],0,tolerance=1e-12)
	expect_equal(ret[3],0,tolerance=1e-12)
})
###################################
###### END: test node tau2
###################################
###################################
###### test node ak
###################################
test_that("ak node functions are consistent.",{
	calc_theta_x=calc_theta_lambda
	calc_Eu_via_theta_x=function(theta){
		theta=t(theta)
		return(calc_Eu_ak_via_theta(theta))
	}
	calc_g_via_theta_x=function(theta){calc_g_ak_via_theta(t(theta))}
	calc_g_x_via_Eu=calc_g_tau2
	ret=test_gamma_root_node(calc_theta_x,calc_Eu_via_theta_x,calc_g_via_theta_x,calc_g_x_via_Eu)
	names(ret)=NULL
	expect_equal(ret[1],0,tolerance=1e-12)
	expect_equal(ret[2],0,tolerance=1e-12)
	expect_equal(ret[3],0,tolerance=1e-12)
})
###################################
###### END: test node ak
###################################
###################################
###### test node alphai
###################################
test_that("alphai node functions are consistent.",{
	ret=test_alphai_node()
	## increasing tolarance somewhat because we sum over multiple values that are close to 0.
	expect_equal(ret[1], 0 + 5e-13)
	expect_equal(ret[2], 0 + 5e-13)
	expect_equal(ret[3], 0 + 5e-13)
})
###################################
###### END: test node alphai
###################################
###################################
###### test node upsilon
###################################
test_that("upsilon node functions are consistent.",{
	ret=test_upsilon_node()
	## increasing tolarance somewhat because we sum over multiple values that are close to 0.
	expect_equal(ret[1], 0 + 1e-14)
	expect_equal(ret[2], 0 + 1e-14)
})
###################################
###### END: test node upsilon
###################################


