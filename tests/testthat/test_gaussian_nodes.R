context("gaussian nodes.")
require(Matrix)
#require(mvtnorm)
#source("bagea/R/bagea_optimized.R")

### testing strategy:
### calculate E_lambda_b[log(P(z|b,lambda))]
### employing all relevant functions 
### calc_theta_z_lambda_via_Eu
get_gaussian_Entropy_from_pars=function(mu,prec){
	myk=dim(prec)[1]
	myldet=determinant(prec,logarithm=TRUE)$modulus
	entropy=myk/2+myk*log(2*pi)/2-0.5*myldet
	return(entropy)
}


test_z_obs_node=function(){
	z_obs=t(t(c(-1,3,2.1,-0.5)))
	PrecBase=solve(toeplitz(c(1,0.8,0.8^2,0.8^3)))
	Eu_lambda=c(0.1,1.2)
	E_y1_list=list()
	E_y1_list[[1]]=t(t(c(0.5,2,-5,0.5)))
	E_y1_list[[2]]=toeplitz(c(1,0.5,0.5^2,0.5^3))+E_y1_list[[1]]%*%t(E_y1_list[[1]])
#	E_y2_list=list()
#	E_y2_list[[1]]=t(t(c(2.5,-4,-1,0.3)))
#	E_y2_list[[2]]=toeplitz(c(1,0.7,0.7^2,0.7^3))+E_y2_list[[1]]%*%t(E_y2_list[[1]])
	n=1000
	ytX=t(z_obs*sqrt(n))### scale
	Sigma=solve(PrecBase)
	XtX=solve(PrecBase)*n
	Eu_b_list=list()
	Eu_b_list[[1]]=PrecBase%*%E_y1_list[[1]]/sqrt(n)
	Eu_b_list[[2]]=PrecBase%*%E_y1_list[[2]]%*%PrecBase/n
	theta_zlambda=calc_theta_z_lambda_via_Eu(myn=n,Eu_b_list=Eu_b_list,ytX=ytX,XtX=XtX,yty=n)
	get_g_zlambda=function(z_obs,XtX,myn){
		k=length(z_obs)
		ldet=determinant(XtX/myn,logarithm=TRUE)$modulus
		ret=-k*log(2*pi)/2-ldet/2
		return(ret)
	}
	get_g_zb=function(z_obs,XtX,myn,Eu_lambda){
		k=length(z_obs)
		ldet=determinant(XtX/myn,logarithm=TRUE)$modulus
		ret1=-k*log(2*pi)/2-ldet/2+k*Eu_lambda[1]/2
		ret2=-Eu_lambda[2]/2*t(z_obs)%*%solve(XtX/myn)%*%z_obs
		ret=ret1+ret2
		return(ret)
	}
	g_zlambda=get_g_zlambda(z_obs,XtX,myn=n)
	ElogPZgLambdaB_1=sum(theta_zlambda*Eu_lambda)+g_zlambda
	theta_z_b=calc_theta_z_b_via_Eu(Eu_lambda=Eu_lambda,ytX=ytX,XtX=XtX)
	g_zb=get_g_zb(z_obs,XtX,myn=n,Eu_lambda)
	ElogPZgLambdaB_2=sum(theta_z_b[[1]]*Eu_b_list[[1]])+sum(theta_z_b[[2]]*Eu_b_list[[2]])+g_zb
	ElogPZgLambdaB_3=calc_L_z_observed(ytX=ytX,XtX=XtX,Eu_b_list=Eu_b_list,Eu_lambda,myn=n,yty=n)
	out1=ElogPZgLambdaB_1-ElogPZgLambdaB_2
	out2=ElogPZgLambdaB_2-ElogPZgLambdaB_3
	######	######	######
	######	######	######
	#### 2. use deterministic parents to test with 
	#### standard density functions.
	Eu_lambda[1]=log(Eu_lambda[2])
	Eu_b_list[[2]]=Eu_b_list[[1]]%*%t(Eu_b_list[[1]])
	first=calc_L_z_observed(ytX=ytX,XtX=XtX,Eu_b_list=Eu_b_list,Eu_lambda,myn=n,yty=n)
	z_obs_mean=Sigma%*%Eu_b_list[[1]]*sqrt(n)
	z_obs_var=Sigma/Eu_lambda[2]
	logPZgLambdaB_1=mvtnorm::dmvnorm(z_obs[,1],mean=z_obs_mean[,1],sigma=z_obs_var,log=TRUE)
	logPZgLambdaB_2=calc_L_z_observed(ytX=ytX,XtX=XtX,Eu_b_list=Eu_b_list,Eu_lambda,myn=n,yty=n)
	out3=logPZgLambdaB_2-logPZgLambdaB_1
	out=c(out1,out2,out3)
	return(out)
}

test_y_obs_node=function(){
	set.seed(11)
	X=cbind(rnorm(100),cbind(rnorm(100),cbind(rnorm(100),rnorm(100))))
	b=t(t(rnorm(4)))
	k=length(b)
	y=X%*%b
	yty=(t(y)%*%y)[1]
	n=length(y)
	ytX=t(y)%*%X
	Eu_b_list=list()
	Eu_b_list[[1]]=b
	Eu_b_list[[2]]=b%*%t(b)+diag(k)
	XtX=t(X)%*%X
	Eu_lambda=cbind(-0.5,1)
	get_g_ylambda=function(myn){
		ret=-myn*log(2*pi)/2
		return(ret)
	}
	get_g_yb=function(myn,Eu_lambda,yty){
		ret=-myn*log(2*pi)/2+myn*Eu_lambda[1]/2-Eu_lambda[2]*yty/2
		return(ret)
	}
	theta_ylambda=calc_theta_y_lambda_via_Eu(myn=n,Eu_b_list=Eu_b_list,ytX=ytX,XtX=XtX,yty=yty)
	ElogPYgLambdaB_1=sum(theta_ylambda*Eu_lambda)+get_g_ylambda(n)
	theta_yb=calc_theta_y_b_via_Eu(Eu_lambda,ytX, XtX)
	g_yb=get_g_yb(myn=n,Eu_lambda,yty)
	ElogPYgLambdaB_2=sum(theta_yb[[1]]*Eu_b_list[[1]])+sum(theta_yb[[2]]*Eu_b_list[[2]])+g_yb
	ElogPYgLambdaB_3=calc_L_y_observed(yty,ytX,XtX,Eu_b_list,Eu_lambda,n=n)
	out1=ElogPYgLambdaB_1-ElogPYgLambdaB_2
	out2=ElogPYgLambdaB_2-ElogPYgLambdaB_3
	######	######	######
	######	######	######
	#### 2. use deterministic parents to test with 
	#### standard density functions.
	Eu_lambda[1]=log(Eu_lambda[2])
	Eu_b_list[[2]]=Eu_b_list[[1]]%*%t(Eu_b_list[[1]])
	z_obs_mean=X%*%Eu_b_list[[1]]
	z_obs_var=diag(n)*1/Eu_lambda[2]
	logPYgLambdaB_1=mvtnorm::dmvnorm(y[,1],mean=z_obs_mean[,1],sigma=z_obs_var,log=TRUE)
	logPYgLambdaB_2=calc_L_y_observed(yty,ytX,XtX,Eu_b_list,Eu_lambda,n=n)
	out3=logPYgLambdaB_2-logPYgLambdaB_1
	out=c(out1,out2,out3)
	return(out)
}

#### test b functions
test_b_node=function(){
	set.seed(11)
	nbs=10
	nomega=5
	nnu=3
	Eu_alphai_1=rgamma(10,shape=1,rate=1)
	Eu_alphai=cbind(log(Eu_alphai_1)-0.1,Eu_alphai_1)
	omega1=t(t(rnorm(nomega)))
	Eu_omega_list=list()
	Eu_omega_list[[1]]=omega1
	Eu_omega_list[[2]]=omega1%*%t(omega1)+diag(nomega)
	nu1=t(t(rnorm(nnu)))
	Eu_nu_list=list()
	Eu_nu_list[[1]]=nu1
	Eu_nu_list[[2]]=nu1%*%t(nu1)+diag(nnu)
	V=matrix(rnorm(nbs*nomega),nbs,nomega)
	F=Matrix(matrix(rnorm(nbs*nnu)>0,nbs,nnu),sparse=TRUE)
	####
	theta_b_list=calc_theta_bi_via_Eu_wshift_nu(Eu_alphai=Eu_alphai,Eu_omega_list=Eu_omega_list,mapping_shift_mat=V,Eu_nu_list=Eu_nu_list,mapping_mat_sub=F)
	Eu_b_list_from_params=calc_Eu_b_via_theta(theta_b_list)
	######	######	######	######
		######	######	######	######
	######	######	######	######
	mySigma=Eu_b_list_from_params[[2]]-Eu_b_list_from_params[[1]]%*%t(Eu_b_list_from_params[[1]])
	theta_b_test=calc_theta_b_from_natural(mu=Eu_b_list_from_params[[1]],Sigma=mySigma)
	out1=sum(abs(theta_b_list[[1]]-theta_b_test[[1]]))
	out1=out1+sum(abs(theta_b_list[[2]]-theta_b_test[[2]]))
	######	######	######	######
		######	######	######	######
	######	######	######	######
#	browser()
	g_b_pa=calc_g_b_via_Eu_wshift_nu(Eu_alphai=Eu_alphai,Eu_omega_list=Eu_omega_list,V=V,F=F,Eu_nu_list=Eu_nu_list)
	E_full_lPb_1=sum(theta_b_list[[1]]*Eu_b_list_from_params[[1]])+sum(theta_b_list[[2]]*Eu_b_list_from_params[[2]])+g_b_pa
	#####
	theta_bi_alphai=calc_theta_bi_alphai_via_Eu_wshift_nu(m=nbs,Eu_b_list=Eu_b_list_from_params,mapping_shift_mat=V,Eu_omega_list=Eu_omega_list,mapping_mat=F,Eu_nu_list=Eu_nu_list)
	calc_g_b_alpha=function(nbs){
		ret=-nbs*log(2*pi)/2
		return(ret)
	}
	g_b_alpha=calc_g_b_alpha(nbs)
	E_full_lPb_2=sum(theta_bi_alphai*Eu_alphai)+g_b_alpha
	#######
	theta_bi_omega=calc_theta_b_omega_via_Eu_nu(Eu_b_list=Eu_b_list_from_params,Eu_alphai=Eu_alphai,mapping_shift_mat=V,mapping_mat=F,Eu_nu_list=Eu_nu_list)
	calc_g_b_omega=function(nbs,Eu_b_list,Eu_alphai){
		ret=-nbs*log(2*pi)/2+sum(Eu_alphai[,1])/2-sum(Eu_alphai[,2]*diag(Eu_b_list[[2]]))/2
		return(ret)
	}
	g_b_omega=calc_g_b_omega(nbs,Eu_b_list_from_params,Eu_alphai)

	E_full_lPb_3=sum(theta_bi_omega[[1]]*Eu_omega_list[[1]])+sum(theta_bi_omega[[2]]*Eu_omega_list[[2]])+g_b_omega
	#######
	calc_g_b_nu=calc_g_b_omega
	g_b_nu=calc_g_b_nu(nbs,Eu_b_list_from_params,Eu_alphai)
	theta_bi_nu=calc_theta_bi_nu_via_Eu_wshift(Eu_b_list=Eu_b_list_from_params,Eu_alphai=Eu_alphai,Eu_omega_list=Eu_omega_list,mapping_shift_mat=V,mapping_mat=F)
	E_full_lPb_4=sum(theta_bi_nu[[1]]*Eu_nu_list[[1]])+sum(theta_bi_nu[[2]]*Eu_nu_list[[2]])+g_b_nu
	out2=E_full_lPb_1-E_full_lPb_2
	out3=E_full_lPb_1-E_full_lPb_3
	out4=E_full_lPb_1-E_full_lPb_4
	#### 2. use deterministic parents to compare Expectation with negative entropy.	#### 
	#### step1 make parent nodes deterministic
	Eu_alphai[,1]=log(Eu_alphai[,2])
#	Eu_nu_list[[1]]=Eu_nu_list[[1]]-Eu_nu_list[[1]]
#	Eu_omega_list[[1]]=Eu_omega_list[[1]]-Eu_omega_list[[1]]
	Eu_nu_list[[2]]=Eu_nu_list[[1]]%*%t(Eu_nu_list[[1]])
	Eu_omega_list[[2]]=Eu_omega_list[[1]]%*%t(Eu_omega_list[[1]])
	#### END step1 make parent nodes deterministic
	#### step2 calculate negEntropy
	mu=(F%*%Eu_nu_list[[1]])*(V%*%Eu_omega_list[[1]])
	prec=diag(Eu_alphai[,2])
	negEntropy=-get_gaussian_Entropy_from_pars(mu,prec)	
	#### step3 compare to results from functions.
	theta_b_list=calc_theta_bi_via_Eu_wshift_nu(Eu_alphai=Eu_alphai,Eu_omega_list=Eu_omega_list,mapping_shift_mat=V,Eu_nu_list=Eu_nu_list,mapping_mat_sub=F)
	Eu_b_list_from_params=calc_Eu_b_via_theta(theta_b_list)
	g_b_pa=calc_g_b_via_Eu_wshift_nu(Eu_alphai=Eu_alphai,Eu_omega_list=Eu_omega_list,V=V,F=F,Eu_nu_list=Eu_nu_list)
	E_full_lPb_5=sum(theta_b_list[[1]]*Eu_b_list_from_params[[1]])+sum(theta_b_list[[2]]*Eu_b_list_from_params[[2]])+g_b_pa
	out5=negEntropy-E_full_lPb_5
	out=c(out1,out2,out3,out4,out5)
	return(out)
}

#### test nu node functions
test_nu_node=function(){
	myc=c(1,1,3,45,3,2)
	p=rep(1,6)
	p[1]=50
	q=length(p)
	theta_nu_list=calc_theta_nu(p=p,q=q,myc=myc)
	Eu_nu_list=calc_Eu_nu_via_theta(theta_nu_list)
	out1=sum(abs(Eu_nu_list[[1]]-myc))
	g1=calc_g_nu_via_theta(theta_nu_list)
	g2=calc_g_nu_via_Eu(p=p,q=q,myc=myc)
	out2=(g1-g2)[1]
	prec=diag(p)
	mu=myc
	negEntropy_dir=-get_gaussian_Entropy_from_pars(mu=mu,prec=prec)
	negEntropy_indir=sum(theta_nu_list[[1]]*Eu_nu_list[[1]])+sum(theta_nu_list[[2]]*Eu_nu_list[[2]])+g1
	out3=negEntropy_dir-negEntropy_indir
	out=c(out1,out2,out3)
	return(out)
}


#### test nu node functions
test_omega_node=function(){
	set.seed(11)	
	d1=rep(c(1:5),10)
	maxd1=max(d1)
	d2=rep(c(1:10),5)
	maxd2=max(d2)
	D_mat=cbind(d1,d2)
	Eu_upsilon_list=list()
	rg=rgamma(maxd1,shape=1)
	Eu_upsilon_list[[1]]=cbind(log(rg)-0.2,rg)
	rg=rgamma(maxd2,shape=2)
	Eu_upsilon_list[[2]]=cbind(log(rg)-0.3,rg)
	theta_omega_list=calc_theta_omegas_via_Eu_upsilon(Eu_upsilon_list,D_mat)
	Eu_omega_list=calc_Eu_omega_via_theta(theta_omega_list)
	g_omega_pa=sum(calc_g_omega_via_Eu_upsilon(Eu_upsilon_list,D_mat))
	E_full_lPb_1=sum(theta_omega_list[[1]]*Eu_omega_list[[1]])+sum(theta_omega_list[[2]]*Eu_omega_list[[2]])+g_omega_pa
	calc_g_omegaupsilon=function(D_mat,Eu_upsilon_list,myotherj){
		nomega=dim(D_mat)[1]
		inds=D_mat[,myotherj]
		ret=sum(Eu_upsilon_list[[myotherj]][inds,1])
		ret=ret/2-log(2*pi)*nomega/2
		return(ret)
	}
	theta_upsilon1=calc_theta_omegas_upsilon_via_Eu(Eu_omega_list,Eu_upsilon_list,D_mat,1)
	g_omegaupsilon1=calc_g_omegaupsilon(D_mat,Eu_upsilon_list,myotherj=2)
	E_full_lPb_2=sum(theta_upsilon1*Eu_upsilon_list[[1]])+g_omegaupsilon1
	theta_upsilon2=calc_theta_omegas_upsilon_via_Eu(Eu_omega_list,Eu_upsilon_list,D_mat,2)
	g_omegaupsilon2=calc_g_omegaupsilon(D_mat,Eu_upsilon_list,myotherj=1)
	E_full_lPb_3=sum(theta_upsilon2*Eu_upsilon_list[[2]])+g_omegaupsilon2
	out1=E_full_lPb_1-E_full_lPb_2
	out2=E_full_lPb_1-E_full_lPb_3
	#### 2. use deterministic parents to compare Expectation with negative entropy.	#### 
	#### step1 make parent nodes deterministic
	Eu_upsilon_list[[1]][,1]=log(Eu_upsilon_list[[1]][,2])
	Eu_upsilon_list[[2]][,1]=log(Eu_upsilon_list[[2]][,2])
		#### step2 calculate negEntropy
	mu=rep(0,dim(D_mat)[1])
	prec=diag(Eu_upsilon_list[[1]][D_mat[,1],2]*Eu_upsilon_list[[2]][D_mat[,2],2])
	negEntropy=-get_gaussian_Entropy_from_pars(mu,prec)	
	#### step3 compare to results from functions.
	theta_omega_list=calc_theta_omegas_via_Eu_upsilon(Eu_upsilon_list,D_mat)
	Eu_omega_list=calc_Eu_omega_via_theta(theta_omega_list)
	g_omega_pa=sum(calc_g_omega_via_Eu_upsilon(Eu_upsilon_list,D_mat))
	E_full_lPb_4=sum(theta_omega_list[[1]]*Eu_omega_list[[1]])+sum(theta_omega_list[[2]]*Eu_omega_list[[2]])+g_omega_pa
	out3=negEntropy-E_full_lPb_4
	out=c(out1,out2,out3)
	return(out)
}

test_that("z observed node functions are consistent.",{
	ret=test_z_obs_node()
	expect_equal(ret[1], 0 + 1e-14)
	expect_equal(ret[2], 0 + 1e-14)
	expect_equal(ret[3], 0 + 1e-14)
})

test_that("y observed node functions are consistent.",{
	ret=test_y_obs_node()
	expect_equal(ret[1], 0 + 1e-14)
	expect_equal(ret[2], 0 + 1e-14)
	expect_equal(ret[3], 0 + 1e-14)
})

test_that("b node functions are consistent.",{
	ret=test_b_node()
	expect_equal(ret[1], 0 + 1e-14)
	expect_equal(ret[2], 0 + 1e-14)
	expect_equal(ret[3], 0 + 1e-14)
	expect_equal(ret[4], 0 + 1e-14)
	expect_equal(ret[5], 0 + 1e-14)
})

test_that("nu node functions are consistent.",{
	ret=test_nu_node()
	expect_equal(ret[1], 0 + 1e-14)
	expect_equal(ret[2], 0 + 1e-14)
	expect_equal(ret[3], 0 + 1e-14)
})

test_that("omega node functions are consistent.",{
	ret=test_omega_node()
	expect_equal(ret[1], 0 + 1e-14)
	expect_equal(ret[2], 0 + 1e-14)
	expect_equal(ret[3], 0 + 1e-14)
})


