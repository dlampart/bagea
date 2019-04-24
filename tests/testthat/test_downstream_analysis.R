context("downstream analysis.")
require(Matrix)


test_that("predict_directed_fast.",{
	set.seed(11)
	n_snps=200
	n_dir_annotation=100
	n_annot=10
	n_indiv=350
	bagea_result=list()
	observed_data=list()
	omega1=matrix(rnorm(n_dir_annotation))
	nu1=matrix(rnorm(n_annot+1))
	bagea_result$Eu_omega_list[[1]]=omega1
	bagea_result$Eu_nu_list[[1]]=nu1
	
	observed_data$settings$KEEP_X_y=TRUE
	Vtmp=matrix(rnorm(n_snps*n_dir_annotation),n_snps,n_dir_annotation)
	Vtmp[abs(Vtmp)<1.5]=0
	V=Matrix(Vtmp,sparse=TRUE)
	###############
	obs_list_el=list()
	obs_list_el$mapping_shift_mat=V
	F=Matrix(matrix(rnorm(n_snps*n_annot),n_snps,n_annot)>1,sparse=TRUE)
	obs_list_el$mapping_mat=F
	F=cbind(TRUE,F)
	X=matrix(rnorm(n_snps*n_indiv),n_indiv,n_snps)
	Sigma=t(X)%*%X/n_indiv
	b_shift=(F%*%nu1)*(V%*%omega1)
	##### init: should have mse_dir=0
	b=b_shift
	y=X%*%b
	yty=sum(y^2)
	tonorm=sqrt(yty/n_indiv)
	X=X/tonorm #(normalize X s.t. yty/n=1)
	y=X%*%b
	ytX=t(y)%*%X
	##### END: init: should have mse_dir=0
	obs_list_el$X=X
	obs_list_el$y=y
	obs_list_el$ytX=ytX
	obs_list_el$n=n_indiv
	obs_list_el$XtX_sqrt=t(X)
	obs_list_el$eig_values=svd(X)$d^2
	###############
	obs_list=list()
	obs_list$gene1=obs_list_el
	observed_data$obs_list=obs_list	
	###############
	ww=predict_directed_fast(observed_data,bagea_result)	
	expect_equal(ww[,mse_dir],0,tolerance=1e-14)
	expect_equal(ww[,mse_dir_approx],0,tolerance=1e-14)
})
