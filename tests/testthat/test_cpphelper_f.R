context("test compiled_functions.")
require(Matrix)
require(Rcpp)

test_getDiagMtCM_optim=function(){
	set.seed(11)
	ind2sub=function(dims=c(0,0),ind){
	    m=dims[1]
	    r = ((ind-1) %% m) + 1
	    c1 = floor((ind-1) / m) + 1
	    return(cbind(r,c1))
	}

	melt_f=function(myM){
		tt2_2=which(myM!=0,2)
		tt2=which(myM!=0)
		new_vals=myM[tt2]
		myMmelted=cbind(tt2_2,new_vals)	
		return(myMmelted)
	}

	myn=500
	mym=700
	nvals=myn*0.05*mym
	totn=myn*mym
	all_vals=ceiling(runif(nvals)*totn)
	tt=ind2sub(c(mym,myn),all_vals)
	vals=rnorm(length(all_vals))
	myM=sparseMatrix(i=tt[,1],j=tt[,2],x=vals,dims=c(mym,myn))

	myMmelted=melt_f(myM)
	Cmat=myM[,1:200]%*%t(myM[,1:200])
	Cmat=as.matrix(Cmat)
	lenM=dim(myM)[2]
	gg=getDiagMtCM_optim(Cmat=Cmat,sparseM=myMmelted,lenM=lenM)
	AA=(Cmat%*%myM)
	AA=t(myM)%*%AA
	gg2=diag(AA)
	return(max(abs((gg-gg2)/gg2)))
}



test_that("getDiagMtCM_optim should .",{

	ret=test_getDiagMtCM_optim()
	expect_equal(ret, 0 + 1e-14)
})

