
get_objsizes=function(envir=parent.frame()){
	all_obj=ls(envir=envir,all.names=TRUE)
	all_sizes=sapply(all_obj,function(x){
		obj=get(x,envir)
		size=pryr::object_size(obj)
	})
	out=data.table(name=all_obj,size=all_sizes)
	setkey(out,size)
	return(out)
}

get_objsizes2=function(envir=parent.frame()){
	return(ls(parent.env(envir)))
}


get_objsizes_round3=function(){
	round=12
	myls_list=list()
	envir=parent.frame()
	for(i in c(1:round)){
		envir=parent.env(envir)
		myls_list[[i]]=ls(envir)
	}
	return(myls_list)
}


get_objsizes_round2=function(){
	round=12
	myls_list=list()
	envir=parent.frame()
	for(i in c(1:round)){
		envir=parent.frame(i)
		all_names=ls(envir)
		all_names=all_names[grepl("^[a-zA-Z]",all_names,perl=TRUE)]
  		all_names=all_names[!grepl("<-",all_names,perl=TRUE)]
 		all_sizes=sapply(all_names,function(x){
			obj=get(x,envir)
			size=pryr::object_size(obj)
		})
		out=data.table(name=all_names,size=all_sizes)
		setkey(out,size)
		myls_list[[i]]=out
	}
	return(myls_list)
}


get_objsizes_round=function(){
	round=12
	myls_list=list()
	envir=parent.frame()
	for(i in c(1:round)){
		all_names=ls(envir)
		all_names=all_names[grepl("^[a-zA-Z0-9_]+$",all_names,perl=TRUE)]
		if(length(all_names)>0){
	  		all_sizes=sapply(all_names,function(x){
				obj=get(x,envir)
				size=pryr::object_size(obj)
			})
			out=data.table(name=all_names,size=all_sizes)
			setkey(out,size)
		}else{
			out=NULL
		}
		myls_list[[i]]=out
		envir=parent.env(envir)
	}
	return(myls_list)
}


get_objsizes3=function(envir=parent.frame()){
	return(ls(parent.env(parent.env(envir))))
}

get_objsizes4=function(envir=parent.frame()){
	return(ls(parent.env(parent.env(parent.env(envir)))))
}

get_objsizes5=function(envir=parent.frame()){
	return(ls(parent.env(parent.env(parent.env(parent.env(envir))))))
}

get_objsizes6=function(envir=parent.frame()){
	return(ls(parent.env(parent.env(parent.env(parent.env(parent.env(envir)))))))
}

get_objsizes7=function(envir=parent.frame()){
	return(ls(parent.env(parent.env(parent.env(parent.env(parent.env(parent.env(envir))))))))
}




get_objsizes_full_stack=function(obj_sizes_list=list(),env = parent.frame()){
	if(identical(env, emptyenv())) {
    	return(obj_sizes_list)
 	}
 	len=length(obj_sizes_list)
 	obj_sizes_list[[len+1]]=get_objsizes(env)
 	get_objsizes_full_stack(obj_sizes_list,parent.env(env))
}


get_objsizes_full_stack_wnamespace=function(obj_sizes_list=list(),env = parent.frame()){
	if(identical(env, emptyenv())) {
    	return(obj_sizes_list)
 	}
 	len=length(obj_sizes_list)
 	appender=list()
 	appender[["package_env"]]=get_objsizes(env)
 	all_names=appender[["package_env"]][,name]
 	all_names=all_names[grepl("^[a-zA-Z]",all_names,perl=TRUE)]
  	all_names=all_names[!grepl("<-",all_names,perl=TRUE)]
 	print("sss")
 	namespace_envs=sapply(all_names,function(myname){environment(eval(parse(text=myname)))})
 	namespace_envs=unlist(unique(namespace_envs))
 	if(length(namespace_envs)!=1){
 		browser()
 		stop("err:should have only one enclosing env for all functions.")
 	}
 	namespace_env=namespace_envs[[1]]
 	if(!is.null(namespace_env)){
 		appender[["namespace_env"]]=get_objsizes(namespace_env)
 		appender[["imports_env"]]=get_objsizes(parent.env(namespace_env))
 	}
 	len=length(obj_sizes_list)
 	obj_sizes_list[[len+1]]=appender
 	get_objsizes_full_stack_wnamespace(obj_sizes_list,parent.env(env))
}




# args: mat: n*m matrix where each col is a phenotype 
# args: vec: n*p matrix where each col is a covariate.
# args: add_intercept: should intercept be added.
# out residual mat of regressing out all covariates for each phenotype col separately.
regress_out_vects=function(mat,vects,add_intercept=FALSE){
    if(dim(mat)[1]!=dim(vects)[1]){
        stop("erramcvsp: dimensions dont agree.")
    }
    if(add_intercept){
        vects=cbind(1,vects)
    }
    myinv=solve(t(vects)%*%vects)
    interim=t(vects)%*%mat
    mybetas=myinv%*%interim
    mypred=vects%*%mybetas
    myResidMat=mat-mypred
    return(myResidMat)
}


my_cov=function(x,use_minus1=FALSE){
	xm=colMeans(x)
	mylen=dim(x)[1]
	ret=t(x)%*%x
	ret=ret/mylen
	ret=ret-(xm%*%t(xm))
	if(use_minus1){
		ret=ret*mylen/(mylen-1)
	}
	return(ret)
}

my_scale=function(x,scale=TRUE,use_minus1=FALSE){
	xm=colMeans(x)
	xm2=colMeans(x^2)
	myvar=xm2-xm^2
	mylen=dim(x)[1]
	if(use_minus1){
		mysd=sqrt(myvar)*sqrt(mylen/(mylen-1))
	}else{
		mysd=sqrt(myvar)
	}
	mMat=rep(1,mylen)%*%t(xm)
	x=x-mMat
	rm("mMat")
	sdMat=rep(1,mylen)%*%t(mysd)
	x=x/sdMat
	return(x)
}

calc_he=function(y,X){
	XXt=X%*%t(X)/dim(X)[2]
	mylen=dim(X)[2]
	XXt=cov(t(X))
	yyt=y%*%t(y)
	logic_mat=upper.tri(XXt)
	yvec=yyt[logic_mat]
	xvec=XXt[logic_mat]
	myb=sum(yvec*xvec)/sum(xvec^2)
	return(myb)
}

calc_he_test=function(){
	mygvar=0.5
	round=300
	gout=rep(0,round)
	gout_tr=rep(0,round)
	X=matrix(rnorm(500*200),500,200)
	for(i in c(1:round)){
		print(i)
		myb=(rnorm(200)/sqrt(dim(X)[2]))*sqrt(mygvar)
		mye=rnorm(500)*sqrt(1-mygvar)
		gvec=X%*%myb
		myy=gvec+mye
		ret=calc_he(myy,X)
		#ret=calc_he(scale(myy),scale(X))
		gout_tr[i]=var(gvec)/var(myy)
		gout[i]=ret
	}
}

is_ascending_integer_vec=function(x){
	problem=sum(diff(sort(unique(x)))!=1)!=0
	if(problem){
		stop("erro2m[ex: D_mat has to contain integer indices starting at one and ascending in steps of one.")
	}
	return(problem)
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
	colnames(my_mat)=names_present_in_all
	return(my_mat)
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

check_range=function(x,varname,lower,upper){
	for(el in x){
		if(el<lower || el>upper){
			stop(paste(varname," argument not within range",sep=""))
		}
	}
}

load_as=function(filepath){
	tt=load(filepath)
	eval(parse(text=paste("out=",tt[1],sep="")))
	return(out)
}
