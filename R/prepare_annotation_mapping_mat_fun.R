#' prepares annotation data oject. 
#' 
#' Maps SNPs to genes in cis and processes annotations. Undirected annotations can be given as bed files. Directed annotations can be given as bed files or as already processed sparse R matrices (see argument \code{SHIFT_RDATA_FILES_PATH_CHR22}). See \code{https://genome.ucsc.edu/FAQ/FAQformat.html#format1} for bed format definition.
#' 
#' @return
#' An R-object of class \code{bagea_annotation_mat}; an \code{R}-list contatining three objects. 
#' \describe{\item{\code{snp_tss_map}}{A data table for annotation of snp-gene symbol pairs.}
#' \item{\code{annotation_mat}}{A 0-1 (sparse) matrix of same column number as snp_tss_map that annotates snp-gene pairs as overlapping a given undirected annotation. Column names are same as the names of the bed-files used (without .bed file format suffix), except for the tss range annotations, which have the format tssrange_(num), (replace (num) with the respective integer range around the tss.}
#' \item{\code{annotation_shift_mat}}{A numerical (sparse) matrix of same column number as snp_tss_map that annotates snp-gene pairs as overlapping a given directed annotation.}
#' \item{\code{settings}}{A list of parameter settings with which prepare_annotation_mapping_mat was run.}}
#' @param RANGE_AROUND_TSS
#' Scalar between 10'000 and 500'000. Window around transcription start site for which snps are considered.
#' @param TSS_CUTOFFS
#' Vector of positive numbers between 100 and value set in \code{RANGE_AROUND_TSS} in ascending order.
#' Defines windows around tss for which separate annotations are produced. Separate annotations are produced for upstream and downstream windows.
#' @param PREPARE_SHIFT_MAT
#' Boolean. Set true if a directed annotation matrix should be built. 
#' @param SHIFT_BED_FILES_PATHS
#' Single character stringor NULL. The path where directed bed files are stored.
#' @param SHIFT_BED_FILES
#' A vector of stringor NULL. Contains the bed-file names that should be processed. Files have to be present in \code{SHIFT_BED_FILES_PATHS}.
#' @param SHIFT_RDATA_FILES_PATH_CHR22
#' Single character string or NULL. Path to a sparse matrix (see package Matrix) with row names rsids and column names directed annotations. The matrix represents the directed annotations for chromosome 22. The filepath has to contain the substring "chr22" and expects that analoguous files are available for chr1 through chr22 in the same path. 
#' @param NCORES
#' Positive integer scalar. Number of cores to be used when loading the data.
#' @param PREPARE_SIMPLE_BED_MAT
#' Boolean. Set true if additional undirected annotations should be considered  
#' in addition to windows around tss'.
#' @param SIMPLE_BED_FILES_PATHS
#' Single character string.
#' The path where bed files are stored.
#' @param SIMPLE_BED_FILES
#' A character vector or NULL. Contains the bed-file names that should be processed. Files have to be present in \code{SIMPLE_BED_FILES_PATHS}.
#' @export
prepare_annotation_mapping_mat=function(RANGE_AROUND_TSS=150000,TSS_CUTOFFS=c(250,500,1000,2000,5000,10000,20000,50000),SHIFT_BED_FILES_PATHS=NULL,SHIFT_BED_FILES=NULL,SHIFT_RDATA_FILES_PATH_CHR22=NULL, PREPARE_SHIFT_MAT=TRUE,NCORES=1, PREPARE_SIMPLE_BED_MAT=FALSE,SIMPLE_BED_FILES_PATHS=NULL,SIMPLE_BED_FILES=NULL){
	print("running prepare_annotation_mapping_mat..")
	## check argumnets
	print("checking arguments..")
	## check SIMPLE_BED_FILES_PATHS
	BAGEA_PATH=Sys.getenv("BAGEA_PATH")
	check_installation()
	bedops_path=paste(BAGEA_PATH,"/bedops_dir/bin/",sep="")
	## check argumnets
	## check PREPARE_SIMPLE_BED_MAT
	if(!file.exists(paste(bedops_path,"/bedops",sep=""))){
		stop("bedops does not seem correctly installed. Have you set BAGEA_PATH and run install_external_data?")
	}
	if(!is.scalar(PREPARE_SIMPLE_BED_MAT) || !is.logical(PREPARE_SIMPLE_BED_MAT)){
		stop("PREPARE_SIMPLE_BED_MAT has to be a single boolean.")
	}
	if(PREPARE_SIMPLE_BED_MAT){
		if(!is.scalar(SIMPLE_BED_FILES_PATHS) || !is.character(SIMPLE_BED_FILES_PATHS)){
			stop("SIMPLE_BED_FILES_PATHS argument not a single character string")
		}
		if(!is.vector(SIMPLE_BED_FILES) || !is.character(SIMPLE_BED_FILES)){
			stop("SIMPLE_BED_FILES argument not a single character vector")
		}
	}
	###
	if(PREPARE_SHIFT_MAT==TRUE){
		shift_bed_counter=(!is.null(SHIFT_BED_FILES_PATHS))+(!is.null(SHIFT_BED_FILES))
		if(is.null(SHIFT_RDATA_FILES_PATH_CHR22) & shift_bed_counter==0){
			SHIFT_RDATA_FILES_PATH_CHR22=paste(BAGEA_PATH,"/Data/DirectedAnnotations/expecto_loading_normed_chr22.RData",sep="")
			print(paste("PREPARE_SHIFT_MAT is set TRUE  but no bed files or shift annotation matrices given. Trying to use the default path",SHIFT_RDATA_FILES_PATH_CHR22))
		}
		if(!is.null(SHIFT_BED_FILES_PATHS) && (!is.scalar(SHIFT_BED_FILES_PATHS) || !is.character(SHIFT_BED_FILES_PATHS))){
				stop("SHIFT_BED_FILES_PATHS argument not NULL or a single character string")
		}
		if(!is.null(SHIFT_BED_FILES) && (!is.vector(SHIFT_BED_FILES) || !is.character(SHIFT_BED_FILES))){
			stop("SHIFT_BED_FILES argument not a single character vector")
		}
		if(shift_bed_counter==1){
			stop("SHIFT_BED_FILES SHIFT_BED_FILES_PATHS should either both NULL or both not NULL.")
		}
		if(!is.null(SHIFT_RDATA_FILES_PATH_CHR22) && (!is.scalar(SHIFT_RDATA_FILES_PATH_CHR22) || !is.character(SHIFT_RDATA_FILES_PATH_CHR22))){
			stop("SHIFT_RDATA_FILES_PATH_CHR22 argument not a single character string")
		}
		if(!is.null(SHIFT_RDATA_FILES_PATH_CHR22)){
			if(!grepl("chr22",SHIFT_RDATA_FILES_PATH_CHR22)){
				stop("SHIFT_RDATA_FILES_PATH_CHR22 has to be either NULL or contain the substring chr22.")
			}
		}
		if(!is.null(SHIFT_RDATA_FILES_PATH_CHR22) && shift_bed_counter==2){
				stop("either SHIFT_RDATA_FILES_PATH_CHR22 or SHIFT_BED_FILES_PATHS has to be NULL.")
		}
	}
	if(!is.numeric(RANGE_AROUND_TSS)){
		stop("RANGE_AROUND_TSS argument not a numeric.")
	}
	if(!is.scalar(RANGE_AROUND_TSS)){
		stop("RANGE_AROUND_TSS should be single scalar.")
	}
	if(RANGE_AROUND_TSS<10000  || RANGE_AROUND_TSS>500000){
		stop("TSS_CUTOFFS has to be between 10'000 and 500'000.")
	}
	## check argumnets
	## check TSS_CUTOFFS
	if(!is.numeric(TSS_CUTOFFS)){
		stop("TSS_CUTOFFS argument not a numeric.")
	}
	if(!is.vector(TSS_CUTOFFS)){
		stop("TSS_CUTOFFS should be a numeric vector.")
	}
	if(length(TSS_CUTOFFS)>1){
		if(min(diff(TSS_CUTOFFS)<0)){
			stop("TSS_CUTOFFS has to be in assending order.")
		}
	}
	if(min(TSS_CUTOFFS)<100  || max(TSS_CUTOFFS)>RANGE_AROUND_TSS){
		stop("all TSS_CUTOFFS have to be between 100 and RANGE_AROUND_TSS argument.")
	}
	## check PREPARE_SIMPLE_BED_MAT
	#############################################
	settings=list()
	settings[["SIMPLE_BED_FILES_PATHS"]]=SIMPLE_BED_FILES_PATHS
	settings[["SIMPLE_BED_FILES"]]=SIMPLE_BED_FILES
	settings[["RANGE_AROUND_TSS"]]=RANGE_AROUND_TSS
	settings[["TSS_CUTOFFS"]]=TSS_CUTOFFS
	settings[["PREPARE_SIMPLE_BED_MAT"]]=PREPARE_SIMPLE_BED_MAT
	settings[["PREPARE_SHIFT_MAT"]]=PREPARE_SHIFT_MAT
	settings[["SHIFT_BED_FILES_PATHS"]]=SHIFT_BED_FILES_PATHS
	settings[["SHIFT_BED_FILES"]]=SHIFT_BED_FILES
	settings[["SHIFT_RDATA_FILES_PATH_CHR22"]]=SHIFT_RDATA_FILES_PATH_CHR22
	#############################################
	range_around_tss_str=sprintf("%d",RANGE_AROUND_TSS)
	#############################################
	#############################################
	## annotate tss
	#############################################
	#############################################
	print("prepare tss ranges ..")
	mytmp=tempfile()
	expand_echomap_to_single_bed_elements_cmd_str="perl -nle '!/\\|$/ && print'| perl -ple 's/;/$1\\n/mg' | perl -ne 'if(/(chr.+\\|)/){$a=$1;print}else{print \"$a$_\"}'"
	snps_in_tss_cmd=paste(bedops_path,"/bedmap --ec --range ",range_around_tss_str," --echo --echo-map ",get_commonsnp_filepath()," ",get_tss_filepath()," |",expand_echomap_to_single_bed_elements_cmd_str," | perl -ple 's/\\|/\\t/' | cut -f1,2,4,8,10,11 > ",mytmp,sep="")

	a=system(snps_in_tss_cmd)		
	if(a!=0){
		stop(paste("couldn't make snps_in_tss_",tss_cutoff,".bed",sep=""))
	}

	snp_tss_map=fread(mytmp)
	setnames(snp_tss_map,c("V1","V2","V3","V4","V5","V6"),c("chrom","snp_pos","snpid","tss_pos","symbol","strand"))
	system(paste("rm ",mytmp,sep=""))
	distance=snp_tss_map[,(snp_pos-tss_pos)]
	minus_strand=snp_tss_map[,strand=="-"]
	distance[minus_strand]=distance[minus_strand]*(-1)
	len_1=dim(snp_tss_map)[1]
	len_2=length(TSS_CUTOFFS)

	indices_pos_list=list()
	for(i in c(1:len_2)){
		print(i)
		indices_pos_list[[i]]=cbind(i,which((distance<TSS_CUTOFFS[i]) & distance>=0))
	}
	combined=do.call("rbind",indices_pos_list)
	tss_pos_mat=sparseMatrix(i=combined[,2],j=combined[,1], dims=c(len_1,len_2))
	colnames(tss_pos_mat)=paste("tssrange_pos",TSS_CUTOFFS,sep="")
	rm("indices_pos_list")
	rm("combined")
	#######
	indices_neg_list=list()
	for(i in c(1:len_2)){
		print(i)
		indices_neg_list[[i]]=cbind(i,which((distance> -TSS_CUTOFFS[i]) & distance<0))
	}
	combined=do.call("rbind",indices_neg_list)
	tss_neg_mat=sparseMatrix(i=combined[,2],j=combined[,1], dims=c(len_1,len_2))
	colnames(tss_neg_mat)=paste("tssrange_neg",TSS_CUTOFFS,sep="")
	rm("indices_neg_list")
	rm("combined")
	tss_mat=cbind(tss_neg_mat,tss_pos_mat)
	tss_mat=tss_mat==TRUE
	rm("tss_neg_mat")
	rm("tss_pos_mat")
	########################################
	########################################
	########################################
	########################################
	## TODO extract into full script.
	########################################
	########################################
	#############################################
	# prepare undirected snp-wise bed files
	#############################################
	if(PREPARE_SIMPLE_BED_MAT){
		print("prepare simple bed files ..")
		bedfiles_to_load=SIMPLE_BED_FILES
		all_bedfiles_present=list.files(SIMPLE_BED_FILES_PATHS)
		nr_of_bedfiles_to_load=length(bedfiles_to_load)
		if(sum(grepl("^tssrange",bedfiles_to_load))){
			print("one bedfile start with 'tssrange'")
			print("this is not allowed")
			stop("abort")
		}
		if(nr_of_bedfiles_to_load==0){
			print(paste("no bedfiles to load. either set PREPARE_SIMPLE_BED_MAT to FALSE or add file names to your simple bed file"))
			stop()
		}
		nr_of_bedfiles_missing=length(setdiff(bedfiles_to_load,all_bedfiles_present))
		if(nr_of_bedfiles_missing>0){
			print("following bedfiles missing:")
			print(setdiff(bedfiles_to_load,all_bedfiles_present))
			stop()
		}
		########################################
		simple_bed_file_list=list()
		for(i in c(1:length(bedfiles_to_load))){
			print(i)
			print(paste("process ",bedfiles_to_load[i],sep=""))
			mytmpfile=tempfile()
			sortcmdstr=paste(bedops_path,"/sort-bed ",SIMPLE_BED_FILES_PATHS,"/",bedfiles_to_load[i], " > ",mytmpfile,sep="")
			a=system(sortcmdstr)
			if(a!=0){
				cmdind=system(paste("rm ",mytmpfile))
				stop("sorting command did not work correctly")
			}
			cmdstr=paste(bedops_path,"/bedops --ec -e ",get_commonsnp_filepath()," ",mytmpfile,sep="")
			annot_dat=fread(cmdstr)
			cmdind=system(paste("rm ",mytmpfile))
			simple_bed_file_list[[i]]=cbind(i,which(is.element(snp_tss_map[,snpid],annot_dat[,V4])))
		}
		combined=do.call("rbind",simple_bed_file_list)
		simple_bed_mat=sparseMatrix(i=combined[,2],j=combined[,1], dims=c(len_1,length(simple_bed_file_list)))
		colnames(simple_bed_mat)=sub("\\.bed$","",bedfiles_to_load)	
		tss_mat_names=colnames(tss_mat)
		simple_bed_mat_names=colnames(simple_bed_mat)
		myintersect=intersect(tss_mat_names,simple_bed_mat_names)
		if(length(myintersect)>0){
			print("the annotation mapping output does not have unique column names, something went wrong.")
			stop("abort")
		}
	}
	#############################################
	# END: prepare undirected snp-wise bed files
	#############################################
	#############################################
	# prepare directed snp-wise bed files
	#############################################
	if(PREPARE_SHIFT_MAT){
		print("prepare shift mat ..")
		if(!is.null(SHIFT_BED_FILES)){
			print("prepare shift bed files ..")
			bedfiles_to_load=SHIFT_BED_FILES
			all_bedfiles_present=list.files(SHIFT_BED_FILES_PATHS)
			nr_of_bedfiles_to_load=length(bedfiles_to_load)
			if(nr_of_bedfiles_to_load==0){
				print(paste("no bedfiles to load. either set PREPARE_SHIFT_MAT to FALSE or add file names to your directed bed file"))
				stop()
			}
			nr_of_bedfiles_missing=length(setdiff(bedfiles_to_load,all_bedfiles_present))
			if(nr_of_bedfiles_missing>0){
				print("following bedfiles missing:")
				print(setdiff(bedfiles_to_load,all_bedfiles_present))
				stop()
			}
			########################################
			#shift_bed_file_list=list()
			shift_bed_file_list=mclapply(c(1:length(bedfiles_to_load)),function(i){
			#for(i in c(1:length(bedfiles_to_load))){
				print(i)
				print(paste("process ",bedfiles_to_load[i],sep=""))
				mytmpfile=tempfile()
				sortcmdstr=paste(bedops_path,"/sort-bed ",SHIFT_BED_FILES_PATHS,"/",bedfiles_to_load[i], " > ",mytmpfile,sep="")
				a=system(sortcmdstr)
				if(a!=0){
					cmdind=system(paste("rm ",mytmpfile))
					stop("sorting command did not work correctly")
				}
	#			cmdstr=paste(bedops_path,"/bedops --ec -e ",mytmpfile," ",get_commonsnp_filepath(),sep="")
				cmdstr=paste(bedops_path,"/bedmap --ec --echo --mean ",get_commonsnp_filepath()," ",mytmpfile," | perl -nle  '!/NAN$/mg && print' | perl -nle  's/\\|/\\t/ && print'",sep="")
				annot_dat=fread(cmd=cmdstr)
				cmdind=system(paste("rm ",mytmpfile))

				indices=which(is.element(snp_tss_map[,snpid],annot_dat[,V4]))
				matcher=match(snp_tss_map[indices,snpid],annot_dat[,V4])
				values=annot_dat[matcher,V7]
				index_value_map=cbind(cbind(i,indices),values)
				return(index_value_map)
			},mc.cores=NCORES)
			#	shift_bed_file_list[[i]]=index_value_map
			#}
			combined=do.call("rbind",shift_bed_file_list)
			shift_mat=sparseMatrix(i=combined[,2],j=combined[,1],x=combined[,3],dims=c(len_1,length(shift_bed_file_list)))
			colnames(shift_mat)=sub("\\.bed$","",bedfiles_to_load)
			tss_mat_names=colnames(tss_mat)
			shift_mat_names=colnames(shift_mat)
			myintersect=intersect(tss_mat_names,shift_mat_names)
			if(length(myintersect)>0){
				print("the annotation mapping output does not have unique column names, something went wrong.")
				stop("abort")
			}
		}else{
			print("use RDATA files to create SHIFT_MAT")
			shift_rdata_files_paths=sapply(c(1:22),function(x){
				gsub("chr22",paste("chr",x,sep=""),SHIFT_RDATA_FILES_PATH_CHR22)
			})
			print("check if files are present.")
			missing_files_l=lapply(shift_rdata_files_paths,function(x){
				if(!file.exists(x)){
					return(x)
				}else{
					return(NULL)
				}
			})
			all_missing_files=do.call("c",missing_files_l)
			if(length(all_missing_files)>0){
				print("expecting the following files but can't find the following")
				lapply(all_missing_files,print)
				stop("abort.")
			}
			myfun=function(mychr,all_filepaths,snp_tss_map){
				print(mychr)
				snp_tss_map_chr=snp_tss_map[chrom==mychr,]
				filepath=all_filepaths[grepl(mychr,all_filepaths)][1]
				ww=load(filepath)
				if(ww!="shift_mat"){
					print(paste("problem: file",filepath,"does not contain an R object of name shift_mat."))
					return(NULL)
				}
				if(sum(is.element("dgCMatrix",class(shift_mat)))!=1){
					print(paste("problem: shift_mat from",filepath,"has to be of class dgCMatrix."))
					return(NULL)
				}
				shift_snpids=rownames(shift_mat)
				ee=is.element(snp_tss_map_chr[,snpid],shift_snpids)
				snp_tss_map_chr=snp_tss_map_chr[ee,]
				mymatch=match(snp_tss_map_chr[,snpid],shift_snpids)
				shift_mat=shift_mat[mymatch,]
				return(shift_mat)
			}
			all_chrs=sort(paste("chr",c(1:22),sep=""))
			if(NCORES==1)
				shift_mat=lapply(all_chrs,all_filepaths=shift_rdata_files_paths,snp_tss_map=snp_tss_map,myfun)
			else{
				shift_mat=mclapply(all_chrs,all_filepaths=shift_rdata_files_paths,myfun,snp_tss_map=snp_tss_map,mc.cores=NCORES)
			}
			print("combining matrices")
			shift_mat=do.call("rbind",shift_mat)
			shift_snpids=rownames(shift_mat)
			tokeep=is.element(snp_tss_map[,snpid],shift_snpids)
			snp_tss_map=snp_tss_map[tokeep,]
			tss_mat=tss_mat[tokeep,]
			if(PREPARE_SIMPLE_BED_MAT){
				simple_bed_mat=simple_bed_mat[tokeep,]
			}
			if(mean(snp_tss_map[,snpid]==rownames(shift_mat))!=1){
				stop("erro2edw:problem.")
			}
		}
	}
	#############################################
	#############################################
	# END: prepare directed snp-wise bed files
	#############################################
	#############################################
	annotation_mat=list()
	annotation_mat[["snp_tss_map"]]=snp_tss_map
	if(PREPARE_SIMPLE_BED_MAT){
		annotation_mat[["annotation_mat"]]=cbind(tss_mat,simple_bed_mat)
		annotation_mat[["annotation_mat"]]=annotation_mat[["annotation_mat"]]==TRUE
	}else{
		annotation_mat[["annotation_mat"]]=tss_mat
	}
	if(PREPARE_SHIFT_MAT){
		annotation_mat[["annotation_shift_mat"]]=shift_mat
	}
	annotation_mat[["settings"]]=settings
	class(annotation_mat)=c("bagea_annotation_mat",class(annotation_mat))
	return(annotation_mat)
}


