get_s3_download_cmdstr=function(filepath="/bagea-data/bagea_data_freeze/annotation_files/ensg2symbol20190410.txt",outpath,use_s3_client=FALSE){
	if(use_s3_client){
		cmdstr_prefix="aws s3 cp s3:/"
		cmdstr=paste(cmdstr_prefix,filepath," ",outpath,sep="")
	}else{
		cmdstr_prefix=paste("wget -O ",outpath," https://s3-us-west-1.amazonaws.com",sep="")
		cmdstr=paste(cmdstr_prefix,filepath,sep="")

	}
	return(cmdstr)
}

my_stop=function(my_path,msg){
	setwd(my_path)
	stop(msg)
}

get_download_cmdstr=function(filepath="/bagea-data/bagea_data_freeze/annotation_files/ensg2symbol20190410.txt",outpath,use_s3_client=FALSE){
	if(use_s3_client){
		cmdstr_prefix="aws s3 cp s3:/"
		cmdstr=paste(cmdstr_prefix,filepath," ",outpath,sep="")
	}else{
		cmdstr=paste("wget -P ",outpath," https://s3-us-west-1.amazonaws.com",filepath,sep="")
	}
	return(cmdstr)
}


get_pickrell_ldblocks_filepath=function(){
	BAGEA_PATH=Sys.getenv("BAGEA_PATH")
	commonsnp_filepath=paste(BAGEA_PATH,"/Data/pickrell_ldblocks.hg19.eur.sorted.bed",sep="")
	return(commonsnp_filepath)
}

get_commonsnp_filepath=function(){
	BAGEA_PATH=Sys.getenv("BAGEA_PATH")
	commonsnp_filepath=paste(BAGEA_PATH,"/Data/common_snps_sorted_hg19.bed",sep="")
	return(commonsnp_filepath)
}

get_tss_filepath=function(){
	BAGEA_PATH=Sys.getenv("BAGEA_PATH")
	tss_filepath=paste(BAGEA_PATH,"/Data/tss_median_hg19_sorted.bed",sep="")
	return(tss_filepath)
}

check_bedops=function(BAGEA_PATH){
	a1=system(paste("which ",BAGEA_PATH,"/bedops_dir/bin/bedops > /dev/null",sep=""))
	a2=system(paste("which ",BAGEA_PATH,"/bedops_dir/bin/bedmap >/dev/null",sep=""))
	a3=system(paste("which ",BAGEA_PATH,"/bedops_dir/bin/sort-bed >/dev/null",sep=""))
	ret=(a1==0 && a2==0 && a3==0)
	return(ret)
}

check_installation=function(){
	BAGEA_PATH=Sys.getenv("BAGEA_PATH")
	if(BAGEA_PATH==""){
		print("environment variable BAGEA_PATH seems not be set.")
		print("set before continuing..")
		stop()
	}
	aa=check_bedops(BAGEA_PATH)
	if(!aa){
		stop("bedops tool not installed:\\n run installation script. Aborting now.")
	}
	if(!file.exists(get_commonsnp_filepath())){
		paste("common snp file does not exist at ")
		paste(get_commonsnp_filepath())
		stop("Abborting now.")
	}
	if(file.size(get_commonsnp_filepath())==0){
		paste("common snp file empty at. ")
		paste(get_commonsnp_filepath())
		stop("Abborting now.")
	}
	if(!file.exists(get_tss_filepath())){
		paste("tss file does not exist at ")
		paste(get_tss_filepath())
		stop("Abborting now.")
	}
	if(file.size(get_tss_filepath())==0){
		paste("tss file empty at. ")
		paste(get_tss_filepath())
		stop("Abborting now.")
	}
}

#' installs external code and data 
#' 
#' Will read the envrionment variable BAGEA_PATH and download necessary code and data to this location.
#' @param proceed_savely:
#' Boolean. if TRUE installation procedure expects BAGEA_PATH location either to not exist yet or not contain either the subdirectory \code{Data} or \code{bedops_dir}.
#' @param download_processed:
#'' Boolean. if TRUE installation procedure pulls already processed files from dedicated s3 bucket rather than processing third party data 
#' @param delete_raw_downloads:
#' Boolean. if TRUE downloaded raw files will be deleted.
#' @export
install_external_data=function(proceed_savely=TRUE,download_processed=TRUE,delete_raw_downloads=FALSE){
	PREPARE_COMMON_SNPS_FILE=TRUE
	PREPARE_TSS_FILE=TRUE
	PREPARE_SNPS_IN_ANNOTATIONS=TRUE
	PREPARE_TSS_FILES=TRUE
	PREPARE_TSS_BED=TRUE
	PREPARE_DIRECTED_ANNOT=TRUE
	PREPARE_LDBLOCK_DATA=FALSE
	INSTALL_EXTERNAL_TOOLS=TRUE
#	INSTALL_EXTERNAL_TOOLS=FALSE
	PREPARE_FOLDERS=TRUE
	R_REPO=NULL#NULL if default repos is used
	BAGEA_PATH=Sys.getenv("BAGEA_PATH")
	if(BAGEA_PATH==""){
		print("environment variable BAGEA_PATH seems not be set.")
		print("set bevore continuing..")
		stop()
	}
	calling_path=getwd()
	if(!dir.exists(BAGEA_PATH)){
		print(paste("making directory",BAGEA_PATH))
		system(paste("mkdir -p ",BAGEA_PATH,sep=""))
	}
	print("changing directory.")
	setwd(BAGEA_PATH)

	if(sum(is.element(c("Data","bedops_dir"),list.files(".")))!=0 && proceed_savely==TRUE){
		print(paste("directory",BAGEA_PATH,"not empty and proceed_savely flag on. empty the directory (remove folder Data and bedops_dir) or set proceed_savely to FALSE"))
		my_stop(calling_path,"Aborting.")
	}

	################################################
	## install external Tools
	################################################
	if(INSTALL_EXTERNAL_TOOLS){
		mysys=Sys.info()[["sysname"]]
		myarch_is64=grepl("64",Sys.info()[["machine"]])
		checkbedops=check_bedops(BAGEA_PATH)
		if(!checkbedops){
		        if(mysys=="Linux" && myarch_is64){
		            system("wget https://github.com/bedops/bedops/releases/download/v2.4.26/bedops_linux_x86_64-v2.4.26.tar.bz2")
					system("mkdir bedops_dir")
		            system("tar jxvf bedops_linux_x86_64-v2.4.26.tar.bz2 -C bedops_dir")
		            system("rm bedops_linux_x86_64-v2.4.26.tar.bz2")
		        }
		        if(mysys=="Linux" && !myarch_is64){
		            system("wget https://github.com/bedops/bedops/releases/download/v2.4.26/bedops_linux_i386-v2.4.26.tar.bz2")
		            system("mkdir bedops_dir")
		            system("tar jxvf bedops_linux_i386-v2.4.26.tar.bz2 -C bedops_dir")
		            system("rm bedops_linux_i386-v2.4.26.tar.bz2")
		        }
		        if(mysys=="Darwin"){
		        	#### TODO
		        	a=system("git version>/dev/null")
		        	if(a!=0){
		        		my_stop(calling_path,"git not installed.")
		        	}
		        	system("git clone https://github.com/bedops/bedops.git")
		        	system("mv bedops bedops_dir")
		        	setwd("bedops_dir")
		        	system("make")
		        	system("make install")
		        	setwd("..")

		        }
			checkbedops=check_bedops(BAGEA_PATH)
			if(!checkbedops){
				my_stop(calling_path,"bedops tools not installed properly and automatic installation didn't work:\\n check http://bedops.readthedocs.io/en/latest/index.html for installation guide. Abborting now.")
			}
			checkperl=system("which perl")
		}
		if(checkperl!=0){
			my_stop(calling_path,"perl not installed: Abborting now.")
		}
		checkgunzip=system("which gunzip")
		if(checkgunzip!=0){
			my_stop(calling_path,"gunzip not installed: Abborting now.")
		}
		checkgzip=system("which gzip")
		if(checkgzip!=0){
			my_stop(calling_path,"gzip not installed: Abborting now.")
		}
		checkcut=system("which cut")
		if(checkcut!=0){
			my_stop(calling_path,"cut not installed: Abborting now.")
		}
		checkwget=system("which wget")
		if(checkwget!=0){
			my_stop(calling_path,"wget not installed: Abborting now.")
		}
	}
	################################################
	## END: install external Tools
	################################################
	################################################
	## prepare folders
	################################################
	if(PREPARE_FOLDERS){
		print("make data folders ..")
		files=list.dirs("Data")
		if(sum(files=="Data")==0){
			print(paste("make folder Data at" ,getwd()))
			a=system("mkdir Data")
			if(a!=0){
				my_stop(calling_path,"couldn't make Data folder")
			}
		}
	}
	################################################
	## END: prepare folders
	################################################

	################################################
	## check tool availability
	################################################
	if(PREPARE_COMMON_SNPS_FILE){
		print("prepare common snps file")
		if(download_processed){
			print("trying donwloading pocessed common-snps file from s3 ..")
			print("checking connection ..")
			checkping=system("which ping")
			mycmd=get_s3_download_cmdstr(filepath="/bagea-data/bagea_data_freeze/annotation_files/common_snps_sorted_hg19.bed",outpath="Data/common_snps_sorted_hg19.bed")
			system(mycmd)
		}else{
			print("trying donwloading common-snps file from ucsc..")
			print("checking connection ..")
			checkping=system("which ping")
			if(checkping!=0){
				stop("ping not installed: trying to download directly.")
			}else{
				ping_res=fread("echo 'ping running'; ping hgdownload.cse.ucsc.edu | head -1", sep="\t",header=F)[,V1]
				if(length(ping_res)==1){
					my_stop(calling_path,"hgdownload.cse.ucsc.edu not available, check your online connection")
				}
			}
			system("wget -P Data/ http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp147Common.txt.gz")
			print("format common SNPs file as bed file..")
			cmdstr=paste("gunzip -cf Data/snp147Common.txt.gz | perl -ne 'if($.>1){print}' | cut -f2,3,5,7,9 | perl -anlF\"\\t\" -e '$a=$F[1]+1;print \"$F[0]\\t$F[1]\\t$a\\t$F[2]\\t$F[3]\\t$F[4]\"'> Data/common_snps_hg19.bed")
			system(cmdstr)
			print("sort bed file..")
			cmdstr=paste(BAGEA_PATH,"/bedops_dir/bin/sort-bed Data/common_snps_hg19.bed > Data/common_snps_sorted_hg19.bed",sep="")
			system(cmdstr)
			system("rm Data/common_snps_hg19.bed")
			if(delete_raw_downloads){		
				system("rm Data/snp147Common.txt.gz")
			}
		}
	}

	if(PREPARE_DIRECTED_ANNOT){
		print("pull directed annotations ..")
		if(!download_processed){
			print("directed annotations only installed if flag 'download_processed' is TRUE")
		}else{
			print("trying donwloading directed annotations from s3 ..")
			print("checking connection ..")
			checkping=system("which ping")
			for(chr in c(1:22)){
				print(paste("get chr",chr))
				cur_download_filepath=paste("/bagea-data/bagea_data_freeze/directed_annotations/expecto_loading_normed_chr",chr,".RData",sep="")
				cur_output_filepath=paste("Data/expecto_loading_normed_chr",chr,".RData",sep="")
				mycmd=get_s3_download_cmdstr(filepath=cur_download_filepath,outpath=cur_output_filepath)
				system(mycmd)
			}
		}
	}

	if(PREPARE_TSS_FILE){
		print("prepare tss file ..")
		if(download_processed){
			print("trying donwloading processed gene annotation file from s3 ..")
			print("checking connection ..")
			checkping=system("which ping")
			mycmd=get_s3_download_cmdstr(filepath="/bagea-data/bagea_data_freeze/annotation_files/symbol2TSSmap_hg19.txt",outpath="Data/symbol2TSSmap_hg19.txt")
			system(mycmd)
			mycmd=get_s3_download_cmdstr(filepath="/bagea-data/bagea_data_freeze/annotation_files/tss_median_hg19_sorted.bed",outpath="Data/tss_median_hg19_sorted.bed")
			system(mycmd)
		}else{
			print("trying donwloading gene annotation file from ucsc..")
			print("checking connection ..")
			checkping=system("which ping")
			if(checkping!=0){
				print("ping not installed: trying to download directly.")
			}else{
				ping_res=fread("echo 'ping running'; ping hgdownload.cse.ucsc.edu | head -1", sep="\t",header=F)[,V1]
				if(length(ping_res)==1){
					my_stop(calling_path,"hgdownload.cse.ucsc.edu not available, check your online connection to server.")
				}
			}
			print("download TSS file..")
			system("wget -P Data/ http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz")
			table_names=c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
			gene_table=fread("gunzip -cf Data/refGene.txt.gz")
			setnames(gene_table,table_names)
			chromosomes=c(paste("chr",c(1:22),sep=""),"chrX")
			gene_table=gene_table[is.element(chrom, chromosomes),]
			gtab=gene_table[,list(maxchrom=max(chrom),minchrom=min(chrom),maxstrand=max(strand),minstrand=min(strand),mintxStart=min(txStart),maxtxEnd=max(txEnd),medtxStart=floor(median(txStart)),medtxEnd=floor(median(txEnd)), name2),by=name2]
			print("remove genes with transcripts on different chromosomes or different strands.")
			gtab=gtab[!(maxchrom!=minchrom | minstrand!=maxstrand),]
			fun=function(strand,st,end){st*(strand=="+") + end*(strand=="-")}
			gtab[,TSSMax:=fun(minstrand,mintxStart,maxtxEnd)]
			gtab[,TSSMedian:=fun(minstrand,medtxStart,medtxEnd)]
			fun=function(strand,st,end){st*(strand=="-") + end*(strand=="+")}
			gtab[,TSEMax:=fun(minstrand,mintxStart,maxtxEnd)]
			gtab[,TSEMedian:=fun(minstrand,medtxStart,medtxEnd)]
			gtab=gtab[,list(name2,minchrom,TSSMax,TSSMedian,TSEMax,TSEMedian,strand=minstrand)]
			setnames(gtab,c("name2","minchrom"),c("symbol","chrom"))
			print(paste("writing tss file"))
			write.table(gtab, file="Data/symbol2TSSmap_hg19.txt", sep="\t", quote=FALSE, col.names=TRUE,row.names=FALSE)
			gtab_bed=gtab[,list(chrom,TSSMedian,TSSMedian+1,symbol,strand)]
			write.table(gtab_bed, file="Data/tss_median_hg19.bed", sep="\t", quote=FALSE, col.names=FALSE,row.names=FALSE)
			a=system(paste(BAGEA_PATH,"/bedops_dir/bin/sort-bed Data/tss_median_hg19.bed > Data/tss_median_hg19_sorted.bed",sep=""))
			if(delete_raw_downloads){		
				system("rm Data/refGene.txt.gz")	
			}
			system("rm Data/tss_median_hg19.bed")
		}
		if(PREPARE_LDBLOCK_DATA){
			system("wget -P Data/ https://data.broadinstitute.org/alkesgroup/SLDP/refpanel/pickrell_ldblocks.hg19.eur.bed")
			system("sort-bed Data/pickrell_ldblocks.hg19.eur.bed > Data/pickrell_ldblocks.hg19.eur.sorted.bed")			
			mytbl=fread("Data/pickrell_ldblocks.hg19.eur.sorted.bed",header=FALSE)
			mytbl[,blocknr:=paste("block",c(1:dim(mytbl)[1]),sep="")]
			write.table(mytbl,file="Data/pickrell_ldblocks.hg19.eur.sorted.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
			system("rm Data/pickrell_ldblocks.hg19.eur.bed")
		}
		if(a!=0){
			my_stop(calling_path,"problem in sort-bed when sorting tss file")
		}

	}
	print("changing directory back.")
	setwd(calling_path)
}
########################################
