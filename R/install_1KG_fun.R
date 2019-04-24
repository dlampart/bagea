#' Downloads and processes 1kg data. 
#'
#' Downloads 1KG data from 1kg ftp site, processes, and installs it in subfolder given by the environmental variable \code{BAGEA_PATH}.
#' @param only_chr22
#' Boolean. Should only chromosome 22 be processed?
#' @param super_pop_selected
#' String. Which superpopulation should be processed.(currently only EUR is allowed.)
#' @param genotype_chr22_filepath
#' String. Gives filepath to genotypes for chr22 on 1KG ftp site.
#' @param indiv_filepath
#' String. Gives filepath to population annotation file on 1KG ftp site.
#' @param relationship_filepath
#' String. Gives filepath to family relationship file on 1KG ftp site.
#' @param maf_cutoff
#' Scalar. gives maf cutoff to be used. has to be at least 0.02.
#' @param keep_ambiguous
#' Boolean. If false, variants where the 1KG REF does not agree with  UCSC REF, are removed, else they are kept.
#' @param proceed_savely
#' Boolean. If true, downloading and processing will only proceed, if there is no files of the relevant super population in the folder ${BAGEA_PATH}/Data/1KG
#' @param download_processed
#' Boolean. If true, preprocessed 1KG files will be downloaded from dedicated aws s3 site. (only available for super population set to EUR and maf_cutoff to 0.02)
#' @export
install_1KG_data=function(
	only_chr22=FALSE,
	super_pop_selected="EUR",
	genotype_chr22_filepath="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
	indiv_filepath="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
	relationship_filepath="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/20140625_related_individuals.txt",
	maf_cutoff=0.02,
	keep_ambiguous=FALSE,
	proceed_savely=TRUE,
	download_processed=TRUE
	){
	if(maf_cutoff<0.02){
		stop("maf_cutoff has to be at least 0.02")
	}
	if(super_pop_selected!="EUR"){
		stop("currently only EUR is the only allowed value for super_pop_selected..")
	}
	mytmpdir=tempdir()
	BAGEA_PATH=Sys.getenv("BAGEA_PATH")
	if(BAGEA_PATH==""){
		print("environment variable BAGEA_PATH seems not be set.")
		print("set bevore continuing..")
		stop()
	}
	if(proceed_savely){
		save_path=paste(BAGEA_PATH,"/Data/1KG/",sep="")
		files_present=list.files(save_path)
		outf=paste("1KG_",tolower(super_pop_selected),("_ambiguous_kept")[keep_ambiguous],"_chr",c(1:22),".txt.gz",sep="")
		superpop_pres=files_present[is.element(files_present,outf)]
		if(length(superpop_pres)>0){
			stop(paste("proceed_savely is activated, but some super_pop files already present, for instance: ",save_path,superpop_pres[1],sep=""))
		}
	}
	calling_path=getwd()
	print("trying to install 1KG data in path given for environment variable BAGEA_PATH which is:")
	print(BAGEA_PATH)
	bagea_data_path=paste(BAGEA_PATH,"/Data",sep="")
	if(!dir.exists(bagea_data_path)){
		stop("${BAGEA_PATH}/Data folder doesn't exist. Make sure you run bagea::install_external_data() first.")
	}
	print("changing directory.")
	setwd(bagea_data_path)	
	if(!file.exists(get_commonsnp_filepath())){
		stop("some necessary files missing. Looks like bagea::install_external_data() was not run successfully for BAGEA_PATH. Make sure bagea::install_external_data() ran successfully first.")
	}
	commonsnp_tbl=fread(get_commonsnp_filepath())
	checkperl=system("which perl")
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
	if(!dir.exists("1KG")){
		system("mkdir 1KG")
	}
	setwd("1KG")
	if(download_processed){
		print("downloading already processed files..")
		if(maf_cutoff!=0.02){
			setwd(calling_path)
			stop("if download_processed==TRUE, then maf_cutoff has to be 0.02")
		}
		if(super_pop_selected!="EUR"){
			setwd(calling_path)
			stop("if download_processed==TRUE, then super_pop_selected has to be EUR")
		}
		chrstorun=c(1:22)
		if(only_chr22){
			chrstorun=22
		}
		lapply(chrstorun,function(x){
			cmdstr=paste("wget https://s3-us-west-1.amazonaws.com/bagea-data/bagea_data_freeze/1KG_EUR/1KG_eur",("_ambiguous_kept")[keep_ambiguous],"_chr",x,".txt.gz",sep="")
			system(cmdstr)
		})
		setwd(calling_path)
		return(TRUE)
	}
	myfilename=sub(".+/","",indiv_filepath)
	if(grepl("^/",myfilename)){
		stop("no double dash at end of filepath allowed.")
	}
	indiv_tmpfile=tempfile(pattern=myfilename,tmpdir=mytmpdir)
	get_indivfile_cmdstr=paste("wget -O ",indiv_tmpfile," ",indiv_filepath,sep="")
	system(get_indivfile_cmdstr)
	indiv_tbl=fread(indiv_tmpfile)

	selected_indivs=indiv_tbl[super_pop==super_pop_selected,sample]
	####################
	myfilename=sub(".+/","",relationship_filepath)
	if(grepl("^/",myfilename)){
		stop("no double dash at end of filepath allowed.")
	}
	relationship_tmpfile=tempfile(pattern=myfilename,tmpdir=mytmpdir)
	get_relationship_cmdstr=paste("wget -O ",relationship_tmpfile," ",relationship_filepath,sep="")
	system(get_relationship_cmdstr)
	relationship_tbl=fread(relationship_tmpfile)
	selected_indivs=setdiff(selected_indivs,relationship_tbl)
	####################
	####################
	chrstorun=c(1:22)
	if(only_chr22){
		chrstorun=22
	}
	for(chr in chrstorun){
		chrstr=paste("chr",chr,sep="")
		print(paste("downloading ",chrstr,"...",sep=""))
		cur_genotype_filepath=gsub("chr22",chrstr,genotype_chr22_filepath)
		myfilename=sub(".+/","",cur_genotype_filepath)
		if(grepl("^/",myfilename)){
			stop("no double dash at end of filepath allowed.")
		}
		geno_tmpfile=tempfile(pattern=myfilename,tmpdir=mytmpdir)
		download_cmdstr=paste("wget -O ",geno_tmpfile," ",cur_genotype_filepath,sep="")
		system(download_cmdstr)
	
		allgeno_tmpfile=tempfile()
		print("remove EUR low frequency SNPs..")
		filter_lowfreq_cmd=paste("gunzip -cf  ",geno_tmpfile," | perl -ne '/(AF=|HG00150)/ && print' | perl -ne '!/EUR_AF=0;/ && print' | perl -ne '!/EUR_AF=0.0(0|1|2)/ && print' |  perl -ne '!/EUR_AF=1/ && print' |  perl -ne '!/EUR_AF=0.9(9|8)/ && print' | perl -ple 's/(0\\|1|1\\|0)/1/mg' | perl -ple 's/0\\|0/0/mg'  | perl -ple 's/1\\|1/2/mg' |  perl -ne '!/\\d\\|\\d/ && print'",sep="")
		#system(filter_lowfreq_cmd)
		all_dat=fread(cmd=filter_lowfreq_cmd,header=TRUE)
		unlink(geno_tmpfile)
#		all_dat=fread(allgeno_tmpfile,header=TRUE)
#		unlink(allgeno_tmpfile)
		euro_dat=all_dat[,which(is.element(names(all_dat),selected_indivs)),with=FALSE]
		snpids=all_dat[,ID]
		######
		######
		if(!keep_ambiguous){
			chrstr=paste("chr",chr,sep="")
			commonsnp_tbl_chr=commonsnp_tbl[V1==chrstr,]
			commonsnp_tbl_chr[,ID:=V4]
			setkey(commonsnp_tbl_chr,ID)
			allele_tbl=all_dat[,list(ID,REF,ALT)]
			setkey(allele_tbl,ID)
			allele_merge_tbl=merge(allele_tbl,commonsnp_tbl_chr)
			flipped_snps=allele_merge_tbl[REF!=V6 & nchar(REF)==nchar(ALT),ID]
			flipped_dels=allele_merge_tbl[REF!=V6 & nchar(REF)>nchar(ALT) & V6=="-",ID]
			flipped_ind=allele_merge_tbl[REF!=V6 & nchar(REF)<nchar(ALT) & V6!="-",ID]
			flipped_indels=c(flipped_ind,flipped_dels)
			print(paste("number of ambigous indels found:", length(flipped_indels)))
			print(paste("number of ambigous snps found:", length(flipped_snps)))
			flipped_ids=c(flipped_snps,c(flipped_ind,flipped_dels))
			not_flipped=!is.element(snpids,flipped_ids)
			euro_dat=euro_dat[not_flipped,]
			snpids=snpids[not_flipped]
		}
		######
		######
		nonan=rowSums(is.na(euro_dat))==0
		euro_dat=euro_dat[nonan,]
		snpids=snpids[nonan]
		mymeans=rowMeans(euro_dat)
		maf_cutoff_double=maf_cutoff*2
		not_too_low=mymeans>maf_cutoff_double & mymeans<(2-maf_cutoff_double)
		mean(not_too_low)
		euro_dat=euro_dat[not_too_low,]
		snpids=snpids[not_too_low]
		euro_dat=data.table(snpid=snpids,euro_dat)
		names(euro_dat)=c("snpid",names(euro_dat)[c(2:dim(euro_dat)[2])])
		outpath=paste("1KG_eur_",("ambiguous_kept_")[keep_ambiguous],chrstr,".txt",sep="")
		write.table(euro_dat,file=outpath,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
		system(paste("gzip ", outpath,sep=""))
	}
	setwd(calling_path)
	return(TRUE)
}
	
