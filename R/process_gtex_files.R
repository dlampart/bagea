#' install GTEx summary statistics
#' 
#' @param GTEX_URLS
#' vector of urls to GTEx summary statistics
#'
#' @export
install_gtex=function(GTEX_URLS=c("https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/all_snp_gene_associations/Cells_Transformed_fibroblasts.allpairs.txt.gz"),NCORES=1){
	GTEX_INFO_URL="https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz"
	ENSG2SYMBOL_S3PATH="/bagea-data/bagea_data_freeze/annotation_files/ensg2symbol20190410.txt"
	CLEAN_PRIOR=TRUE
	#require(data.table)
	basename_gtex_infofile=basename(GTEX_INFO_URL)
	basename_gtex_ensg2symbolfile=basename(ENSG2SYMBOL_S3PATH)
	BAGEA_PATH=Sys.getenv("BAGEA_PATH")
	gtex_data_path=gsub("//","/",paste(BAGEA_PATH,"/Data/gtex_summary/",sep=""))
	if(!file.exists(gtex_data_path)){
		cmdstr=paste("mkdir ",gtex_data_path,sep="")
		system(cmdstr)
	}
	print("getting gtex info file")
	gtex_info_outpath=paste(gtex_data_path,basename_gtex_infofile,sep="")
	if(file.exists(gtex_info_outpath) && CLEAN_PRIOR){
		system(paste("rm ",gtex_info_outpath,sep=""))
	}
	if(!file.exists(gtex_info_outpath)){
		cmdstr=paste("wget -P ",gtex_data_path," ",GTEX_INFO_URL,sep="")
		system(cmdstr)
	}
	ensg2symbol_outpath=paste(gtex_data_path,basename_gtex_ensg2symbolfile,sep="")
	if(file.exists(ensg2symbol_outpath) && CLEAN_PRIOR){
		system(paste("rm ",ensg2symbol_outpath,sep=""))
	}
	if(!file.exists(ensg2symbol_outpath)){
		mycmdstr=get_download_cmdstr(filepath=ENSG2SYMBOL_S3PATH,outpath=gtex_data_path)
		system(mycmdstr)
	}
	ensg2symbol_tbl=fread(ensg2symbol_outpath)
	######################
	######################
	## 
	## 
	######################
	######################
	tss_tbl=fread(get_tss_filepath())
	##### subset ensg2symbol_tbl to the genes that are mappable.
	ensg2symbol_tbl=ensg2symbol_tbl[is.element(`HGNC symbol`,tss_tbl[,V4]),]
	ensg2symbol_tbl[,ensg:=`Gene stable ID`]
	setkey(ensg2symbol_tbl,ensg)
	commonsnp_tbl=fread(get_commonsnp_filepath())
	gtex_info_tbl=fread(cmd=paste("gunzip -cf ",gtex_info_outpath,sep=""))
	lapply(GTEX_URLS,function(gtex_url){
		basename_gtex_file=basename(gtex_url)
		print("trying to download the GTEX file..")
		print(gtex_url)
		gtex_outpath=paste(gtex_data_path,"/",basename_gtex_file,sep="")
		gtex_outpath_unzip=sub(".gz$","",gtex_outpath)
		gtex_outpath_stub=sub(".allpairs.txt.gz$","",gtex_outpath)
		if(file.exists(gtex_outpath) && CLEAN_PRIOR){
			system(paste("rm ",gtex_outpath,sep=""))
		}
		cmdstr=paste("wget -P ",gtex_data_path," ",gtex_url,sep="")
		system(cmdstr)
		system(paste("gunzip ",gtex_outpath,sep=""))
		gtex_out_tbl=fread(gtex_outpath_unzip)
		setkey(gtex_out_tbl,variant_id)		
		mclapply(c(22:1),function(mychr){
			#### subset to chr
			chrstr=paste("chr",mychr,sep="")
			print(chrstr)
			commonsnp_tbl_chr=commonsnp_tbl[V1==chrstr,]
			commonsnp_tbl_chr[,snps:=V4]
			commonsnp_tbl_chr[,V4:=NULL]
			commonsnp_tbl_chr[,refUCSC:=V6]
			commonsnp_tbl_chr[,V6:=NULL]

			gtex_info_tbl_chr=gtex_info_tbl[chr==mychr,]
			#### merge gtex data with its info (to map snpids)
			gtex_info_tbl_chr=gtex_info_tbl_chr["."!=rs_id_dbSNP147_GRCh37p13,]
			setkey(gtex_info_tbl_chr,variant_id)
			gtex_tbl_chr=merge(gtex_out_tbl,gtex_info_tbl_chr)
			#### subset to  chr400KB
			gtex_tbl_chr=gtex_tbl_chr[abs(tss_distance)<400000,]
			#### merge symbol table
			gtex_tbl_chr[,ensg:=sub("\\.\\d+","",gene_id)]
			gtex_tbl_chr=gtex_tbl_chr[is.element(ensg,ensg2symbol_tbl[,ensg]),]
			setkey(gtex_tbl_chr,ensg)
			setkey(ensg2symbol_tbl,ensg)
			full_tbl=merge(gtex_tbl_chr,ensg2symbol_tbl)
			full_tbl=full_tbl[,list(chr=chr,ref=ref,snps=rs_id_dbSNP147_GRCh37p13,symbol=`Gene name`,zscore=slope/slope_se)]
			#### merge to commonsnps table to be able to check directionality
			setkey(full_tbl,snps)
			setkey(commonsnp_tbl_chr,snps)
			full_tbl=merge(full_tbl,commonsnp_tbl_chr)
			full_tbl=full_tbl[ref==refUCSC,]
			#### write table out
			out_tbl=full_tbl[,list(gene=symbol,snps,statistic=zscore)]
			outfile_path=paste(gtex_outpath_stub,"_",chrstr,".txt",sep="")
			write.table(out_tbl,file=outfile_path,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
		},mc.cores=NCORES)
		if(file.exists(gtex_outpath_unzip) && CLEAN_PRIOR){
			system(paste("rm ",gtex_outpath_unzip,sep=""))
		}
	})
}