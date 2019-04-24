#########################################################
#########################################################
# install relevant R packages.
#########################################################
#########################################################
install.packages("devtools")
install.packages("Rcpp")
install.packages("roxygen2")
install.packages("data.table")
install.packages("pryr")
#########################################################
#########################################################
# END: install relevant R packages.
#########################################################
#########################################################
#########################################################
#########################################################
# compile R package
##
## make sure that your working directory is the bagea packag path
#########################################################
#########################################################
devtools::use_rcpp()
dummy=tryCatch({devtools::document()},
	error=function(e){return(FALSE)}
)
Rcpp::compileAttributes()
devtools::load_all()
devtools::document()
devtools::install()
#########################################################
#########################################################
# END: compile R package
#########################################################
#########################################################
#########################################################
#########################################################
# download 
#########################################################
#########################################################
## environment variable BAGEA_PATH has to be set and point 
## to an empty folder
## for instance create an empty folder:
##> mkdir  "/home/BAGEA_DATA"
## append the following line to you ~/.bashrc file:
##> export BAGEA_PATH="/home/BAGEA_DATA"
require(bagea)
##install_external_data(proceed_savely=TRUE,download_processed=FALSE)
install_external_data(proceed_savely=TRUE)
install_1KG_data()
#########################################################
