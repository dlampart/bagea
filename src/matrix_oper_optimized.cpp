#include <Rcpp.h>
using namespace Rcpp;


/* calculates diag(t(M)%*%C%*%M), where M is sparse.
 It's optimized for when M has a large number of columns. M is supplied as sparseM,
 where the first two columns are indexes which can be derived in R by which(M,2).
 The third column is the values in M. lenM. is the  number of columns of M.

*/
// [[Rcpp::export]]
NumericVector getDiagMtCM_optim(NumericMatrix Cmat, NumericMatrix sparseM, int lenM) {
  int toprow =sparseM.nrow();
  NumericVector out(lenM);
  int snp_ind = 0;
  double modif = 0.0; 	// multiplier that multiplies off-diagonal elements by 2.
  int lagging_i = 0;		// lagging annot index in sparse melted format
  int lagging_snp_ind = sparseM(0,1)-1;		// lagging snp index in sparse melted format
  for(int i = 1; i<(toprow+1);++i){		// iterate over sparse matrix. Add one to write last snp 
  	if(i<toprow){
  		snp_ind = sparseM(i,1)-1;		// nonlagging snp index in sparse melted format
  	}
	if(snp_ind!=lagging_snp_ind || i==(toprow)){		// new snp starting
		int stretch_len = i - lagging_i;		// number of active annotations for current snp
		for(int k = 0; k < stretch_len; ++k){		//iterate over active annotations
			for(int l = k; l < stretch_len; ++l){		//iterate over active annotations
				modif = 2.0;
				if(k==l){
					modif = 1.0;
				}
 	 			int sparse_ind_k = lagging_i + k;		// l index in melted sparse
 	 			int sparse_ind_l = lagging_i + l;		// k index in melted sparse
 	 			int k_v = (int) sparseM(sparse_ind_k,0)-1;		// annot index k
 	 			int l_v = (int) sparseM(sparse_ind_l,0)-1;		// annot index l
				out(lagging_snp_ind) = out(lagging_snp_ind) + Cmat(k_v,l_v) * sparseM(sparse_ind_l,2) * sparseM(sparse_ind_k,2) * modif;
////				out(lagging_snp_ind)= out(lagging_snp_ind) + Cmat(k_v,l_v);
////				out(lagging_snp_ind)=out(lagging_snp_ind) + sparseM(sparse_ind_l,2) * sparseM(sparse_ind_l,2);
			}
		}
		lagging_snp_ind = snp_ind; // set new lagging snps.
		lagging_i = i;  
	}
   }
  return out;
}

// [[Rcpp::export]]
NumericVector getDiagMtCM(NumericMatrix Cmat, NumericMatrix sparseM, int lenM) {
  int toprow =sparseM.nrow();
  NumericVector out(lenM);
  int snp_ind = 0;
  int lagging_i = 0;		// lagging annot index in sparse melted format
  int lagging_snp_ind = 0;		// lagging snp index in sparse melted format
  for(int i = 1; i<(toprow+1);++i){		// iterate over sparse matrix. Add one to write last snp 
  	if(i<toprow){
  		snp_ind=sparseM(i,1)-1;		// nonlagging snp index in sparse melted format
  	}
	if(snp_ind!=lagging_snp_ind || i==(toprow)){		// new snp starting
		int stretch_len = i - lagging_i;		// number of active annotations for current snp
		for(int k = 0; k < stretch_len; ++k){		//iterate over active annotations
			for(int l = 0; l < stretch_len; ++l){		//iterate over active annotations
 	 			int sparse_ind_k = lagging_i + k;		// l index in melted sparse
 	 			int sparse_ind_l = lagging_i + l;		// k index in melted sparse
 	 			int k_v = (int) sparseM(sparse_ind_k,0)-1;		// annot index k
 	 			int l_v = (int) sparseM(sparse_ind_l,0)-1;		// annot index l
				out(lagging_snp_ind)=out(lagging_snp_ind) + Cmat(k_v,l_v) * sparseM(sparse_ind_l,2) * sparseM(sparse_ind_k,2);
////				out(lagging_snp_ind)= out(lagging_snp_ind) + Cmat(k_v,l_v);
////				out(lagging_snp_ind)=out(lagging_snp_ind) + sparseM(sparse_ind_l,2) * sparseM(sparse_ind_l,2);
			}
		}
		lagging_snp_ind = snp_ind; // set new lagging snps.
		lagging_i = i;  
	}
   }
  return out;
}