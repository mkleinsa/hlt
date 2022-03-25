#ifndef standardize_lambda_H
#define standardize_lambda_H

void standardize_lambda(Rcpp::NumericMatrix & post, int start_idx, int end_idx,
                        int nT, int n, Rcpp::NumericMatrix & corr_theta);

#endif