#ifndef gena_mating_H
#define gena_mating_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

NumericVector mating_tournament(const NumericVector pop_ind,
                                const int candidates_n,
                                const NumericVector fitness,
                                const int parents_n);

#endif
