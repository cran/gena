#include <RcppArmadillo.h>
#include "mating.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export(rng = false)]]
NumericVector mating_tournament(const NumericVector pop_ind,
                                const int candidates_n,
                                const NumericVector fitness,
                                const int parents_n)
{
  // The number of chromosomes in population
  const int pop_n = pop_ind.size();
  const NumericVector pop_ind_adj = pop_ind - 1; 
  
  // Indexes of chromosomes selected to be parents
  NumericVector parents_ind(pop_n);
  
  // Begin the tournamenet
  for (int i = 0; i < parents_n; i++)
  {
    // select candidates
    NumericVector candidates = sample(pop_ind_adj, candidates_n, false);
    // get their fitness
    NumericVector candidates_fitness = fitness[candidates];
    // select the fittest candidate
    NumericVector best_ind = which_max(candidates_fitness);
    // index of the best candidate in the whole population
    NumericVector best_ind_pop = candidates[best_ind];
    // assign the index of the best candidate to be
    // the index of the i-th parent
    parents_ind[i] = best_ind_pop[0] + 1;
  }
  
  // Return the result
  return(parents_ind);
}
