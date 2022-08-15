#' Constraints
#' @description Impose constraints on chromosomes.
#' @param population numeric matrix which rows are chromosomes i.e. vectors of 
#' parameters values.
#' @param method method used to impose constraints.
#' @param par additional parameters to be passed depending on the \code{method}.
#' @param iter iteration number of the genetic algorithm.
#' @details If \code{method = "bounds"} then chromosomes will be bounded
#' between \code{par$lower} and \code{par$upper}.
#' @return The function returns a list with the following elements:
#' \itemize{
#' \item \code{population} - matrix which rows are chromosomes after
#' constraints have been imposed.
#' \item \code{constr.ind} - matrix of logical values which (i, j)-th
#' elements equals \code{TRUE} (\code{FALSE} otherwise) if j-th jene of
#' i-th chromosome is a subject to constraint.
#' }
#' @examples 
#' # Randomly initialize population
#' set.seed(123)
#' population <- gena.population(pop.n = 10,
#'                               lower = c(-5, -5), 
#'                               upper = c(5, 5))
#'                            
#' # Impose lower and upper bounds constraints
#' pop.constr <- gena.constr(population, 
#'                           method = "bounds",
#'                           par = list(lower = c(-1, 2),
#'                                      upper = c(1, 5)))
#' print(pop.constr)
#' 
gena.constr <- function(population,
                        method = "bounds",
                        par,
                        iter)              
{
  # Get the number of genes
  
  genes.n <- ncol(population)
  pop.n <- nrow(population)
  
  # Transform lower and upper bounds into matrices
  
  lower.mat <- matrix(par$lower, 
                      nrow = pop.n, ncol = genes.n, byrow=TRUE)
  upper.mat <- matrix(par$upper, 
                      nrow = pop.n, ncol = genes.n, byrow=TRUE)
  
  # Prepare the output list
  
  return.list <- list(population = NA,
                      constr.ind = NA)

  if (method == "bounds")
  {
      constr.ind.lower <- population < lower.mat
      constr.ind.upper <- population > upper.mat
      constr.ind <- constr.ind.lower | constr.ind.upper

      population[constr.ind.lower] <- lower.mat[constr.ind.lower]
      population[constr.ind.upper] <- upper.mat[constr.ind.upper]
      
      return.list$population <- population
      return.list$constr.ind <- constr.ind
  }

  return(return.list)
}

# Assign default parameters for
# constraint algorithm depending
# on the "method"
gena.constr.validate <- function(method, par)
{
  # Validate the "method"
  
  methods <- c("bounds")                                   # the list of all
                                                           # available methods
  
  if (!(method %in% methods))
  {
    stop(paste0("Incorrect constr.method argument. ",      # if the user has provided
                "It should be one of: ",                   # incorrect argument
                paste(methods, collapse = ", "),
                "\n"))
  }
}
