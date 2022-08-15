#' Hybridization
#' @description Hybridization method (algorithm) to be used in the
#' genetic algorithm.
#' @param population numeric matrix which rows are chromosomes i.e. vectors of 
#' parameters values.
#' @param fitness numeric vector which \code{i}-th element is the value of 
#' \code{fn} at point \code{population[i, ]}.
#' @param hybrid.n positive integer representing the number of hybrids.
#' @param method hybridization method to improve chromosomes via local search.
#' @param par additional parameters to be passed depending on the \code{method}.
#' @param opt.par parameters of the local optimization function
#' to be used for hybridization algorithm (including \code{fn} and \code{gr}).
#' @param info logical; if \code{TRUE} then some optimization related 
#' information will be printed each iteration.
#' @param iter iteration number of the genetic algorithm.
#' @param ... additional parameters to be passed to 
#' \code{fn} and \code{gr} functions.
#' @details This function uses \code{\link[gena]{gena.mating}} function to 
#' select hybrids. Therefore \code{method} and \code{par} arguments will
#' be passed to this function. If some chromosomes selected to become hybrids
#' are duplicated then these duplicates will not be subject to local
#' optimization i.e. the number of hybrids will be decreased by the number
#' of duplicates (actual number of hybrids during some iterations may be 
#' lower than \code{hybrid.n}).
#' 
#' Currently \code{\link[stats]{optim}} is the only available local
#' optimizer. Therefore \code{opt.par} is a list containing parameters
#' that should be passed to \code{\link[stats]{optim}}.
#' 
#' For more information on hybridization
#' please see El-mihoub et. al. (2006).
#' 
#' @return The function returns a list with the following elements:
#' \itemize{
#' \item \code{population} - matrix which rows are chromosomes including hybrids.
#' \item \code{fitness} - vector which i-th element is the fitness of the
#' i-th chromosome.
#' \item \code{hybrids.ind} - vector of indexes of chromosomes selected for
#' hybridization.
#' \item \code{counts} a two-element integer vector giving the number of
#' calls to \code{fn} and \code{gr} respectively.
#' }
#' 
#' @references T. El-mihoub, A. Hopgood, L. Nolle, B. Alan (2006). 
#' Hybrid Genetic Algorithms: A Review.
#' \emph{Engineering Letters}, 13 (3), 124-137.
#' 
#' @examples 
#' # Consider the following fitness function
#' fn <- function(x)
#' {
#'   val <- x[1] * x[2] - x[1] ^ 2 - x[2] ^ 2
#' }
#' 
#' # Also let's provide it's gradient (optional)
#' gr <- function(x)
#' {
#'   val <- c(x[2] - 2 * x[1],
#'            x[1] - 2 * x[2])
#' }
#' 
#' # Randomly initialize the population
#' set.seed(123)
#' n_population <- 10
#' population <- gena.population(pop.n = n_population,
#'                               lower = c(-5, -5), 
#'                               upper = c(5, 5))
#'
#' # Calculate fitness of each chromosome
#' fitness <- rep(NA, n_population)
#' for(i in 1:n_population)
#' {
#'   fitness[i] <- fn(population[i, ])
#' }
#' 
#' # Perform hybridization
#' hybrids <- gena.hybrid(population = population,
#'                        fitness = fitness,
#'                        opt.par = list(fn = fn,
#'                                       gr = gr,
#'                                       method = "BFGS",
#'                                       control = list(fnscale = -1,
#'                                                      abstol = 1e-10,
#'                                                      reltol = 1e-10,
#'                                                      maxit = 1000)),
#'                        hybrid.n = 2,
#'                        method = "rank",
#'                        par = 0.8)
#' print(hybrids)
#' 

# Make hybrids
gena.hybrid <- function(population,
                        fitness,
                        hybrid.n = 1,
                        method,
                        par,
                        opt.par,
                        info = FALSE,
                        iter = NULL,
                        ...)          
{      
  # ---
  # The algorithm brief description:
  # 1. Select the chromosomes to become the hybrids
  # 2. Perform local optimization on chromosomes
  #    to make them the hybrids
  # ---
  
  # Prepare some values
  
  n_pop <- nrow(population)                                # number of chromosomes
  genes.n <- ncol(population)                              # number of genes
  
  dots <- list(...)                                        # fn specific parameters
  all.args <- append(opt.par, dots)                        # all parameters to be
                                                           # passed to optim
  hybrids <- matrix(NA, nrow = n_pop, ncol = genes.n)      # matrix to store hybrids
  hybrids_fitness <- rep(NA, n_pop)                        # vector to store fitnesses 
                                                           # of hybrids
  
  # Get the name of the first argument of the function
  
  fn.par <- names(formals(all.args$fn))[1]
  
  # To select hybrids use the same function
  # as for the parents
  
  hybrids_list <- gena.mating(population = population,     # select the chromosomes
                              fitness = fitness,           # to become a hybrids
                              parents.n = hybrid.n,
                              method = method,
                              par = par,
                              iter = iter,
                              self = TRUE)
  hybrids.ind <- unique(hybrids_list$ind)                  # indexes of chromosomes that
                                                           # will become a hybrids
  # Perform the optimization
  
  counter <- 1
  counts <- c(0, 0)
  for (i in hybrids.ind)
  {
    hybrid_message <- paste0("hybrid ", counter," fitness: before = ",
                             fitness[i],
                             ", after = ")
    
    all.args$par <- population[i, ]
    opt.result <- NULL
    tryCatch(
    {
      opt.result <- do.call(optim, args = all.args)
      if (!is.null(opt.result[["counts"]]))
      {
        optim.counts <- opt.result$counts
        optim.counts[is.na(optim.counts)] <- 0
        counts <- counts + optim.counts
      }
      if (fitness[i] < opt.result$value)
      {
        population[i, ] <- opt.result$par
        fitness[i] <- opt.result$value
      }
    })
    hybrid_message <- paste0(hybrid_message, fitness[i], "\n")
    if (info)
    {
      cat(hybrid_message)
    }
    counter <- counter + 1
  }
  
  return_list <- list(population = population,
                      fitness = fitness,
                      hybrids.ind = hybrids.ind,
                      counts = counts)
  
  return(return_list)
  
  # ---
  # Ideas for the future development
  # ---
  # 1. Add argument type - string representing a type of the selection mechanism.
  #   If \code{hybrid.type = "Lamarck"} (default) then hybridization changes
  #   both genes and fitness of the corresponding chromosome.
  #   If \code{hybrid.type = "Baldwin"} then only the fitness changes.
  #   If this argument is double between \code{0} and \code{1} then this 
  #   number represents probability of the Lamarck search.
  # ---
}

# References
# Best BHGA: https://www.hindawi.com/journals/jam/2013/103591/