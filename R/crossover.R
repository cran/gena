#' Crossover
#' @description Crossover method (algorithm) to be used in the
#' genetic algorithm.
#' @param parents numeric matrix which rows are parents i.e. vectors of 
#' parameters values.
#' @param fitness numeric vector which \code{i}-th element is the value of 
#' \code{fn} at point \code{population[i, ]}.
#' @param prob probability of crossover.
#' @param method crossover method to be used for making children.
#' @param par additional parameters to be passed depending on the \code{method}.
#' @param iter iteration number of the genetic algorithm.
#' @details Denote \code{parents} by \eqn{C^{parent}} which \code{i}-th row 
#' \code{parents[i, ]} is a chromosome \eqn{c_{i}^{parent}} i.e. the vector of 
#' parameter values of the function being optimized \eqn{f(.)} that is
#' provided via \code{fn} argument of \code{\link[gena]{gena}}.
#' The elements of chromosome \eqn{c_{ij}^{parent}} are genes 
#' representing parameters values.
#' 
#' Crossover algorithm determines the way parents produce children. 
#' During crossover each of randomly selected pairs of parents 
#' \eqn{c_{i}^{parent}}, \eqn{c_{i + 1}^{parent}}
#' produce two children 
#' \eqn{c_{i}^{child}}, \eqn{c_{i + 1}^{child}}, 
#' where \eqn{i} is odd. Each pair of parents is selected with 
#' probability \code{prob}. If pair of parents have not been selected
#' for crossover then corresponding children and parents are coincide i.e. 
#' \eqn{c_{i}^{child}=c_{i}^{parent}} and
#' \eqn{c_{i+1}^{child}=c_{i+1}^{parent}}.
#' 
#' Argument \code{method} determines particular crossover algorithm to 
#' be applied. Denote by \eqn{\tau} the vector of parameters used by the 
#' algorithm. Note that \eqn{\tau} corresponds to \code{par}.
#' 
#' If \code{method = "split"} then each gene of the first child will
#' be equiprobably picked from the first or from the second parent. So 
#' \eqn{c_{ij}^{child}} may be equal to \eqn{c_{ij}^{parent}}
#' or \eqn{c_{(i+1)j}^{parent}} with equal probability. The second
#' child is the reversal of the first one in a sense that if the first child 
#' gets particular gene of the first (second) parent then the second child gets
#' this gene from the first (second) parent i.e. if
#' \eqn{c_{ij}^{child}=c_{ij}^{parent}} then 
#' \eqn{c_{(i+1)j}^{child}=c_{(i+1)j}^{parent}}; if 
#' \eqn{c_{ij}^{child}=c_{(i+1)j}^{parent}} then 
#' \eqn{c_{(i+1)j}^{child}=c_{ij}^{parent}}.
#' 
#' If \code{method = "arithmetic"} then:
#' \deqn{c_{i}^{child}=\tau_{1}c_{i}^{parent}+
#' \left(1-\tau_{1}\right)c_{i+1}^{parent}}
#' \deqn{c_{i+1}^{child}=\left(1-\tau_{1}\right)c_{i}^{parent}+
#' \tau_{1}c_{i+1}^{parent}}
#' where \eqn{\tau_{1}} is \code{par[1]}. By default \code{par[1] = 0.5}.
#' 
#' If \code{method = "local"} then the procedure is the same as 
#' for "arithmetic" method but \eqn{\tau_{1}} is a uniform random
#' value between 0 and 1.
#' 
#' If \code{method = "flat"} then \eqn{c_{ij}^{child}} is a uniform
#' random number between \eqn{c_{ij}^{parent}} and 
#' \eqn{c_{(i+1)j}^{parent}}. 
#' Similarly for the second child \eqn{c_{(i+1)j}^{child}}.
#' 
#' For more information on crossover algorithms
#' please see Kora, Yadlapalli (2017).
#' 
#' @return The function returns a matrix which rows are children.
#' 
#' @references P. Kora, P. Yadlapalli. (2017). 
#' Crossover Operators in Genetic Algorithms: A Review.
#' \emph{International Journal of Computer Applications}, 162 (10), 34-36,
#' <doi:10.5120/ijca2017913370>.
#' @examples 
#' # Randomly initialize the parents
#' set.seed(123)
#' parents.n <- 10
#' parents <- gena.population(pop.n = parents.n,
#'                            lower = c(-5, -5), 
#'                            upper = c(5, 5))
#'                            
#' # Perform the crossover
#' children <- gena.crossover(parents = parents,
#'                            prob = 0.6,
#'                            method = "local")
#' print(children)

# Crossover the parents
gena.crossover <- function(parents,                        # the parents
                           fitness = NULL,                 # fitness values of chromosomes
                           prob = 0.8,                     # crossover probability
                           method = "local",               # crossover type
                           par = NULL,                     # crossover parameters
                           iter = NULL)                    # iteration of the
                                                           # genetic algorithm)                            # additional parameters                
{
  # ---
  # The algorithm brief description:
  # 1. Randomly pick the parents for the crossover
  # 2. Perform the crossover for each of randomly 
  #    selected pair of parents
  # ---
  
  # Prepare some variables
  
  parents.n <- nrow(parents)                               # number of parents
  genes.n <- ncol(parents)                                 # number of genes
  
  children <- parents                                      # matrix to store children
                                                           # initially contains
                                                           # the parents

  # Perform the crossover for each chromosome
  
  random_values <- runif(parents.n / 2, 0, 1)              # values to control for
                                                           # weather crossover should
                                                           # take place for a given
                                                           # pair of chromosomes
  crossover_ind <- which(random_values <= prob) * 2        # select the pairs
                                                           # for a crossover
  
  for (i in crossover_ind)
  {
    ind <- c(i - 1, i)                                     # indexes of the fist 
                                                           # and the second parents 
    
    if (method == "split")
    {
      split_ind <- rbinom(genes.n, 1, 0.5)                 # indexes to get from
                                                           # the other parent
      children[ind[1], split_ind] <- parents[ind[2],       
                                             split_ind] 
      children[ind[2], split_ind] <- parents[ind[1], 
                                             split_ind]
    }
    
    if (method == "arithmetic")
    {
      children[ind[1], ] <- parents[ind[1], ] * par[1] +      # mix the parents with
        parents[ind[2], ] * (1 - par[1])                      # a given weight
      children[ind[2], ] <- parents[ind[2], ] * par[1] +
        parents[ind[1], ] * (1 - par[1])
    }
    
    if (method == "local")                                    # local crossover
    {
      weight <- runif(n = genes.n, min = 0, max = 1)          # random weight
      children[ind[1], ] <- parents[ind[1], ] * weight +      # mix the parents with
                            parents[ind[2], ] * (1 - weight)  # a random weight
      children[ind[2], ] <- parents[ind[2], ] * weight +
                            parents[ind[1], ] * (1 - weight)
    }
    
    if (method == "flat")                                          # flat crossover
    {
      for (j in 1:genes.n)                                         # for each gene
      {                                                            # and for both
        for (t in 1:2)                                             # children
        {
          children[ind[t], j] <- runif(n = 1,                      # randomly pick    
                                       min(parents[ind[1], j],     # the value between   
                                           parents[ind[2], j]),    # the genes of
                                       max(parents[ind[1], j],     # the parents
                                           parents[ind[2], j]))
        }
      }
    }
  }

  return(children)                                         # return the children
}

# Assign default parameters for
# crossover algorithm depending
# on the "method"
gena.crossover.validate <- function(method, par)
{
  # Validate the "method"
  
  methods <- c("split", "arithmetic", "local", "flat")     # the list of all
                                                           # available methods
  
  if (!(method %in% methods))
  {
    stop(paste0("Incorrect crossover.method argument. ",   # if the user has provided
                "It should be one of: ",                   # incorrect argument
                paste(methods, collapse = ", "),
                "\n"))
  }
  
  # Assign default parameters
  
  if (method == "arithmetic")
  {
    if (!is.null(par))
    {
      if ((length(par) != 1) | (!is.numeric(par)))
      {
        stop(paste0("Incorrect crossover.par agrument. Please, insure that ",
                    "(length(crossover.par) == 1) and ",
                    "is.numeric(crossover.par)",
                    "\n"))
      }
    } else {
      par <- 0.5
    }
  }
  
  return(par)
}