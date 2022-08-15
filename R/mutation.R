#' Mutation
#' @description Mutation method (algorithm) to be used in the
#' genetic algorithm.
#' @param children numeric matrix which rows are children i.e. vectors of 
#' parameters values.
#' @param lower lower bound of the search space.
#' @param upper upper bound of the search space.
#' @param prob probability of mutation for a child.
#' @param prob.genes numeric vector or numeric value representing the
#' probability of mutation of a child's gene. See 'Details'.
#' @param method mutation method to be used for transforming genes of children.
#' @param par additional parameters to be passed depending on the \code{method}.
#' @param iter iteration number of the genetic algorithm.
#' 
#' @details Denote \code{children} by \eqn{C^{child}} which \code{i}-th row 
#' \code{children[i, ]} is a chromosome \eqn{c_{i}^{child}} i.e. the vector of 
#' parameter values of the function being optimized \eqn{f(.)} that is
#' provided via \code{fn} argument of \code{\link[gena]{gena}}.
#' The elements of chromosome \eqn{c_{ij}^{child}} are genes 
#' representing parameters values.
#' 
#' Mutation algorithm determines random transformation of children's genes. 
#' Each child may be selected for mutation with probability \code{prob}.
#' If \eqn{i}-th child is selected for mutation and \code{prob.genes} is a 
#' vector then \eqn{j}-th gene of this child
#' is transformed with probability \code{prob.genes[j]}. If \code{prob.genes}
#' is a constant then this probability is the same for all genes.
#' 
#' Argument \code{method} determines particular mutation algorithm to 
#' be applied. Denote by \eqn{\tau} the vector of parameters used by the 
#' algorithm. Note that \eqn{\tau} corresponds to \code{par}.
#' Also let's denote by \eqn{c_{ij}^{mutant}} the value of 
#' gene \eqn{c_{ij}^{child}} after mutation.
#' 
#' If \code{method = "constant"} then \eqn{c_{ij}^{mutant}}
#' is a uniform random variable between \code{lower[j]} and \code{upper[j]}.
#' 
#' If \code{method = "normal"} then \eqn{c_{ij}^{mutant}}
#' equals to the sum of \eqn{c_{ij}^{child}} and normal random variable
#' with zero mean and standard deviation \code{par[j]}. 
#' By default \code{par} is identity vector of length \code{ncol(children)}
#' so \code{par[j] = 1} for all \code{j}.
#' 
#' If \code{method = "percent"} then \eqn{c_{ij}^{mutant}} is generated 
#' from \eqn{c_{ij}^{child}}  by equiprobably increasing or decreasing it 
#' by \eqn{q} percent,
#' where \eqn{q} is a uniform random variable between \eqn{0} and \code{par[j]}. 
#' Note that \code{par} may also be a constant then all
#' genes have the same maximum possible percentage change.
#' By default \code{par = 20}.
#' 
#' For more information on mutation algorithms
#' please see Patil, Bhende (2014).
#' 
#' @return The function returns a matrix which rows are children
#' (after mutation has been applied to some of them).
#' 
#' @references S. Patil, M. Bhende. (2014). 
#' Comparison and Analysis of Different Mutation Strategies to improve the 
#' Performance of Genetic Algorithm.
#' \emph{International Journal of Computer Science and 
#' Information Technologies}, 5 (3), 4669-4673.
#' 
#' @examples 
#' # Randomly initialize some children
#' set.seed(123)
#' children.n <- 10
#' children <- gena.population(pop.n = children.n,
#'                             lower = c(-5, -5), 
#'                             upper = c(5, 5))
#'                            
#' # Perform the mutation
#' mutants <- gena.mutation(children = children,
#'                          prob = 0.6,
#'                          prob.genes = c(0.7, 0.8),
#'                          par = 30,
#'                          method = "percent")
#' print(mutants)

# Mutation operator
gena.mutation <- function(children,                        # children
                          lower,                           # lower bound of search space
                          upper,                           # upper bound of search space
                          prob = 0.2,                      # mutation probability for children
                          prob.genes = 1 / nrow(children), # mutation probabilities for genes
                          method = "constant",             # mutation type
                          par = 1,                         # mutation parameters
                          iter = NULL)                     # iteration of the
                                                           # genetic algorithm
{
  # ---
  # The algorithm brief description:
  # 1. Randomly pick the children for mutation
  # 2. Perform the mutation for each gene of each children
  # ---
  
  # Prepare some variables
  
  children.n <- nrow(children)                             # number of children
  genes.n <- ncol(children)                                # number of genes
  
  # Select children for mutation

  random_values <- runif(children.n, 0, 1)                 # values to control for
                                                           # weather mutation should
                                                           # take place for a given
                                                           # child
  mutation_ind <- which(random_values <= prob)             # select the children
                                                           # to mutate
  
  for (i in mutation_ind)
  {
    random_values_genes <- runif(genes.n, 0, 1)
    genes_ind <- which(prob.genes <= random_values_genes)
    
    if (method == "constant")
    {
      for (j in genes_ind)
      {
        children[i, j] <- runif(1, lower[j], upper[j])
      }
    }
    
    if (method == "normal")
    {
      for (j in genes_ind)
      {
        children[i, j] <- children[i, j] + rnorm(1, mean = 0, sd = par[j])
        #children[i, j] <- min(max(children[i, j], lower[j]), upper[j])
      }
    }
    
    if (method == "percent")
    {
      par_adj <- par / 100
      children_adj <- abs(children[i, ]) * par_adj
      mutation_lower <- children[i, ] - children_adj
      mutation_upper <- children[i, ] + children_adj
      for (j in genes_ind)
      {
        children[i, j] <- runif(1, mutation_lower[j], mutation_upper[j])
        #children[i, j] <- min(max(children[i, j], lower[j]), upper[j])
      }
    }
  }
  
  return(children)                                         # return the children
}

# Assign default parameters for
# mutation algorithm depending
# on the "method"
gena.mutation.validate <- function(method, par, genes.n)
{
  # Validate the "method"
  
  methods <- c("constant", "normal", "percent")
  
  if (!(method %in% methods))
  {
    stop(paste0("Incorrect mutation.method argument. ",    # if the user has provided
                "It should be one of: ",                   # incorrect argument
                paste(methods, collapse = ", "),
                "\n"))
  }
  
  # Assign default parameters
  
  if (method == "normal")
  {
    if (is.null(par))
    {
      par <- 1
    } else {
      if (any((par <= 0)))
      {
        stop("par should be a vector with posive values")
      }
    }
    
    if (length(par) != genes.n)
    {
      if (length(par) == 1)
      {
        par <- rep(par, genes.n)
      } else {
        stop("length(par) should be equal to length(lower)")
      }
    }
  }
  
  if (method == "percent")
  {
    if (is.null(par))
    {
      par <- rep(20, genes.n)
    } else {
      if (length(par) == 1)
      {
        par <- rep(par, genes.n)
      } else {
        if (length(par) > genes.n)
        {
          stop(paste0("Please, insure that 'par' has the same size as the ",
                      "number of genes i.e. 'length(par)==ncol(lower)'"))
        }
      }
      
      if (any((par < 0)))
      {
        stop("par should not have non-negative values")
      }
    }
    
    if (length(par) != genes.n)
    {
      if (length(par) == 1)
      {
        par <- rep(par, genes.n)
      } else {
        stop("length(par) should be equal to length(lower)")
      }
    }
  }
  
  return(par)
}