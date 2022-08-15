#' Population
#' @description Initialize the population of chromosomes.
#' @param pop.n positive integer representing the number of chromosomes 
#' in population.
#' @param lower numeric vector which i-th element determines the minimum
#' possible value for i-th gene.
#' @param upper numeric vector which i-th element determines the maximum
#' possible value for i-th gene.
#' @param pop.initial numeric matrix which rows are initial chromosomes
#' suggested by user.
#' @param method string representing the initialization method to be used.
#' For a list of possible values see Details.
#' @details If \code{"method = uniform"} then i-th gene of each chromosome is randomly 
#' (uniformly) chosen between \code{lower[i]} and \code{upper[i]} bounds. If
#' \code{"method = normal"} then i-th gene is generated from a truncated
#' normal distribution with mean \code{(upper[i] + lower[i]) / 2} and
#' standard deviation \code{(upper[i] - lower[i]) / 6} where \code{lower[i]}
#' and \code{upper[i]} are lower and upper truncation bounds correspondingly.
#' If \code{"method = hypersphere"} then population is simulated uniformly
#' from the hypersphere with center \code{upper - lower} and radius
#' \code{sqrt(sum((upper - lower) ^ 2))} via \code{\link[gena]{rhypersphere}}
#' function setting \code{type = "inside"}.
#' @returns This function returns a matrix which rows are chromosomes.
#' @references B. Kazimipour, X. Li, A. Qin (2014). 
#' A review of population initialization techniques for evolutionary algorithms.
#' \emph{2014 IEEE Congress on Evolutionary Computation}, 2585-2592,
#' <doi:10.1109/CEC.2014.6900618>.
#' @examples
#' set.seed(123)
#' gena.population(pop.n = 10,
#'                 lower = c(-1, -2, -3),
#'                 upper = c(1, 0, -1),
#'                 pop.initial = rbind(c(0, -1, -2),
#'                                     c(0.1, -1.2, -2.3)),
#'                 method = "normal")
                
gena.population <- function(pop.n,                         # the size of the population
                            lower,                         # lower bounds for the parameters
                            upper,                         # upper bounds for the parameters
                            pop.initial = NULL,            # initial chromosomes suggested by user
                            method = "uniform")            # population initialization algorithm
{
  # Validation
  
  methods <- c("uniform", "normal", "hypersphere")         # the list of all
                                                           # available methods
    
  if (!(method %in% methods))
  {
    stop(paste0("Incorrect population.method argument. ",  # if the user has provided
                "It should be one of: ",                   # incorrect argument
                paste(methods, collapse = ", "),
                "\n"))
  }

  # Prepare some values
  
  genes.n <- length(lower)                                 # the number of parameters
  pop.n.0 <- NULL                                          # the number of population members
                                                           # suggested by the user
  pop.n.1 <- pop.n                                         # the number of population members
                                                           # to be created    
  # Deal with chromosomes proposed
  # by the user
  
  if (!is.null(pop.initial))                               # if the user has provided any
  {                                                        # initial chromosomes
    pop.n.0 <- nrow(pop.initial)                           # the number of chromosomes 
                                                           # provided by the user
    if (pop.n.0 == pop.n)                                  # if the user has provided
    {                                                      # the entire population
      return(pop.initial)                                  # simply return the user
    }                                                      # provided population
    
    pop.n.1 <- pop.n - pop.n.0                             # the number of chromosomes
                                                           # to be created
  }

  # Initialize the population
  
  population <- matrix(NA,                                 # population matrix
                       nrow = pop.n.1, 
                       ncol = genes.n)   
  
  if (method == "uniform")                                 # create the chromosomes using
  {                                                        # the uniform distribution
    for (i in 1:genes.n)
    {
      population[, i] <- runif(pop.n.1, 
                               lower[i], upper[i])
    }
  }
  
  if (method == "normal")                                  # create the chromosomes using
  {                                                        # normal distribution
    mu <- (upper + lower) / 2  
    sigma <- (upper - mu) / 3
    p_lower <- pnorm(lower)
    p_upper <- pnorm(upper)
    p_diff <- p_upper - p_lower
    for (i in 1:genes.n)
    {
      u <- runif(pop.n.1)
      population[, i] <- qnorm(p_lower[i] + 
                               u * p_diff[i]) *
                         sigma[i] + mu[i]
    }
  }
  
  if (method == "hypersphere")
  {
    population <- rhypersphere(n = pop.n.1, dim = genes.n,
                               radius = sqrt(sum((upper - lower) ^ 2)),
                               center = upper - lower,
                               type = "inside")
  }
  
  population <- rbind(pop.initial, population)             # combine the chromosomes provided
                                                           # by the user and the chromosomes
                                                           # that were initialized
                                                         
  return(population)                                       # return the population of chromosomes
}
