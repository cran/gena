#' Velocity
#' @description Calculates (updates) velocities of the particles.
#' @param population numeric matrix which rows are particles i.e. vectors of 
#' parameters values.
#' @param method string representing method to be used for velocities 
#' calculation. See 'Details' for additional information.
#' @param velocity matrix which i-th row is a velocity of the i-th particle. 
#' @param par additional parameters to be passed depending on the \code{method}.
#' @param best.pn numeric matrix which i-th row is a best personal position 
#' known by the i-th particle.
#' @param best.nh numeric matrix which i-th row is a best personal position 
#' in a neighbourhood of the i-th particle.
#' @param best.pn.fitness numeric vector which i-th row is the value of
#' a fitness function at point \code{best.pn[i, ]}.
#' @param best.nh.fitness numeric vector which i-th row is the value of
#' a fitness function at point \code{best.nh[i, ]}.
#' @param iter iteration number of the genetic algorithm.
#' @details If \code{method = "classic"} then classical velocity formula
#' is used:
#' \deqn{v_{i,j,(t+1)}=w\times v_{i,j,t} + 
#'                      c_{1}\times u_{1,i,j} \times b^{pn}_{i,j,t} + 
#'                      c_{2}\times u_{2,i,j} \times b^{nh}_{i,j,t}}
#'                      
#'  where \eqn{v_{i, j, t}} is a velocity of the \eqn{i}-th particle 
#'  respect to the \eqn{j}-th component at time \eqn{t}. Random variables
#'  \eqn{u_{1,i,j}} and \eqn{u_{2,i,j}} are i.i.d. respect to all indexes and
#'  follow standard uniform distribution \eqn{U(0, 1)}. 
#'  Variable \eqn{b^{pn}_{i,j,t}} is \eqn{j}-th component of the best known
#'  particle's (personal) position up to time period \eqn{t}. 
#'  Similarly \eqn{b^{nh}_{i,j,t}} is \eqn{j}-th component of the best of best 
#'  known particle's position in a neighbourhood of the \eqn{i}-th particle.
#'  Hyperparameters \eqn{w}, \eqn{c_{1}} and \eqn{c_{2}} may be provided
#'  via \code{par} argument as a list with elements \code{par$w}, \code{par$c1}
#'  and \code{par$c2} correspondingly. 
#'  
#'  If \code{method = "hypersphere"} then rotation invariant formula from
#'  sections 3.4.2 and 3.4.3 of M. Clerc (2012) is used with arguments
#'  identical to the classical method. To simulate a random variate from
#'  the hypersphere function \code{\link[gena]{rhypersphere}} is used
#'  setting \code{type = "non-uniform"}.
#'  
#'  In accordance with M. Clerc (2012) 
#'  default values are \code{par$w = 1/(2 * log(2))},
#'  \code{par$c1 = 0.5 + log(2)} and \code{par$c2 = 0.5 + log(2)}.
#'  
#' @return  This function returns a matrix which i-th row represents
#' updated velocity of the i-th particle.
#'  
#' @references Maurice Clerc (2012). 
#' Standard Particle Swarm Optimisation.
#' \emph{HAL archieve}.
#'              
pso.velocity <- function(population,
                         method = "hypersphere",
                         par = list(w = 1 / (2 * log(2)), 
                                    c1 = 0.5 + log(2),
                                    c2 = 0.5 + log(2)),
                         velocity,
                         best.pn,
                         best.nh,
                         best.pn.fitness,
                         best.nh.fitness,
                         iter = 1)              
{
  # Dimensional variables
  
  pop.n <- nrow(population)
  genes.n <- ncol(population)
  total.n <- genes.n * pop.n
  
  # Determine best particles in a neighbourhood
  
  is.best <- best.pn.fitness == best.nh.fitness
  
  # Matrix of updated velocities
  
  velocity.new <- NULL
  
  # Update velocities according to method
  
  if (method == "classic")
  {
    # Simulate from uniform distribution
    
    unif.1 <- matrix(runif(total.n), nrow = pop.n, ncol = genes.n)
    unif.2 <- matrix(runif(total.n), nrow = pop.n, ncol = genes.n)
    unif.2[is.best] <- 0
    
    # Calculate the velocity
    
    velocity.new <- par$w * velocity +
                    par$c1 * unif.1 * (best.pn - population) +
                    par$c2 * unif.2 * (best.nh - population)
  }
  
  if (method == "hypersphere")
  {
    # Simulate uniform random variates from the hypersphere
    
    adj <- matrix(NA, nrow = pop.n, ncol = genes.n)
    adj[is.best] <- (par$c1 / 2) * (best.pn[is.best, ] - population[is.best, ])
    adj[!is.best] <- (par$c1 / 3) * (best.pn[!is.best, ] - population[!is.best, ]) +
                     (par$c2 / 3) * (best.nh[!is.best, ] - population[!is.best, ])
    center <- population + adj
    radius <- sqrt(rowSums(adj ^ 2))
    x <- rhypersphere(n = pop.n, 
                      dim = genes.n, 
                      radius = radius,
                      center = center,
                      type = "non-uniform")

    # Calculate the velocity
    
    velocity.new <- par$w * velocity + x - population
  }
  
  return(velocity.new)
}

# Assign default parameters for
# velocity algorithm depending
# on the "method"
pso.velocity.validate <- function(method, par)
{
  # Validate the "method"
  
  methods <- c("classic", "hypersphere")                   # the list of all
                                                           # available methods
  if (!(method %in% methods))
  {
    stop(paste0("Incorrect nh.method argument. ",          # if the user has provided
                "It should be one of: ",                   # incorrect argument
                paste(methods, collapse = ", "),
                "\n"))
  }
  
  if (method %in% c("classic", "hypersphere"))
  {
    if (is.null(par))
    {
      par <- list(w = 1 / (2 * log(2)), 
                  c1 = 0.5 + log(2),
                  c2 = 0.5 + log(2))
    }
    else
    {
      
      if (!is.list(par))
      {
        par.vec <- par
        par.vec.n <- length(par.vec)
        if ((par.vec.n != 2) & (par.vec.n != 3))
        {
          stop("If 'nh.par' is a vector then it should contain 2 or 3 elements.")
        }
        par <- list(w = par.vec[1], 
                    c1 = par.vec[2])
        if (length(par.vec) == 2)
        {
          par$c2 <- par$c1
        }
      }
      else
      {
        if (!all(c("w", "c1", "c2") %in% names(par)))
        {
          stop("Argument 'nh.par' should be a list containing 'w', 'c1' and 'c2'.")
        }
      }
    }
  }
  
  return(par)
}