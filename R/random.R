#' Hypersphere
#' @description Simulates uniform random variates from the hypersphere.
#' @param n number of observations to simulate.
#' @param dim dimensions of hypersphere.
#' @param radius radius of hypersphere.
#' @param center center of hypersphere.
#' @param type character; if \code{"boundary"} (default) then random variates 
#' are simulated from the hypersphere. If \code{"inside"} random
#' variates are points lying inside the hypersphere. If \code{"non-uniform"} 
#' then random variates are non-uniform and simulated from the inner part
#' of the hypersphere simply by making radius a uniform random variable
#' between \code{0} and \code{radius}.
#' @return The function returns a vector of random variates.
#' @examples 
#' set.seed(123)
#' # Get 5 random uniform variates from 3D hypersphere
#' # of radius 10 centered at (2, 3, 1)
#' rhypersphere(n = 5, dim = 3, radius = 10, center = c(2, 3, 1))
#' 
rhypersphere <- function(n, dim = 2,
                         radius = 1, center = rep(0, dim),
                         type = "boundary")              
{
  # Generate independent standard normal variates
  
  x <- matrix(rnorm(n * dim), nrow = n, ncol = dim)
  
  # Transform these variates
  
  y <- sweep(x, 1, sqrt(rowSums(x ^ 2)), "/")
  
  # Select points inside sphere or on the boundary
  
  if (type == "inside")
  {
    radius <- (radius * runif(n)) ^ (1 / dim)
  }
  
  if (type == "non-uniform")
  {
    radius <- radius * runif(n)
  }
  
  # Adjust for radius
  
  if (length(radius) == 1)
  {
    y <- y * radius
  }
  else
  {
    y <- sweep(y, 1, radius, "*")
  }

  # Adjust for center
  
  if (!is.matrix(center))
  {
    y <- sweep(y, 2, center, "+")
  }
  else
  {
    y <- y + center
  }
  
  # Return the results
  
  return(y)
}