#' Neighbourhood
#' @description Constructs a neighbourhood of each particle using 
#' particular topology.
#' @param pop.n integer representing the size of the population.
#' @param method string representing the topology to be used for construction
#' of the neighbourhood. See 'Details' for additional information.
#' @param par additional parameters to be passed depending on the \code{method}.
#' @param iter iteration number of the genetic algorithm.
#' @details If \code{method = "ring"} then each particle will have 
#' \code{par[1]} neighbours. By default \code{par[1] = 3}.
#' See section 3.2.1 of M. Clerc (2012) for 
#' additional details.
#' If \code{method = "wheel"} then there is a single (randomly selected) 
#' particle which informs (and informed by) other particles while there is 
#' no direct communication between other particles. 
#' If \code{method = "random"} then each particle randomly informs other
#' \code{par[1]} particles and itself. Note that duplicates are possible so 
#' sometimes each particle may inform less then \code{par[1]} particles. 
#' By default \code{par[1] = 3}.
#' See section 3.2.2 of M. Clerc (2012) for more details.
#' If \code{method = "star"} then all particles are fully informed
#' by each other.
#' If \code{method = "random2"} then each particle will be self-informed
#' and informed by the j-th particle with probability \code{par[1]}
#' (value between 0 and 1). By default \code{par[1] = 0.1}.
#' @return This function returns a list which i-th element is a vector of
#' particles' indexes which inform i-th particle i.e. neighbourhood of the
#' i-th particle.
#' 
#' @references Maurice Clerc (2012). 
#' Standard Particle Swarm Optimisation.
#' \emph{HAL archieve}.
#' 
#' @examples 
#' # Prepare random number generator
#' set.seed(123)
#' 
#' # Ring topology with 5 neighbours
#' pso.nh(pop.n = 10, method = "ring", par = 5)
#' 
#' # Wheel topology
#' pso.nh(pop.n = 10, method = "wheel")
#' 
#' # Star topology
#' pso.nh(pop.n = 10, method = "star")
#' 
#' # Random topology where each particle 
#' # randomly informs 3 other particles
#' pso.nh(pop.n = 10, method = "random", par = 3)
#' 
#' # Random2 topology wehere each particle could
#' # be informed by the other with probability 0.2
#' pso.nh(pop.n = 10, method = "random2", par = 0.2)
#' 
pso.nh <- function(pop.n = 40,
                   method = "ring",
                   par = 3,
                   iter = 1)              
{
  # Prepare list to store neighbourhood
  # of each particle
  
  nh <- vector(mode = "list", length = pop.n)
  
  if (method == "ring")
  {
    nh.adj <- floor((par[1] - 1) / 2)
    is.odd <- (par[1] %% 2 == 0)
    nh.elements <- c((pop.n - nh.adj + 1):pop.n,
                     1:pop.n,
                     1:(nh.adj + is.odd))
    nh.ind <- 0:(par[1] - 1)
    for (i in 1:pop.n)
    {
      nh[[i]] <- nh.elements[nh.ind + i]
    }
  }
  
  if (method == "wheel")
  {
    wheel.particle = sample(1:pop.n, size = 1)
    nh.mat <- cbind(wheel.particle, 1:pop.n)
    nh <- as.list(as.data.frame(t(nh.mat)))
  }
  
  if (method == "star")
  {
    for (i in 1:pop.n)
    {
      nh[[i]] <- 1:pop.n
    }
  }
  
  if (method == "random")
  {
    pop.ind <- 1:pop.n
    for (i in pop.ind)
    {
      nh[[i]] <- i
    }
    for (i in pop.ind)
    {
      informs <- sample(x = pop.ind[-i], 
                        size = par[1], 
                        replace = TRUE)
      for (j in informs)
      {
        nh[[j]] <- c(nh[[j]], i)
      }
    }
    for (i in pop.ind)
    {
      nh[[i]] <- unique(nh[[i]])
    }
  }
  
  if (method == "random2")
  {
    pop.ind <- 1:pop.n
    for (i in pop.ind)
    {
      informs.ind <- as.logical(rbinom(n = pop.n - 1, size = 1, prob = par[1]))
      nh[[i]] <- c(i, pop.ind[-i][informs.ind])
    }
  }
  
  return(nh)
}

# Assign default parameters for
# neighbourhood algorithm depending
# on the "method"
pso.nh.validate <- function(method, par, pop.n)
{
  # Validate the "method"
  
  methods <- c("ring", "random", "random2",                # the list of all
               "star", "wheel")                            # available methods
                                                           
  if (!(method %in% methods))
  {
    stop(paste0("Incorrect 'nh.method' argument. ",        # if the user has provided
                "It should be one of: ",                   # incorrect argument
                paste(methods, collapse = ", "),
                "\n"))
  }
  
  if (method %in% c("ring"))
  {
    if (!is.null(par))
    {
      if ((par[1] > pop.n) | (par[1] < 3))
      {
        stop(paste0("Incorrect 'nh.par' agrument. Please, insure that ",
                    "'nh.par' is not less then 3 and not greater",
                    "then 'pop.n'. \n"))
      }
    }
  }
  
  if (method %in% c("random"))
  {
    if (!is.null(par))
    {
      if (!is.numeric(par))
      {
        stop(paste0("Incorrect 'nh.par' agrument. Please, insure that ",
                    "'nh.par' is numeric. \n"))
      }
    } else {
      par <- 3
    }
  }
  
  if (method %in% c("random2"))
  {
    if (!is.null(par))
    {
      if ((par[1] < 0) | (par[1] > 1))
      {
        stop(paste0("Incorrect 'nh.par' agrument. Please, insure that ",
                    "'nh.par' is numeric value between 0 and 1. \n"))
      }
    }
    else
    {
      par <- 0.1
    }
  }
  
  return(par)
}