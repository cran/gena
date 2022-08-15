#' Particle Swarm Optimization
#' @description This function allows to use particle swarm algorithm for
#' numeric global optimization of real-valued functions.
#' @param fn function to be maximized i.e. fitness function.
#' @param gr gradient of the \code{fn}.
#' @param lower lower bound of the search space.
#' @param upper upper bound of the search space.
#' @param pop.initial numeric matrix which rows are particles to be
#' included into the initial population. Numeric vector will be coerced to
#' single row matrix.
#' @param pop.n integer representing the size of the population.
#' @param pop.method the algorithm to be applied for a creation of 
#' the initial population. See 'Details' for additional information.
#' @param nh.method string representing the method (topology) to be used for
#' the creation of neighbourhoods. See 'Details' for additional information.
#' @param nh.par parameters of the topology algorithm.
#' @param nh.adaptive logical; if \code{TRUE} (default) then neighbourhoods
#' change every time when the best known (to the swarm) fitnesses value have
#' not increased. Neighbourhoods are updated according to the topology
#' defined via \code{nh.method} argument.
#' @param velocity.method string representing the method to be used for
#' the update of velocities.
#' @param velocity.par parameters of the velocity formula.
#' @param hybrid.method hybrids selection algorithm i.e. mechanism 
#' determining which particles should be subject to local optimization.
#' See 'Details' for additional information.
#' @param hybrid.par parameters of the hybridization algorithm.
#' @param hybrid.prob probability of generating the hybrids each iteration.
#' @param hybrid.opt.par parameters of the local optimization function
#' to be used for hybridization algorithm (including \code{fn} and \code{gr}).
#' @param hybrid.n number of hybrids that appear if hybridization
#' should take place during the iteration.
#' @param constr.method the algorithm to be applied for imposing constraints
#' on the particles. See 'Details' for additional information.
#' @param constr.par parameters of the constraint algorithm.
#' @param random.order logical; if \code{TRUE} (default) then particles
#' related routine will be implemented in a random order.
#' @param maxiter maximum number of iterations of the algorithm.
#' @param is.max logical; if \code{TRUE} (default) then fitness function
#' will be maximized. Otherwise it will be minimized.
#' @param info logical; if \code{TRUE} (default) then some optimization related 
#' information will be printed each iteration.
#' @param ... additional parameters to be passed to 
#' \code{fn} and \code{gr} functions.
#' 
#' @details Default arguments have been set in accordance with SPSO 2011
#' algorithm proposed by M. Clerc (2012).
#' 
#' To find information on particular methods available via
#' \code{pop.method}, \code{nh.method}, \code{velocity.method},
#' \code{hybrid.method} and \code{constr.method}
#' arguments please see 'Details' section of 
#' \code{\link[gena]{gena.population}}, \code{\link[gena]{pso.nh}},
#' \code{\link[gena]{pso.velocity}}, \code{\link[gena]{gena.hybrid}} 
#' and \code{\link[gena]{gena.constr}} correspondingly.
#' 
#' It is possible to provide manually implemented functions for population
#' initialization, neighbourhoods creation, velocity updated, hybridization
#' and constraints in a similar way as for \code{\link[gena]{gena}}.
#' 
#' By default function does not impose any constraints upon the parameters.
#' If \code{constr.method = "bounds"} then \code{lower} and \code{upper} 
#' constraints will be imposed. Lower bounds should be strictly smaller
#' then upper bounds.
#' 
#' Currently the only available termination condition is \code{maxiter}. We
#' are going to provide some additional termination conditions during
#' future updates.
#' 
#' Infinite values in \code{lower} and \code{upper} are substituted with
#' \code{-(.Machine$double.xmax * 0.9)} and \code{.Machine$double.xmax * 0.9}
#' correspondingly.
#' 
#' By default if \code{gr} is provided then BFGS algorithm will be used inside
#' \code{\link[stats]{optim}} during hybridization.
#' Otherwise \code{Nelder-Mead} will be used. 
#' Manual values for \code{\link[stats]{optim}} arguments may be provided 
#' (as a list) through \code{hybrid.opt.par} argument.
#' 
#' For more information on particle swarm optimization 
#' please see M. Clerc (2012).
#' 
#' @return This function returns an object of class \code{pso} that is a list
#' containing the following elements:
#' \itemize{
#' \item \code{par} - particle (solution) with the highest fitness
#' (objective function) value.
#' \item \code{value} - value of \code{fn} at \code{par}.
#' \item \code{population} - matrix of particles (solutions) of the 
#' last iteration of the algorithm.
#' \item \code{counts} - a two-element integer vector giving the number of
#' calls to \code{fn} and \code{gr} respectively.
#' \item \code{is.max} - identical to \code{is.max} input argument.
#' \item \code{fitness.history} - vector which i-th element is fitness
#' of the best particle in i-th iteration.
#' \item \code{iter} - last iteration number.
#' }
#' 
#' @references M. Clerc (2012). 
#' Standard Particle Swarm Optimisation.
#' \emph{HAL archieve}.
#' 
#' @examples
#' ## Consider Ackley function
#' \donttest{
#' fn <- function(par, a = 20, b = 0.2)
#' {
#'   val <- a * exp(-b * sqrt(0.5 * (par[1] ^ 2 + par[2] ^ 2))) +
#'          exp(0.5 * (cos(2 * pi * par[1]) + cos(2 * pi * par[2]))) -
#'          exp(1) - a
#'   return(val)
#' }
#' 
#' # Maximize this function using particle swarm algorithm
#' 
#' set.seed(123)
#' lower <- c(-5, -100)
#' upper <- c(100, 5)
#' opt <- pso(fn = fn, 
#'            lower = lower, upper = upper,
#'            a = 20, b = 0.2)
#' print(opt$par)
#' }
#' 
#' ## Consider Bukin function number 6
#' 
#' fn <- function(x, a = 20, b = 0.2)
#' {
#'   val <- 100 * sqrt(abs(x[2] - 0.01 * x[1] ^ 2)) + 0.01 * abs(x[1] + 10)
#'   return(val)
#' }
#' 
#' # Minimize this function using initially provided
#' # position for one of the particles
#' set.seed(777)
#' lower <- c(-15, -3)
#' upper <- c(-5, 3)
#' opt <- pso(fn = fn, 
#'            pop.init = c(8, 2),
#'            lower = lower, upper = upper,
#'            is.max = FALSE)
#' print(opt$par)
#' 
pso <-  function(fn, 
                 gr = NULL, 
                 lower, 
                 upper,
                 pop.n = 40,
                 pop.initial = NULL,
                 pop.method = "uniform",
                 nh.method = "random",
                 nh.par = 3,
                 nh.adaptive = TRUE,
                 velocity.method = "hypersphere",
                 velocity.par = list(w = 1 / (2 * log(2)), 
                                     c1 = 0.5 + log(2),
                                     c2 = 0.5 + log(2)),
                 hybrid.method = "rank",
                 hybrid.par = 2,
                 hybrid.prob = 0,
                 hybrid.opt.par = NULL,
                 hybrid.n = 1,
                 constr.method = NULL,
                 constr.par = NULL,
                 random.order = TRUE,
                 maxiter = 100,
                 is.max = TRUE,
                 info = TRUE,
                 ...)
{
  # Get the parameters to be passed to fitness
  
  dots <- list(...)
  
  # Initialize some variables
  
  genes.n <- length(lower)
  
  # Perform initial validation
  
    # fn
  
  if (!is.function(fn))
  {
    stop("Please, insure that 'fn' is a function.\n")
  }
  
    # gr
  if(!is.null(gr))
  {
    if (!is.function(gr))
    {
      stop("Please, insure that 'gr' is a function.\n")
    }
  }
  
    # lower and upper
  
  if (length(lower) != length(upper))
  {
    stop("Please, insure that 'lower' and 'upper' are of the same length.\n")
  }
  
  if (any(lower > upper))
  {
    stop(paste0("Please, insure that all elements of the 'lower' are ",
                "smaller than the corresponding elements of the 'upper.'",
                "\n"))
  }
  
  lower[is.infinite(lower)] <- -(.Machine$double.xmax * 0.9)
  upper[is.infinite(upper)] <- .Machine$double.xmax * 0.9
  
    # pop.initial
  
  if (!is.null(pop.initial))
  {
    if (!is.numeric(pop.initial))
    {
      stop("Please, insure that 'pop.initial' is a numeric matrix.\n")
    }
    
    if (!is.matrix(pop.initial))
    {
      pop.initial <- matrix(pop.initial, nrow = 1)
    }
    
    if (ncol(pop.initial) != genes.n)
    {
      stop(paste0("Please, insure that the number of columns of 'pop.initial' ",
                  "equals to the number of estimated parameters i.e. genes.",
                  "\n"))
    }
  }
  
    # pop.n
  
  if (pop.n <= 1)
  {
    stop("Please, insure that 'pop.n' is an integer greater than one.")
  }
  
  if ((pop.n %% 2) != 0) 
  {
    pop.n <- pop.n + 1
    warning("Since 'pop.n' is not even it has been incremented by one.\n")
  }
  
  # Initialize some additional variables  
  
  fitness <- rep(NA, pop.n)
  
  # Continue validation
  
    # hybrid.method and hybrid.par
  
  if (!is.function(hybrid.method))
  {
    hybrid.par <- gena.mating.validate(method = hybrid.method,
                                       par = hybrid.par,
                                       parents.n = pop.n)
  }
  else
  {
    gena.hybrid <- hybrid.method
  }
  
    # hybrid.prob
  
  if (is.numeric(hybrid.prob))
  {
    if ((hybrid.prob < 0) | (hybrid.prob > 1))
    {
      stop("Please, insure that 'hybrid.prob' is between 0 and 1.\n")
    }
  } else {
    stop("Please, insure that 'hybrid.prob' is a numeric value.\n")
  }
  
    # hybrid.n
  
  if (hybrid.n < 0)
  {
    stop("Please, insure that 'hybrid.n' is a non-negative integer.\n")
  }
  
  if (hybrid.n > pop.n)
  {
    stop("Please, insure that 'hybrid.n' is not greater than 'pop.n'.")
  }
  
    # maxiter
  
  if (maxiter < 1)
  {
    stop("Please, insure that 'maxiter' is a positive integer.")
  }
  
    # constr.method and constr.par
  
  if(!is.null(constr.method))
  {
    if (!(is.function(constr.method) | is.character(constr.method)))
    {
      stop("Please, insure that 'constr.method' is a function or character.\n")
    }
    if (is.character(constr.method))
    {
      if ((constr.method == "bounds") & is.null(constr.par))
      {
        constr.par$lower <- lower
        constr.par$upper <- upper
      }
      gena.constr.validate(method = constr.method, par = constr.par)
    }
    if (is.function(constr.method))
    {
      gena.constr <- constr.method
    }
  }
  
    # is.max
  
  if (!is.logical(info))
  {
    if ((info == 0) | (info == 1))
    {
      info <- as.logical(info)
    } else {
      stop("Please, insure that 'is.max' is a logical.\n")
    }
  }
  
    # info
  
  if (!is.logical(info))
  {
    if ((info == 0) | (info == 1))
    {
      info <- as.logical(info)
    } else {
      stop("Please, insure that 'info' is a logical.\n")
    }
  }
  
    # random.order
  
  if (!is.logical(random.order))
  {
    if ((random.order == 0) | (random.order == 1))
    {
      random.order <- as.logical(random.order)
    } else {
      stop("Please, insure that 'random.order' is a logical.\n")
    }
  }
  
    # nh.adaptive
  
  if (!is.logical(nh.adaptive))
  {
    if ((nh.adaptive == 0) | (nh.adaptive == 1))
    {
      nh.adaptive <- as.logical(nh.adaptive)
    } else {
      stop("Please, insure that 'nh.adaptive' is a logical.\n")
    }
  }
  
    # nh.method and nh.par
  
  if (!is.function(nh.method))
  {
    nh.par <- pso.nh.validate(method = nh.method,
                              par = nh.par,
                              pop.n = pop.n)
  }
  else
  {
    pso.nh <- nh.method
  }
  
    # velocity.method and velocity.par
  
  if (!is.function(velocity.method))
  {
    velocity.par <- pso.velocity.validate(method = velocity.method,
                                          par = velocity.par)
  }
  else
  {
    pso.velocity <- velocity.method
  }
  
  # Value to reverse signs of fitnesses if need
  
  is.max.val <- 2 * is.max - 1
  
  # Deal with hybrid optimizer default parameters
  
  hybrid.opt.par.default = list(fn = fn,
                                gr = gr,
                                method = ifelse(is.null(gr), 
                                                "Nelder-Mead", 
                                                "BFGS"),
                                control = list(fnscale = -1 * is.max.val,
                                               abstol = 1e-10,
                                               reltol = 1e-10,
                                               maxit = 1000000000))
  
  if (!is.null(hybrid.opt.par))
  {
    if (!is.null(hybrid.opt.par$control))
    {
      hybrid.opt.par.default$control <- modifyList(hybrid.opt.par.default$control,
                                                   hybrid.opt.par$control)
    }
    hybrid.opt.par <- modifyList(hybrid.opt.par.default, 
                                 hybrid.opt.par)
  } else {
    hybrid.opt.par <- hybrid.opt.par.default
  }
  
  # Get the name of the first argument of the function
  
  fn.par <- names(formals(fn))[1]
  
  # Number of function and gradient evaluations
  
  counts <- c(0, 0)
  names(counts) <- c("function", "gradient")
  
  # -------------------------------------
  # Perform preliminary iteration of PSO
  # -------------------------------------
  
  # Vector of population members indexes
  
  pop.ind <- 1:pop.n
  
  # Define neighbourhoods of particles
  
  nh <- pso.nh(pop.n = pop.n,
               method = nh.method,
               par = nh.par,
               iter = 0)

  # Initialize the population
  
  if (is.function(pop.method))
  {
    gena.population <- pop.method
  }
  
  population <- gena.population(pop.n = pop.n,
                                lower = lower, 
                                upper = upper,
                                method = pop.method,
                                pop.initial = pop.initial)
  
  # Calculate initial fitnesses
  
  fitness <- rep(NA, pop.n)
  for(i in 1:pop.n)
  {
    dots[[fn.par]] <- population[i, ] 
    fitness[i] <- do.call(what = fn, args = dots)
  }
  fitness[is.na(fitness)] <- -Inf
  
  # Adjust for minimization if need
  
  if (!is.max)
  {
    fitness <- -fitness
  }

  # Assign best personal positions
  
  best.pn <- population
  best.pn.fitness <- fitness
  
  # Assign best positions in the neighbourhood
  
  best.nh <- matrix(NA, nrow = pop.n, ncol = genes.n)
  best.nh.fitness <- rep(NA, pop.n)
  for (j in 1:pop.n)
  {
    best.nh.ind <- which.max(best.pn.fitness[nh[[j]]])
    best.nh[j, ] <- population[nh[[j]], , drop = FALSE][best.nh.ind, ]
    best.nh.fitness[j] <- best.pn.fitness[nh[[j]]][best.nh.ind]
  }

  # Calculate initial velocities of the particles
  
  velocity <- matrix(NA, nrow = pop.n, ncol = genes.n)
  for (j in 1:genes.n)
  {
    runif.tmp <- runif(pop.n)
    velocity[, j] <- runif.tmp * lower[j] + (1 - runif.tmp) * upper[j] - 
                     population[, j]
  }

  # Store information on the best particle
  
  best.final.ind <- which.max(best.pn.fitness)
  best.final <- best.pn[best.final.ind, ]
  best.final.fitness <- best.pn.fitness[best.final.ind]
  
  # Historical best fitness values
  
  fitness.history <- NULL

  # -------------------------------------
  # Perform the particle swarm algorithm
  # -------------------------------------
  
  for(i in 1:maxiter)
  {
    # Get the parameters to be passed to fitness
    
    dots <- list(...)
    
    # Randomize the order of particles
    
    if (random.order)
    {
      pop.ind <- sample(x = pop.ind, size = pop.n, replace = FALSE)
    }
    
    # Impose the constraints on the parameters
    
    constr.ind <- NULL
    if (!is.null(constr.method))
    {
      if (is.function(constr.method))
      {
        constr.par$population <- population
        constr.par$iter <- i
        constr.list <- do.call(what = constr.method, 
                               args = constr.par)
        population <- constr.list$population
        constr.ind <- constr.list$constr.ind
      }
      else
      {
        constr.list <- gena.constr(population = population, 
                                   iter = i,
                                   method = constr.method,
                                   par = constr.par)
        population <- constr.list$population
        constr.ind <-constr.list$constr.ind
      }
    }
    if (!is.null(constr.ind))
    {
      velocity[constr.ind] <- -0.5 * velocity[constr.ind]
    }
    
    # Update velocities and position
    
    velocity <- pso.velocity(population = population,
                             velocity = velocity,
                             method = velocity.method,
                             par = velocity.par,
                             best.pn = best.pn,
                             best.nh = best.nh,
                             best.pn.fitness = best.pn.fitness, 
                             best.nh.fitness = best.nh.fitness,
                             iter = i)
    population <- population + velocity
    
    # Estimate the fitness of each particle
    
    for(j in pop.ind)
    {
      dots[[fn.par]] <- population[j, ] 
      fitness[j] <- do.call(what = fn, args = dots)
    }
    fitness[is.na(fitness) | is.nan(fitness)] <- -Inf
    
    # Make the hybrids
    
    hybrid.value <- runif(1, 0, 1)
    
    if (((hybrid.value <= hybrid.prob) |
         (i == 1) | (i == maxiter)) & 
        (hybrid.n >= 1) & (hybrid.prob > 0))
    {
      hybrids.list <- gena.hybrid(population = population,
                                  fitness = fitness,
                                  hybrid.n = hybrid.n,
                                  method = hybrid.method,
                                  par = hybrid.par,
                                  opt.par = hybrid.opt.par,
                                  info = info,
                                  ... = ...)
      population <- hybrids.list$population
      fitness <- hybrids.list$fitness
      counts <- hybrids.list$counts + counts
      # velocity[hybrids.list$hybrids.ind, ] <- rep(0, genes.n)
    }
    
    # Reverse signs of fitnesses for minimization task
    
    if (!is.max)
    {
      fitness <- -fitness
    }
    
    # Update best local position
  
    is_better <- fitness > best.pn.fitness
    best.pn[is_better, ] <- population[is_better, ]
    best.pn.fitness[is_better] <- fitness[is_better]

    # Update best global position

    for (j in pop.ind)
    {
      best.nh.ind <- which.max(best.pn.fitness[nh[[j]]])
      candidate.fitness <- best.pn.fitness[nh[[j]]][best.nh.ind]
      if (candidate.fitness > best.nh.fitness[j])
      {
        best.nh[j, ] <- best.pn[nh[[j]], , drop = FALSE][best.nh.ind, ]
        best.nh.fitness[j] <- candidate.fitness
      }
    }
    
    # Update the best particle
    
    best.old.fitness <- best.final.fitness
    best.final.ind <- which.max(best.pn.fitness)
    best.final <- best.pn[best.final.ind, ]
    best.final.fitness <- best.pn.fitness[best.final.ind]
    
    # Store fitness of the best particle
    
    fitness.history <- c(fitness.history, best.final.fitness)

    # Change the topology if need along with best
    # known solution in the neighbourhood
    
    if (nh.adaptive & (best.old.fitness >= best.final.fitness))
    {
      # Update the neighbourhood
      
      nh <- pso.nh(pop.n = pop.n,
                   method = nh.method,
                   par = nh.par,
                   iter = i)
      
      # Update best known solutions in the neighbourhood
      
      for (j in pop.ind)
      {
        best.nh.ind <- which.max(best.pn.fitness[nh[[j]]])
        best.nh[j, ] <- best.pn[nh[[j]], , drop = FALSE][best.nh.ind, ]
        best.nh.fitness[j] <- best.pn.fitness[nh[[j]]][best.nh.ind]
      }
    }
    
    if (info)
    {
      cat(paste0("Iter = ", i, ", ",
                 "mean = ", is.max.val * signif(mean(fitness), 5), ", ",
                 "median = ", is.max.val * signif(median(fitness), 5), ", ",
                 "best = ", is.max.val * signif(best.final.fitness, 5), "\n"))
    }

    # Return the results if need
    
    if (i == maxiter)
    {
      counts[1] <- counts[1] + i * pop.n
      return.list <- list(population = population,
                          par = is.max.val * best.final,
                          value = is.max.val * best.final.fitness,
                          counts = counts,
                          fitness.history = is.max.val * fitness.history,
                          is.max = is.max,
                          iter = i)
      class(return.list) <- "pso"
      return(return.list)
    }
    
  }
  
  return(NULL)
}

#' Print method for "pso" object
#' @param x Object of class "pso"
#' @param ... further arguments (currently ignored)
#' @returns This function does not return anything.
print.pso <- function(x, ...)
{
  cat("We have optimized, optimized and finally optimized! \n")
  cat("---\n")
  cat(paste0(ifelse(x$is.max, "The largest ", "The smallest "),
             "value of the fitness function found is ", 
             signif(x$value, 5), " at point: \n"))
  cat(signif(x$par, 5))
  cat("\n")
  cat("---\n")
  cat(paste0("There was ", x$counts[1], " calls to function and ",
             x$counts[2], " calls to gradient. \n"))
  cat("---\n")
  cat("This object is a list containing the following elements:\n")
  cat(names(x))
}

#' Summarizing pso Fits
#' @param object Object of class "pso"
#' @param ... further arguments (currently ignored)
#' @return This function returns the same list as \code{\link[gena]{pso}} 
#' function changing its class to "summary.pso".
summary.pso <- function(object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  class(object) <- "summary.pso"
  return(object)
}

#' Summary for "pso" object
#' @param x Object of class "pso"
#' @param ... further arguments (currently ignored)
#' @returns This function returns \code{x} input argument.
print.summary.pso <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  print.pso(x)
  return(x)
}

#' Plot best found fitnesses during genetic algorithm
#' @param x Object of class "pso"
#' @param y this parameter currently ignored
#' @param ... further arguments (currently ignored)
#' @returns This function does not return anything.
plot.pso <- function (x, y = NULL, ...) 
{
  plot(x = 1:x$iter, y = x$fitness.history, 
       main = "Historical values of fitnesses", 
       xlab = "Iteration", 
       ylab = paste0(ifelse(x$is.max, "The largest ", "The smallest "),
                     "fitness value"))
}