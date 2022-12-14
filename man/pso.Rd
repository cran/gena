% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pso.R
\name{pso}
\alias{pso}
\title{Particle Swarm Optimization}
\usage{
pso(
  fn,
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
  velocity.par = list(w = 1/(2 * log(2)), c1 = 0.5 + log(2), c2 = 0.5 + log(2)),
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
  ...
)
}
\arguments{
\item{fn}{function to be maximized i.e. fitness function.}

\item{gr}{gradient of the \code{fn}.}

\item{lower}{lower bound of the search space.}

\item{upper}{upper bound of the search space.}

\item{pop.n}{integer representing the size of the population.}

\item{pop.initial}{numeric matrix which rows are particles to be
included into the initial population. Numeric vector will be coerced to
single row matrix.}

\item{pop.method}{the algorithm to be applied for a creation of 
the initial population. See 'Details' for additional information.}

\item{nh.method}{string representing the method (topology) to be used for
the creation of neighbourhoods. See 'Details' for additional information.}

\item{nh.par}{parameters of the topology algorithm.}

\item{nh.adaptive}{logical; if \code{TRUE} (default) then neighbourhoods
change every time when the best known (to the swarm) fitnesses value have
not increased. Neighbourhoods are updated according to the topology
defined via \code{nh.method} argument.}

\item{velocity.method}{string representing the method to be used for
the update of velocities.}

\item{velocity.par}{parameters of the velocity formula.}

\item{hybrid.method}{hybrids selection algorithm i.e. mechanism 
determining which particles should be subject to local optimization.
See 'Details' for additional information.}

\item{hybrid.par}{parameters of the hybridization algorithm.}

\item{hybrid.prob}{probability of generating the hybrids each iteration.}

\item{hybrid.opt.par}{parameters of the local optimization function
to be used for hybridization algorithm (including \code{fn} and \code{gr}).}

\item{hybrid.n}{number of hybrids that appear if hybridization
should take place during the iteration.}

\item{constr.method}{the algorithm to be applied for imposing constraints
on the particles. See 'Details' for additional information.}

\item{constr.par}{parameters of the constraint algorithm.}

\item{random.order}{logical; if \code{TRUE} (default) then particles
related routine will be implemented in a random order.}

\item{maxiter}{maximum number of iterations of the algorithm.}

\item{is.max}{logical; if \code{TRUE} (default) then fitness function
will be maximized. Otherwise it will be minimized.}

\item{info}{logical; if \code{TRUE} (default) then some optimization related 
information will be printed each iteration.}

\item{...}{additional parameters to be passed to 
\code{fn} and \code{gr} functions.}
}
\value{
This function returns an object of class \code{pso} that is a list
containing the following elements:
\itemize{
\item \code{par} - particle (solution) with the highest fitness
(objective function) value.
\item \code{value} - value of \code{fn} at \code{par}.
\item \code{population} - matrix of particles (solutions) of the 
last iteration of the algorithm.
\item \code{counts} - a two-element integer vector giving the number of
calls to \code{fn} and \code{gr} respectively.
\item \code{is.max} - identical to \code{is.max} input argument.
\item \code{fitness.history} - vector which i-th element is fitness
of the best particle in i-th iteration.
\item \code{iter} - last iteration number.
}
}
\description{
This function allows to use particle swarm algorithm for
numeric global optimization of real-valued functions.
}
\details{
Default arguments have been set in accordance with SPSO 2011
algorithm proposed by M. Clerc (2012).

To find information on particular methods available via
\code{pop.method}, \code{nh.method}, \code{velocity.method},
\code{hybrid.method} and \code{constr.method}
arguments please see 'Details' section of 
\code{\link[gena]{gena.population}}, \code{\link[gena]{pso.nh}},
\code{\link[gena]{pso.velocity}}, \code{\link[gena]{gena.hybrid}} 
and \code{\link[gena]{gena.constr}} correspondingly.

It is possible to provide manually implemented functions for population
initialization, neighbourhoods creation, velocity updated, hybridization
and constraints in a similar way as for \code{\link[gena]{gena}}.

By default function does not impose any constraints upon the parameters.
If \code{constr.method = "bounds"} then \code{lower} and \code{upper} 
constraints will be imposed. Lower bounds should be strictly smaller
then upper bounds.

Currently the only available termination condition is \code{maxiter}. We
are going to provide some additional termination conditions during
future updates.

Infinite values in \code{lower} and \code{upper} are substituted with
\code{-(.Machine$double.xmax * 0.9)} and \code{.Machine$double.xmax * 0.9}
correspondingly.

By default if \code{gr} is provided then BFGS algorithm will be used inside
\code{\link[stats]{optim}} during hybridization.
Otherwise \code{Nelder-Mead} will be used. 
Manual values for \code{\link[stats]{optim}} arguments may be provided 
(as a list) through \code{hybrid.opt.par} argument.

For more information on particle swarm optimization 
please see M. Clerc (2012).
}
\examples{
## Consider Ackley function
\donttest{
fn <- function(par, a = 20, b = 0.2)
{
  val <- a * exp(-b * sqrt(0.5 * (par[1] ^ 2 + par[2] ^ 2))) +
         exp(0.5 * (cos(2 * pi * par[1]) + cos(2 * pi * par[2]))) -
         exp(1) - a
  return(val)
}

# Maximize this function using particle swarm algorithm

set.seed(123)
lower <- c(-5, -100)
upper <- c(100, 5)
opt <- pso(fn = fn, 
           lower = lower, upper = upper,
           a = 20, b = 0.2)
print(opt$par)
}

## Consider Bukin function number 6

fn <- function(x, a = 20, b = 0.2)
{
  val <- 100 * sqrt(abs(x[2] - 0.01 * x[1] ^ 2)) + 0.01 * abs(x[1] + 10)
  return(val)
}

# Minimize this function using initially provided
# position for one of the particles
set.seed(777)
lower <- c(-15, -3)
upper <- c(-5, 3)
opt <- pso(fn = fn, 
           pop.init = c(8, 2),
           lower = lower, upper = upper,
           is.max = FALSE)
print(opt$par)

}
\references{
M. Clerc (2012). 
Standard Particle Swarm Optimisation.
\emph{HAL archieve}.
}
