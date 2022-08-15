#' Numeric Differentiation
#' @name genaDiff
#' @description Numeric estimation of the gradient and Hessian.
#' @param fn function for which gradient or Hessian should be calculated.
#' @param gr gradient function of \code{fn}.
#' @param par point (parameters' value) at which \code{fn} should be 
#' differentiated.
#' @param eps numeric vector representing increment of the \code{par}. 
#' So \code{eps[i]} represents increment of \code{par[i]}. If \code{eps} is
#' a constant then all increments are the same.
#' @param method numeric differentiation method: "central-difference" or
#' "forward-difference".
#' @param fn.args list containing arguments of \code{fn} except \code{par}.
#' @param gr.args list containing arguments of \code{gr} except \code{par}.
#' @details It is possible to substantially improve numeric Hessian accuracy 
#' by using analytical gradient \code{gr}. If both \code{fn} and \code{gr}
#' are provided then only \code{gr} will be used. If only \code{fn} is provided 
#' for \code{gena.hessian} then \code{eps} will be transformed to 
#' \code{sqrt(eps)} for numeric stability purposes.
#'
#' @return Function \code{gena.grad} returns a vector that is a gradient of 
#' \code{fn} at point \code{par} calculated via \code{method} numeric 
#' differentiation approach using increment \code{eps}.
#' 
#' Function \code{gena.hessian} returns a matrix that is a Hessian of 
#' \code{fn} at point \code{par}.
#' @examples 
#' # Consider the following function
#' fn <- function(par, a = 1, b = 2)
#' {
#'   val <- par[1] * par[2] - a * par[1] ^ 2 - b * par[2] ^ 2
#' }
#' 
#' # Calculate the gradient at point (2, 5) respect to 'par' 
#' # when 'a = 1' and 'b = 1'
#' par <- c(2, 5)
#' fn.args = list(a = 1, b = 1)
#' gena.grad(fn = fn, par = par, fn.args = fn.args)
#' 
#' # Calculate Hessian at the same point
#' gena.hessian(fn = fn, par = par, fn.args = fn.args)
#' 
#' # Repeat calculation of the Hessian using analytical gradient
#' gr <- function(par, a = 1, b = 2)
#' {
#'   val <- c(par[2] - 2 * a * par[1],
#'            par[1] - 2 * b * par[2])
#' }
#' gena.hessian(gr = gr, par = par, gr.args = fn.args)
#' 
#' @rdname genaDiff
#' @export
gena.grad <- function(fn, par, 
                      eps = sqrt(.Machine$double.eps) * abs(par),
                      method = "central-difference",
                      fn.args = NULL)
{
  n_par <- length(par)
  
  if(length(eps) == 1)
  {
    eps <- rep(eps, n_par)
  }
    
  par_plus <- par
  par_minus <- par
    
  val <- rep(0, n_par)

  for(i in 1:n_par)
  {
    if (method == "central-difference")
    {
      par_plus[i] <- par[i] + eps[i]
    }

    par_minus[i] <- par[i] - eps[i]
      
    fn.args$par <- par_plus
    fn_plus <- do.call(fn, fn.args)
   
    fn.args$par <- par_minus
    fn_minus <- do.call(fn, fn.args)
      
    par_plus <- par
    par_minus <- par
       
    if (method == "central-difference")
    {
      val[i] <- (fn_plus - fn_minus) / (2 * eps[i])
    }
    if (method == "forward-difference")
    {
      val[i] <- (fn_plus - fn_minus) / eps[i]
    }
  }
  
  return(val)
}

#' @rdname genaDiff
#' @export
gena.hessian <- function(fn = NULL, gr = NULL, 
                         par, eps = sqrt(.Machine$double.eps) * abs(par), 
                         fn.args = NULL, gr.args = NULL)
{
  n_par <- length(par)
  
  if (is.null(gr))
  {
    eps <- sqrt(eps)
  }
  
  if(length(eps) == 1)
  {
    eps <- rep(eps, n_par)
  }
  
  par_plus <- par
  par_minus <- par
  
  val <- matrix(0, n_par, n_par)
  
  # With analytical gradient
  if (!is.null(gr))
  {
    for(i in 1:n_par)
    {
      par_plus[i] <- par[i] + eps[i]
      par_minus[i] <- par[i] - eps[i]
        
      gr.args$par <- par_plus
      gr_plus <- do.call(gr, gr.args)
        
      gr.args$par <- par_minus
      gr_minus <- do.call(gr, gr.args)
        
      par_plus <- par
      par_minus <- par
        
      val[i, ] <- (gr_plus - gr_minus) / (2 * eps[i])
    }
    
    for (i in 1:n_par)
    {
      for (j in 1:i)
      {
        if (i != j)
        {
          val[i, j] <- (val[i, j] + val[j, i]) / 2
          val[j, i] <- val[i, j]
        }
      }
    }
  }
  # Without analytical gradient
  else
  {
    fn.args$par <- par
    fn_0 <- do.call(fn, fn.args)
    
    fn_plus <- rep(NA, n_par)
    for (i in 1:n_par)
    {
      par_plus[i] <- par[i] + eps[i]
      fn.args$par <- par_plus
      fn_plus[i] <- do.call(fn, fn.args)
      par_plus[i] <- par[i]
    }
    
    par_ij <- par
    for (i in 1:n_par)
    {
      par_ij[i] <- par[i] + eps[i]
      for (j in 1:i)
      {
        par_ij[j] <- par_ij[j] + eps[j]
        fn.args$par <- par_ij
        
        fn_ij <- do.call(fn, fn.args)
        val[i, j] <- (fn_ij - fn_plus[i] - fn_plus[j] + fn_0) / 
                     (eps[i] * eps[j])
        val[j, i] <- val[i, j]
          
        par_ij[j] <- par[j]
      }
      par_ij[i] <- par[i]
    }
  }
  
  return(val)
}

gena.lineSearch <- function(fn, gr = NULL,
                            fn.val, gr.val,
                            par,
                            dir, lr = 1,
                            c1 = 1e-4, c2 = 0.9,
                            mult = 0.7,
                            maxiter = 100,
                            fn.args = list(), gr.args = list())
{
  # Variable to control whether learning rate is increasing
  is_decrease = TRUE
  
  # Learning rate values
  lr.vec <- rep(NA, maxiter)
  lr.vec[1] <- lr
  
  # Parameters values
  par.mat <- matrix(NA, ncol = length(par), nrow = maxiter)
  
  # Function and gradient calls
  fn.calls <- 0
  gr.calls <- 0
  
  # Non-monotone line search
  # M <- min(length(fn.val), 5)
  # fn.val <- mean(tail(fn.val, M))
  
  # Wolfe's conditions
  armajilo <- FALSE
  curvature <- FALSE
  
  for (i in 1:maxiter)
  {
    # Calculate new point and function
    # value at this point
    par.mat[i, ] <- par + lr.vec[i] * dir
    fn.args$par <- par.mat[i, ]
    fn.val.new <- do.call(fn, fn.args)
    fn.calls <- fn.calls + 1
    
    # Check Armijo condition
    gr.dir <- sum(gr.val * dir)
    armajilo <- fn.val.new <= (fn.val + c1 * lr.vec[i] * gr.dir)
    if (armajilo)
    {
      # Check strong Curvature condition
      gr.args$par <- par.mat[i, ]
      gr.val.new <- NULL
      if(!is.null(gr))
      {
        gr.val.new <- matrix(do.call(gr, gr.args), ncol = 1)
      }
      else
      {
        gr.val.new <- gena.grad(fn = fn, par = par.mat[i, ], 
                                fn.args = fn.args)
      }
      gr.calls <- gr.calls + 1
      gr.div.new <- sum(gr.val.new * dir)
      curvature <- gr.div.new <= (c2 * gr.dir)
      if (curvature)
      {
        if (is_decrease)
        {
          if (i == 1)
          {
            is_decrease <- FALSE
          }
          else
          {
            return.list <- list(lr = lr.vec[i],
                                par = par.mat[i, ],
                                fail = FALSE,
                                fn.calls = fn.calls,
                                gr.calls = gr.calls)
            return(return.list)
          }
        }
      }
    }
    if(!(armajilo & curvature) & !is_decrease)
    {
      return.list <- list(lr = lr.vec[i - 1],
                          par = par.mat[i - 1, ],
                          fail = FALSE,
                          fn.calls = fn.calls,
                          gr.calls = gr.calls)
      return(return.list)
    }
    if (is_decrease)
    {
      lr.vec[i + 1] <- lr.vec[i] * mult
    }
    else
    {
      lr.vec[i + 1] <- lr.vec[i] / mult
    }
  }
  
  # If line search is not successful
  return.list <- list(lr = 0,
                      par = par,
                      fail = TRUE,
                      fn.calls = fn.calls,
                      gr.calls = gr.calls)

  return(return.list)
}

gena.local <- function(par,                                   # initial point
                       fn,                                    # function to optimize
                       gr = NULL,                             # gradient of the function
                       maxiter = 1e+6,                        # the number of iterations
                       reltol = 1e-10,                        # relative tolerance
                       abstol = 1e-10,                        # absolute tolerance
                       lr = 1024,                             # learning rate
                       type = "BFGS",                         # optimization type
                       ...)                                   # additional arguments for fn and gr
{
  dots <- list(...)
  iter <- NULL
  dir <- NULL
  
  fn.calls <- 0
  gr.calls <- 0
  
  reltol.cond <- FALSE
  maxiter.cond <- FALSE
  abstol.cond <- FALSE

  # Get the number of estimated parameters
  n_x <- length(par)

  # Initialize the matrix to store the
  # points for each iteration
  x <- matrix(NA,
              nrow = maxiter,
              ncol = length(par))
  x[1, ] <- par

  # Initialize variable to store function 
  # value and gradient information
  fn.val <- rep(NA, n_x)
  gr.val <- matrix(NA, nrow = maxiter, ncol = n_x)
  
  # Initialize variable to store Hessians and their inverse
  hessian <- array(dim = c(maxiter + 1, n_x, n_x))
  hessian.inv <- hessian
  hessian[1, , ] <- diag(rep(1, n_x)) #solve(gena.hessian(fn = fn, gr = gr, par = par))
  hessian.inv[1, , ] <- hessian[1, , ]
  
  # Estimate gradient and function value
  # at the initial point
  dots$par <- par
  fn.val[1] <- do.call(fn, dots)
  if(is.null(gr))
  {
    gr.val[1, ] <- gena.grad(fn, x[1, ], fn.args = dots)
  } else {
    dots$par <- x[1, ]
    gr.val[1, ] <- do.call(gr, dots)
  }
  gr.calls <- gr.calls + 1

  # Start the algorithm
  for (i in 1:maxiter)
  {
    # Determine the direction
    dir <- -hessian.inv[i, , ] %*% matrix(gr.val[i, ], ncol = 1)

    # Perform the line search
    ls.result <- gena.lineSearch(fn = fn, gr = gr,
                                 fn.val = fn.val[i], 
                                 gr.val = gr.val[i, ],
                                 par = x[i, ],
                                 dir = dir, lr = lr,
                                 fn.args = dots,
                                 gr.args = dots)
    # if (ls.result$fail)
    # {
    #   print(123)
    #   ls.result <- gena.lineSearch(fn = fn, gr = gr,
    #                                fn.val = fn.val[i], 
    #                                gr.val = gr.val[i, ],
    #                                par = x[i, ],
    #                                dir = -dir, lr = lr,
    #                                fn.args = dots,
    #                                gr.args = dots)
    #   hessian.inv[i, , ] <- diag(rep(1, n_x))
    # }
    x[i + 1, ] <- ls.result$par
    lr <- ls.result$lr
    fn.calls <- fn.calls + ls.result$fn.calls
    gr.calls <- gr.calls + ls.result$gr.calls

    # Estimate the function value at new point
    dots$par <- x[i + 1, ]
    fn.val[i + 1] <- do.call(fn, dots)
    fn.calls <- fn.calls + 1
    print(c(fn.val[i], fn.val[i + 1] , lr))
    
    # Estimate the gradient and new point
    if(is.null(gr))
    {
      gr.val[i + 1, ] <- gena.grad(fn, x[i + 1, ], fn.args = dots)
    } else {
      dots$par <- x[i + 1, ]
      gr.val[i + 1, ] <- do.call(gr, dots)
    }

    # Calculate inverse of the Hessian
    y <- matrix(gr.val[i + 1, ] - gr.val[i, ], ncol = 1)
    s <- matrix(lr * dir, ncol = 1)
    
    hessian.inv[i + 1, , ] <- hessian.inv[i, , ] +
                              (s %*% t(s)) * 
                              (as.vector(t(s) %*% y + 
                                         t(y) %*% hessian.inv[i, , ] %*% y) /
                              as.vector((t(s) %*% y) ^ 2)) -
                              (hessian.inv[i, , ] %*% y %*% t(s) + 
                               s %*% t(y) %*% hessian.inv[i, , ]) / 
                              as.vector(t(s) %*% y)
    
    # I <- diag(rep(1, n_x))
    # hessian.inv[i + 1, ,] <- (I - s %*% t(y) / as.vector(t(y) %*% s)) %*% 
    #                          hessian.inv[i, , ] %*%
    #                          (I - y %*% t(s) / as.vector(t(y) %*% s)) +
    #                          (s %*% t(s)) / as.vector(t(y) %*% s)

    # Print the information
    cat(paste0("iter-", i, ": ", fn.val[i], "\n"))

    # Check whether termination conditions
    # are satisfied
    abstol.cond <- abs(fn.val[i] - fn.val[i + 1]) < abstol
    if (fn.val[i] != 0)
    {
      reltol.cond <- all(abs((fn.val[i] - fn.val[i + 1]) /
                              fn.val[i]) < reltol)
    }
    else
    {
      reltol.cond <- FALSE
    }
    maxiter.cond <- i == maxiter
    if (abstol.cond | reltol.cond | maxiter.cond)
    {
      iter <- i
      break
    }
  }

  # Get the final point
  return_list <- list(points = x[1:iter],
                      solution = x[iter, ],
                      value = fn.val[iter],
                      fn_grad = gr.val[iter, ],
                      hessian = solve(hessian.inv[iter, , ]),
                      hessian.inv = hessian.inv[iter, , ],
                      fn.calls = fn.calls,
                      gr.calls = gr.calls,
                      dir = dir)

  return(return_list)
}