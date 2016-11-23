#' @name greta-samplers
#' @title sample model variables
#' @description After defining a greta model in R, draw samples of the random
#'   variables of interest
#' @param ... nodes to sample values from, probably parameters of a
#'   model. Observed nodes cannot be sampled from.
#' @param method the method used to sample values. Currently only \code{hmc} is
#'   implemented
#' @param n_samples the number of samples to draw (after any warm-up, but before
#'   thinning)
#' @param thin the thinning rate; every \code{thin} samples is retained, the
#'   rest are discarded
#' @param warmup the number of samples to spend warming up the sampler. During
#'   this phase the sampler moves toward the highest density area and may tune
#'   sampler hyperparameters.
#' @param verbose whether to print progress information to the console
#' @param control an optional named list of hyperparameters and options to
#'   control behaviour of the sampler
#' @export
#' @examples
#' # define a simple model
#' mu = free()
#' sigma = lognormal(1, 0.1)
#' x = observed(rnorm(10))
#' x %~% normal(mu, sigma)
#'
#' draws <- samples(mu, sigma,
#'                 n_samples = 100,
#'                 warmup = 10)
samples <- function (...,
                    method = c('hmc', 'nuts'),
                    n_samples = 1000,
                    thin = 1,
                    warmup = 100,
                    verbose = TRUE,
                    control = list()) {

  method <- match.arg(method)

  # nodes required
  target_nodes <- list(...)

  # find variable names to label samples
  names <- substitute(list(...))[-1]
  names <- vapply(names, deparse, '')
  names(target_nodes) <- names

  # check they're not data nodes, provide a useful error message if they are
  type <- vapply(target_nodes, member, 'type', FUN.VALUE = '')
  bad <- type == 'data'
  if (any(bad)) {
    is_are <- ifelse(sum(bad) == 1, 'is an observed node', 'are observed nodes')
    bad_nodes <- paste(names[bad], collapse = ', ')
    msg <- sprintf('%s %s, observed nodes cannot be sampled',
                   bad_nodes,
                   is_are)
    stop (msg)
  }

  # get the dag containing the target nodes
  dag <- dag_class$new(target_nodes)

  if (verbose)
    message('compiling model')

  # define the TF graph
  dag$define_tf()

  # random starting locations
  init <- dag$example_parameters()
  init[] <- rnorm(length(init))

  # get default control options
  con <- switch(method,
                hmc = list(Lmin = 10,
                           Lmax = 20,
                           epsilon = 0.005),
                nuts = list())

  # update them with user overrides
  con[names(control)] <- control

  # fetch the algorithm
  method <- switch(method,
                   hmc = hmc,
                   nuts = nuts)

  # if warmup is required, do that now and update init
  if (warmup > 0) {

    if (verbose)
      message('warming up')

    # run it
    warmup_draws <- method(dag = dag,
                           init = init,
                           n_samples = warmup,
                           thin = thin,
                           verbose = verbose,
                           guess_epsilon = FALSE,
                           control = con)

    # use the last draw of the full parameter vector as the init
    init <- attr(warmup_draws, 'last_x')
    con$epsilon <- attr(warmup_draws, 'epsilon')

    if (verbose)
      message('sampling')

  }

  # run the sampler
  draws <- method(dag = dag,
                  init = init,
                  n_samples = n_samples,
                  thin = thin,
                  verbose = verbose,
                  guess_epsilon = FALSE,
                  control = con)

  # coerce to data.frame, but keep the sample density
  draws_df <- data.frame(draws)
  attr(draws_df, 'density') <- attr(draws, 'density')
  attr(draws_df, 'last_x') <- attr(draws, 'last_x')
  draws_df

}

# run NUTS HMC sampler
nuts <- function (dag, init, n_samples, thin, verbose, control = list()) {
  stop ('not yet implemented')
}


hmc <- function (dag,
                 init,
                 n_samples,
                 thin,
                 verbose,
                 guess_epsilon = FALSE,
                 control = list(Lmin = 10,
                                Lmax = 20,
                                epsilon = 0.005)) {

  # unpack options
  Lmin <- control$Lmin
  Lmax <- control$Lmax
  epsilon <- control$epsilon

  # set initial location, log joint density and gradients
  x <- init
  dag$send_parameters(x)
  grad <- dag$gradients()
  logprob <- dag$log_density()

  if (guess_epsilon)
    epsilon <- guess_epsilon(dag, init, logprob, grad)


  # set up trace store (grab values of target variables from graph to get
  # dimension and names)
  init_trace <- dag$trace_values()
  n_target <- length(init_trace)
  trace <- matrix(NA,
                  nrow = n_samples %/% thin,
                  ncol = n_target)
  colnames(trace) <- names(init_trace)

  # set up log joint density store
  ljd <- rep(NA, n_samples)

  # get free parameter dimension
  npar <- length(x)

  accept_count <- 0

  # set up progress bar
  if (verbose)
    pb <- txtProgressBar(max = n_samples, style = 3)

  # loop through iterations
  for (i in 1:n_samples) {

    # copy old state
    x_old <- x
    logprob_old <- logprob
    grad_old <- grad
    p_old <- rnorm(npar)

    # start leapfrog steps
    reject <- FALSE
    p <- p_old + 0.5 * epsilon * grad
    n_steps <- sample(Lmin:Lmax, 1)

    for (l in seq_len(n_steps)) {

      # do a step
      new_step <- leapfrog(list(x = x_old, p = p), dag, epsilon)

      # if it went awry, quit now
      if (is.null(new_step)) {

        reject <- TRUE
        break()

      } else {

        # otherwise, update the location an momemntum and continue
        x <- new_step$x
        p <- new_step$p

      }

      p <- p + epsilon * grad

    }

    p <- p - 0.5 * epsilon * grad

    # if the step was bad, reject it out of hand
    if (reject) {

      if (verbose)
        message ('proposal rejected due to numerical instability')

      x <- x_old
      logprob <- logprob_old
      grad <- grad_old

    } else {

      # otherwise do the Metropolis accept/reject step

      # acceptance ratio
      log_accept_ratio  <- log_accept_prob(logprob, p, logprob_old, p_old)
      log_u = log(runif(1))

      if (log_u < log_accept_ratio) {

        # on acceptance, iterate the counter and leave the parameters in the dag
        # to be put in the trace
        accept_count <- accept_count + 1

      } else {

        # on rejection, reset all the parameters and push old parameters to the
        # graph for the trace
        x <- x_old
        logprob <- logprob_old
        grad <- grad_old

      }

    }

    # either way, store density and location of target parameters straight from the graph
    # reset dag parameters for extracting the trace
    if (i %% thin == 0) {
      dag$send_parameters(x)
      trace[i / thin, ] <- dag$trace_values()
      ljd[i / thin] <- dag$log_density()
    }

    if (verbose)
      setTxtProgressBar(pb, i)

  }

  if (verbose)
    close(pb)

  attr(trace, 'density') <- -ljd
  attr(trace, 'last_x') <- x
  attr(trace, 'epsilon') <- epsilon
  trace

}

# do a single leapfrog step and return x and p in a list
# if the step went bad, return a NULL instead
leapfrog <- function(list, dag, epsilon) {
  x <- list$x
  p <- list$p

  # move to new location
  x <- x + epsilon * p

  # evaluate
  dag$send_parameters(x)
  logprob <- dag$log_density()
  grad <- dag$gradients()

  # check gradients and density are finite
  if (any(!is.finite(grad)) | !is.finite(logprob)) {

    return (NULL)

  } else {

    # otherwise update and return
    p <- p + epsilon * grad
    return (list(x = x, p = p, logprob = logprob, grad = grad))

  }
}

# log acceptance probability
log_accept_prob <- function(logprob, p, logprob_old, p_old) {

  # inner products
  p_prod <- (t(p) %*% p)[1, 1]
  p_prod_old <- (t(p_old) %*% p_old)[1, 1]

  # log ratio
  logprob - 0.5 * p_prod - logprob_old + 0.5 * p_prod_old
}

# make a reasonable guess at epsilon from the starting location
guess_epsilon <- function(dag, x_old, logprob_old, grad_old) {

  # initial conditions
  epsilon <- 1
  p_old <- rnorm(length(x_old))

  # do a leapfrog step
  new_step <- leapfrog(list(x = x_old, p = p_old),
                 dag, epsilon)

  # see how far to step in that direction
  k <- 1
  while(!is.finite(new_step$logprob) | !all(is.finite(new_step$grad))) {
    k <- k * 0.5
    new_step <- leapfrog(list(x = x, p = p_old), dag, epsilon * k)
  }

  # shrink epsilon to halfway there
  epsilon <- epsilon * 0.5 * k

  # get required elements
  logprob <- new_step$logprob
  p <- new_step$p

  # now, get acceptance probability and direction toward 0.5 acceptance
  accept_prob <- exp(log_accept_prob(logprob, p, logprob_old, p_old))
  a <- ifelse(accept_prob > 0.5, 1, -1)

  while ((accept_prob ^ a)  > (2 ^ -a) ) {
    epsilon <- epsilon * 2 ^ a
    new_step <- leapfrog(list(x = x_old, p = p_old), dag, epsilon)
    log_accept <- log_accept_prob(new_step$logprob, new_step$p,
                                  logprob_old, p_old)
    accept_prob <- exp(log_accept)
  }

  return (epsilon)

}

#
# # leapfrog once
# x <- x_old + epsilon * p
#
# # send parameters
# dag$send_parameters(x)
# logprob <- dag$log_density()
# grad <- dag$gradients()
#
# # check gradients are finite
# if (any(!is.finite(grad))) {
#   reject <- TRUE
#   break()
# }
#
# p <- p + epsilon * grad
#
#
#
#
# def find_reasonable_epsilon(theta0, grad0, logp0, f):
#   """ Heuristic for choosing an initial value of epsilon """
# epsilon = 1.
# r0 = np.random.normal(0., 1., len(theta0))
#
# # Figure out what direction we should be moving epsilon.
# _, rprime, gradprime, logpprime = leapfrog(theta0, r0, grad0, epsilon, f)
# # brutal! This trick make sure the step is not huge leading to infinite
# # values of the likelihood. This could also help to make sure theta stays
# # within the prior domain (if any)
# k = 1.
# while np.isinf(logpprime) or np.isinf(gradprime).any():
#   k *= 0.5
# _, rprime, _, logpprime = leapfrog(theta0, r0, grad0, epsilon * k, f)
#
# epsilon = 0.5 * k * epsilon
#
# acceptprob = np.exp(logpprime - logp0 - 0.5 * (np.dot(rprime, rprime.T) - np.dot(r0, r0.T)))
#
# a = 2. * float((acceptprob > 0.5)) - 1.
# # Keep moving epsilon in that direction until acceptprob crosses 0.5.
# while ( (acceptprob ** a) > (2. ** (-a))):
#   epsilon = epsilon * (2. ** a)
# _, rprime, _, logpprime = leapfrog(theta0, r0, grad0, epsilon, f)
# acceptprob = np.exp(logpprime - logp0 - 0.5 * ( np.dot(rprime, rprime.T) - np.dot(r0, r0.T)))
#
# print "find_reasonable_epsilon=", epsilon
