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

  # fetch the algorithm
  method <- switch(method,
                   hmc = hmc,
                   nuts = nuts6)

  # get user and default control options
  control_user <- control
  control <- eval(as.list(args(method))$control)

  # update existing control with user overrides
  control[names(control_user)] <- control_user

  # if the user didn't specify an epsilon, take a guess
  if (! 'epsilon' %in% names(control_user))
    control$epsilon <- guess_epsilon(dag, init)

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
                           tune_epsilon = TRUE,
                           control = control)

    # use the last draw of the full parameter vector as the init, and grab epsilon
    init <- attr(warmup_draws, 'last_x')
    control <- attr(warmup_draws, 'control')

    if (verbose)
      message('sampling')

  }

  # run the sampler
  draws <- method(dag = dag,
                  init = init,
                  n_samples = n_samples,
                  thin = thin,
                  verbose = verbose,
                  tune_epsilon = FALSE,
                  control = control)

  # coerce to data.frame, but keep the sample density
  draws_df <- data.frame(draws)
  attr(draws_df, 'density') <- attr(draws, 'density')
  attr(draws_df, 'last_x') <- attr(draws, 'last_x')
  attr(draws_df, 'control') <- attr(draws, 'control')
  draws_df

}

# ~~~~~~~~~~~~~
# functions used by samples() (may be used by samplers too)

# make a reasonable guess at epsilon from the starting location
guess_epsilon <- function(dag, x_old) {

  dag$send_parameters(x_old)
  grad_old <- dag$gradients()
  logprob_old <- dag$log_density()

  # initial conditions (start with random inertia)
  epsilon <- 1
  p_old <- rnorm(length(x_old))

  # do a leapfrog step
  new_step <- leapfrog(list(x = x_old, p = p_old),
                 dag, epsilon)

  # if the initial epsilon leads to numerical errors, row back until it's stable
  k <- 1
  while(!is.finite(new_step$logprob) | !all(is.finite(new_step$grad))) {
    k <- k * 0.5
    new_step <- leapfrog(list(x = x, p = p_old), dag, epsilon * k)
  }

  # shrink epsilon to halfway there
  epsilon <- epsilon * 0.5 * k

  # now, get acceptance probability and direction toward 0.5 acceptance
  accept_prob <- accept_prob_fun(new_step$logprob, new_step$p, logprob_old, p_old)
  a <- ifelse(accept_prob > 0.5, 1, -1)

  while ((accept_prob ^ a)  > (2 ^ -a) ) {
    epsilon <- epsilon * 2 ^ a
    new_step <- leapfrog(list(x = x_old, p = p_old), dag, epsilon)
    accept_prob <- accept_prob_fun(new_step$logprob, new_step$p,
                                   logprob_old, p_old)
  }

  return (epsilon)

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
accept_prob_fun <- function(logprob, p, logprob_old, p_old, log = FALSE) {

  # inner products
  p_prod <- (t(p) %*% p)[1, 1]
  p_prod_old <- (t(p_old) %*% p_old)[1, 1]

  # log ratio
  res <- (logprob - 0.5 * p_prod) - (logprob_old - 0.5 * p_prod_old)

  # unlog if required
  if (!log)
    res <- exp(res)

  res
}
