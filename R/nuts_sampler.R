# NUTS6 sampler, inspired by R code by Cole Monnahan
# (https://github.com/colemonnahan/gradmcmc/blob/master/algorithms/mcmc.R) and
# Python code by Morgan Fouesneau
# (https://github.com/mfouesneau/NUTS/blob/master/nuts.py)

nuts6 <- function (dag,
                   init,
                   n_samples,
                   thin,
                   verbose,
                   tune_epsilon = FALSE,
                   estimate_mass_matrix = FALSE,
                   control = list(max_doublings = 4,
                                  epsilon = 0.05,
                                  mass_cholesky = NULL,
                                  init_proposal = NULL,
                                  target_acceptance = 0.5,
                                  accept_group = 100,
                                  t0 = 10,
                                  gamma = 0.05,
                                  kappa = 0.75,
                                  init_prop = 0.15,
                                  term_prop = 0.3)) {

  # unpack control options
  unpack(control)

  # store epsilon as epsilon_mean, so that epsilon can be jittered on each
  # iteration
  epsilon_mean <- epsilon

  # set initial location, log joint density and gradients
  x <- init

  # initial epsilon tuning parameters
  if (tune_epsilon) {

    # keep track of tuner progress
    epsilon_trace <- c(epsilon, rep(NA, n_samples))
    epsilon_bar_trace <- epsilon_trace
    Hbar_trace <- c(0, rep(NA, length = n_samples))
    nstep_trace <- rep(NA, n_samples)

    mu <- log(10 * epsilon)

  }

  # if we're estimating the mass matrix
  if (estimate_mass_matrix) {

    # set up to trace the second half of the values
    trace_x <- matrix(NA,
                      nrow = floor(n_samples / 2),
                      ncol = length(x))

  }

  # make sure to untransform by the provided mass matrix before estimating!

  # set up trace store (grab values of target variables from graph to get
  # dimension and names)
  x_projected <- (mass_cholesky %*% x)[, 1]
  dag$send_parameters(x_projected)
  init_trace <- dag$trace_values()
  trace <- matrix(NA,
                  nrow = n_samples %/% thin,
                  ncol = length(init_trace))
  colnames(trace) <- names(init_trace)

  # set up log joint density store
  ljd <- rep(NA, n_samples)

  # track acceptance
  accept_trace <- rep(0, n_samples)

  # get free parameter dimension
  npar <- length(x)

  # set up progress bar
  if (verbose)
    pb <- txtProgressBar(max = n_samples, style = 3)

  # loop through iterations
  for (i in 1:n_samples) {

    # copy old state
    x_minus <- x_plus <- x_old <- x
    p_minus <- p_plus <- p_old <- p <- rnorm(npar)

    # draw slice variable
    u_max <- exp(calc_H(x, p, dag, mass_cholesky))
    u <- runif(1, 0, u_max)

    # counters
    j <- 0
    n <- 1
    continue <- TRUE

    while (continue) {

      # pick a direction
      v <- sample(c(1, -1), 1)

      if (v == 1) {

        # move right
        new_step <- build_tree(x_plus,
                               p_plus,
                               u = u,
                               v = v,
                               j = j,
                               epsilon = epsilon,
                               x_old = x_old,
                               p_old = p_old,
                               mass_cholesky = mass_cholesky,
                               dag = dag)

        x_plus <- new_step$x_plus
        p_plus <- new_step$p_plus

      } else {

        # move right
        new_step <- build_tree(x_minus,
                               p_minus,
                               u = u,
                               v = v,
                               j = j,
                               epsilon = epsilon,
                               x_old = x_old,
                               p_old = p_old,
                               mass_cholesky = mass_cholesky,
                               dag = dag)

        x_minus <- new_step$x_minus
        p_minus <- new_step$p_minus

      }

      # add premature rejection here, checked in build_tree!

      ## test whether to accept this state
      if (is.na(new_step$continue) | is.nan(new_step$continue))
        new_step$continue <- FALSE

      # if it's valid, do accept/reject
      if (new_step$continue && runif(1, 0, 1) <= (new_step$n / n))
        x <- new_step$x_new

      n <- n + new_step$n

      # see if it did a U-turn
      u_turn <- turned(x_plus, x_minus, p_plus, p_minus)
      continue <- new_step$continue & !u_turn

      ## Stop trajectory if there are any problems, probably happens
      ## when jumping way too far into the tails and the model isn't
      ## defined
      if(is.na(continue) | is.nan(continue))
        continue <- FALSE

      j <- j + 1

      ## Stop doubling if too many or it's diverged enough
      if (continue & j > max_doublings)
        break

    }

    # tune epsilon using the dual averaging algorithm
    if (tune_epsilon) {

      nstep_trace[i] <- j - 1

      Hbar_trace [i + 1] <- (1 - 1 / (i + t0)) * Hbar_trace[i] +
        (target_acceptance - new_step$alpha / new_step$nalpha) / (i + t0)

      ## If logalpha not defined, skip this updating step and use
      ## the last one.
      if (is.nan(Hbar_trace[i + 1]))
        Hbar_trace[i + 1] <- abs(Hbar_trace[i])

      epsilon_trace[i + 1] <- exp(mu - sqrt(i) * Hbar_trace[i + 1] / gamma)
      epsilon_bar_trace[i + 1] <- exp(i ^ (-kappa) *
                                        log(epsilon_trace[i + 1]) +
                                        (1 - i ^ (-kappa)) *
                                        log(epsilon_bar_trace[i]))
      epsilon <- epsilon_trace[i + 1]

    } else {

      # otherwise, jitter epsilon
      epsilon <- epsilon_mean * rnorm(1, 0.9, 1.1)

    }

    # if we're estimating the mass matrix, and in the second half of the chain,
    # store the trace
    if (estimate_mass_matrix & i > floor(n_samples / 2)) {
      idx <- i - floor(n_samples / 2)
      trace_x[idx, ] <- x
    }

    # store the parameter values
    if (i %% thin == 0) {
      dag$send_parameters(x)
      trace[i / thin, ] <- dag$trace_values()
      ljd[i / thin] <- dag$log_density()
    }


    if (verbose)
      setTxtProgressBar(pb, i)

  }

  if (verbose) {

    close(pb)

    acc <- mean(accept_trace)
    message(sprintf('acceptance rate: %s',
                    prettyNum(round(acc, 3))))

  }

  # store the tuned epsilon
  if (tune_epsilon)
    control$epsilon <- epsilon_bar_trace[n_samples + 1]

  # get the Cholesky of the regularised, estimated mass matrix
  if (estimate_mass_matrix) {

    # back-transform the trace using the current mass matrix
    trace_x <- t(solve(mass_cholesky, t(trace_x)))

    # get the covariance
    covar <- cov(trace_x)

    # regularize it
    n_est = nrow(trace_x)
    mass_matrix_new <- (n_est / (n_est + 5)) *
      covar + 1e-3 * (5 / (n_est + 5)) * diag(length(x))

    # cholesky it
    mass_cholesky_new <- chol(mass_matrix_new)

    # store it
    control$mass_cholesky <- mass_cholesky_new
  }

  attr(trace, 'density') <- -ljd
  attr(trace, 'last_x') <- x
  attr(trace, 'control') <- control
  trace

}

calc_H <- function(x, p, dag, mass_cholesky) {
  x_projected <- (mass_cholesky %*% x)[, 1]
  dag$send_parameters(x_projected)
  dag$log_density() - 0.5 * sum(p ^ 2)
}

turned <- function (x_plus, x_minus, p_plus, p_minus) {
  diff <- (x_plus - x_minus)
  minus_turned <- diff %*% p_minus >= 0
  plus_turned <- diff %*% p_plus >= 0
  (minus_turned | plus_turned)[1, 1]
}


build_tree <- function(x,
                       p,
                       u,
                       v,
                       j,
                       epsilon,
                       x_old,
                       p_old,
                       dag,
                       mass_cholesky,
                       delta_max=1000) {

  # if at root of tree
  if (j == 0) {

    # send parameters once, and get gradients
    x_projected <- (mass_cholesky %*% x)[, 1]
    dag$send_parameters(x_projected)
    grad <- dag$gradients()

    # take one step in direction v
    epsilon <- v * epsilon
    p <- p + 0.5 * epsilon * grad
    x <- x + epsilon * p
    p <- p + 0.5 * epsilon * grad

    # check the trajectory
    H <- calc_H(x, p, dag, mass_cholesky)
    continue <- (H - log(u) + delta_max) > 0

    # switch this to premature reject!
    # check it's valid
    if (is.na(continue) | is.nan(continue))
      continue <- FALSE

    n <- ifelse(log(u) <= H, 1, 0)

    # get acceptance probability
    H_diff <- calc_H(x, p, dag, mass_cholesky) -
      calc_H(x_old, p_old, dag, mass_cholesky)
    alpha <- min(1, exp(H_diff))

    # combine step info
    tree <- list(x_minus = x,
                 x_plus = x,
                 x_new = x,
                 p_minus = p,
                 p_plus = p,
                 continue = continue,
                 n = n,
                 alpha = alpha,
                 nalpha = 1)

  } else {
    # if recursing, build subtrees

    ## recursion - build left and right subtrees
    first_tree <- build_tree(x,
                             p,
                             u = u,
                             v = v,
                             j = j - 1,
                             epsilon = epsilon,
                             x_old = x_old,
                             p_old = p_old,
                             mass_cholesky = mass_cholesky,
                             dag = dag)

    # unpack the results into this environment
    unpack(first_tree)

    # check for bad step
    if(is.na(continue) | is.nan(continue))
      continue <- FALSE

    # try the opposite direction
    if (continue) {
      if (v == -1) {

        second_tree <- build_tree(x_minus,
                                  p_minus,
                                  u = u,
                                  v = v,
                                  j = j-1,
                                  epsilon = epsilon,
                                  x_old = x_old,
                                  p_old = p_old,
                                  mass_cholesky = mass_cholesky,
                                  dag = dag)

        x_minus <- second_tree$x_minus
        p_minus <- second_tree$p_minus

      } else {

        second_tree <- build_tree(x_plus,
                                  p_plus,
                                  u = u,
                                  v = v,
                                  j = j-1,
                                  epsilon = epsilon,
                                  x_old = x_old,
                                  p_old = p_old,
                                  mass_cholesky = mass_cholesky,
                                  dag = dag)

        x_plus <- second_tree$x_plus
        p_plus <- second_tree$p_plus

      }

      # from CM:
      # # This isn't in the paper but if both slice variables failed,
      # # then you get 0/0. So I skip this test. Likewise if model
      # # throwing errors, don't keep that x.
      n <- first_tree$n + second_tree$n

      if (!is.finite(n))
        n <- 0

      if (n != 0){
        # choose whether to use the second tree (rather than the first)
        if (runif(1, 0, 1) <= second_tree$n / n)
          x_new <- second_tree$x_new
      }

      # check whether to continue
      u_turn <- turned(x_plus, x_minus, p_plus, p_minus)
      continue <- first_tree$continue & second_tree$continue & !u_turn

    }

    tree <- list(x_minus = x_minus,
                 x_plus = x_plus,
                 x_new = x_new,
                 p_minus = p_minus,
                 p_plus = p_plus,
                 continue = continue,
                 n = n,
                 alpha = alpha,
                 nalpha = 1)
  }

  tree

}


