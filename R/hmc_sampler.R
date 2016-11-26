# basic HMC sampler with capacity to tune epsilon inspired by the HMC code by
# James Hensman and Alex G.G. Matthews in the GPFlow python package
# (https://github.com/GPflow/GPflow/blob/master/GPflow/hmc.py)

hmc <- function (dag,
                 init,
                 n_samples,
                 thin,
                 verbose,
                 tune_epsilon = FALSE,
                 control = list(Lmin = 10,
                                Lmax = 20,
                                epsilon = 0.005,
                                tune_control = list(target_acceptance = 0.651,
                                                    accept_group = 100,
                                                    gamma = 0.1,
                                                    kappa = 0.75))) {

  # unpack control options
  unpack(control)

  # set initial location, log joint density and gradients
  x <- init
  dag$send_parameters(x)
  grad <- dag$gradients()
  logprob <- dag$log_density()

  # initial epsilon tuning parameters
  if (tune_epsilon) {

    # unpack tuning options
    unpack(tune_control)

    # keep track of progress
    epsilon_trace <- rep(NA, n_samples)

  }

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
      new_step <- leapfrog(list(x = x, p = p), dag, epsilon)

      # if it went awry, quit now
      if (is.null(new_step)) {

        reject <- TRUE
        break()

      } else {

        # otherwise, update the location an momemntum and continue
        x <- new_step$x
        p <- new_step$p

      }

    }

    p <- p - 0.5 * epsilon * grad

    # if the step was bad, reject it out of hand
    if (reject) {

      if (verbose)
        message ('proposal rejected due to numerical instability')

      x <- x_old

    } else {

      # otherwise do the Metropolis accept/reject step

      logprob <- new_step$logprob
      p <- new_step$p

      # acceptance ratio
      log_accept_ratio  <- accept_prob_fun(logprob, p, logprob_old, p_old, log = TRUE)
      log_u = log(runif(1))

      if (log_u < log_accept_ratio) {

        # on acceptance, store the success and leave the parameters in the dag
        # to be put in the trace
        accept_trace[i] <- 1

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

    # optionally tune epsilon
    if (tune_epsilon) {

      # acceptance rate over the last accept_group runs
      start <- max(1, i - accept_group)
      end <- i
      accept_rate <- mean(accept_trace[start:end], na.rm = TRUE)

      # decrease the adaptation rate as we go
      adapt_rate <- min(1, gamma * i ^ (-kappa))

      # shift epsilon in the right direction, making sure it never goes negative
      epsilon <- epsilon + pmax(-(epsilon + sqrt(.Machine$double.eps)),
                                adapt_rate * (accept_rate - target_acceptance))

      # keep track of epsilon
      epsilon_trace[i] <- epsilon

    }

  }

  if (verbose) {

    close(pb)

    acc <- mean(accept_trace)
    message(sprintf('acceptance rate: %s',
                    prettyNum(round(acc, 3))))

  }

  # store the tuned epsilon as the mean of the last half
  if (tune_epsilon) {
    start <- floor(n_samples/2)
    end <- n_samples
    control$epsilon <- mean(epsilon_trace[start:end], na.rm = TRUE)
  }

  attr(trace, 'density') <- -ljd
  attr(trace, 'last_x') <- x
  attr(trace, 'control') <- control
  trace

}

