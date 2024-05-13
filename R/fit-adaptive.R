fit_adaptive <- function(
    args, epochs, max_iter = 5, stepsize = 2, alpha = 0.1, ws = NULL,
    modFUN = "BoxCoxDA", start_xi = TRUE, lr = 0.05, cb = \() NULL,
    indep_over_unif = 1, ...
) {
  if (start_xi) {
    ### Build model
    mod <- do.call(modFUN, args)
    ### Warmstart if supplied
    if (!is.null(ws))
      set_weights(mod$model, ws)
    ### Initialize lambda s.t. weighting is equal
    iPIT <- predict(mod, type = "cdf")
    cmv <- mean((iPIT - ecdf(iPIT)(iPIT))^2)
    Z <- dare:::.rm_int(model.matrix(args$anchor, data = args$data))
    hsic <- dHSIC::dhsic.test(iPIT, Z, method = "gamma")$statistic
    xi_start <- indep_over_unif * (cmv / hsic)
    args$xi <- xi_start
    cat("\nInitializing xi =", round(xi_start, 3), "\n")
  }
  for (iter in seq_len(max_iter)) {
    ### Build model
    mod <- do.call(modFUN, c(args, list(
      tf_seed = iter, optimizer = optimizer_adam(lr))))
    ### Warmstart if supplied
    if (!is.null(ws))
      set_weights(mod$model, ws)
    ### Fit DIVE
    fit(mod, epochs = epochs, callbacks = cb(), ...)
    ### Check if uniform and independent
    iPIT <- predict(mod, type = "cdf")
    unif <- ks.test(iPIT, "punif")$p.value
    Z <- dare:::.rm_int(model.matrix(args$anchor, data = args$data))
    indep <- dHSIC::dhsic.test(iPIT, Z, method = "gamma")$p.value
    ### Return or update lambda and restart
    mod$xi <- args$xi
    mod$p.unif <- unif
    mod$p.indep <- indep
    if (min(unif, indep) > alpha) {
      message("\nSolution found at level alpha. Returning fitted model.\n")
      return(mod)
    }
    else {
      args$xi <- ifelse(indep < unif, args$xi * (1 + stepsize / iter),
                        args$xi / (1 + stepsize / iter))
      cat("\nUpdating xi from", mod$xi, "to", args$xi, "\n")
    }
    if (iter == max_iter) {
      message("\nNo solution for which uniformity and independence is not
          rejected at level alpha. Returning last fitted model.\n")
      return(mod)
    }
  }
}
