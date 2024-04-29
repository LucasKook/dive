fit_adaptive <- function(
    args, epochs, max_iter = 5, stepsize = 2, alpha = 0.1, ws = NULL,
    modFUN = "BoxCoxDA", ...
) {
  for (iter in seq_len(max_iter)) {
    mod <- do.call(modFUN, c(args, list(tf_seed = iter)))
    if (!is.null(ws))
      set_weights(mod$model, ws)
    fit(mod, epochs = epochs, ...)
    iPIT <- predict(mod, type = "cdf")
    unif <- ks.test(iPIT, "punif")$p.value
    Z <- dare:::.rm_int(model.matrix(args$anchor, data = args$data))
    indep <- dHSIC::dhsic.test(iPIT, Z, method = "gamma")$p.value
    mod$xi <- args$xi
    mod$p.unif <- unif
    mod$p.indep <- indep
    if (min(unif, indep) > alpha)
      return(mod)
    else
      args$xi <- ifelse(indep < unif, args$xi * (1 + stepsize),
                        args$xi / (1 + stepsize))
  }
  message("No solution for which uniformity and independence is not
          rejected at level alpha.")
  return(do.call(modFUN, args))
}
