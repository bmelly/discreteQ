uncond_cdfs_dr <- function(ys, ye, xe, we,  x, w, method, cl, estim.glm, par.estim) {
  if (is.null(cl)) {
    sapply(
      ys,
      uncond_cdfs_dr_int,
      ye = ye,
      xe = xe,
      we = we,
      x = x,
      w = w,
      method = method,
      estim.glm = estim.glm,
      par.estim
    )
  } else{
    # parSapply(cl, ys, function(level)
    #   uncond_cdfs_dr_int(
    #     level,
    #     ye = ye,
    #     xe = xe,
    #     we = we,
    #     x = x,
    #     w = w,
    #     method = method
    #   ))
    i <- NULL
    c(unlist(foreach::`%dopar%`(foreach::foreach(i = 1:length(ys)),{uncond_cdfs_dr_int(ys[i], ye, xe, we, x, w, method, estim.glm, par.estim)})))
  }
}

uncond_cdfs_dr_int <-
  function(ys, ye, xe, we, x, w, method, estim.glm, par.estim) {
    if (method != "sample") {
      if (method == "logit" |
          method == "probit" |
          method == "cauchit" |
          method == "cloglog" |
          method == "log") {
        suppressWarnings(fit  <-
                           do.call(
                             estim.glm,
                             c(list(xe,
                             y = 1*(ye <= ys),
                             weights = we,
                             family = stats::binomial(link = method)),
                             par.estim)
                           )
                         )
        if (fit$converge == FALSE) {
          warning(paste("The glm algorithm did not converge at the threshold ys=",ys))
        }
        fit <- fit$coef
      } else if (method == "lpm") {
        fit  <- stats::lm.wfit(xe, (ye <= ys), w = we)$coef
      } else
        stop("The selected method has not yet been implemented.")
      F <- x %*% fit
      if (method == "logit" |
          method == "probit" |
          method == "cauchit" | method == "cloglog") {
        F <- stats::binomial(method)$linkinv(F)
      }
      stats::weighted.mean(F, w)
    } else {
      stats::weighted.mean((ye <= ys), w = we)
    }
  }

uncond_cdfs_po <- function(ys, ye, xe, we,  x, w) {
  fit  <- stats::glm(ye ~ xe - 1, weights = we, family = stats::poisson)$coef
  lambda <- exp(x %*% fit)
  sapply(ys, function(l)
    stats::weighted.mean(stats::ppois(l, lambda = lambda), w))
}

# Distribution regression with Poisson link
# Maximum likelihood function
objective <- function(beta, y, binary, x, w = 1) {
  lambda <- exp(x %*% beta)
  prob <- pmin(pmax(stats::ppois(y, lambda), 10 ^ -15), 1 - 10 ^ -15)
  - sum(w * (binary * log(prob) + (1 - binary) * log(1 - prob)))
}

uncond_cdfs_drp <- function(ys, ye, xe, we,  x, w, cl) {
  start <-
    stats::glm(ye ~ xe - 1, weight = we, family = stats::poisson)$coef
  if (is.null(cl)) {
    sapply(
      ys,
      uncond_cdfs_drp_int,
      ye = ye,
      xe = xe,
      we = we,
      x = x,
      w = w,
      start = start
    )
  } else{
    # parSapply(cl, ys, function(level)
    #   uncond_cdfs_drp_int(
    #     level,
    #     ye = ye,
    #     xe = xe,
    #     we = we,
    #     x = x,
    #     w = w,
    #     start = start
    #   ))
    i <- NULL
    c(unlist(foreach::`%dopar%`(foreach::foreach(i = 1:length(ys)),{uncond_cdfs_drp_int(ys[i], ye, xe, we, x, w, start)})))
  }
}

uncond_cdfs_drp_int <-
  function(ys, ye, xe, we, x, w, start) {
    fit  <- stats::optim(
      start,
      objective,
      y = ys,
      binary = (ye <= ys),
      x = xe,
      w = we
    )$par
    lambda <- exp(x %*% fit)
    stats::weighted.mean(stats::ppois(ys, lambda = lambda), w)
}
