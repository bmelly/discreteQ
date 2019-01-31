#Definition of the main function
dq_univariate <-
  function(y,
           q.range = c(0.05, 0.95),
           w = NULL,
           bsrep = 200,
           alpha = 0.05,
           ys,
           old.res = NULL,
           return.boot = FALSE,
           list_of_seeds = NULL,
           return.seeds = FALSE) {
    if (is.null(old.res)) {
      n <- length(y)
      if (is.null(w)) w <- rep(1, n)
      yu <- unique(y)
      if (is.null(ys)) {
        if (length(yu) < 1000){
          ys <- sort(yu)
        } else {
          ys <-
            sort(unique(stats::quantile(
              y, seq(1 / 1000, 999 / 1000, 1 / 1000), type = 1
            )))
        }
      } else if (length(ys) == 1) {
        ys <-
          sort(unique(stats::quantile(y, seq(
            1 / (ys + 1), ys / (ys + 1), 1 / ys
          ), type = 1)))
      } else {
        ys <- unique(sort(ys[ys %in% yu]))
      }
      F <-  sapply(ys, function(ys) stats::weighted.mean((y <= ys), w))
      F.func <- stats::stepfun(ys, c(0, F))
      Q.func  <-  stats::stepfun(F, c(ys, max(ys)), right = TRUE)
      if(is.null(list_of_seeds)){
        list_of_seeds <- rngtools::RNGseq(bsrep, simplify=FALSE)
      }
      F.b <- sapply(
                1:bsrep,
                function(i){
                      rngtools::RNGseed(list_of_seeds[[i]])
                      bw  <-  w * stats::rexp(n)
                      sapply(ys, function(ys) stats::weighted.mean((y <= ys), bw))
                  }
              )
    } else {
      F.func <- old.res$F
      ys <- old.res$ys
      F <- F.func(ys)
      F.b <- old.res$F.b
      Q.func <- old.res$Q
    }
    delta <- F.b - F
    se <- apply( F.b, 1, function(x) stats::IQR(x) / 1.349)
    select <-
      (F.b >= q.range[1]) * (rbind(0, F.b[1:(nrow(F.b) - 1), ]) < q.range[2])
    zs <- apply(rbind(abs(delta * select) / se),  2, max, na.rm = TRUE)

    crt.q <- stats::quantile(zs, 1 - alpha)
    #bands for QFs
    ub.F.i <- sort(F + crt.q * se)
    lb.F.i <- sort(F - crt.q * se)
    ub.F.i <- ifelse(ub.F.i <= 1, ub.F.i, 1)
    lb.F.i <- ifelse(lb.F.i >= 0, lb.F.i, 0)
    lb.F.func <- stats::stepfun(ys, c(0, lb.F.i))
    ub.F.func <- stats::stepfun(ys, c(0, ub.F.i))
    ub.Q.func <- stats::stepfun(lb.F.i, c(ys, max(ys)), right = FALSE)
    lb.Q.func <- stats::stepfun(ub.F.i, c(ys, max(ys)), right = TRUE)
    res <-
      list(
        Q = Q.func,
        ub.Q = ub.Q.func,
        lb.Q = lb.Q.func,
        F = F.func,
        lb.F = lb.F.func,
        ub.F = ub.F.func,
        q.range = q.range,
        ys = ys,
        bsrep = bsrep
    )
    if(return.boot) res$F.b <- F.b
    if(return.seeds) res$seeds <- list_of_seeds
    res
  }

#summary
dq_summary.univariate <- function(object, taus = 0.05) {
  if (length(taus) == 1 & taus<1) {
    taus <- seq(object$q.range[1], object$q.range[2], taus)
  } else if (length(taus) == 1 & taus<1){
    taus <- seq(object$q.range[1], object$q.range[2], length.out = taus)
  } else if (max(taus) > object$q.range[2] |
             min(taus) < object$q.range[1]) {
    stop(
      "The quantiles specified by the argument `tau' must be within the range defined when calling cb.univariate()."
    )
  }
  tab <-
    cbind(taus, object$Q(taus), object$lb.Q(taus), object$ub.Q(taus))
  colnames(tab) <- c("quantile", "QF", "lower bound", "upper bound")
  tab
}

#plot
dq_plot.univariate <-
  function(object,
           xlim = NULL,
           ylim = NULL,
           main = NULL,
           xlab = NULL,
           ylab = NULL,
           add = FALSE,
           col.l = "dark blue",
           col.b = "light blue",
           shift = NULL,
           lty.l = 1,
           lwd.l = 1,
           lty.b = 1,
           lwd.b = 5,
           support,
           ...) {
    if (is.null(shift))
      shift <- 0
    if (is.null(xlim)) {
      xlim <- object$q.range
    } else if (max(xlim) > object$q.range[2] |
               min(xlim) < object$q.range[1]) {
      stop(
        "The quantiles specified by the argument `xlim' must be within the range defined when calling cb.univariate()."
      )
    }
    if (is.null(main))
      main <- "Quantile function and uniform bands"
    if (is.null(xlab))
      xlab <- "Probability"
    if (is.null(ylab))
      ylab <- "Quantile function"
    if (is.null(ylim))
      ylim <- c(object$lb.Q(xlim[1]), object$ub.Q(xlim[2])) + shift
    kx <-
      sort(unique(c(
        knots(object$Q), knots(object$lb.Q), knots(object$ub.Q)
      )))
    kx <- c(xlim[1], kx[kx >= xlim[1] & kx <= xlim[2]], xlim[2])
    if (add == FALSE){
      graphics::plot(
        NA,
        xlim = xlim,
        ylab = ylab,
        xlab = xlab,
        ylim = ylim,
        main = main,
        ...
      )
    }
    all.y <- object$ys
    if (is.null(support)){
      if(length(all.y) > 200){
        support <- "continuous"
      } else {
      support <- "empirical"
      }
    } else if (support[1] != "empirical" &
             support[1] != "continuous"){
      all.y <- sort(unique(c(all.y, support)))
    }

    if (support[1] != "continuous"){
      for (i in 2:length(kx)){
        for (j in all.y[all.y >= object$lb.Q(kx[i]) &
                        all.y <= object$ub.Q(kx[i] - .Machine$double.eps)]){
          graphics::segments(
            kx[i - 1],
            j + shift,
            kx[i],
            j + shift,
            col = col.b,
            lty = lty.b,
            lwd = lwd.b,
            lend = 1
          )
        }
      }
    } else{
      for (i in 2:length(kx)){
        graphics::polygon(
          c(kx[i - 1], kx[i - 1], kx[i], kx[i]),
          c(
            object$lb.Q(kx[i]) + shift,
            object$ub.Q(kx[i] - .Machine$double.eps) + shift,
            object$ub.Q(kx[i] - .Machine$double.eps) + shift,
            object$lb.Q(kx[i]) + shift
          ),
          col = col.b,
          border = col.b
        )
      }
    }
    for (i in 2:length(kx)){
      graphics::segments(
        kx[i - 1],
        object$Q(kx[i]) + shift,
        kx[i],
        object$Q(kx[i]) + shift,
        col = col.l,
        lty = lty.l,
        lwd = lwd.l
      )
    }
  }
