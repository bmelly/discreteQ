#Definition of the main function
dq_decomposition <-
  function(y,
           x,
           g,
           w = NULL,
           q.range = c(0.05, 0.95),
           method = "logit",
           bsrep = 200,
           alpha = 0.05,
           ys,
           cl = NULL,
           cluster = NULL,
           old.res = NULL,
           return.boot = FALSE,
           list_of_seeds = NULL,
           return.seeds = FALSE,
           estim.glm = fastglm::fastglm,
           par.estim = NULL) {
    if (is.null(old.res)) {
      x <- as.matrix(x)
      y0 <- y[g == 0]
      y1 <- y[g == 1]
      x0 <- x[g == 0, , drop = FALSE]
      x1 <- x[g == 1, , drop = FALSE]
      n0 <- length(y0)
      n1 <- length(y1)
      n <- length(y)
      if (is.null(w) == 1) w = rep(1, n)
      w0 <- w[g == 0]
      w1 <- w[g == 1]
      yu0 <- unique(y[g == 0])
      yu1 <- unique(y[g == 1])
      if (is.null(ys)) {
        if (length(yu0) < 100) {
          ys0 <- sort(yu0)
        } else {
          ys0 <-
            sort(unique(stats::quantile(y[g == 0], seq(
              0.01, 0.99, 0.01
            ), type = 1)))
        }
        if (length(yu1) < 100) {
          ys1 <- sort(yu1)
        } else {
          ys1 <-
            sort(unique(stats::quantile(y[g == 1], seq(
              0.01, 0.99, 0.01
            ), type = 1)))
        }
      } else if (length(ys) == 1) {
        if (ys > length(yu0)) {
          ys0 <- sort(yu0)
        } else {
          ys0 <-
            sort(unique(stats::quantile(y[g == 0],
                                        seq(
                                          1 / (ys + 1), ys / (ys + 1), 1 / ys
                                        ),
                                        type = 1)))
        }
        if (ys > length(yu1)) {
          ys1 <- sort(yu1)
        } else {
          ys1 <- sort(unique(stats::quantile(y[g == 1],
                                             seq(
                                               1 / (ys + 1), ys / (ys + 1), 1 / ys
                                             ),
                                             type = 1)))
        }
      } else{
        ys0 <- unique(sort(ys[ys %in% yu0]))
        ys1 <- unique(sort(ys[ys %in% yu1]))
      }
      max0 <- max(ys0)
      max1 <- max(ys1)
      if (!is.null(cl)) {
        parallel::clusterExport(
          cl,
          c(
            "boot_decomp",
            "uncond_cdfs_po",
            "uncond_cdfs_dr",
            "uncond_cdfs_dr_int",
            "uncond_cdfs_drp",
            "uncond_cdfs_drp_int",
            "objective"
          ),
          envir=environment()
        )
        doParallel::registerDoParallel(cl)
      }

      if (method == "poisson") {
        Fc <- uncond_cdfs_po(ys1, y1, x1, w1, x0, w0)
      } else if (method == "drp") {
        Fc <- uncond_cdfs_drp(ys1, y1, x1, w1, x0, w0, cl)
      } else{
        Fc <- uncond_cdfs_dr(ys1, y1, x1, w1, x0, w0, method, cl, estim.glm, par.estim)
      }

      Fc <- sort(Fc)
      F0 <- sapply(ys0, function(y) stats::weighted.mean((y0 <= y), w0))
      F1 <- sapply(ys1, function(y) stats::weighted.mean((y1 <= y), w1))

      F0.func <- stats::stepfun(ys0, c(0, F0))
      F1.func <- stats::stepfun(ys1, c(0, F1))
      Fc.func <- stats::stepfun(ys1, c(0, Fc))
      Q0.func <- stats::stepfun(F0, c(ys0, max0), right = TRUE)
      Q1.func <- stats::stepfun(F1, c(ys1, max1), right = TRUE)
      Qc.func <- stats::stepfun(Fc, c(ys1, max1), right = TRUE)

      if (!is.null(cluster)) {
        cluster <- c(cluster[g == 0], cluster[g == 1])
      }
      if(is.null(list_of_seeds)){
        list_of_seeds <- rngtools::RNGseq(bsrep, simplify=FALSE)
      }
      if (is.null(cl)) {
        F.b <-
          sapply(1:bsrep, function(i)
            boot_decomp(list_of_seeds[[i]], ys0, ys1, y0, y1, x0, x1, w0, w1, n0, n1, method, cluster, estim.glm, par.estim))
      } else{
        i <- NULL
        F.b <-  foreach::`%dopar%`(foreach::foreach(i = 1:bsrep),{boot_decomp(list_of_seeds[[i]], ys0, ys1, y0, y1, x0, x1, w0, w1, n0, n1, method, cluster, estim.glm, par.estim)})
        F.b <- matrix(unlist(F.b),ncol=bsrep)
      }
    } else {
      F0.func <- old.res$F0
      F1.func <- old.res$F1
      Fc.func <- old.res$Fc
      ys0 <- old.res$ys0
      ys1 <- old.res$ys1
      F0 <- F0.func(ys0)
      F1 <- F1.func(ys1)
      Fc <- Fc.func(ys1)
      F.b <- old.res$F.b
      Q0.func <- old.res$Q0
      Q1.func <- old.res$Q1
      Qc.func <- old.res$Qc
      max0 <- max(ys0)
      max1 <- max(ys1)
    }
    F0.b  <-  F.b[1:length(ys0), ]
    F1.b  <-  F.b[(length(ys0) + 1):(length(ys0) + length(ys1)), ]
    Fc.b  <-
      F.b[(length(ys0) + length(ys1) + 1):(length(ys0) + 2*length(ys1)), ]
    delta.0 <- F0.b - F0
    se.0 <- pmax(apply(F0.b, 1, function(x) stats::IQR(x) / 1.349), 1e-06)
    delta.1 <- F1.b - F1
    se.1 <- pmax(apply(F1.b, 1, function(x) stats::IQR(x) / 1.349), 1e-06)
    delta.c <- Fc.b - Fc
    se.c <- pmax(apply(Fc.b, 1, function(x) stats::IQR(x) / 1.349), 1e-06)
    select.0 <-
      (F0.b >= q.range[1]) * (rbind(0, F0.b[1:(nrow(F0.b) - 1),]) < q.range[2])
    select.1 <-
      (F1.b >= q.range[1]) * (rbind(0, F1.b[1:(nrow(F1.b) - 1),]) < q.range[2])
    select.c <-
      (Fc.b >= q.range[1]) * (rbind(0, Fc.b[1:(nrow(Fc.b) - 1),]) < q.range[2])
    zs.j <-
      apply(rbind(
        abs(delta.1 * select.1) / se.1,
        abs(delta.0 * select.0) / se.0,
        abs(delta.c * select.c) / se.c
      ),  2, max, na.rm = TRUE)
    crt.j <- stats::quantile(zs.j, 1 - alpha)
    ub.F0j.i <- sort(F0 + crt.j * se.0)
    lb.F0j.i <- sort(F0 - crt.j * se.0)
    ub.F0j.i <- ifelse(ub.F0j.i <= 1, ub.F0j.i, 1)
    lb.F0j.i <- ifelse(lb.F0j.i >= 0, lb.F0j.i, 0)
    lb.F0.func <- stats::stepfun(ys0, c(0, lb.F0j.i))
    ub.F0.func <- stats::stepfun(ys0, c(0, ub.F0j.i))
    ub.F1j.i <- sort(F1 + crt.j * se.1)
    lb.F1j.i <- sort(F1 - crt.j * se.1)
    ub.F1j.i <- ifelse(ub.F1j.i <= 1, ub.F1j.i, 1)
    lb.F1j.i <- ifelse(lb.F1j.i >= 0, lb.F1j.i, 0)
    lb.F1.func <- stats::stepfun(ys1, c(0, lb.F1j.i))
    ub.F1.func <- stats::stepfun(ys1, c(0, ub.F1j.i))
    ub.Fcj.i <- sort(Fc + crt.j * se.c)
    lb.Fcj.i <- sort(Fc - crt.j * se.c)
    ub.Fcj.i <- ifelse(ub.Fcj.i <= 1, ub.Fcj.i, 1)
    lb.Fcj.i <- ifelse(lb.Fcj.i >= 0, lb.Fcj.i, 0)
    lb.Fc.func <- stats::stepfun(ys1, c(0, lb.Fcj.i))
    ub.Fc.func <- stats::stepfun(ys1, c(0, ub.Fcj.i))
    ub.Q0j.func <-
      stats::stepfun(lb.F0j.i, c(ys0, max0), right = FALSE)
    lb.Q0j.func <-
      stats::stepfun(ub.F0j.i, c(ys0, max0), right = TRUE)
    ub.Q1j.func <-
      stats::stepfun(lb.F1j.i, c(ys1, max1), right = FALSE)
    lb.Q1j.func <-
      stats::stepfun(ub.F1j.i, c(ys1, max1), right = TRUE)
    ub.Qcj.func <-
      stats::stepfun(lb.Fcj.i, c(ys1, max0), right = FALSE)
    lb.Qcj.func <-
      stats::stepfun(ub.Fcj.i, c(ys1, max0), right = TRUE)

    knots.new <- sort(unique(c(stats::knots(Q1.func), stats::knots(Q0.func))))
    observed <-
      stats::stepfun(knots.new, c(
        Q1.func(knots.new - .Machine$double.eps) - Q0.func(knots.new - .Machine$double.eps),
        max1 - max0
      ), right = TRUE)
    knots.new <-
      sort(unique(c(
        stats::knots(lb.Q1j.func), stats::knots(ub.Q0j.func)
      )))
    lb.observed <-
      stats::stepfun(knots.new, c(
        lb.Q1j.func(knots.new - .Machine$double.eps) - ub.Q0j.func(knots.new - .Machine$double.eps),
        max1 - max0
      ), right = TRUE)
    knots.new <-
      sort(unique(c(
        stats::knots(ub.Q1j.func), stats::knots(lb.Q0j.func)
      )))
    ub.observed <-
      stats::stepfun(knots.new, c(
        ub.Q1j.func(knots.new - .Machine$double.eps) - lb.Q0j.func(knots.new - .Machine$double.eps),
        max1 - max0
      ), right = TRUE)

    knots.new <- sort(unique(c(stats::knots(Q1.func), stats::knots(Qc.func))))
    composition <-
      stats::stepfun(knots.new, c(
        Q1.func(knots.new - .Machine$double.eps) - Qc.func(knots.new - .Machine$double.eps),
        max1 - max0
      ), right = TRUE)
    knots.new <-
      sort(unique(c(
        stats::knots(lb.Q1j.func), stats::knots(ub.Qcj.func)
      )))
    lb.composition <-
      stats::stepfun(knots.new, c(
        lb.Q1j.func(knots.new - .Machine$double.eps) - ub.Qcj.func(knots.new - .Machine$double.eps),
        max1 - max0
      ), right = TRUE)
    knots.new <-
      sort(unique(c(
        stats::knots(ub.Q1j.func), stats::knots(lb.Qcj.func)
      )))
    ub.composition <-
      stats::stepfun(knots.new, c(
        ub.Q1j.func(knots.new - .Machine$double.eps) - lb.Qcj.func(knots.new - .Machine$double.eps),
        max1 - max0
      ), right = TRUE)
    knots.new <- sort(unique(c(stats::knots(Qc.func), stats::knots(Q0.func))))

    unexplained <-
      stats::stepfun(knots.new, c(
        Qc.func(knots.new - .Machine$double.eps) - Q0.func(knots.new - .Machine$double.eps),
        0
      ), right = TRUE)
    knots.new <-
      sort(unique(c(
        stats::knots(lb.Qcj.func), stats::knots(ub.Q0j.func)
      )))
    lb.unexplained <-
      stats::stepfun(knots.new, c(
        lb.Qcj.func(knots.new - .Machine$double.eps) - ub.Q0j.func(knots.new - .Machine$double.eps),
        0
      ), right = TRUE)
    knots.new <-
      sort(unique(c(
        stats::knots(ub.Qcj.func), stats::knots(lb.Q0j.func)
      )))
    ub.unexplained <-
      stats::stepfun(knots.new, c(
        ub.Qcj.func(knots.new - .Machine$double.eps) - lb.Q0j.func(knots.new - .Machine$double.eps),
        0
      ), right = TRUE)

    res <-
      list(
        observed = observed,
        lb.observed = lb.observed,
        ub.observed = ub.observed,
        composition = composition,
        lb.composition = lb.composition,
        ub.composition = ub.composition,
        unexplained = unexplained,
        lb.unexplained = lb.unexplained,
        ub.unexplained = ub.unexplained,
        Q0 = Q0.func,
        lb.Q0 = lb.Q0j.func,
        ub.Q0 = ub.Q0j.func,
        Q1 = Q1.func,
        lb.Q1 = lb.Q1j.func,
        ub.Q1 = ub.Q1j.func,
        Qc = Qc.func,
        lb.Qc = lb.Qcj.func,
        ub.Qc = ub.Qcj.func,
        F0 = F0.func,
        F1 = F1.func,
        Fc = Fc.func,
        lb.F0 = lb.F0.func,
        lb.F1 = lb.F1.func,
        lb.Fc = lb.Fc.func,
        ub.F0 = ub.F0.func,
        ub.F1 = ub.F1.func,
        ub.Fc = ub.Fc.func,
        ys0 = ys0,
        ys1 = ys1,
        q.range = q.range,
        bsrep = bsrep
      )
    if(return.boot) res$F.b <- F.b
    if(return.seeds) res$seeds <- list_of_seeds
    res
  }

#summary
dq_summary.decomposition <-
  function(object,
           taus = 0.05,
           which = "unexplained") {
    if (length(taus) == 1)
      taus <- seq(object$q.range[1], object$q.range[2], taus)
    else if (max(taus) > object$q.range[2] |
             min(taus) < object$q.range[1])
      stop(
        "The quantiles specified by the argument `tau' must be within the range defined when calling cb.univariate()."
      )
    if (which == "observed") {
      tab <-
        cbind(
          taus,
          object$observed(taus),
          object$lb.observed(taus),
          object$ub.observed(taus)
        )
      colnames(tab) <-
        c("quantile", "Observed: Q1-Q0", "lower bound", "upper bound")
    } else if (which == "unexplained") {
      tab <-
        cbind(
          taus,
          object$unexplained(taus),
          object$lb.unexplained(taus),
          object$ub.unexplained(taus)
        )
      colnames(tab) <-
        c("quantile",
          "Unexplained: Qc-Q0",
          "lower bound",
          "upper bound")
    } else if (which == "composition") {
      tab <-
        cbind(
          taus,
          object$composition(taus),
          object$lb.composition(taus),
          object$ub.composition(taus)
        )
      colnames(tab) <-
        c("quantile",
          "Composition: Q1-Qc",
          "lower bound",
          "upper bound")
    } else if (which == "Q0") {
      tab <-
        cbind(taus,
              object$Q0(taus),
              object$lb.Q0(taus),
              object$ub.Q0(taus))
      colnames(tab) <-
        c("quantile", "Q0", "lower bound", "upper bound")
    } else if (which == "Q1") {
      tab <-
        cbind(taus,
              object$Q1(taus),
              object$lb.Q1(taus),
              object$ub.Q1(taus))
      colnames(tab) <-
        c("quantile", "Q1", "lower bound", "upper bound")
    } else if (which == "Qc") {
      tab <-
        cbind(taus,
              object$Qc(taus),
              object$lb.Qc(taus),
              object$ub.Qc(taus))
      colnames(tab) <-
        c("quantile",
          "Counterfactual QF",
          "lower bound",
          "upper bound")
    }
    tab
  }

#plot
dq_plot.decomposition <-
  function(object,
           which = "decomposition",
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
    if (is.null(which))
      which <- "decomposition"
    if (which != "F0" & which != "F1" & which != "Fc") {
      if (is.null(xlim)) {
        xlim <- object$q.range
      } else if (max(xlim) > object$q.range[2] |
                 min(xlim) < object$q.range[1]) {
        stop(
          "The quantiles specified by the argument `xlim' must be within the range defined when calling cb.univariate()."
        )
      }
      if (is.null(xlab))
        xlab <- "Probability"
    }
    if (which == "decomposition") {
      keep.mfrow <- par()$mfrow
      on.exit(par(mfrow = keep.mfrow), add = TRUE)
      par(mfrow = c(2, 2))
      if (is.null(ylim)) {
        ylim1 <-
          c(
            min(
              object$lb.Q0(xlim[1]),
              object$lb.Q1(xlim[1]),
              object$lb.Qc(xlim[1])
            ),
            max(
              object$ub.Q0(xlim[2]),
              object$ub.Q1(xlim[2]),
              object$ub.Qc(xlim[2])
            )
          )
      } else {
        ylim1 <- ylim
      }
      if (is.null(shift)) {
        shift1 <- min(diff(object$ys0), diff(object$ys1)) / 5
        shift2 = 0
      }  else
        shift1 <- shift2 <- shift
      if (length(col.l) != 3 | length(col.b) != 3) {
        col.l <- c(col.l, "dark green", "black")
        col.b <- c(col.b, "light green", "grey")
      }
      if (is.null(main)) {
        main <-
          c(
            "Quantile functions",
            "Observed difference (Q1-Q0)",
            "Composition effect (Q1-Qc)",
            "Unexplained difference (Qc-Q0)"
          )
      } else if (length(main) != 4) {
        stop("For the decomposition plots, the argument main must be a vector of length 4.")
      }
      plot_decomposition_internal(
        object,
        "Q0",
        xlim,
        ylim1,
        main[1],
        xlab,
        ylab,
        add = FALSE,
        col.l[1],
        col.b[1],
        shift = 0,
        lty.l,
        lwd.l,
        lty.b,
        lwd.b,
        support,
        ...
      )
      plot_decomposition_internal(
        object,
        "Qc",
        xlim,
        ylim1,
        NULL,
        xlab,
        ylab,
        add = TRUE,
        col.l[2],
        col.b[2],
        shift = shift1,
        lty.l,
        lwd.l,
        lty.b,
        lwd.b,
        support
      )
      plot_decomposition_internal(
        object,
        "Q1",
        xlim,
        ylim1,
        NULL,
        xlab,
        ylab,
        add = TRUE,
        col.l[3],
        col.b[3],
        shift = -shift1,
        lty.l,
        lwd.l,
        lty.b,
        lwd.b,
        support
      )
      graphics::legend(
        "topleft",
        legend = c("Q0", "Qc", "Q1"),
        col = col.b,
        lty = lty.l,
        lwd = lwd.b,
        bty = "n"
      )
      plot_decomposition_internal(
        object,
        "observed",
        xlim,
        ylim,
        main[2],
        xlab,
        ylab,
        add = FALSE,
        col.l[1],
        col.b[1],
        shift2,
        lty.l,
        lwd.l,
        lty.b,
        lwd.b,
        support,
        ...
      )
      plot_decomposition_internal(
        object,
        "composition",
        xlim,
        ylim,
        main[3],
        xlab,
        ylab,
        add = FALSE,
        col.l[1],
        col.b[1],
        shift2,
        lty.l,
        lwd.l,
        lty.b,
        lwd.b,
        support,
        ...
      )
      plot_decomposition_internal(
        object,
        "unexplained",
        xlim,
        ylim,
        main[4],
        xlab,
        ylab,
        add = FALSE,
        col.l[1],
        col.b[1],
        shift2,
        lty.l,
        lwd.l,
        lty.b,
        lwd.b,
        support,
        ...
      )
    } else if (which == "observed" |
               which == "composition" |
               which == "unexplained" | which == "Q0" |
               which == "Q1" | which == "Qc") {
      if (is.null(shift))
        shift <- 0
      plot_decomposition_internal(
        object,
        which,
        xlim,
        ylim,
        main,
        xlab,
        ylab,
        add,
        col.l,
        col.b,
        shift,
        lty.l,
        lwd.l,
        lty.b,
        lwd.b,
        support,
        ...
      )
    } else if (which == "F0" | which == "F1" | which == "Fc") {
      plot_decomposition_internal(
        object,
        which,
        xlim,
        ylim,
        main,
        xlab,
        ylab,
        add,
        col.l,
        col.b,
        shift,
        lty.l,
        lwd.l,
        lty.b,
        lwd.b,
        support,
        ...
      )
    }
  }

plot_decomposition_internal <-
  function(object,
           which = "QTE",
           xlim = NULL,
           ylim = NULL,
           main = NULL,
           xlab = NULL,
           ylab = NULL,
           add = FALSE,
           col.l = "dark blue",
           col.b = "light blue",
           shift = 0,
           lty.l = 1,
           lwd.l = 1,
           lty.b = 1,
           lwd.b = 5,
           support,
           ...) {
    if (which != "F0" & which != "F1" & which != "Fc") {
      if (which == "observed") {
        temp <- object$observed
        templ <- object$lb.observed
        tempu <- object$ub.observed
        allv <- expand.grid(object$ys1, object$ys0)
        allv <- sort(unique(allv[, 1] - allv[, 2]))
        if (is.null(main))
          main <- "Observed difference (Q1-Q0)"
        if (is.null(ylab))
          ylab <- "Quantile difference"
      } else if (which == "unexplained") {
        temp <- object$unexplained
        templ <- object$lb.unexplained
        tempu <- object$ub.unexplained
        allv <- expand.grid(object$ys1, object$ys0)
        allv <- sort(unique(allv[, 1] - allv[, 2]))
        if (is.null(main))
          main <- "Unexplained difference (Qc-Q0)"
        if (is.null(ylab))
          ylab <- "Quantile difference"
      } else if (which == "composition") {
        temp <- object$composition
        templ <- object$lb.composition
        tempu <- object$ub.composition
        allv <- expand.grid(object$ys0, object$ys0)
        allv <- sort(unique(allv[, 1] - allv[, 2]))
        if (is.null(main))
          main <- "Composition effect (Q1-Qc)"
        if (is.null(ylab))
          ylab <- "Quantile difference"
      } else if (which == "Q0") {
        temp <- object$Q0
        templ <- object$lb.Q0
        tempu <- object$ub.Q0
        allv <- object$ys0
        if (is.null(main))
          main <- "QF for Group 0"
        if (is.null(ylab))
          ylab <- "Quantile function"
      } else if (which == "Q1") {
        temp <- object$Q1
        templ <- object$lb.Q1
        tempu <- object$ub.Q1
        allv <- object$ys1
        if (is.null(main))
          main <- "QF for Group 1"
        if (is.null(ylab))
          ylab <- "Quantile function"
      } else if (which == "Qc") {
        temp <- object$Qc
        templ <- object$lb.Qc
        tempu <- object$ub.Qc
        allv <- object$ys0
        if (is.null(main))
          main <- "Counterfactual QF"
        if (is.null(ylab))
          ylab <- "Quantile function"
      }
      kx <- sort(unique(c(stats::knots(templ), stats::knots(tempu))))
      kx <- c(xlim[1], kx[kx >= xlim[1] & kx <= xlim[2]], xlim[2])
      allv <- allv[allv >= min(templ(kx)) & allv <= max(tempu(kx))]
      if (is.null(support))
        if (length(allv) > 200)
          support <- "continuous"
      else
        support <- "empirical"
      else if (support[1] != "empirical" &
               support[1] != "continuous")
        allv <- sort(unique(c(allv, support)))
      if (is.null(ylim))
        ylim <- c(min(templ(kx)), max(tempu(kx))) + shift
      if (add == FALSE)
        graphics::plot(
          NA,
          xlim = xlim,
          ylab = ylab,
          xlab = xlab,
          ylim = ylim,
          main = main,
          ...
        )
      if (support[1] != "continuous")
        for (i in 2:length(kx))
          for (j in allv[allv >= templ(kx[i]) &
                         allv <= tempu(kx[i] - .Machine$double.eps)])
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
      else
        for (i in 2:length(kx))
          graphics::polygon(
            c(kx[i - 1], kx[i - 1], kx[i], kx[i]),
            c(
              templ(kx[i]) + shift,
              tempu(kx[i] - .Machine$double.eps) + shift,
              tempu(kx[i] - .Machine$double.eps) + shift,
              templ(kx[i]) + shift
            ),
            col = col.b,
            border = col.b
          )
      kx <- sort(stats::knots(temp))
      kx <- c(xlim[1], kx[kx >= xlim[1] & kx <= xlim[2]], xlim[2])
      for (i in 2:length(kx))
        graphics::segments(
          kx[i - 1],
          temp(kx[i]) + shift,
          kx[i],
          temp(kx[i]) + shift,
          col = col.l,
          lty = lty.l,
          lwd = lwd.l
        )
    } else{
      if (which == "F0") {
        if (is.null(xlim)) {
          xlim <- range(object$ys0)
        } else if (max(xlim) > max(object$ys0) |
                   min(xlim) < min(object$ys0)) {
          stop(
            "The values specified by the argument `xlim' must be within the range of the observed values."
          )
        }
        temp <- object$F0
        templ <- object$lb.F0
        tempu <- object$ub.F0
        allv <-
          object$ys0[object$ys0 >= min(xlim) & object$ys0 <= max(xlim)]
        if (is.null(main))
          main <- "DF for Group 0"
      } else if (which == "F1") {
        if (is.null(xlim)) {
          xlim <- range(object$ys1)
        } else if (max(xlim) > max(object$ys1) |
                   min(xlim) < min(object$ys1)) {
          stop(
            "The values specified by the argument `xlim' must be within the range of the observed values."
          )
        }
        temp <- object$F1
        templ <- object$lb.F1
        tempu <- object$ub.F1
        allv <-
          object$ys1[object$ys1 >= min(xlim) & object$ys1 <= max(xlim)]
        if (is.null(main))
          main <- "DF for Group 1"
      } else if (which == "Fc") {
        if (is.null(xlim)) {
          xlim <- range(object$ys0)
        } else if (max(xlim) > max(object$ys0) |
                   min(xlim) < min(object$ys0)) {
          stop(
            "The values specified by the argument `xlim' must be within the range of the observed values."
          )
        }
        temp <- object$Fc
        templ <- object$lb.Fc
        tempu <- object$ub.Fc
        allv <-
          object$ys0[object$ys0 >= min(xlim) & object$ys0 <= max(xlim)]
        if (is.null(main))
          main <- "DF of the counterfactual outcome"
      }
      if (is.null(xlab))
        xlab <- "y"
      if (is.null(ylab))
        ylab <- "Distribution function"
      if (is.null(ylim))
        ylim <- c(min(templ(allv)), max(tempu(allv)))
      if (add == FALSE) {
        graphics::plot(
          temp(allv),
          xlim = xlim,
          ylab = ylab,
          xlab = xlab,
          ylim = ylim,
          main = main,
          sub = " ",
          type = "n",
          ...
        )
      }
      for (v in 1:(length(allv) - 1)) {
        graphics::polygon(
          c(allv[v], allv[v + 1], allv[v + 1], allv[v]),
          c(templ(allv[v]),
            templ(allv[v]),
            tempu(allv[v]),
            tempu(allv[v])),
          col = col.b,
          border = NA
        )
      }
      graphics::lines(
        temp,
        verticals = FALSE,
        do.points = FALSE,
        col = col.l,
        lty = lty.l,
        lend = "butt"
      )
    }
  }

boot_decomp <-
  function(seed,
           ys0,
           ys1,
           y0,
           y1,
           x0,
           x1,
           w0,
           w1,
           n0,
           n1,
           method,
           cluster,
           estim.glm,
           par.estim) {
    rngtools::RNGseed(seed)
    if (is.null(cluster)) {
      bw0 <- stats::rexp(n0) * w0
      bw1 <- stats::rexp(n1) * w1
      bw <- c(bw0, bw1)
    } else{
      cluster_id <- unique(cluster)
      nc <- length(cluster_id)
      bwc <- stats::rexp(nc)
      bw <-
        plyr::join(
          x = data.frame(cluster = cluster),
          y = data.frame(cluster = cluster_id, bw = bwc),
          by = "cluster",
          type = "left"
        )[, "bw"] * c(w0, w1)
      bw0 <- bw[1:n0]
      bw1 <- bw[(n0 + 1):(n0 + n1)]
    }
    if (method == "poisson") {
      c(
        sapply(ys0, function(y)
          stats::weighted.mean((y0 <= y), bw0)),
        sapply(ys1, function(y)
          stats::weighted.mean((y1 <= y), bw1)),
        uncond_cdfs_po(ys1, y1, x1, bw1, x0, bw0)
      )
    } else if (method == "drp") {
      c(
        sapply(ys0, function(y)
          stats::weighted.mean((y0 <= y), bw0)),
        sapply(ys1, function(y)
          stats::weighted.mean((y1 <= y), bw1)),
        uncond_cdfs_drp(ys1, y1, x1, bw1, x0, bw0, cl = NULL)
      )
    } else{
      c(
        sapply(ys0, function(y)
          stats::weighted.mean((y0 <= y), bw0)),
        sapply(ys1, function(y)
          stats::weighted.mean((y1 <= y), bw1)),
        uncond_cdfs_dr(ys1, y1, x1, bw1, x0, bw0, method, cl = NULL, estim.glm, par.estim)
      )
    }
  }
