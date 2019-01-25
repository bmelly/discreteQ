#Definition of the main function
dq_qte <-
  function(y,
           d,
           x = NULL,
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
           return.seeds = FALSE) {
    if (!is.null(cl)){
      parallel::clusterExport(cl, c("boot_cdfs", "uncond_cdfs_po", "uncond_cdfs_dr", "uncond_cdfs_dr_int", "uncond_cdfs_drp", "uncond_cdfs_drp_int", "objective"))
      doParallel::registerDoParallel(cl)
    }
    if(is.null(old.res)){
      y0 <- y[d == 0]
      y1 <- y[d == 1]
      n0 <- length(y0)
      n1 <- length(y1)
      n <- length(y)
      if (is.null(w)) w <- rep(1, n)
      if (sq <- is.null(x)) {
        method <- "sample"
        x0 <- x1 <- x <- NULL
      } else {
        x <- as.matrix(x)
        x0 <- x[d == 0, , drop = FALSE]
        x1 <- x[d == 1, , drop = FALSE]
        x <- rbind(x0, x1)
      }
      if (!is.null(cluster)){
        cluster <- c(cluster[d == 0], cluster[d == 1])
      }
      w0 <- w[d == 0]
      w1 <- w[d == 1]
      w <- c(w0, w1)
      yu0 <- unique(y[d == 0])
      yu1 <- unique(y[d == 1])
      if (is.null(ys)) {
        if (length(yu0) < 100)
          ys0 <-
            sort(yu0)
        else
          ys0 <-
            sort(unique(quantile(y[d == 0], seq(
              1 / 100, 99 / 100, 1 / 100
            ), type = 1)))
        if (length(yu1) < 100)
          ys1 <-
            sort(yu1)
        else
          ys1 <-
            sort(unique(quantile(y[d == 1], seq(
              1 / 100, 99 / 100, 1 / 100
            ), type = 1)))
      } else if (length(ys) == 1) {
        ys0 <-
          sort(unique(quantile(y[d == 0], seq(
            1 / (ys + 1), ys / (ys + 1), 1 / ys
          ), type = 1)))
        ys1 <-
          sort(unique(quantile(y[d == 1], seq(
            1 / (ys + 1), ys / (ys + 1), 1 / ys
          ), type = 1)))
      } else{
        ys0 <- unique(sort(ys[ys %in% yu0]))
        ys1 <- unique(sort(ys[ys %in% yu1]))
      }
      max0 <- max(ys0)
      max1 <- max(ys1)

      if (method=="poisson"){
        F0 <- uncond_cdfs_po(ys0, y0, x0, w0, x, w)
        F1 <- uncond_cdfs_po(ys1, y1, x1, w1, x, w)
      } else if (method=="drp"){
        F0 <- uncond_cdfs_drp(ys0, y0, x0, w0, x, w, cl)
        F1 <- uncond_cdfs_drp(ys1, y1, x1, w1, x, w, cl)
      } else{
        F0 <- uncond_cdfs_dr(ys0, y0, x0, w0, x, w, method, cl)
        F1 <- uncond_cdfs_dr(ys1, y1, x1, w1, x, w, method, cl)
      }
      F0 <- sort(F0)
      F1 <- sort(F1)
      F0.func <- stats::stepfun(ys0, c(0, F0))
      F1.func <- stats::stepfun(ys1, c(0, F1))
      Q0.func <- stats::stepfun(F0, c(ys0, max0), right = TRUE)
      Q1.func <- stats::stepfun(F1, c(ys1, max1), right = TRUE)
      if(is.null(list_of_seeds)){
        list_of_seeds <- rngtools::RNGseq(bsrep, simplify=FALSE)
      }
      if(is.null(cl)) {
        F.b <- sapply(1:bsrep, function(i) boot_cdfs(list_of_seeds[[i]], ys0, ys1, y0, y1, x0, x1, w0, w1, n0, n1, method, cluster))
      } else{
        F.b <- foreach::`%dopar%`(foreach::foreach(i = 1:bsrep),{boot_cdfs(list_of_seeds[[i]], ys0, ys1, y0, y1, x0, x1, w0, w1, n0, n1, method, cluster)})
        F.b <- matrix(unlist(F.b),ncol=bsrep)
        #F.b <- parSapply(cl, 1:bsrep, function(iter) boot_cdfs(iter, ys0, ys1, y0, y1, x0, x1, w0, w1, n0, n1, method, cluster))
      }
    } else {
      F0.func <- old.res$F0
      F1.func <- old.res$F1
      ys0 <- old.res$ys0
      ys1 <- old.res$ys1
      F0 <- F0.func(ys0)
      F1 <- F1.func(ys1)
      F.b <- old.res$F.b
      q.range <- old.res$q.range
      Q0.func <- old.res$Q0
      Q1.func <- old.res$Q1
      max0 <- max(ys0)
      max1 <- max(ys1)
    }
    F0.b <- F.b[1:length(ys0),]
    F1.b <- F.b[(length(ys0)+1):(length(ys0)+length(ys1)),]
    delta.0 <- F0.b - F0
    se.0 <- apply(F0.b, 1, function(x) IQR(x) / 1.349)
    delta.1 <- F1.b - F1
    se.1 <- apply(F1.b, 1, function(x) IQR(x) / 1.349
    )
    select.1 <-
      (F1.b >= q.range[1]) * (rbind(0, F1.b[1:(nrow(F1.b) - 1), ]) < q.range[2])
    select.0 <-
      (F0.b >= q.range[1]) * (rbind(0, F0.b[1:(nrow(F0.b) - 1), ]) < q.range[2])
    zs.j <-
      apply(rbind(abs(delta.1 * select.1) / se.1, abs(delta.0 * select.0) / se.0),  2, max, na.rm = TRUE)

    crt.j <- quantile(zs.j, 1 - alpha)
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

    ub.Q0j.func <- stats::stepfun(lb.F0j.i, c(ys0, max0), right = FALSE)
    lb.Q0j.func <- stats::stepfun(ub.F0j.i, c(ys0, max0), right = TRUE)
    ub.Q1j.func <- stats::stepfun(lb.F1j.i, c(ys1, max1), right = FALSE)
    lb.Q1j.func <- stats::stepfun(ub.F1j.i, c(ys1, max1), right = TRUE)

    knots.new <- sort(unique(c(knots(Q1.func), knots(Q0.func))))
    QTE.func <-
      stats::stepfun(knots.new, c(
        Q1.func(knots.new - .Machine$double.eps) - Q0.func(knots.new - .Machine$double.eps),
        max1 - max0
      ), right = TRUE)
    knots.new <-
      sort(unique(c(
        knots(lb.Q1j.func), knots(ub.Q0j.func)
      )))
    lb.QTE.func <-
      stats::stepfun(knots.new, c(
        lb.Q1j.func(knots.new - .Machine$double.eps) - ub.Q0j.func(knots.new - .Machine$double.eps),
        max1 - max0
      ), right = TRUE)
    knots.new <-
      sort(unique(c(
        knots(ub.Q1j.func), knots(lb.Q0j.func)
      )))
    ub.QTE.func <-
      stats::stepfun(knots.new, c(
        ub.Q1j.func(knots.new - .Machine$double.eps) - lb.Q0j.func(knots.new - .Machine$double.eps),
        max1 - max0
      ), right = TRUE)

  res <-
    list(
      QTE = QTE.func,
      lb.QTE = lb.QTE.func,
      ub.QTE = ub.QTE.func,
      Q0 = Q0.func,
      lb.Q0 = lb.Q0j.func,
      ub.Q0 = ub.Q0j.func,
      Q1 = Q1.func,
      lb.Q1 = lb.Q1j.func,
      ub.Q1 = ub.Q1j.func,
      q.range = q.range,
      ys0 = ys0,
      ys1 = ys1,
      F0 = F0.func,
      F1 = F1.func,
      lb.F0 = lb.F0.func,
      lb.F1 = lb.F1.func,
      ub.F0 = ub.F0.func,
      ub.F1 = ub.F1.func
    )
  if(return.boot) res$F.b <- F.b
  if(return.seeds) res$seeds <- list_of_seeds
  res
}

#summary
dq_summary.qte <- function(object,
                           taus = 0.05,
                           which = "QTE") {
  if (length(taus) == 1) {
    taus <- seq(object$q.range[1], object$q.range[2], taus)
  } else if (max(taus) > object$q.range[2] |
             min(taus) < object$q.range[1]) {
    stop(
      "The quantiles specified by the argument `tau' must be within the range defined when calling cb.univariate()."
    )
  }
  if (which == "QTE") {
    tab <-
      cbind(taus,
            object$QTE(taus),
            object$lb.QTE(taus),
            object$ub.QTE(taus))
    colnames(tab) <-
      c("quantile", "QTE", "lower bound", "upper bound")
  } else if (which == "Q0") {
    tab <-
      cbind(taus,
            object$Q0(taus),
            object$lb.Q0(taus),
            object$ub.Q0(taus))
    colnames(tab) <-
      c("quantile", "QF of Y0", "lower bound", "upper bound")
  } else if (which == "Q1") {
    tab <-
      cbind(taus,
            object$Q1(taus),
            object$lb.Q1(taus),
            object$ub.Q1(taus))
    colnames(tab) <-
      c("quantile", "QF of Y1", "lower bound", "upper bound")
  } else if (which == "F0") {
    tab <-
      cbind(object$ys0,
            object$F0(object$ys0),
            object$lb.F0(object$ys0),
            object$ub.F0(object$ys0))
    colnames(tab) <-
      c("y", "DF of Y0", "lower bound", "upper bound")
  } else if (which == "F1") {
    tab <-
      cbind(object$ys1,
            object$F1(object$ys1),
            object$lb.F1(object$ys1),
            object$ub.F1(object$ys1))
    colnames(tab) <-
      c("y", "DF of Y0", "lower bound", "upper bound")
  }
  tab
}

#plot
dq_plot.qte <-
  function(object,
           which = NULL,
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
      which <- "QTE"
    if (is.null(shift))
      shift <- 0
    if (which!="F0" & which!="F1"){
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
      if (which == "QTE") {
        temp <- object$QTE
        templ <- object$lb.QTE
        tempu <- object$ub.QTE
        allv <- expand.grid(object$ys1, object$ys0)
        allv <- sort(unique(allv[, 1] - allv[, 2]))
        if (is.null(main))
          main <- "QTE and uniform bands"
        if (is.null(ylab))
          ylab <- "Quantile treatment effect"
      }
      if (which == "Q0") {
        temp <- object$Q0
        templ <- object$lb.Q0
        tempu <- object$ub.Q0
        allv <- object$ys0
        if (is.null(main))
          main <- "QF for the control outcome and uniform bands"
        if (is.null(ylab))
          ylab <- "Quantile function"
      }
      if (which == "Q1") {
        temp <- object$Q1
        templ <- object$lb.Q1
        tempu <- object$ub.Q1
        allv <- object$ys1
        if (is.null(main))
          main <- "QF for the treated outcome and uniform bands"
        if (is.null(ylab))
          ylab <- "Quantile function"
      }
      kx <- sort(unique(c(knots(templ), knots(tempu))))
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
        plot(
          NA,
          xlim = xlim,
          ylab = ylab,
          xlab = xlab,
          ylim = ylim,
          main = main,
          ...
        )
      if (support[1] != "continuous") {
        for (i in 2:length(kx))
          for (j in allv[allv >= templ(kx[i]) &
                         allv <= tempu(kx[i] - .Machine$double.eps)])
            segments(
              kx[i - 1],
              j + shift,
              kx[i],
              j + shift,
              col = col.b,
              lty = lty.b,
              lwd = lwd.b,
              lend = "butt"
            )
      } else {
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
      }
      kx <- sort(knots(temp))
      kx <- c(xlim[1], kx[kx >= xlim[1] & kx <= xlim[2]], xlim[2])
      for (i in 2:length(kx))
        segments(
          kx[i - 1],
          temp(kx[i]) + shift,
          kx[i],
          temp(kx[i]) + shift,
          col = col.l,
          lty = lty.l,
          lwd = lwd.l,
          lend="butt"
        )
    }
    else{
      if (which=="F0"){
        if (is.null(xlim)){
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
        allv <- object$ys0[object$ys0>=min(xlim) & object$ys0<=max(xlim)]
        if (is.null(main))
          main <- "DF for the control outcome and uniform bands"
      } else if (which=="F1") {
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
        allv <- object$ys1[object$ys1>=min(xlim) & object$ys1<=max(xlim)]
        if (is.null(main))
          main <- "DF for the treated outcome and uniform bands"
      }
      if (is.null(xlab))
        xlab <- "y"
      if (is.null(ylab))
        ylab <- "Distribution function"
      if (is.null(ylim))
        ylim <- c(min(templ(allv)), max(tempu(allv)))
      if (add == FALSE){
        plot(
          temp(allv),
          xlim = xlim,
          ylab = ylab,
          xlab = xlab,
          ylim = ylim,
          main = main,
          sub = " ",
          type="n",
          ...
        )
      }
      for (v in 1:(length(allv)-1)){
        polygon(
          c(allv[v], allv[v + 1], allv[v + 1], allv[v]),
          c(
            templ(allv[v]),
            templ(allv[v]),
            tempu(allv[v]),
            tempu(allv[v])
          ),
          col = col.b,
          border = NA
        )
      }
        lines(
          temp,
          verticals = FALSE,
          do.points = FALSE,
          col = col.l,
          lty = lty.l,
          lend = "butt"
        )
    }
  }

uncond_cdfs_dr_int <-
  function(ys, ye, xe, we, x, w, method) {
    if (method != "sample") {
      if (method == "logit" |
          method == "probit" |
          method == "cauchit" | method == "cloglog") {
        suppressWarnings(fit  <-
                           glm.fit(
                             xe,
                             (ye <= ys),
                             weights = we,
                             family = binomial(link = method)
                           )$coef)
      } else if (method == "lpm") {
        fit  <- lm.wfit(xe, (ye <= ys), w = we)$coef
      } else
        stop("The selected method has not yet been implemented.")
      F <- x %*% fit
      if (method == "logit" |
          method == "probit" |
          method == "cauchit" | method == "cloglog") {
        F <- binomial(method)$linkinv(F)
      }
      stats::weighted.mean(F, w)
    } else {
      stats::weighted.mean((ye <= ys), w = we)
    }
  }

uncond_cdfs_dr <- function(ys, ye, xe, we,  x, w, method, cl){
  if (is.null(cl)) {
    sapply(
      ys,
      uncond_cdfs_dr_int,
      ye = ye,
      xe = xe,
      we = we,
      x = x,
      w = w,
      method = method
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
     c(unlist(foreach::`%dopar%`(foreach::foreach(i = 1:length(ys)),{uncond_cdfs_dr_int(ys[i], ye, xe, we, x, w, method)})))
  }
}

uncond_cdfs_po <- function(ys, ye, xe, we,  x, w) {
  fit  <- glm(ye~xe-1, weights = we, family = poisson)$coef
  lambda <- exp(x %*% fit)
  sapply(ys, function(l) stats::weighted.mean(ppois(l, lambda = lambda), w))
}

boot_cdfs <- function(seed, ys0, ys1, y0, y1, x0, x1, w0, w1, n0, n1, method, cluster) {
  rngtools::RNGseed(seed)
  if (is.null(cluster)) {
    bw0 <- rexp(n0) * w0
    bw1 <- rexp(n1) * w1
    bw <- c(bw0, bw1)
  } else{
    cluster_id <- unique(cluster)
    nc <- length(cluster_id)
    bwc <- rexp(nc)
    bw <- plyr::join(x = data.frame(cluster=cluster), y = data.frame(cluster = cluster_id, bw=bwc), by = "cluster", type = "left")[,"bw"] * c(w0, w1)
    bw0 <- bw[1:n0]
    bw1 <- bw[(n0+1):(n0+n1)]
  }
  if (method == "poisson") {
    c(uncond_cdfs_po(ys0, y0, x0, bw0, rbind(x0, x1), bw),
      uncond_cdfs_po(ys1, y1, x1, bw1, rbind(x0, x1), bw))
  } else if (method == "drp") {
    c(uncond_cdfs_drp(ys0, y0, x0, bw0, rbind(x0, x1), bw, cl=NULL),
      uncond_cdfs_drp(ys1, y1, x1, bw1, rbind(x0, x1), bw, cl=NULL))
  } else{
    c(
      uncond_cdfs_dr(ys0, y0, x0, bw0, rbind(x0, x1), bw, method, cl=NULL),
      uncond_cdfs_dr(ys1, y1, x1, bw1, rbind(x0, x1), bw, method, cl=NULL)
    )
  }
}

# Distribution regression with Poisson link
# Maximum likelihood function
objective <- function(beta, y, binary, x, w = 1) {
  lambda <- exp(x %*% beta)
  prob <- pmin(pmax(ppois(y, lambda), 10 ^ -15), 1 - 10 ^ -15)
  - sum(w * (binary * log(prob) + (1 - binary) * log(1 - prob)))
}

uncond_cdfs_drp <- function(ys, ye, xe, we,  x, w, cl) {
  start <-
    glm(ye ~ xe - 1, weight = we, family = poisson)$coef
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
    c(unlist(foreach::`%dopar%`(foreach::foreach(i = 1:length(ys)),{uncond_cdfs_drp_int(ys[i], ye, xe, we, x, w, method)})))
  }
}

uncond_cdfs_drp_int <-
  function(ys, ye, xe, we, x, w, start) {
    fit  <- optim(
      start,
      objective,
      y = ys,
      binary = (ye <= ys),
      x = xe,
      w = we
    )$par
    lambda <- exp(x %*% fit)
    stats::weighted.mean(ppois(ys, lambda = lambda), w)
  }
