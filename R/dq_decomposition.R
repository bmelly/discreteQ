#Definition of the main function
dq_decomposition <- function(y, x, g, w=NULL, q.range=c(0.05,0.95), method="logit", bsrep=200, alpha=0.05){
  x <- as.matrix(x)
  y0 <- y[g==0]
  y1 <- y[g==1]
  x0 <- x[g==0, , drop=FALSE]
  x1 <- x[g==1, , drop=FALSE]
  n0 <- length(y0)
  n1 <- length(y1)
  n <- length(y)
  if(is.null(w)==1) w=rep(1,n)
  w0 <- w[g==0]
  w1 <- w[g==1]
  ys <- sort(unique(y))
  ns <- length(ys)

  F.dr <- sapply(ys, decomp_cdfs_dr, y0=y0, y1=y1, x0=x0, x1=x1, w0=w0, w1=w1, method=method)
  F0 <- sort(F.dr[1,])
  F1 <- sort(F.dr[2,])
  Fc <- sort(F.dr[3,])

  Q0.func  <-  left.inv(ys, F0)
  Q1.func  <-  left.inv(ys, F1)
  Qc.func  <-  left.inv(ys, Fc)

  F0.b  <-  F1.b  <-  Fc.b  <-  matrix(NA, length(ys), bsrep)
  for(r in 1:bsrep){
    bw0 <- rexp(n0)*w0
    bw1 <- rexp(n1)*w1
    F.dr <- sapply(ys, decomp_cdfs_dr, y0=y0, y1=y1, x0=x0, x1=x1, w0=bw0, w1=bw1, method=method)
    F0.b[,r] <- F.dr[1,]
    F1.b[,r] <- F.dr[2,]
    Fc.b[,r] <- F.dr[3,]
  }

  delta.0 <- F0.b - F0
  se.0 <- apply(F0.b, 1, FUN=function(x) IQR(x)/1.349)
  delta.1 <- F1.b - F1
  se.1 <- apply(F1.b, 1, FUN=function(x) IQR(x)/1.349)
  delta.c <- Fc.b - Fc
  se.c <- apply(Fc.b, 1, FUN=function(x) IQR(x)/1.349)
  select.0 <- (F0.b>=q.range[1])*(rbind(0,F0.b[1:(nrow(F0.b)-1),])<q.range[2])
  select.1 <- (F1.b>=q.range[1])*(rbind(0,F1.b[1:(nrow(F1.b)-1),])<q.range[2])
  select.c <- (Fc.b>=q.range[1])*(rbind(0,Fc.b[1:(nrow(Fc.b)-1),])<q.range[2])
  zs.j <- apply(rbind(abs(delta.1*select.1) / se.1, abs(delta.0*select.0) / se.0, abs(delta.c*select.c) / se.c),  2, max, na.rm = TRUE)

  crt.j <- quantile(zs.j, 1-alpha)
  ub.F0j.i <- sort(F0 + crt.j*se.0)
  lb.F0j.i <- sort(F0 - crt.j*se.0)
  ub.F0j.i <- ifelse(ub.F0j.i <= 1, ub.F0j.i, 1)
  lb.F0j.i <- ifelse(lb.F0j.i >= 0, lb.F0j.i, 0)

  ub.F1j.i <- sort(F1 + crt.j*se.1)
  lb.F1j.i <- sort(F1 - crt.j*se.1)
  ub.F1j.i <- ifelse(ub.F1j.i <= 1, ub.F1j.i, 1)
  lb.F1j.i <- ifelse(lb.F1j.i >= 0, lb.F1j.i, 0)

  ub.Fcj.i <- sort(Fc + crt.j*se.c)
  lb.Fcj.i <- sort(Fc - crt.j*se.c)
  ub.Fcj.i <- ifelse(ub.Fcj.i <= 1, ub.Fcj.i, 1)
  lb.Fcj.i <- ifelse(lb.Fcj.i >= 0, lb.Fcj.i, 0)

  ub.Q0j.func <- left.inv(ys, lb.F0j.i)
  lb.Q0j.func <- left.inv(ys, ub.F0j.i)
  ub.Q1j.func <- left.inv(ys, lb.F1j.i)
  lb.Q1j.func <- left.inv(ys, ub.F1j.i)
  ub.Qcj.func <- left.inv(ys, lb.Fcj.i)
  lb.Qcj.func <- left.inv(ys, ub.Fcj.i)

  knots.new <- sort(unique(c(knots(Q1.func), knots(Q0.func))))
  observed <- stats::stepfun(knots.new, c(Q1.func(knots.new) - Q0.func(knots.new), max(ys)), right=TRUE)
  knots.new <- sort(unique(c(knots(lb.Q1j.func), knots(ub.Q0j.func))))
  lb.observed <- stats::stepfun(knots.new, c(lb.Q1j.func(knots.new) - ub.Q0j.func(knots.new), max(ys)), right=TRUE)
  knots.new <- sort(unique(c(knots(ub.Q1j.func), knots(lb.Q0j.func))))
  ub.observed <- stats::stepfun(knots.new, c(ub.Q1j.func(knots.new) - lb.Q0j.func(knots.new), max(ys)), right=TRUE)

  knots.new <- sort(unique(c(knots(Q1.func), knots(Qc.func))))
  unexplained <- stats::stepfun(knots.new, c(Q1.func(knots.new) - Qc.func(knots.new), max(ys)), right=TRUE)
  knots.new <- sort(unique(c(knots(lb.Q1j.func), knots(ub.Qcj.func))))
  lb.unexplained <- stats::stepfun(knots.new, c(lb.Q1j.func(knots.new) - ub.Qcj.func(knots.new), max(ys)), right=TRUE)
  knots.new <- sort(unique(c(knots(ub.Q1j.func), knots(lb.Qcj.func))))
  ub.unexplained <- stats::stepfun(knots.new, c(ub.Q1j.func(knots.new) - lb.Qcj.func(knots.new), max(ys)), right=TRUE)

  knots.new <- sort(unique(c(knots(Qc.func), knots(Q0.func))))
  composition <- stats::stepfun(knots.new, c(Qc.func(knots.new) - Q0.func(knots.new), max(ys)), right=TRUE)
  knots.new <- sort(unique(c(knots(lb.Qcj.func), knots(ub.Q0j.func))))
  lb.composition <- stats::stepfun(knots.new, c(lb.Qcj.func(knots.new) - ub.Q0j.func(knots.new), max(ys)), right=TRUE)
  knots.new <- sort(unique(c(knots(ub.Qcj.func), knots(lb.Q0j.func))))
  ub.composition <- stats::stepfun(knots.new, c(ub.Qcj.func(knots.new) - lb.Q0j.func(knots.new), max(ys)), right=TRUE)

  list(observed=observed, lb.observed=lb.observed, ub.observed=ub.observed, composition=composition, lb.composition=lb.composition, ub.composition=ub.composition, unexplained=unexplained, lb.unexplained=lb.unexplained, ub.unexplained=ub.unexplained, Q0=Q0.func, lb.Q0=lb.Q0j.func, ub.Q0=ub.Q0j.func, Q1=Q1.func, lb.Q1=lb.Q1j.func, ub.Q1=ub.Q1j.func, Qc=Qc.func, lb.Qc=lb.Qcj.func, ub.Qc=ub.Qcj.func, ys=ys, q.range=q.range)
}

#summary
dq_summary.decomposition <- function(object, taus=0.05, which="unexplained"){
  if(length(taus)==1) taus <- seq(object$q.range[1], object$q.range[2], taus)
  else if(max(taus)>object$q.range[2] | min(taus)<object$q.range[1]) stop("The quantiles specified by the argument `tau' must be within the range defined when calling cb.univariate().")
  if(which=="observed"){
    tab <- cbind(taus, object$observed(taus), object$lb.observed(taus), object$ub.observed(taus))
    colnames(tab) <- c("quantile", "Observed: Q1-Q0", "lower bound", "upper bound")
  } else if(which=="unexplained"){
    tab <- cbind(taus, object$unexplained(taus), object$lb.unexplained(taus), object$ub.unexplained(taus))
    colnames(tab) <- c("quantile", "Unexplained: Q1-Qc", "lower bound", "upper bound")
  } else if(which=="composition"){
    tab <- cbind(taus, object$composition(taus), object$lb.composition(taus), object$ub.composition(taus))
    colnames(tab) <- c("quantile", "Composition: Qc-Q0", "lower bound", "upper bound")
  } else if(which=="Q0"){
    tab <- cbind(taus, object$Q0(taus), object$lb.Q0(taus), object$ub.Q0(taus))
    colnames(tab) <- c("quantile", "Q0", "lower bound", "upper bound")
  } else if(which=="Q1"){
    tab <- cbind(taus, object$Q1(taus), object$lb.Q1(taus), object$ub.Q1(taus))
    colnames(tab) <- c("quantile", "Q1", "lower bound", "upper bound")
  } else if(which=="Qc"){
    tab <- cbind(taus, object$Qc(taus), object$lb.Qc(taus), object$ub.Qc(taus))
    colnames(tab) <- c("quantile", "Counterfactual QF", "lower bound", "upper bound")
  }
  tab
}

#plot
dq_plot.decomposition <- function(object, which="decomposition", xlim=NULL, ylim=NULL, main=NULL, xlab=NULL, ylab=NULL, add=FALSE, col.l="dark blue", col.b="light blue", shift=NULL, lty.l=1, lwd.l=1, lty.b=1, lwd.b=5){
  if (is.null(which)) which <- "decomposition"
  if (is.null(xlim)) xlim <- object$q.range
  else if (max(xlim)>object$q.range[2] | min(xlim)<object$q.range[1]) stop("The quantiles specified by the argument `xlim' must be within the range defined when calling cb.univariate().")
  if (is.null(xlab)) xlab <- "Probability"
  if (which=="decomposition") {
    keep.mfrow <- par()$mfrow
    on.exit(par(mfrow=keep.mfrow), add = TRUE)
    par(mfrow=c(2,2))
    if(is.null(ylim)){
      ylim1 <- c(min(object$lb.Q0(xlim[1]), object$lb.Q1(xlim[1]), object$lb.Qc(xlim[1])), max(object$ub.Q0(xlim[2]), object$ub.Q1(xlim[2]), object$ub.Qc(xlim[2])))
    } else {
        ylim1 <- ylim
    }
    if(is.null(shift)){
      shift1 <- min(diff(object$ys))/5
      shift2=0
    }  else shift1 <- shift2 <- shift
    if(length(col.l)!=3 | length(col.b)!=3){
      col.l <- c(col.l, "dark green", "black")
      col.b <- c(col.b, "light green", "grey")
    }
    if(is.null(main)){
      main <- c("Quantile functions", "Observed difference (QF1-QF0)", "Composition effect (QFc-QF0)", "Unexplained difference (QF1-QFc)")
    } else if(length(main)!=4){
      stop("For the decomposition plots, the argument main must be a vector of length 4.")
    }
    plot_decomposition_internal(object, "Q0", xlim, ylim1, main[1], xlab, ylab, add=FALSE, col.l[1], col.b[1], shift=0, lty.l, lwd.l, lty.b, lwd.b)
    plot_decomposition_internal(object, "Qc", xlim, ylim1, NULL, xlab, ylab, add=TRUE, col.l[2], col.b[2], shift=shift1, lty.l, lwd.l, lty.b, lwd.b)
    plot_decomposition_internal(object, "Q1", xlim, ylim1, NULL, xlab, ylab, add=TRUE, col.l[3], col.b[3], shift=-shift1, lty.l, lwd.l, lty.b, lwd.b)
    legend("topleft", legend=c("Q0","Qc","Q1"), col=col.b, lty=lty.l, lwd=lwd.b, bty="n")
    plot_decomposition_internal(object, "observed", xlim, ylim, main[2], xlab, ylab, add=FALSE, col.l[1], col.b[1], shift2, lty.l, lwd.l, lty.b, lwd.b)
    plot_decomposition_internal(object, "composition", xlim, ylim, main[3], xlab, ylab, add=FALSE, col.l[1], col.b[1], shift2, lty.l, lwd.l, lty.b, lwd.b)
    plot_decomposition_internal(object, "unexplained", xlim, ylim, main[4], xlab, ylab, add=FALSE, col.l[1], col.b[1], shift2, lty.l, lwd.l, lty.b, lwd.b)
  } else if(which=="observed" | which=="composition" | which=="unexplained" | which=="Q0" | which=="Q1" | which=="Qc"){
    if(is.null(shift)) shift <- 0
    plot_decomposition_internal(object, which, xlim, ylim, main, xlab, ylab, add, col.l, col.b, shift, lty.l, lwd.l, lty.b, lwd.b)
  }
}

plot_decomposition_internal <- function(object, which="QTE", xlim=NULL, ylim=NULL, main=NULL, xlab=NULL, ylab=NULL, add=FALSE, col.l="dark blue", col.b="light blue", shift=0, lty.l=1, lwd.l=1, lty.b=1, lwd.b=5){
  if(which=="observed"){
    temp <- object$observed
    templ <- object$lb.observed
    tempu <- object$ub.observed
    allv <- expand.grid(object$ys, object$ys)
    allv <- sort(unique(allv[,1]-allv[,2]))
    if(is.null(main)) main <- "Observed difference (QF1-QF0)"
    if(is.null(ylab)) ylab <- "Quantile difference"
  } else if(which=="unexplained"){
    temp <- object$unexplained
    templ <- object$lb.unexplained
    tempu <- object$ub.unexplained
    allv <- expand.grid(object$ys, object$ys)
    allv <- sort(unique(allv[,1]-allv[,2]))
    if(is.null(main)) main <- "Unexplained difference (QF1-QFc)"
    if(is.null(ylab)) ylab <- "Quantile difference"
  } else if(which=="composition"){
    temp <- object$composition
    templ <- object$lb.composition
    tempu <- object$ub.composition
    allv <- expand.grid(object$ys, object$ys)
    allv <- sort(unique(allv[,1]-allv[,2]))
    if(is.null(main)) main <- "Composition effect (QFc-QF0)"
    if(is.null(ylab)) ylab <- "Quantile difference"
  } else if(which=="Q0"){
    temp <- object$Q0
    templ <- object$lb.Q0
    tempu <- object$ub.Q0
    allv <- object$ys
    if(is.null(main)) main <- "QF for Group 0"
    if(is.null(ylab)) ylab <- "Quantile function"
  } else if(which=="Q1"){
    temp <- object$Q1
    templ <- object$lb.Q1
    tempu <- object$ub.Q1
    allv <- object$ys
    if(is.null(main)) main <- "QF for Group 1"
    if(is.null(ylab)) ylab <- "Quantile function"
  } else if(which=="Qc"){
    temp <- object$Qc
    templ <- object$lb.Qc
    tempu <- object$ub.Qc
    allv <- object$ys
    if(is.null(main)) main <- "Counterfactual QF"
    if(is.null(ylab)) ylab <- "Quantile function"
  }
  kx <- sort(unique(c(knots(temp), knots(templ), knots(tempu))))
  kx <- c(xlim[1],kx[kx>=xlim[1] & kx<=xlim[2]],xlim[2])
  if(is.null(ylim)) ylim <- c(min(templ(kx)), max(tempu(kx)))+shift
  if(add==FALSE) plot(NA, xlim=xlim, ylab=ylab, xlab=xlab, ylim=ylim, main=main)
  for(i in 2:length(kx)) for(j in allv[allv>=templ(kx[i]) & allv<=tempu(kx[i])]) segments(kx[i-1],j+shift,kx[i],j+shift,col=col.b, lty = lty.b, lwd=lwd.b, lend=1)
  for(i in 2:length(kx)) segments(kx[i-1],temp(kx[i])+shift,kx[i],temp(kx[i])+shift,col=col.l, lty = lty.l, lwd=lwd.l)
}

decomp_cdfs_dr <- function(ys, y0, y1, x0, x1, w0, w1, method){
  if(method=="logit" | method=="probit" | method=="cauchit" | method=="cloglog"){
    fit0  <- glm.fit(x0, (y0<=ys), weights=w0, family = binomial(link = method))$coef
    fit1  <- glm.fit(x1, (y1<=ys), weights=w1, family = binomial(link = method))$coef
  } else if(method=="lpm"){
    fit0  <- lm.wfit(x0, (y0<=ys), w=w0)$coef
    fit1  <- lm.wfit(x1, (y1<=ys), w=w1)$coef
  } else stop("The selected method has not yet been implemented.")
  F0 <- x0%*%fit0
  F1 <- x1%*%fit1
  Fc <- x1%*%fit0
  if(method=="logit" | method=="probit" | method=="cauchit" | method=="cloglog"){
    F0 <- binomial(method)$linkinv(F0)
    F1 <- binomial(method)$linkinv(F1)
    Fc <- binomial(method)$linkinv(Fc)
  }
  F0  <- weighted.mean(F0, w0)
  F1  <- weighted.mean(F1, w1)
  Fc  <- weighted.mean(Fc, w1)
  c(F0,F1,Fc)
}