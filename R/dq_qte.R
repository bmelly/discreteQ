#Definition of the main function
dq_qte <- function(y, d, x=NULL, w=NULL, q.range=c(0.05,0.95), method="logit", bsrep=200, alpha=0.05){
  y0 <- y[d==0]
  y1 <- y[d==1]
  n0 <- length(y0)
  n1 <- length(y1)
  n <- length(y)
  if(is.null(w)) w <- rep(1,n)
  if (sq <- is.null(x)){
    method <- "sample"
  } else {
    x <- as.matrix(x)
    x0 <- x[d==0, , drop=FALSE]
    x1 <- x[d==1, , drop=FALSE]
  }
  w0 <- w[d==0]
  w1 <- w[d==1]
  ys <- sort(unique(y))
  ns <- length(ys)

  F.dr <- sapply(ys, uncond_cdfs_dr, y0=y0, y1=y1, x0=x0, x1=x1, w0=w0, w1=w1, x=x, w=w, method=method)

  F0 <- sort(F.dr[1,])
  F1 <- sort(F.dr[2,])

  Q0.func  <-  left.inv(ys, F0)
  Q1.func  <-  left.inv(ys, F1)
  if(!sq) x <- as.matrix(rbind(x0,x1))
  w <- c(w0,w1)
  F0.b  <-  F1.b  <-  matrix(NA, length(ys), bsrep)
  for(r in 1:bsrep){
    bw0 <- rexp(n0)*w0
    bw1 <- rexp(n1)*w1
    F.dr <- sapply(ys, uncond_cdfs_dr, y0=y0, y1=y1, x0=x0, x1=x1, w0=bw0, w1=bw1, x=x, w=w, method=method)
    F0.b[,r] <- F.dr[1,]
    F1.b[,r] <- F.dr[2,]
  }

  delta.0 <- F0.b - F0
  se.0 <- apply(F0.b, 1, FUN=function(x) IQR(x)/1.349)
  delta.1 <- F1.b - F1
  se.1 <- apply(F1.b, 1, FUN=function(x) IQR(x)/1.349)
  select.1 <- (F1.b>=q.range[1])*(rbind(0,F1.b[1:(nrow(F1.b)-1),])<q.range[2])
  select.0 <- (F0.b>=q.range[1])*(rbind(0,F0.b[1:(nrow(F0.b)-1),])<q.range[2])
  zs.j <- apply(rbind(abs(delta.1*select.1) / se.1, abs(delta.0*select.0) / se.0),  2, max, na.rm = TRUE)

  crt.j <- quantile(zs.j, 1-alpha)
  ub.F0j.i <- sort(F0 + crt.j*se.0)
  lb.F0j.i <- sort(F0 - crt.j*se.0)
  ub.F0j.i <- ifelse(ub.F0j.i <= 1, ub.F0j.i, 1)
  lb.F0j.i <- ifelse(lb.F0j.i >= 0, lb.F0j.i, 0)

  ub.F1j.i <- sort(F1 + crt.j*se.1)
  lb.F1j.i <- sort(F1 - crt.j*se.1)
  ub.F1j.i <- ifelse(ub.F1j.i <= 1, ub.F1j.i, 1)
  lb.F1j.i <- ifelse(lb.F1j.i >= 0, lb.F1j.i, 0)

  ub.Q0j.func <- left.inv(ys, lb.F0j.i)
  lb.Q0j.func <- left.inv(ys, ub.F0j.i)
  ub.Q1j.func <- left.inv(ys, lb.F1j.i)
  lb.Q1j.func <- left.inv(ys, ub.F1j.i)

  knots.new <- sort(unique(c(knots(Q1.func), knots(Q0.func))))
  QTE.func <- stats::stepfun(knots.new, c(Q1.func(knots.new) - Q0.func(knots.new), max(ys)), right=TRUE)
  knots.new <- sort(unique(c(knots(lb.Q1j.func), knots(ub.Q0j.func))))
  lb.QTE.func <- stats::stepfun(knots.new, c(lb.Q1j.func(knots.new) - ub.Q0j.func(knots.new), max(ys)), right=TRUE)
  knots.new <- sort(unique(c(knots(ub.Q1j.func), knots(lb.Q0j.func))))
  ub.QTE.func <- stats::stepfun(knots.new, c(ub.Q1j.func(knots.new) - lb.Q0j.func(knots.new), max(ys)), right=TRUE)

  list(QTE=QTE.func, lb.QTE=lb.QTE.func, ub.QTE=ub.QTE.func, Q0=Q0.func, lb.Q0=lb.Q0j.func, ub.Q0=ub.Q0j.func, Q1=Q1.func, lb.Q1=lb.Q1j.func, ub.Q1=ub.Q1j.func, q.range=q.range, ys=ys)
}

#summary
dq_summary.qte <- function(object, taus=0.05, which="QTE"){
  if(length(taus)==1){
    taus <- seq(object$q.range[1], object$q.range[2], taus)
  } else if(max(taus)>object$q.range[2] | min(taus)<object$q.range[1]){
    stop("The quantiles specified by the argument `tau' must be within the range defined when calling cb.univariate().")
  }
  if(which=="QTE"){
    tab <- cbind(taus, object$QTE(taus), object$lb.QTE(taus), object$ub.QTE(taus))
    colnames(tab) <- c("quantile", "QTE", "lower bound", "upper bound")
  } else if(which=="Q0"){
    tab <- cbind(taus, object$Q0(taus), object$lb.Q0(taus), object$ub.Q0(taus))
    colnames(tab) <- c("quantile", "QF of Y0", "lower bound", "upper bound")
  } else if(which=="Q1"){
    tab <- cbind(taus, object$Q1(taus), object$lb.Q1(taus), object$ub.Q1(taus))
    colnames(tab) <- c("quantile", "QF of Y1", "lower bound", "upper bound")
  }
  tab
}

#plot
dq_plot.qte <- function(object, which=NULL, xlim=NULL, ylim=NULL, main=NULL, xlab=NULL, ylab=NULL, add=FALSE, col.l="dark blue", col.b="light blue", shift=NULL, lty.l=1, lwd.l=1, lty.b=1, lwd.b=5){
  if (is.null(which)) which <- "QTE"
  if (is.null(shift)) shift <- 0
  if (is.null(xlim)) xlim <- object$q.range
  else if (max(xlim)>object$q.range[2] | min(xlim)<object$q.range[1]) stop("The quantiles specified by the argument `xlim' must be within the range defined when calling cb.univariate().")
  if (is.null(xlab)) xlab <- "Probability"
  if (which=="QTE"){
    temp <- object$QTE
    templ <- object$lb.QTE
    tempu <- object$ub.QTE
    allv <- expand.grid(object$ys, object$ys)
    allv <- sort(unique(allv[,1]-allv[,2]))
    if(is.null(main)) main <- "QTE and uniform bands"
    if(is.null(ylab)) ylab <- "Quantile treatment effect"
  }
  if(which=="Q0"){
    temp <- object$Q0
    templ <- object$lb.Q0
    tempu <- object$ub.Q0
    allv <- object$ys
    if(is.null(main)) main <- "QF for the control outcome and uniform bands"
    if(is.null(ylab)) ylab <- "Quantile function"
  }
  if(which=="Q1"){
    temp <- object$Q1
    templ <- object$lb.Q1
    tempu <- object$ub.Q1
    allv <- object$ys
    if(is.null(main)) main <- "QF for the treated outcome and uniform bands"
    if(is.null(ylab)) ylab <- "Quantile function"
  }
  kx <- sort(unique(c(knots(temp), knots(templ), knots(tempu))))
  kx <- c(xlim[1],kx[kx>=xlim[1] & kx<=xlim[2]],xlim[2])
  if(is.null(ylim)) ylim <- c(min(templ(kx)), max(tempu(kx)))+shift
  if(add==FALSE) plot(NA, xlim=xlim, ylab=ylab, xlab=xlab, ylim=ylim, main=main)
  for(i in 2:length(kx)) for(j in allv[allv>=templ(kx[i]) & allv<=tempu(kx[i])]) segments(kx[i-1],j+shift,kx[i],j+shift,col=col.b, lty = lty.b, lwd=lwd.b, lend=1)
  for(i in 2:length(kx)) segments(kx[i-1],temp(kx[i])+shift,kx[i],temp(kx[i])+shift,col=col.l, lty = lty.l, lwd=lwd.l)
}

uncond_cdfs_dr <- function(ys, y0, y1, x0, x1, w0, w1, x, w, method){
  if (method!="sample") {
    if (method=="logit" | method=="probit" | method=="cauchit" | method=="cloglog") {
      fit0  <- glm.fit(x0, (y0<=ys), weights=w0, family = binomial(link = method))$coef
      fit1  <- glm.fit(x1, (y1<=ys), weights=w1, family = binomial(link = method))$coef
    } else if (method=="lpm") {
      fit0  <- lm.wfit(x0, (y0<=ys), w=w0)$coef
      fit1  <- lm.wfit(x1, (y1<=ys), w=w1)$coef
    } else stop("The selected method has not yet been implemented.")
    F0 <- x%*%fit0
    F1 <- x%*%fit1
    if(method=="logit" | method=="probit" | method=="cauchit" | method=="cloglog"){
      F0 <- binomial(method)$linkinv(F0)
      F1 <- binomial(method)$linkinv(F1)
    }
    F0 <- weighted.mean(F0,w)
    F1 <- weighted.mean(F1,w)
  } else {
    F0 <- weighted.mean((y0<=ys), w=w0)
    F1 <- weighted.mean((y1<=ys), w=w1)
  }
  c(F0,F1)
}
