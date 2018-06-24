
####################################################################################################

# Monte Carlo Simulation

# Paper: Generic Inference on Quantile and Quantile Effect Functions for Discrete Outcomes

# Authors: V. Chernozhukov, I. Fernandez-Val, B. Melly, K. Wuethrich

####################################################################################################


library(doParallel)
library(doRNG)
library(xtable)

### Preliminaries
rm(list = ls())

###############
# Functions
###############

cdf=  function(ys, Fs){
  ys= sort(ys)
  Fs= sort(Fs)
  F= stepfun(ys, c(0,Fs))
  return(F)
}

left.inv= function(ys, Fs) {
  ys= sort(ys)
  Fs= sort(Fs)  
  iF= stepfun(Fs, c(ys,max(ys)), right=TRUE)
  return(iF)
}


wecdf = function(y,outcome,w){
  Ff    = weighted.mean((outcome <= y),w); 
  return(Ff)
}

ecdf = function(y,outcome){
  Ff    = mean(outcome <= y); 
  return(Ff)
}

ci <-  function(y0, y1, ycdf.0, ycdf.1, bsrep, alphas,  taus, taus.grid, all, ys.all=1){

  n.0 <- length(y0)
  n.1 <- length(y1)
  if(ys.all[1]==1 & length(ys.all)==1){
    ys.0 <- sort(unique(y0))
    ys.1 <- sort(unique(y1))
  } else ys.0 <- ys.1 <- ys.all

  ### Empirical CDF
  F0  <-  sapply(ys.0, ecdf, outcome=y0)
  F1  <-  sapply(ys.1, ecdf, outcome=y1)

  Q0.func  <-  left.inv(ys.0, F0)
  Q1.func  <-  left.inv(ys.1, F1)
  qte.func <- function(q) Q1.func(q)-Q0.func(q)
 
  Q0  <-  Q0.func(taus.grid)
  Q1  <-  Q1.func(taus.grid)
  qte  <-  qte.func(taus.grid)
  
  ### Boostrapped critical values
  
  F0.b  <-  matrix(NA, length(ys.0), bsrep)
  F1.b  <-  matrix(NA, length(ys.1), bsrep)
  Q1.b  <-  Q0.b <-  qte.b <-  matrix(NA, length(taus.grid), bsrep)
  for(r in 1:bsrep){
    w.0  <-  rexp(n.0)
    w.1  <-  rexp(n.1)
    
    F0.b[,r]  <-  sapply(ys.0, wecdf, outcome=y0, w.0)
    F1.b[,r]  <-  sapply(ys.1, wecdf, outcome=y1, w.1)
    
    if(all){
      Q0.b.r.func <-  left.inv(ys.0, F0.b[,r])
      Q1.b.r.func <-  left.inv(ys.1, F1.b[,r])
  
      Q0.b[,r]  <-  Q0.b.r.func(taus.grid)
      Q1.b[,r]  <-  Q1.b.r.func(taus.grid)
      qte.b[,r]  <-  Q1.b[,r] - Q0.b[,r]
    }
  }

  ### Our Procedure
  delta.0    = F0.b - F0;
  se.0 = apply(F0.b, 1, FUN=function(x) IQR(x)/1.349)
  delta.1    = F1.b - F1;
  se.1 = apply(F1.b, 1, FUN=function(x) IQR(x)/1.349)
  
  # Algorithm 1
  select.1 <- (F1.b>=min(taus))*(rbind(0,F1.b[1:(nrow(F1.b)-1),])<max(taus))
  select.0 <- (F0.b>=min(taus))*(rbind(0,F0.b[1:(nrow(F0.b)-1),])<max(taus))
  
  #crt for quantile range
  zs.0 = apply(rbind(abs(delta.0*select.0) / se.0),  2, max, na.rm = TRUE)
  zs.1 = apply(rbind(abs(delta.1*select.1) / se.1),  2, max, na.rm = TRUE)
  zs.j = apply(rbind(abs(delta.1*select.1) / se.1, abs(delta.0*select.0) / se.0),  2, max, na.rm = TRUE)

  ub.F0 <- lb.F0 <- ub.F0j <- lb.F0j <- ub.F0jj <- lb.F0jj <- matrix(NA, length(ycdf.0), length(alphas))
  ub.F1 <- lb.F1 <- ub.F1j <- lb.F1j <- ub.F1jj <- lb.F1jj <- matrix(NA, length(ycdf.1), length(alphas))
  ub.Q0 <- lb.Q0 <- ub.Q1 <- lb.Q1 <- ub.Q0j <- lb.Q0j <- ub.Q1j <- lb.Q1j <- ub.Q0jj <- lb.Q0jj <- ub.Q1jj <- lb.Q1jj <- matrix(NA, length(taus), length(alphas))
  temp.taus <- ub.qte <- lb.qte <- list()
  mean.length <- max.length <- rep(NA,length(alphas))
  
  for(i in 1:length(alphas)){
    alpha <- alphas[i]
    crt.q.0 = quantile(zs.0, 1-alpha)
    crt.q.1 = quantile(zs.1, 1-alpha)
  
    #bands for QFs
    ub.F0.i = sort(F0 + crt.q.0*se.0)
    lb.F0.i = sort(F0 - crt.q.0*se.0)
    ub.F0.i = ifelse(ub.F0.i <= 1, ub.F0.i, 1)
    lb.F0.i = ifelse(lb.F0.i >= 0, lb.F0.i, 0)
  
    ub.F1.i = sort(F1 + crt.q.1*se.1)
    lb.F1.i = sort(F1 - crt.q.1*se.1)
    ub.F1.i = ifelse(ub.F1.i <= 1, ub.F1.i, 1)
    lb.F.i1 = ifelse(lb.F1.i >= 0, lb.F1.i, 0)
    
    ub.Q0.func  = left.inv(ys.0, lb.F0.i)
    lb.Q0.func  = left.inv(ys.0, ub.F0.i)
    ub.Q1.func  = left.inv(ys.1, lb.F1.i)
    lb.Q1.func  = left.inv(ys.1, ub.F1.i)
    
    ub.Q0[,i] = ub.Q0.func(taus)
    lb.Q0[,i] = lb.Q0.func(taus)
    ub.Q1[,i]  = ub.Q1.func(taus)
    lb.Q1[,i]  = lb.Q1.func(taus)
  
    ub.F0.func = cdf(ys.0,ub.F0.i)
    lb.F0.func = cdf(ys.0,lb.F0.i)
    ub.F1.func = cdf(ys.1,ub.F1.i)
    lb.F1.func = cdf(ys.1,lb.F1.i)
    
    ub.F0[,i]  = ub.F0.func(ycdf.0)
    lb.F0[,i]  = lb.F0.func(ycdf.0)
    ub.F1[,i]  = ub.F1.func(ycdf.1)
    lb.F1[,i]  = lb.F1.func(ycdf.1)
    
    # Algorithm 2
    #joint bands for QFs and QE
    crt.j <- quantile(zs.j, 1-alpha)
  
    ub.F0j.i = sort(F0 + crt.j*se.0)
    lb.F0j.i = sort(F0 - crt.j*se.0)
    ub.F0j.i = ifelse(ub.F0j.i <= 1, ub.F0j.i, 1)
    lb.F0j.i = ifelse(lb.F0j.i >= 0, lb.F0j.i, 0)
   
    ub.F1j.i = sort(F1 + crt.j*se.1)
    lb.F1j.i = sort(F1 - crt.j*se.1)
    ub.F1j.i = ifelse(ub.F1j.i <= 1, ub.F1j.i, 1)
    lb.F1j.i = ifelse(lb.F1j.i >= 0, lb.F1j.i, 0)
    
    ub.Q0j.func  = left.inv(ys.0, lb.F0j.i)
    lb.Q0j.func  = left.inv(ys.0, ub.F0j.i)
    ub.Q1j.func  = left.inv(ys.1, lb.F1j.i)
    lb.Q1j.func  = left.inv(ys.1, ub.F1j.i)
    
    ub.Q0j[,i] = ub.Q0j.func(taus)
    lb.Q0j[,i] = lb.Q0j.func(taus)
    ub.Q1j[,i]  = ub.Q1j.func(taus)
    lb.Q1j[,i]  = lb.Q1j.func(taus)
    temp.taus[[i]] <- unique(sort(pmax(pmin(c(taus, ub.F0j.i+1e-6, lb.F0j.i+1e-6, ub.F1j.i+1e-6, lb.F1j.i+1e-6, ub.F0j.i-1e-6, lb.F0j.i-1e-6, ub.F1j.i-1e-6, lb.F1j.i-1e-6),max(taus)),min(taus))))
    lb.qte[[i]] = lb.Q1j.func(temp.taus[[i]])-ub.Q0j.func(temp.taus[[i]])
    ub.qte[[i]] = ub.Q1j.func(temp.taus[[i]])-lb.Q0j.func(temp.taus[[i]])
  
    ub.F0j.func = cdf(ys.0,ub.F0j.i)
    lb.F0j.func = cdf(ys.0,lb.F0j.i)
    ub.F1j.func = cdf(ys.1,ub.F1j.i)
    lb.F1j.func = cdf(ys.1,lb.F1j.i)
    
    ub.F0j[,i]  = ub.F0j.func(ycdf.0)
    lb.F0j[,i]  = lb.F0j.func(ycdf.0)
    ub.F1j[,i]  = ub.F1j.func(ycdf.1)
    lb.F1j[,i]  = lb.F1j.func(ycdf.1)
    #length of the bands
    mean.length[i] <- mean(ub.Q1j.func(taus.grid)-lb.Q0j.func(taus.grid)-lb.Q1j.func(taus.grid)+ub.Q0j.func(taus.grid))
    max.length[i] <- max(ub.Q1j.func(taus.grid)-lb.Q0j.func(taus.grid)-lb.Q1j.func(taus.grid)+ub.Q0j.func(taus.grid))
   }
  
  if(all==0){
    return(list(ub.F0=ub.F0, lb.F0=lb.F0, ub.F1=ub.F1, lb.F1=lb.F1, ub.Q0=ub.Q0, lb.Q0=lb.Q0, ub.Q1=ub.Q1, lb.Q1=lb.Q1, ub.F0j=ub.F0j, lb.F0j=lb.F0j, ub.F1j=ub.F1j, lb.F1j=lb.F1j, ub.Q0j=ub.Q0j, lb.Q0j=lb.Q0j, ub.Q1j=ub.Q1j, lb.Q1j=lb.Q1j, ub.qte=ub.qte, lb.qte=lb.qte, temp.taus=temp.taus, mean.length=mean.length, max.length=max.length))
  } else{
    
    ### Bootstrapping the quantiles: equal width CI
    delta.Q0 = Q0.b - Q0
    delta.Q1 = Q1.b - Q1
    delta.qte = qte.b - qte
  
    zs.Q0 = apply(rbind(abs(delta.Q0)),2,max,na.rm = TRUE)
    zs.Q1 = apply(rbind(abs(delta.Q1)),2,max,na.rm = TRUE)
    zs.qte = apply(rbind(abs(delta.qte)),2,max,na.rm = TRUE)
    
    se.0 = pmax(apply(Q0.b,1,FUN=function(x) IQR(x)/1.349), 1e-6)
    se.1 = pmax(apply(Q1.b,1,FUN=function(x) IQR(x)/1.349), 1e-6)
    se.qte = pmax(apply(qte.b,1,FUN=function(x) IQR(x)/1.349), 1e-6)
    
    zs.Q0.w = apply(rbind(abs(delta.Q0) / se.0),  2, max, na.rm = TRUE)
    zs.Q1.w = apply(rbind(abs(delta.Q1) / se.1),  2, max, na.rm = TRUE)
    zs.qte.w = apply(rbind(abs(delta.qte) /se.qte),  2, max, na.rm = TRUE)
    
    ub.Q0.ew <- lb.Q0.ew <- ub.Q1.ew <- lb.Q1.ew <- ub.qte.ew <- lb.qte.ew <- matrix(NA, length(taus.grid), length(alphas))
    ub.Q0.w <- lb.Q0.w <- ub.Q1.w <- lb.Q1.w <- ub.qte.w <- lb.qte.w <- matrix(NA, length(taus.grid), length(alphas))
    ub.Q0.s <- lb.Q0.s <- ub.Q1.s <- lb.Q1.s <- ub.qte.s <- lb.qte.s <- matrix(NA, length(taus.grid), length(alphas))
    ub.Q0.ss <- lb.Q0.ss <- ub.Q1.ss <- lb.Q1.ss <- ub.qte.ss <- lb.qte.ss <- matrix(NA, length(taus.grid), length(alphas))
    
    for(i in 1:length(alphas)){
      alpha <- alphas[i]
      crt.Q0 = quantile(zs.Q0, 1-alpha)
      crt.Q1 = quantile(zs.Q1, 1-alpha)
      crt.qte = quantile(zs.qte, 1-alpha)
    
      lb.Q0.ew[,i]= Q0 - crt.Q0
      ub.Q0.ew[,i] = Q0 + crt.Q0
      
      lb.Q1.ew[,i] = Q1 - crt.Q1
      ub.Q1.ew[,i] = Q1 + crt.Q1
      
      lb.qte.ew[,i] = qte - crt.qte
      ub.qte.ew[,i] = qte + crt.qte
    
      ### Bootstrapping the quantiles: weighted sup t 
      
      #crt
      crt.Q0 = quantile(zs.Q0.w, 1-alpha)
      crt.Q1 = quantile(zs.Q1.w, 1-alpha)
      crt.qte = quantile(zs.qte.w, 1-alpha)
      
      lb.Q0.w[,i] = Q0 - crt.Q0*se.0
      ub.Q0.w[,i] = Q0 + crt.Q0*se.0
      
      lb.Q1.w[,i] = Q1 - crt.Q1*se.1
      ub.Q1.w[,i] = Q1 + crt.Q1*se.1
      
      lb.qte.w[,i] = qte - crt.qte*se.qte
      ub.qte.w[,i] = qte + crt.qte*se.qte
    }
    
    #Adding noise and bootstrapping
    y.smooth.0 <- y0 + runif(n.0)
    y.smooth.1 <- y1 + runif(n.1)
    if(ys.all[1]==1 & length(ys.all)==1){
      ys.0 <- sort(unique(y.smooth.0))
      ys.1 <- sort(unique(y.smooth.1))
    } else ys.0 <- ys.1 <- ys.all
    F0  <-  sapply(ys.0, ecdf, outcome=y.smooth.0)
    F1  <-  sapply(ys.1, ecdf, outcome=y.smooth.1)
    Q0.func  <-  left.inv(ys.0, F0)
    Q1.func  <-  left.inv(ys.1, F1)
    Q.smooth.0  <-  Q0.func(taus.grid)
    Q.smooth.1  <-  Q1.func(taus.grid)
    qte.smooth  <-  Q.smooth.1 - Q.smooth.0
    Q1.b  <-  Q0.b  <-  qte.b  <-  matrix(NA, length(taus.grid), bsrep)
    for(r in 1:bsrep){
      w.0 <-  rexp(n.0)
      w.1 <-  rexp(n.1)
      F0.b.r <-  sapply(ys.0, wecdf, outcome=y.smooth.0, w.0)
      F1.b.r <-  sapply(ys.1, wecdf, outcome=y.smooth.1, w.1)
      Q0.b.r.func <-  left.inv(ys.0, F0.b.r)
      Q1.b.r.func <-  left.inv(ys.1, F1.b.r)
      Q0.b[,r] <-  Q0.b.r.func(taus.grid)
      Q1.b[,r] <-  Q1.b.r.func(taus.grid)
      qte.b[,r] <-  Q1.b[,r] - Q0.b[,r]
    }
    se.0 <-  apply(Q0.b, 1, FUN=function(x) IQR(x)/1.349)
    se.1 <-  apply(Q1.b, 1, FUN=function(x) IQR(x)/1.349)
    se.qte <-  apply(qte.b, 1, FUN=function(x) IQR(x)/1.349)
    delta.Q0 <-  Q0.b - Q.smooth.0
    delta.Q1 <-  Q1.b - Q.smooth.1
    delta.qte <-  qte.b - qte.smooth
    zs.Q0 <-  apply(abs(delta.Q0)/se.0,  2, max, na.rm = TRUE)
    zs.Q1 <-  apply(abs(delta.Q1)/se.1,  2, max, na.rm = TRUE)
    zs.qte <-  apply(abs(delta.qte)/se.qte,  2, max, na.rm = TRUE)
    for(i in 1:length(alphas)){
      alpha <- alphas[i]
      crt.Q0 <-  quantile(zs.Q0, 1-alpha)
      crt.Q1 <-  quantile(zs.Q1, 1-alpha)
      crt.qte <-  quantile(zs.qte, 1-alpha)
      lb.Q0.s[,i] <-  Q0 - crt.Q0*se.0
      ub.Q0.s[,i] <-  Q0 + crt.Q0*se.0
      lb.Q1.s[,i] <-  Q1 - crt.Q1*se.1
      ub.Q1.s[,i] <-  Q1 + crt.Q1*se.1
      lb.qte.s[,i] <-  qte - crt.qte*se.qte
      ub.qte.s[,i] <-  qte + crt.qte*se.qte
      lb.Q0.ss[,i] <-  Q.smooth.0 - crt.Q0*se.0
      ub.Q0.ss[,i] <-  Q.smooth.0 + crt.Q0*se.0
      lb.Q1.ss[,i] <-  Q.smooth.1 - crt.Q1*se.1
      ub.Q1.ss[,i] <-  Q.smooth.1 + crt.Q1*se.1
      lb.qte.ss[,i] <-  qte.smooth - crt.qte*se.qte
      ub.qte.ss[,i] <-  qte.smooth + crt.qte*se.qte
    }

    # Results
    return(list(ub.F0=ub.F0, lb.F0=lb.F0, ub.F1=ub.F1, lb.F1=lb.F1, ub.Q0=ub.Q0, lb.Q0=lb.Q0, ub.Q1=ub.Q1, lb.Q1=lb.Q1, ub.F0j=ub.F0j, lb.F0j=lb.F0j, ub.F1j=ub.F1j, lb.F1j=lb.F1j, ub.Q0j=ub.Q0j, lb.Q0j=lb.Q0j, ub.Q1j=ub.Q1j, lb.Q1j=lb.Q1j, ub.qte=ub.qte, lb.qte=lb.qte, lb.Q1.ew=lb.Q1.ew, ub.Q1.ew=ub.Q1.ew, lb.Q0.ew=lb.Q0.ew, ub.Q0.ew=ub.Q0.ew, lb.qte.ew=lb.qte.ew, ub.qte.ew=ub.qte.ew, lb.Q1.w=lb.Q1.w, ub.Q1.w=ub.Q1.w, lb.Q0.w=lb.Q0.w, ub.Q0.w=ub.Q0.w, lb.qte.w=lb.qte.w, ub.qte.w=ub.qte.w, lb.Q1.s=lb.Q1.s, ub.Q1.s=ub.Q1.s, lb.Q0.s=lb.Q0.s, ub.Q0.s=ub.Q0.s, lb.qte.s=lb.qte.s, ub.qte.s=ub.qte.s, lb.Q1.ss=lb.Q1.ss, ub.Q1.ss=ub.Q1.ss, lb.Q0.ss=lb.Q0.ss, ub.Q0.ss=ub.Q0.ss, lb.qte.ss=lb.qte.ss, ub.qte.ss=ub.qte.ss, temp.taus=temp.taus, mean.length=mean.length, max.length=max.length))
  }
}

evaluate <- function(ci.s, taus, tau.grid, f0.t, f1.t, q0.t.func, q1.t.func, all){
  
  tau.l <- min(taus)
  tau.u <- max(taus)
  q0.t <- q0.t.func(taus)
  q1.t <- q1.t.func(taus)  
  
  if(all==0) results <- matrix(NA,ncol(ci.s$lb.F0), 15) else results <- matrix(NA,ncol(ci.s$lb.F0), 39)
    #cov.F1 <- cov.Q0 <- cov.Q1 <- cov.F0j <- cov.F1j <- cov.Q0j <- cov.Q1j <- cov.qte <- cov.joint <- cov.joint.F <- cov.joint.Q <- cov.Q0.ew <- cov.Q1.ew <- cov.qte.ew <- cov.Q0.w <- cov.Q1.w <- cov.qte.w <- cov.Q0.s <- cov.Q1.s <- cov.qte.s <- cov.Q0.ss <- cov.Q1.ss <- cov.qte.ss <- 
  for(i in 1:ncol(ci.s$lb.F0)){
    qte.t <- q1.t.func(ci.s$temp.taus[[i]])-q0.t.func(ci.s$temp.taus[[i]])
    # Algorithm 1
    results[i,1] <- !any((f0.t[f0.t<=tau.u]<ci.s$lb.F0[f0.t<=tau.u,i])*(tau.l<ci.s$lb.F0[f0.t<=tau.u,i]),
                  (f0.t[f0.t>=tau.l]>ci.s$ub.F0[f0.t>=tau.l,i])*(tau.u>ci.s$ub.F0[f0.t>=tau.l,i]))
    results[i,2] <- !any((f1.t[f1.t<=tau.u]<ci.s$lb.F1[f1.t<=tau.u,i])*(tau.l<ci.s$lb.F1[f1.t<=tau.u,i]),
                  (f1.t[f1.t>=tau.l]>ci.s$ub.F1[f1.t>=tau.l,i])*(tau.u>ci.s$ub.F1[f1.t>=tau.l,i]))
    results[i,3] <- 1-((sum(q0.t+1e-12<ci.s$lb.Q0[,i])+sum(q0.t-1e-12>ci.s$ub.Q0[,i]))>0)
    results[i,4] <- 1-((sum(q1.t+1e-12<ci.s$lb.Q1[,i])+sum(q1.t-1e-12>ci.s$ub.Q1[,i]))>0)

    # Algorithm 2  
    results[i,5] <- !any((f0.t[f0.t<=tau.u]<ci.s$lb.F0j[f0.t<=tau.u,i])*(tau.l<ci.s$lb.F0j[f0.t<=tau.u,i]),
                   (f0.t[f0.t>=tau.l]>ci.s$ub.F0j[f0.t>=tau.l,i])*(tau.u>ci.s$ub.F0j[f0.t>=tau.l,i]))
    results[i,6] <- !any((f1.t[f1.t<=tau.u]<ci.s$lb.F1j[f1.t<=tau.u,i])*(tau.l<ci.s$lb.F1j[f1.t<=tau.u,i]),
                   (f1.t[f1.t>=tau.l]>ci.s$ub.F1j[f1.t>=tau.l,i])*(tau.u>ci.s$ub.F1j[f1.t>=tau.l,i]))
    results[i,7] <- 1-((sum(q0.t+1e-12<ci.s$lb.Q0j[,i])+sum(q0.t-1e-12>ci.s$ub.Q0j[,i]))>0)
    results[i,8] <- 1-((sum(q1.t+1e-12<ci.s$lb.Q1j[,i])+sum(q1.t-1e-12>ci.s$ub.Q1j[,i]))>0)
    results[i,9] <- 1-((sum(qte.t+1e-12<ci.s$lb.qte[[i]])+sum(qte.t-1e-12>ci.s$ub.qte[[i]]))>0)
    results[i,10] <- results[i,5]*results[i,6]
    results[i,11] <- results[i,7]*results[i,8]
    results[i,12] <- results[i,5]*results[i,6]*results[i,7]*results[i,8]*results[i,9]
  }
  
  if(all==0){
    #power
    for(i in 1:ncol(ci.s$lb.F0)){
      results[i,13] <- ((sum(0<ci.s$lb.qte[[i]])+sum(0>ci.s$ub.qte[[i]]))>0)
      results[i,14] <- ci.s$max.length[i]
      results[i,15] <- ci.s$mean.length[i]
    } 

    
    colnames(results) <- c("cov.F0", "cov.F1", "cov.Q0", "cov.Q1", "cov.F0.j", "cov.F1.j", "cov.Q0.j", "cov.Q1.j", "cov.qte", "cov.all.F", "cov.all.Q", "cov.all", "power", "max length", "mean length")
    return(results)

  } else {
    #Bootstrapping the quantiles: equal width CI
    q0.t <- q0.t.func(tau.grid)
    q1.t <- q1.t.func(tau.grid)  
    qte.t <- q1.t-q0.t
    
    for(i in 1:ncol(ci.s$lb.F0)){
      results[i,13] <- 1-((sum(q0.t<ci.s$lb.Q0.ew[,i])+sum(q0.t>ci.s$ub.Q0.ew[,i]))>0)
      results[i,14] <- 1-((sum(q1.t<ci.s$lb.Q1.ew[,i])+sum(q1.t>ci.s$ub.Q1.ew[,i]))>0)
      results[i,15] <- 1-((sum(qte.t<ci.s$lb.qte.ew[,i])+sum(qte.t>ci.s$ub.qte.ew[,i]))>0)
    
    #Bootstrapping the quantiles: weights
      results[i,16] <- 1-((sum(q0.t<ci.s$lb.Q0.w[,i])+sum(q0.t>ci.s$ub.Q0.w[,i]))>0)
      results[i,17] <- 1-((sum(q1.t<ci.s$lb.Q1.w[,i])+sum(q1.t>ci.s$ub.Q1.w[,i]))>0)
      results[i,18] <- 1-((sum(qte.t<ci.s$lb.qte.w[,i])+sum(qte.t>ci.s$ub.qte.w[,i]))>0)
    
    #smoothing centered around smoothed QF
      results[i,19] <- 1-((sum(q0.t<ci.s$lb.Q0.ss[,i])+sum(q0.t>ci.s$ub.Q0.ss[,i]))>0)
      results[i,20] <- 1-((sum(q1.t<ci.s$lb.Q1.ss[,i])+sum(q1.t>ci.s$ub.Q1.ss[,i]))>0)
      results[i,21] <- 1-((sum(qte.t<ci.s$lb.qte.ss[,i])+sum(qte.t>ci.s$ub.qte.ss[,i]))>0)
  
    #smoothing centered around unsmoothed QF
      results[i,22] <- 1-((sum(q0.t<ci.s$lb.Q0.s[,i])+sum(q0.t>ci.s$ub.Q0.s[,i]))>0)
      results[i,23] <- 1-((sum(q1.t<ci.s$lb.Q1.s[,i])+sum(q1.t>ci.s$ub.Q1.s[,i]))>0)
      results[i,24] <- 1-((sum(qte.t<ci.s$lb.qte.s[,i])+sum(qte.t>ci.s$ub.qte.s[,i]))>0)
  
    # Width of the CB
      results[i,25] <- ci.s$max.length[i]
      results[i,26] <- max(abs(ci.s$ub.qte.ew[,i]-ci.s$lb.qte.ew[,i]))
      results[i,27] <- max(abs(ci.s$ub.qte.w[,i]-ci.s$lb.qte.w[,i]))
      results[i,28] <- max(abs(ci.s$ub.qte.s[,i]-ci.s$lb.qte.s[,i]))
      results[i,29] <- max(abs(ci.s$ub.qte.ss[,i]-ci.s$lb.qte.ss[,i]))
      results[i,30] <- ci.s$mean.length[i]
      results[i,31] <- mean(abs(ci.s$ub.qte.ew[,i]-ci.s$lb.qte.ew[,i]))
      results[i,32] <- mean(abs(ci.s$ub.qte.w[,i]-ci.s$lb.qte.w[,i]))
      results[i,33] <- mean(abs(ci.s$ub.qte.s[,i]-ci.s$lb.qte.s[,i]))
      results[i,34] <- mean(abs(ci.s$ub.qte.ss[,i]-ci.s$lb.qte.ss[,i]))
    
    # Power to reject 0 effect
      results[i,35] <- ((sum(0<ci.s$lb.qte[[i]])+sum(0>ci.s$ub.qte[[i]]))>0)
      results[i,36] <- ((sum(0<ci.s$lb.qte.ew[,i])+sum(0>ci.s$ub.qte.ew[,i]))>0)
      results[i,37] <- ((sum(0<ci.s$lb.qte.w[,i])+sum(0>ci.s$ub.qte.w[,i]))>0)
      results[i,38] <- ((sum(0<ci.s$lb.qte.s[,i])+sum(0>ci.s$ub.qte.s[,i]))>0)
      results[i,39] <- ((sum(0<ci.s$lb.qte.ss[,i])+sum(0>ci.s$ub.qte.ss[,i]))>0)
    }
    # Results
    colnames(results) <- c("cov.F0", "cov.F1", "cov.Q0", "cov.Q1", "cov.F0.j", "cov.F1.j", "cov.Q0.j", "cov.Q1.j", "cov.qte", "cov.all.F", "cov.all.Q", "cov.all", "cov.Q0.ew", "cov.Q1.ew", "cov.qte.ew", "cov.Q0.w", "cov.Q1.w", "cov.qte.w", "cov.Q0.s", "cov.Q1.s", "cov.qte.s", "cov.Q0.ss", "cov.Q1.ss", "cov.qte.ss", "max.length", "max.length.ew","max.length.w", "max.length.s", "max.length.ss", "mean.length", "mean.length.ew", "mean.length.w", "mean.length.s", "mean.length.ss", "power", "power.ew", "power.w", "power.s", "power.ss")
    return(results)
  }
  
}

sim.poisson = function(n, lambda0, lambda1, bsrep, alpha, ys.0, ys.1, taus, taus.grid, f0.t, f1.t, q0.t, q1.t, all=1){
  
  # DGP
  y0 = rpois(n,lambda0)
  y1 = rpois(n,lambda1)

  # Calculate the bands
  ci.s = ci(y0, y1, ys.0, ys.1, bsrep, alpha, taus, taus.grid, all)
  
  #evaluate the bands
  evaluate(ci.s, taus, taus.grid, f0.t, f1.t, q0.t, q1.t, all)
  
}

sim.ordered = function(n, mu1, nc, bsrep, alpha, ys.0, ys.1, taus, taus.grid, f0.t, f1.t, q0.t, q1.t, all=1){
  
  # DGP
  cut <- seq(from=qnorm(0.5/nc), to=qnorm(1-0.5/nc), length.out=nc)
  y0 = rnorm(n)
  y0 <- sapply(1:n, function(i) sum(cut<y0[i]))
  y1 = rnorm(n, mu1)
  y1 <- sapply(1:n, function(i) sum(cut<y1[i]))
  
  # Calculate the bands
  ci.s = ci(y0, y1, ys.0, ys.1, bsrep, alpha, taus, taus.grid, all)
  #evaluate the bands
  evaluate(ci.s, taus, taus.grid, f0.t, f1.t, q0.t, q1.t, all)
  
}


###############
# Simulations
###############

registerDoParallel(cores = 7)

bsrep = 399
alpha = c(0.01,0.05,0.1)
taus.grid = seq(0.1,0.9,0.01)

#Poisson
parameter.mat <- cbind(3, rep(c(3,2.75,2.5),3), rep(c(400, 1600, 6400), each=3), 5005, c(rep(1,6), rep(0,3)))
results.poisson = c()
for (l in 1:dim(parameter.mat)[1]){
  set.seed(123456)
  #true values
  #values at which we evaluate the cdf
  ys.0 <- unique(qpois(taus.grid, parameter.mat[l,1]))
  if(min(ys.0)>0) ys.0 <- c(min(ys.0)-1, ys.0)
  ys.1 <- unique(qpois(taus.grid, parameter.mat[l,2]))
  if(min(ys.1)>0) ys.1 <- c(min(ys.1)-1, ys.1)
  #true cdf
  f0.t = ppois(ys.0, parameter.mat[l,1])
  f1.t = ppois(ys.1, parameter.mat[l,2])
  #quantiles at which we evaluate the qfs (chosen such that the whole true function in the range defined by taus is tested)
  tau <- sort(unique(pmax(pmin(c(min(taus.grid),f0.t-1e-6, f0.t+1e-6, f1.t-1e-6, f1.t+1e-6), max(taus.grid)), min(taus.grid))))
  #true QFs
  q0.t <- function(q) qpois(q, parameter.mat[l,1])
  q1.t <- function(q) qpois(q, parameter.mat[l,2])  
  #simulations
  results.par = foreach(i=1:parameter.mat[l,4]) %dorng% {
    sim.poisson(parameter.mat[l,3], parameter.mat[l,1], parameter.mat[l,2], bsrep, alpha, ys.0, ys.1, tau, taus.grid, f0.t, f1.t, q0.t, q1.t, all=parameter.mat[l,5])
  }
  if(parameter.mat[l,5]==0){
    results.poisson[[l]] = matrix(rowMeans(matrix(unlist(results.par), 45)),length(alpha))
  } else{
    results.poisson[[l]] = matrix(rowMeans(matrix(unlist(results.par), 117)),length(alpha))
  }
  print(results.poisson[[l]])
}

# Ordered normal distribution
parameter.mat = cbind(rep(c(0, 0.2, 0.4),3), 5, rep(c(400, 1600, 6400), each=3), 5005, c(rep(1,6),rep(0,3)))
taus.grid <- seq(0.01,0.99,0.01)
results.ordered <- c()
for (l in 1:dim(parameter.mat)[1]){
  set.seed(123456)
  #true values
  #values at which we evaluate the cdf
  nc=parameter.mat[l,2]
  cut <- c(seq(from=qnorm(0.5/nc), to=qnorm(1-0.5/nc), length.out=nc), Inf)
  temp.0 <- pnorm(cut)
  ys.0 <- unique(sapply(1:length(taus.grid), function(i) (0:nc)[sum(temp.0<=taus.grid[i])+1]))
  if(min(ys.0)>0) ys.0 <- c(min(ys.0)-1, ys.0)
  temp.1 <- pnorm(cut, parameter.mat[l,1])
  ys.1 <- unique(sapply(1:length(taus.grid), function(i) (0:nc)[sum(temp.1<=taus.grid[i])+1]))
  if(min(ys.1)>0) ys.1 <- c(min(ys.1)-1, ys.1)
  #true cdf
  f0.t = pnorm(cut[ys.0+1])
  f1.t = pnorm(cut[ys.1+1], parameter.mat[l,1])
  #quantiles at which we evaluate the qfs (chosen such that the whole true function in the range defined by taus is tested)
  tau <- sort(unique(pmax(pmin(c(min(taus.grid),f0.t-1e-6, f0.t+1e-6, f1.t-1e-6, f1.t+1e-6), max(taus.grid)), min(taus.grid))))
  #true QFs
  q0.t <- left.inv(ys.0, f0.t)
  q1.t <- left.inv(ys.1, f1.t)
  results.par = foreach(i=1:parameter.mat[l,4]) %dorng% {
    sim.ordered(parameter.mat[l,3], parameter.mat[l,1], parameter.mat[l,2], bsrep, alpha, ys.0, ys.1, tau, taus.grid, f0.t, f1.t, q0.t, q1.t, all=parameter.mat[l,5])
  }
  if(parameter.mat[l,5]==0){
    results.ordered[[l]] = matrix(rowMeans(matrix(unlist(results.par), 45)),length(alpha))
  } else{
    results.ordered[[l]] = matrix(rowMeans(matrix(unlist(results.par), 117)),length(alpha))
  }
  print(results.ordered[[l]])
}

#Tables
table1 <- rbind(
  cbind(rep(c(400, 1600, 6400), each=3), rep(c(0.99,0.95,0.9),3), rbind(results.poisson[[1]][,c(3,4,12,9,35)], results.poisson[[4]][,c(3,4,12,9,35)], results.poisson[[7]][,c(3,4,12,9,13)])),
  cbind(rep(c(400, 1600, 6400), each=3), rep(c(0.99,0.95,0.9),3), rbind(results.poisson[[2]][,c(3,4,12,9,35)], results.poisson[[5]][,c(3,4,12,9,35)], results.poisson[[8]][,c(3,4,12,9,13)])),
  cbind(rep(c(400, 1600, 6400), each=3), rep(c(0.99,0.95,0.9),3), rbind(results.poisson[[3]][,c(3,4,12,9,35)], results.poisson[[6]][,c(3,4,12,9,35)], results.poisson[[9]][,c(3,4,12,9,13)]))
)
colnames(table1) <- c("n", "p", "cov.Q0", "cov.Q1", "cov.all", "cov.qte", "Reject 0")
print(xtable(table1, caption="Poisson"), include.rownames = FALSE)

table2 <- rbind(
  cbind(rep(c(400, 1600, 6400), each=3), rep(c(0.99,0.95,0.9),3), rbind(results.ordered[[1]][,c(3,4,12,9,35)], results.ordered[[4]][,c(3,4,12,9,35)], results.ordered[[7]][,c(3,4,12,9,13)])),
  cbind(rep(c(400, 1600, 6400), each=3), rep(c(0.99,0.95,0.9),3), rbind(results.ordered[[2]][,c(3,4,12,9,35)], results.ordered[[5]][,c(3,4,12,9,35)], results.ordered[[8]][,c(3,4,12,9,13)])),
  cbind(rep(c(400, 1600, 6400), each=3), rep(c(0.99,0.95,0.9),3), rbind(results.ordered[[3]][,c(3,4,12,9,35)], results.ordered[[6]][,c(3,4,12,9,35)], results.ordered[[9]][,c(3,4,12,9,13)]))
)
colnames(table2) <- c("n", "p", "cov.Q0", "cov.Q1", "cov.all", "cov.qte", "Reject 0")
print(xtable(table2, caption="Ordered"), include.rownames = FALSE)

table3 <- rbind(
  cbind(rep(c(400, 1600), each=3), rep(c(0.99,0.95,0.9),2), rbind(results.poisson[[1]][,c(9, 15, 21, 24, 30, 31,  33:34)], results.poisson[[4]][,c(9, 15, 21, 24, 30, 31, 33:34)])),
  cbind(rep(c(400, 1600), each=3), rep(c(0.99,0.95,0.9),2), rbind(results.poisson[[2]][,c(9, 15, 21, 24, 30, 31,  33:34)], results.poisson[[5]][,c(9, 15, 21, 24, 30, 31, 33:34)])),
  cbind(rep(c(400, 1600), each=3), rep(c(0.99,0.95,0.9),2), rbind(results.poisson[[3]][,c(9, 15, 21, 24, 30, 31,  33:34)], results.poisson[[6]][,c(9, 15, 21, 24, 30, 31, 33:34)]))
)
colnames(table3) <- c("n", "p", "cov", "cov.ew", "cov.j1", "cov.j2", "l", "l.ew", "l.j1", "l.j2")
print(xtable(table3, digits=c(0,0,2,2,2,2,2,2,2,2,2), caption="Comparison Poisson"), include.rownames = FALSE)

table4 <- rbind(
  cbind(rep(c(400, 1600), each=3), rep(c(0.99,0.95,0.9),2), rbind(results.ordered[[1]][,c(9, 15, 21, 24, 30, 31,  33:34)], results.ordered[[4]][,c(9, 15, 21, 24, 30, 31, 33:34)])),
  cbind(rep(c(400, 1600), each=3), rep(c(0.99,0.95,0.9),2), rbind(results.ordered[[2]][,c(9, 15, 21, 24, 30, 31,  33:34)], results.ordered[[5]][,c(9, 15, 21, 24, 30, 31, 33:34)])),
  cbind(rep(c(400, 1600), each=3), rep(c(0.99,0.95,0.9),2), rbind(results.ordered[[3]][,c(9, 15, 21, 24, 30, 31,  33:34)], results.ordered[[6]][,c(9, 15, 21, 24, 30, 31, 33:34)]))
)
colnames(table4) <- c("n", "p", "cov", "cov.ew", "cov.j1", "cov.j2", "l", "l.ew", "l.j1", "l.j2")
print(xtable(table4, digits=c(0,0,2,2,2,2,2,2,2,2,2), caption="Comparison Ordered"), include.rownames = FALSE)

  