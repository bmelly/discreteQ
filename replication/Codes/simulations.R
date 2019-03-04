####################################################################################################

# Monte Carlo Simulations

# Paper: Generic Inference on Quantile and Quantile Effect Functions for Discrete Outcomes

# Authors: V. Chernozhukov, I. Fernandez-Val, B. Melly, K. Wuethrich

####################################################################################################

##############################
#Definition of the functions
##############################

##Function that calculate the confidence bands using different methods
ci <-
  function(y0, y1,
           ycdf.0, ycdf.1,
           bsrep,
           alphas,
           taus,
           taus.grid,
           all) {
#Number of observations    
    n.0 <- length(y0)
    n.1 <- length(y1)
#Values at which the CDF will be estimated
    ys0 <- sort(unique(y0))
    ys1 <- sort(unique(y1))
    ys <- sort(unique(c(ys0, ys1)))
#Definition of the matrices that will contain the bands    
    ub.F0 <-
      lb.F0 <-
      ub.F0j <- lb.F0j <- matrix(NA, length(ycdf.0), length(alphas))
    ub.F1 <-
      lb.F1 <-
      ub.F1j <- lb.F1j <- matrix(NA, length(ycdf.1), length(alphas))
    ub.Q0 <-
      lb.Q0 <-
      ub.Q1 <-
      lb.Q1 <-
      ub.Q0j <-
      lb.Q0j <-
      ub.Q1j <- lb.Q1j <-  matrix(NA, length(taus), length(alphas))
    temp.taus <- ub.qte <- lb.qte <- list()
    mean.length <- max.length <- rep(NA, length(alphas))
#Bands proposed in the paper
#Loop over the confidence levels    
    for (i in 1:length(alphas)) {
      if (i == 1) {
        joint.res <-
          keep <-
          keep0.res <-
          keep1.res <-
          discreteQ(
            c(y0, y1),
            c(rep(0, n.0), rep(1, n.1)),
            q.range = range(taus.grid),
            bsrep = bsrep,
            alpha = alphas[i],
            ys = ys,
            return.boot = TRUE
          )
        keep0.res$F <- keep0.res$F0
        keep0.res$Q <- keep0.res$Q0
        keep0.res$ys <- keep0.res$ys0
        keep0.res$F.b <- keep0.res$F.b[1:length(ys0), ]
        keep1.res$F <- keep1.res$F1
        keep1.res$Q <- keep1.res$Q1
        keep1.res$ys <- keep1.res$ys1
        keep1.res$F.b <-
          keep1.res$F.b[(length(ys0) + 1):(length(ys0) + length(ys1)), ]
      } else{
        joint.res <-
          discreteQ(
            c(y0, y1),
            c(rep(0, n.0), rep(1, n.1)),
            q.range = range(taus.grid),
            bsrep = bsrep,
            alpha = alphas[i],
            ys = ys,
            old.res = keep
          )
      }
#Bounds for a single outcome      
      univariate0.res <-
        discreteQ(
          y0,
          q.range = range(taus.grid),
          bsrep = bsrep,
          alpha = alphas[i],
          ys = ys0,
          old.res = keep0.res
        )
      univariate1.res <-
        discreteQ(
          y1,
          q.range = range(taus.grid),
          bsrep = bsrep,
          alpha = alphas[i],
          ys = ys1,
          old.res = keep1.res
        )
      
      ub.Q0[, i] <-  univariate0.res$ub.Q(taus)
      lb.Q0[, i] <- univariate0.res$lb.Q(taus)
      ub.Q1[, i] <- univariate1.res$ub.Q(taus)
      lb.Q1[, i] <- univariate1.res$lb.Q(taus)
      ub.F0[, i] <- univariate0.res$ub.F(ycdf.0)
      lb.F0[, i] <- univariate0.res$lb.F(ycdf.0)
      ub.F1[, i] <- univariate1.res$ub.F(ycdf.1)
      lb.F1[, i] <- univariate1.res$lb.F(ycdf.1)
#Bands that cover both the treated and control outcomes      
      ub.Q0j[, i] <- joint.res$ub.Q0(taus)
      lb.Q0j[, i] <- joint.res$lb.Q0(taus)
      ub.Q1j[, i] <- joint.res$ub.Q1(taus)
      lb.Q1j[, i] <- joint.res$lb.Q1(taus)
      ub.F0j[, i] <- joint.res$ub.F0(ycdf.0)
      lb.F0j[, i] <- joint.res$lb.F0(ycdf.0)
      ub.F1j[, i] <- joint.res$ub.F1(ycdf.1)
      lb.F1j[, i] <- joint.res$lb.F1(ycdf.1)
#Values at which the bands for the QE function jump      
      temp.taus[[i]] <-
        unique(sort(pmax(pmin(
          c(
            taus,
            environment(joint.res$lb.F0)$y + 1e-6,
            environment(joint.res$lb.F0)$y - 1e-6,
            environment(joint.res$ub.F0)$y + 1e-6,
            environment(joint.res$ub.F0)$y - 1e-6,
            environment(joint.res$lb.F1)$y + 1e-6,
            environment(joint.res$lb.F1)$y - 1e-6,
            environment(joint.res$ub.F1)$y + 1e-6,
            environment(joint.res$ub.F1)$y - 1e-6
          ),
          max(taus)
        ), min(taus))))
#Bands for the QE function      
      lb.qte[[i]] <- joint.res$lb.QTE(temp.taus[[i]])
      ub.qte[[i]] <- joint.res$ub.QTE(temp.taus[[i]])
      #length of the bands
      mean.length[i] <-
        mean(joint.res$ub.QTE(taus.grid) - joint.res$lb.QTE(taus.grid))
      max.length[i] <-
        max(joint.res$ub.QTE(taus.grid) - joint.res$lb.QTE(taus.grid))
    }
    
    if (all == 0) {
      return(
        list(
          ub.F0 = ub.F0,
          lb.F0 = lb.F0,
          ub.F1 = ub.F1,
          lb.F1 = lb.F1,
          ub.Q0 = ub.Q0,
          lb.Q0 = lb.Q0,
          ub.Q1 = ub.Q1,
          lb.Q1 = lb.Q1,
          ub.F0j = ub.F0j,
          lb.F0j = lb.F0j,
          ub.F1j = ub.F1j,
          lb.F1j = lb.F1j,
          ub.Q0j = ub.Q0j,
          lb.Q0j = lb.Q0j,
          ub.Q1j = ub.Q1j,
          lb.Q1j = lb.Q1j,
          ub.qte = ub.qte,
          lb.qte = lb.qte,
          temp.taus = temp.taus,
          mean.length = mean.length,
          max.length = max.length
        )
      )
    } else{
### Bootstrapping the quantiles: equal width CI
      Q1.b  <-
        Q0.b <-  qte.b <-  matrix(NA, length(taus.grid), bsrep)
      Q0 <- keep0.res$Q(taus.grid)
      Q1 <- keep1.res$Q(taus.grid)
      qte <- Q1 - Q0
      
      for (r in 1:bsrep) {
        Q0.b.r.func <-
          stepfun(keep0.res$F.b[, r], c(keep0.res$ys, max(keep0.res$ys)), right =
                    TRUE)
        Q1.b.r.func <-
          stepfun(keep1.res$F.b[, r], c(keep1.res$ys, max(keep1.res$ys)), right =
                    TRUE)
        Q0.b[, r]  <-  Q0.b.r.func(taus.grid)
        Q1.b[, r]  <-  Q1.b.r.func(taus.grid)
        qte.b[, r]  <-  Q1.b[, r] - Q0.b[, r]
      }
      
      delta.qte <- qte.b - qte
      
      zs.qte <- apply(rbind(abs(delta.qte)), 2, max, na.rm = TRUE)
      
      se.qte <- pmax(apply(
        qte.b,
        1,
        FUN = function(x)
          IQR(x) / 1.349
      ), 1e-6)
      
      zs.qte.w <- apply(rbind(abs(delta.qte) / se.qte),  2, max, na.rm = TRUE)
      
      ub.qte.ew <-
        lb.qte.ew <- matrix(NA, length(taus.grid), length(alphas))
      ub.qte.s <-
        lb.qte.s <- matrix(NA, length(taus.grid), length(alphas))
      ub.qte.ss <-
        lb.qte.ss <- matrix(NA, length(taus.grid), length(alphas))
      
      for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        crt.qte <- quantile(zs.qte, 1 - alpha)
        
        lb.qte.ew[, i] <- qte - crt.qte
        ub.qte.ew[, i] <- qte + crt.qte
      }
      
##Adding noise and bootstrapping
      y.smooth.0 <- y0 + runif(n.0)
      y.smooth.1 <- y1 + runif(n.1)
      Q.smooth.0  <-  quantile(y.smooth.0, taus.grid, type = 1)
      Q.smooth.1  <-  quantile(y.smooth.1, taus.grid, type = 1)
      qte.smooth  <-  Q.smooth.1 - Q.smooth.0
      Q1.b  <-
        Q0.b  <-  qte.b  <-  matrix(NA, length(taus.grid), bsrep)
      for (r in 1:bsrep) {
        w.0 <-  rexp(n.0)
        w.1 <-  rexp(n.1)
        Q0.b[, r] <-
          wtd.quantile(y.smooth.0,
                       type = "i/(n+1)",
                       weights = w.0,
                       probs = taus.grid)
        Q1.b[, r] <-
          wtd.quantile(y.smooth.1,
                       type = "i/(n+1)",
                       weights = w.1,
                       probs = taus.grid)
        qte.b[, r] <-  Q1.b[, r] - Q0.b[, r]
      }
      se.qte <-  apply(
        qte.b,
        1,
        FUN = function(x)
          IQR(x) / 1.349
      )
      delta.qte <-  qte.b - qte.smooth
      zs.qte <-  apply(abs(delta.qte) / se.qte,  2, max, na.rm = TRUE)
      for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        crt.qte <-  quantile(zs.qte, 1 - alpha)
        lb.qte.s[, i] <-  qte - crt.qte * se.qte
        ub.qte.s[, i] <-  qte + crt.qte * se.qte
        lb.qte.ss[, i] <-  qte.smooth - crt.qte * se.qte
        ub.qte.ss[, i] <-  qte.smooth + crt.qte * se.qte
      }
## Results
      return(
        list(
          ub.F0 = ub.F0,
          lb.F0 = lb.F0,
          ub.F1 = ub.F1,
          lb.F1 = lb.F1,
          ub.Q0 = ub.Q0,
          lb.Q0 = lb.Q0,
          ub.Q1 = ub.Q1,
          lb.Q1 = lb.Q1,
          ub.F0j = ub.F0j,
          lb.F0j = lb.F0j,
          ub.F1j = ub.F1j,
          lb.F1j = lb.F1j,
          ub.Q0j = ub.Q0j,
          lb.Q0j = lb.Q0j,
          ub.Q1j = ub.Q1j,
          lb.Q1j = lb.Q1j,
          ub.qte = ub.qte,
          lb.qte = lb.qte,
          lb.qte.ew = lb.qte.ew,
          ub.qte.ew = ub.qte.ew,
          lb.qte.s = lb.qte.s,
          ub.qte.s = ub.qte.s,
          lb.qte.ss = lb.qte.ss,
          ub.qte.ss = ub.qte.ss,
          temp.taus = temp.taus,
          mean.length = mean.length,
          max.length = max.length
        )
      )
    }
  }

#Function that evaluates the different confidence bands
evaluate <-
  function(ci.s,
           taus,
           tau.grid,
           f0.t,
           f1.t,
           q0.t.func,
           q1.t.func,
           all) {
#Range of quantiles    
    tau.l <- min(taus)
    tau.u <- max(taus)
#True values
    q0.t <- q0.t.func(taus)
    q1.t <- q1.t.func(taus)
#Definition of the matrix that will contain the results
    if (all == 0){
      results <- matrix(NA, ncol(ci.s$lb.F0), 15)
    } else {
      results <- matrix(NA, ncol(ci.s$lb.F0), 39)
    }
#Loop over the different confidence levels    
    for (i in 1:ncol(ci.s$lb.F0)) {
      qte.t <-
        q1.t.func(ci.s$temp.taus[[i]]) - q0.t.func(ci.s$temp.taus[[i]])
#Results for the bands for a single function
      results[i, 1] <-
        !any((f0.t[f0.t <= tau.u] < ci.s$lb.F0[f0.t <= tau.u, i]) * (tau.l < ci.s$lb.F0[f0.t <=
                                                                                          tau.u, i]),
             (f0.t[f0.t >= tau.l] > ci.s$ub.F0[f0.t >= tau.l, i]) *
               (tau.u > ci.s$ub.F0[f0.t >= tau.l, i])
        )
      results[i, 2] <-
        !any((f1.t[f1.t <= tau.u] < ci.s$lb.F1[f1.t <= tau.u, i]) * (tau.l < ci.s$lb.F1[f1.t <=
                                                                                          tau.u, i]),
             (f1.t[f1.t >= tau.l] > ci.s$ub.F1[f1.t >= tau.l, i]) *
               (tau.u > ci.s$ub.F1[f1.t >= tau.l, i])
        )
      results[i, 3] <-
        1 - ((
          sum(q0.t + 1e-12 < ci.s$lb.Q0[, i]) + sum(q0.t - 1e-12 > ci.s$ub.Q0[, i])
        ) > 0)
      results[i, 4] <-
        1 - ((
          sum(q1.t + 1e-12 < ci.s$lb.Q1[, i]) + sum(q1.t - 1e-12 > ci.s$ub.Q1[, i])
        ) > 0)
#Results for the bands that cover both outcomes
      results[i, 5] <-
        !any((f0.t[f0.t <= tau.u] < ci.s$lb.F0j[f0.t <= tau.u, i]) * (tau.l < ci.s$lb.F0j[f0.t <=
                                                                                            tau.u, i]),
             (f0.t[f0.t >= tau.l] > ci.s$ub.F0j[f0.t >= tau.l, i]) *
               (tau.u > ci.s$ub.F0j[f0.t >= tau.l, i])
        )
      results[i, 6] <-
        !any((f1.t[f1.t <= tau.u] < ci.s$lb.F1j[f1.t <= tau.u, i]) * (tau.l < ci.s$lb.F1j[f1.t <=
                                                                                            tau.u, i]),
             (f1.t[f1.t >= tau.l] > ci.s$ub.F1j[f1.t >= tau.l, i]) *
               (tau.u > ci.s$ub.F1j[f1.t >= tau.l, i])
        )
      results[i, 7] <-
        1 - ((
          sum(q0.t + 1e-12 < ci.s$lb.Q0j[, i]) + sum(q0.t - 1e-12 > ci.s$ub.Q0j[, i])
        ) > 0)
      results[i, 8] <-
        1 - ((
          sum(q1.t + 1e-12 < ci.s$lb.Q1j[, i]) + sum(q1.t - 1e-12 > ci.s$ub.Q1j[, i])
        ) > 0)
      results[i, 9] <-
        1 - ((
          sum(qte.t + 1e-12 < ci.s$lb.qte[[i]]) + sum(qte.t - 1e-12 > ci.s$ub.qte[[i]])
        ) > 0)
      results[i, 10] <- results[i, 5] * results[i, 6]
      results[i, 11] <- results[i, 7] * results[i, 8]
      results[i, 12] <-
        results[i, 5] * results[i, 6] * results[i, 7] * results[i, 8] * results[i, 9]
    }
    
    if (all == 0) {
#Power of the suggested bands to reject 0 QE
      for (i in 1:ncol(ci.s$lb.F0)) {
        results[i, 13] <-
          ((sum(0 < ci.s$lb.qte[[i]]) + sum(0 > ci.s$ub.qte[[i]])) > 0)
        results[i, 14] <- ci.s$max.length[i]
        results[i, 15] <- ci.s$mean.length[i]
      }
      colnames(results) <-
        c(
          "cov.F0",
          "cov.F1",
          "cov.Q0",
          "cov.Q1",
          "cov.F0.j",
          "cov.F1.j",
          "cov.Q0.j",
          "cov.Q1.j",
          "cov.qte",
          "cov.all.F",
          "cov.all.Q",
          "cov.all",
          "power",
          "max length",
          "mean length"
        )
      return(results)
    } else {
#Bootstrapping the quantiles: equal width CI
      q0.t <- q0.t.func(tau.grid)
      q1.t <- q1.t.func(tau.grid)
      qte.t <- q1.t - q0.t
      
      for (i in 1:ncol(ci.s$lb.F0)) {
        results[i, 13] <-
          1 - ((sum(q0.t < ci.s$lb.Q0.ew[, i]) + sum(q0.t > ci.s$ub.Q0.ew[, i])) >
                 0)
        results[i, 14] <-
          1 - ((sum(q1.t < ci.s$lb.Q1.ew[, i]) + sum(q1.t > ci.s$ub.Q1.ew[, i])) >
                 0)
        results[i, 15] <-
          1 - ((
            sum(qte.t < ci.s$lb.qte.ew[, i]) + sum(qte.t > ci.s$ub.qte.ew[, i])
          ) > 0)
        
#smoothing centered around smoothed QF
        results[i, 19] <-
          1 - ((sum(q0.t < ci.s$lb.Q0.ss[, i]) + sum(q0.t > ci.s$ub.Q0.ss[, i])) >
                 0)
        results[i, 20] <-
          1 - ((sum(q1.t < ci.s$lb.Q1.ss[, i]) + sum(q1.t > ci.s$ub.Q1.ss[, i])) >
                 0)
        results[i, 21] <-
          1 - ((
            sum(qte.t < ci.s$lb.qte.ss[, i]) + sum(qte.t > ci.s$ub.qte.ss[, i])
          ) > 0)
        
#smoothing centered around unsmoothed QF
        results[i, 22] <-
          1 - ((sum(q0.t < ci.s$lb.Q0.s[, i]) + sum(q0.t > ci.s$ub.Q0.s[, i])) >
                 0)
        results[i, 23] <-
          1 - ((sum(q1.t < ci.s$lb.Q1.s[, i]) + sum(q1.t > ci.s$ub.Q1.s[, i])) >
                 0)
        results[i, 24] <-
          1 - ((sum(qte.t < ci.s$lb.qte.s[, i]) + sum(qte.t > ci.s$ub.qte.s[, i])) >
                 0)
        
# Width of the confidence bands
        results[i, 25] <- ci.s$max.length[i]
        results[i, 26] <-
          max(abs(ci.s$ub.qte.ew[, i] - ci.s$lb.qte.ew[, i]))
        results[i, 28] <-
          max(abs(ci.s$ub.qte.s[, i] - ci.s$lb.qte.s[, i]))
        results[i, 29] <-
          max(abs(ci.s$ub.qte.ss[, i] - ci.s$lb.qte.ss[, i]))
        results[i, 30] <- ci.s$mean.length[i]
        results[i, 31] <-
          mean(abs(ci.s$ub.qte.ew[, i] - ci.s$lb.qte.ew[, i]))
        results[i, 33] <-
          mean(abs(ci.s$ub.qte.s[, i] - ci.s$lb.qte.s[, i]))
        results[i, 34] <-
          mean(abs(ci.s$ub.qte.ss[, i] - ci.s$lb.qte.ss[, i]))
        
#Power to reject 0 QE
        results[i, 35] <-
          ((sum(0 < ci.s$lb.qte[[i]]) + sum(0 > ci.s$ub.qte[[i]])) > 0)
        results[i, 36] <-
          ((sum(0 < ci.s$lb.qte.ew[, i]) + sum(0 > ci.s$ub.qte.ew[, i])) > 0)
        results[i, 38] <-
          ((sum(0 < ci.s$lb.qte.s[, i]) + sum(0 > ci.s$ub.qte.s[, i])) > 0)
        results[i, 39] <-
          ((sum(0 < ci.s$lb.qte.ss[, i]) + sum(0 > ci.s$ub.qte.ss[, i])) > 0)
      }
      # Results
      colnames(results) <-
        c(
          "cov.F0",
          "cov.F1",
          "cov.Q0",
          "cov.Q1",
          "cov.F0.j",
          "cov.F1.j",
          "cov.Q0.j",
          "cov.Q1.j",
          "cov.qte",
          "cov.all.F",
          "cov.all.Q",
          "cov.all",
          "cov.Q0.ew",
          "cov.Q1.ew",
          "cov.qte.ew",
          "cov.Q0.w",
          "cov.Q1.w",
          "cov.qte.w",
          "cov.Q0.s",
          "cov.Q1.s",
          "cov.qte.s",
          "cov.Q0.ss",
          "cov.Q1.ss",
          "cov.qte.ss",
          "max.length",
          "max.length.ew",
          "max.length.w",
          "max.length.s",
          "max.length.ss",
          "mean.length",
          "mean.length.ew",
          "mean.length.w",
          "mean.length.s",
          "mean.length.ss",
          "power",
          "power.ew",
          "power.w",
          "power.s",
          "power.ss"
        )
      return(results)
    }
    
  }

#Function that performs one Poisson simulation and evaluate the bands
sim.poisson <- function(n,
                       lambda0, lambda1,
                       bsrep,
                       alpha,
                       ys.0, ys.1,
                       taus,
                       taus.grid,
                       f0.t, f1.t,
                       q0.t, q1.t,
                       all = 1) {
  # DGP
  y0 <- rpois(n, lambda0)
  y1 <- rpois(n, lambda1)
  
  # Calculate the bands
  ci.s <- ci(y0, y1, ys.0, ys.1, bsrep, alpha, taus, taus.grid, all)
  
  #evaluate the bands
  evaluate(ci.s, taus, taus.grid, f0.t, f1.t, q0.t, q1.t, all)
  
}

#Function that performs one ordered simulation and evaluate the bands
sim.ordered <- function(n,
                       mu1,
                       nc,
                       bsrep,
                       alpha,
                       ys.0, ys.1,
                       taus,
                       taus.grid,
                       f0.t, f1.t,
                       q0.t, q1.t,
                       all = 1) {
  # DGP
  cut <- seq(
    from = qnorm(0.5 / nc),
    to = qnorm(1 - 0.5 / nc),
    length.out = nc
  )
  y0 <- rnorm(n)
  y0 <- sapply(1:n, function(i)
    sum(cut < y0[i]))
  y1 <- rnorm(n, mu1)
  y1 <- sapply(1:n, function(i)
    sum(cut < y1[i]))
  
  # Calculate the bands
  ci.s <- ci(y0, y1, ys.0, ys.1, bsrep, alpha, taus, taus.grid, all)
  #evaluate the bands
  evaluate(ci.s, taus, taus.grid, f0.t, f1.t, q0.t, q1.t, all)
  
}

###############
# Simulations
###############

if(!is.null(cl)){
  registerDoParallel(cl)
  clusterEvalQ(cl, library("Hmisc"))
  clusterEvalQ(cl, library("discreteQ"))
}

bsrep <- 399
alpha <- c(0.01, 0.05, 0.1)
taus.grid <- seq(0.1, 0.9, 0.01)

#Poisson
parameter.mat <-
  cbind(3, 
        rep(c(3, 2.75, 2.5), 3), 
        rep(c(400, 1600, 6400), each = 3), 
        5005, 
        c(rep(1, 6), rep(0, 3)))
results.poisson <- c()
for (l in 1:dim(parameter.mat)[1]) {
  set.seed(123456)
  #true values
  #values at which we evaluate the cdf
  ys.0 <- unique(qpois(taus.grid, parameter.mat[l, 1]))
  if (min(ys.0) > 0)
    ys.0 <- c(min(ys.0) - 1, ys.0)
  ys.1 <- unique(qpois(taus.grid, parameter.mat[l, 2]))
  if (min(ys.1) > 0)
    ys.1 <- c(min(ys.1) - 1, ys.1)
  #true cdf
  f0.t <- ppois(ys.0, parameter.mat[l, 1])
  f1.t <- ppois(ys.1, parameter.mat[l, 2])
  #quantiles at which we evaluate the qfs (chosen such that the whole true function in the range defined by taus is tested)
  tau <-
    sort(unique(pmax(pmin(
      c(min(taus.grid), f0.t - 1e-6, f0.t + 1e-6, f1.t - 1e-6, f1.t + 1e-6),
      max(taus.grid)
    ), min(taus.grid))))
  #true QFs
  q0.t <- function(q)
    qpois(q, parameter.mat[l, 1])
  q1.t <- function(q)
    qpois(q, parameter.mat[l, 2])
  #simulations
  results.par <- foreach(i = 1:parameter.mat[l, 4]) %dorng% {
    sim.poisson(
      parameter.mat[l, 3],
      parameter.mat[l, 1],
      parameter.mat[l, 2],
      bsrep,
      alpha,
      ys.0,
      ys.1,
      tau,
      taus.grid,
      f0.t,
      f1.t,
      q0.t,
      q1.t,
      all = parameter.mat[l, 5]
    )
  }
  if (parameter.mat[l, 5] == 0) {
    results.poisson[[l]] <- matrix(rowMeans(matrix(unlist(results.par), 45)), length(alpha))
  } else{
    results.poisson[[l]] <- matrix(rowMeans(matrix(unlist(results.par), 117)), length(alpha))
  }
}

# Ordered normal distribution
parameter.mat <- cbind(rep(c(0, 0.2, 0.4), 3), 
                       5, 
                       rep(c(400, 1600, 6400), each = 3),
                       5005, 
                       c(rep(1, 6), rep(0, 3)))
taus.grid <- seq(0.01, 0.99, 0.01)
results.ordered <- c()
for (l in 1:dim(parameter.mat)[1]) {
  set.seed(123456)
  #true values
  #values at which we evaluate the cdf
  nc <- parameter.mat[l, 2]
  cut <-
    c(seq(
      from = qnorm(0.5 / nc),
      to = qnorm(1 - 0.5 / nc),
      length.out = nc
    ), Inf)
  temp.0 <- pnorm(cut)
  ys.0 <-
    unique(sapply(1:length(taus.grid), function(i)
      (0:nc)[sum(temp.0 <= taus.grid[i]) + 1]))
  if (min(ys.0) > 0)
    ys.0 <- c(min(ys.0) - 1, ys.0)
  temp.1 <- pnorm(cut, parameter.mat[l, 1])
  ys.1 <-
    unique(sapply(1:length(taus.grid), function(i)
      (0:nc)[sum(temp.1 <= taus.grid[i]) + 1]))
  if (min(ys.1) > 0)
    ys.1 <- c(min(ys.1) - 1, ys.1)
  #true cdf
  f0.t <- pnorm(cut[ys.0 + 1])
  f1.t <- pnorm(cut[ys.1 + 1], parameter.mat[l, 1])
  #quantiles at which we evaluate the qfs (chosen such that the whole true function in the range defined by taus is tested)
  tau <-
    sort(unique(pmax(pmin(
      c(min(taus.grid), f0.t - 1e-6, f0.t + 1e-6, f1.t - 1e-6, f1.t + 1e-6),
      max(taus.grid)
    ), min(taus.grid))))
  #true QFs
  q0.t <- stepfun(f0.t, c(ys.0, max(ys.0)), right = TRUE)
  q1.t <- stepfun(f1.t, c(ys.1, max(ys.1)), right = TRUE)
  results.par <- foreach(i = 1:parameter.mat[l, 4]) %dorng% {
    sim.ordered(
      parameter.mat[l, 3],
      parameter.mat[l, 1],
      parameter.mat[l, 2],
      bsrep,
      alpha,
      ys.0,
      ys.1,
      tau,
      taus.grid,
      f0.t,
      f1.t,
      q0.t,
      q1.t,
      all = parameter.mat[l, 5]
    )
  }
  if (parameter.mat[l, 5] == 0) {
    results.ordered[[l]] <- matrix(rowMeans(matrix(unlist(results.par), 45)), length(alpha))
  } else{
    results.ordered[[l]] <- matrix(rowMeans(matrix(unlist(results.par), 117)), length(alpha))
  }
}

#Tables
table1 <- rbind(cbind(
  rep(c(400, 1600, 6400), each = 3),
  rep(c(0.99, 0.95, 0.9), 3),
  rbind(results.poisson[[1]][, c(3, 4, 12, 9, 35)], results.poisson[[4]][, c(3, 4, 12, 9, 35)], results.poisson[[7]][, c(3, 4, 12, 9, 13)])
),
cbind(
  rep(c(400, 1600, 6400), each = 3),
  rep(c(0.99, 0.95, 0.9), 3),
  rbind(results.poisson[[2]][, c(3, 4, 12, 9, 35)], results.poisson[[5]][, c(3, 4, 12, 9, 35)], results.poisson[[8]][, c(3, 4, 12, 9, 13)])
),
cbind(
  rep(c(400, 1600, 6400), each = 3),
  rep(c(0.99, 0.95, 0.9), 3),
  rbind(results.poisson[[3]][, c(3, 4, 12, 9, 35)], results.poisson[[6]][, c(3, 4, 12, 9, 35)], results.poisson[[9]][, c(3, 4, 12, 9, 13)])
))
colnames(table1) <-
  c("n", "p", "cov.Q0", "cov.Q1", "cov.all", "cov.qte", "Reject 0")
print(xtable(table1, caption = "Poisson", digits = c(0, 0, rep(2, 6))),
      include.rownames = FALSE, file = paste0(code.dir, "/Results/Table1.txt"))

table2 <- rbind(cbind(
  rep(c(400, 1600, 6400), each = 3),
  rep(c(0.99, 0.95, 0.9), 3),
  rbind(results.ordered[[1]][, c(3, 4, 12, 9, 35)], results.ordered[[4]][, c(3, 4, 12, 9, 35)], results.ordered[[7]][, c(3, 4, 12, 9, 13)])
),
cbind(
  rep(c(400, 1600, 6400), each = 3),
  rep(c(0.99, 0.95, 0.9), 3),
  rbind(results.ordered[[2]][, c(3, 4, 12, 9, 35)], results.ordered[[5]][, c(3, 4, 12, 9, 35)], results.ordered[[8]][, c(3, 4, 12, 9, 13)])
),
cbind(
  rep(c(400, 1600, 6400), each = 3),
  rep(c(0.99, 0.95, 0.9), 3),
  rbind(results.ordered[[3]][, c(3, 4, 12, 9, 35)], results.ordered[[6]][, c(3, 4, 12, 9, 35)], results.ordered[[9]][, c(3, 4, 12, 9, 13)])
))
colnames(table2) <-
  c("n", "p", "cov.Q0", "cov.Q1", "cov.all", "cov.qte", "Reject 0")
print(xtable(table2, caption = "Ordered", digits = c(0, 0, rep(2, 6))),
      include.rownames = FALSE, file = paste0(code.dir, "/Results/Table2.txt"))

table3 <- rbind(cbind(
  rep(c(400, 1600), each = 3),
  rep(c(0.99, 0.95, 0.9), 2),
  rbind(results.poisson[[1]][, c(9, 15, 21, 24, 30, 31,  33:34)], results.poisson[[4]][, c(9, 15, 21, 24, 30, 31, 33:34)])
),
cbind(
  rep(c(400, 1600), each = 3),
  rep(c(0.99, 0.95, 0.9), 2),
  rbind(results.poisson[[2]][, c(9, 15, 21, 24, 30, 31,  33:34)], results.poisson[[5]][, c(9, 15, 21, 24, 30, 31, 33:34)])
),
cbind(
  rep(c(400, 1600), each = 3),
  rep(c(0.99, 0.95, 0.9), 2),
  rbind(results.poisson[[3]][, c(9, 15, 21, 24, 30, 31,  33:34)], results.poisson[[6]][, c(9, 15, 21, 24, 30, 31, 33:34)])
))
colnames(table3) <-
  c("n",
    "p",
    "cov",
    "cov.ew",
    "cov.j1",
    "cov.j2",
    "l",
    "l.ew",
    "l.j1",
    "l.j2")
print(xtable(
  table3,
  digits = c(0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2),
  caption = "Comparison Poisson"
),
include.rownames = FALSE, file = paste0(code.dir, "/Results/Table3.txt"))

table4 <- rbind(cbind(
  rep(c(400, 1600), each = 3),
  rep(c(0.99, 0.95, 0.9), 2),
  rbind(results.ordered[[1]][, c(9, 15, 21, 24, 30, 31,  33:34)], results.ordered[[4]][, c(9, 15, 21, 24, 30, 31, 33:34)])
),
cbind(
  rep(c(400, 1600), each = 3),
  rep(c(0.99, 0.95, 0.9), 2),
  rbind(results.ordered[[2]][, c(9, 15, 21, 24, 30, 31,  33:34)], results.ordered[[5]][, c(9, 15, 21, 24, 30, 31, 33:34)])
),
cbind(
  rep(c(400, 1600), each = 3),
  rep(c(0.99, 0.95, 0.9), 2),
  rbind(results.ordered[[3]][, c(9, 15, 21, 24, 30, 31,  33:34)], results.ordered[[6]][, c(9, 15, 21, 24, 30, 31, 33:34)])
))
colnames(table4) <-
  c("n",
    "p",
    "cov",
    "cov.ew",
    "cov.j1",
    "cov.j2",
    "l",
    "l.ew",
    "l.j1",
    "l.j2")
print(xtable(
  table4,
  digits = c(0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2),
  caption = "Comparison Ordered"
),
include.rownames = FALSE, file = paste0(code.dir, "/Results/Table4.txt"))
