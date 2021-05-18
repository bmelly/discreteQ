#' Uniform inference on quantile and quantile effect functions for discrete outcomes
#'
#' The function \code{discreteQ} provides uniform confidence bands for the unconditional
#' quantile function, the quantile treatment effect function or the
#' decomposition of the observed difference between the quantile function of an
#' outcome for two groups. This function implements the algorithms suggested in Chernozhukov et al. (2019).
#' See also the vignette available with \code{vignette("discreteQ", package="discreteQ")}.
#'
#' The function \code{discreteQ} can be used in three different ways:
#' \enumerate{
#' \item
#' First, if no treatment variable is specified in the argument \code{d}, then the the
#' command will provide uniform bands for the unconditional quantile and distribution functions of
#' \code{y}.
#' \item
#' Second, if a treatment variable is specified in the argument \code{d} but
#' \code{decomposition=FALSE}, then the command will provide uniform bands that
#' cover the quantile and distribution functions of the treated and control outcomes and the
#' quantile treatment effect function (the difference between both quantile functions).
#' \item Third, if a treatment variable is
#' specified in the argument \code{d} and \code{decomposition=TRUE}, then the command provides
#' uniform bands that cover both distribution functions, both quantile functions as well as the
#' decomposition of their difference into an explained (by the regressors
#' \code{x}) and unexplained component.
#' }
#'
#' The output of the function is a long list of step functions (see below for the details). We recommend to use
#' the \code{\link{plot.discreteQ}} and \code{\link{summary.discreteQ}}
#' functions to analyze the results. For further details see the vignette available with \code{vignette("discreteQ", package="discreteQ")}.
#'
#' @param y outcome (vector of length n).
#' @param d treatment or group variable (binary vector of length n).
#' @param x matrix of regressors (n x p matrix). \code{x} must include a constant if appropriate (i.e. the constant is NOT added automatically). The function \code{\link[stats]{model.matrix}} can be used to create the matrix \code{x} if factor variables or interaction terms are present.
#' @param w sampling weights (vector of length n).
#' @param decomposition logical, indicating if the decomposition of the observed
#'   difference between the control and treated quantile functions should be
#'   performed. By default, inference about the quantile treatment effect
#'   function is performed.
#' @param q.range vector of length 2 that provides the lowest and highest
#'   quantile indexes. The uniform bands will cover the whole QF and QE
#'   function in this range. Default is c(0.05,0.95).
#' @param method link function for the distribution regression model. Possible
#'   values: "logit" (the default), "probit", "cloglog", "lpm" (linear
#'   probability model, i.e. OLS estimation of binary regressions), "cauchit", "drp" (incomplete gamma link function suggested in Chernozhukov, Fernandez-Val, Melly and Wütrich (2019). The maximum likelihood estimator of the fully parametric Poisson regression model is used if \code{method = "poisson"}.
#'   This argument is relevant only if there are regressors in \code{x}.
#' @param bsrep number of bootstrap replications. Default: 200.
#' @param alpha confidence level. Default: 0.05.
#' @param ys specifies the thresholds at which the cumulative distribution
#'   function will be estimated. This argument can be specified either as a
#'   scalar that will be interpreted as the number of thresholds or as a vector
#'   that will contain the values of the thresholds. By default, the cdf is
#'   estimated at all distinct observed values of the outcome in the sample if
#'   there are less than 100 unique values and at the 99 empirical percentiles
#'   of the outcomes if there are more than 100 distinct values.
#' @param cl a cluster object as returned by the function \code{\link[parallel]{makeCluster}}. Parallel computing is
#'   not used if this argument is not specified.
#' @param cluster vector that specifies to which group each observation belongs.
#'   The cluster weighted bootstrap is used if this argument is specified. Otherwise,
#'   simple random sampling is assumed.
#' @param old.res a discreteQ object (obtained with the argument return.boot set
#'   to TRUE). This argument allows for instance to change the size \code{alpha} or the quantile range \code{q.range}
#'   without recomputing the estimates.
#' @param return.boot logical scalar. The results of the bootstrap are return in
#'   the matrix F.b when this argument is set to TRUE.
#' @param list_of_seeds list of seeds for L'Ecuyer RNG. The length of this list
#'   must be the same as the value of the argument \code{bsrep}.
#' @param return.seeds logical scalar. The list of seeds is returned by the
#'   function if this argument is set to TRUE.
#' @param estim.glm function used to estimate the binary regressions if \code{method} is "logit", "probit", "cloglog" or "cauchit". The default is the function fastglm() from the package fastglm. Tested alternatives: glm.fit, glm2::glm.fit2, speedglm::speedglm.wfit.
#' @param par.estim arguments to be passed to the function selected by \code{estim.glm}. For instance, the arguments \code{method}, \code{tol} and \code{maxit} of \code{fastglm} can be set.
#' @return \code{discreteQ} returns an object of \code{\link{class}} "\code{discreteQ}". There are methods available for plotting ("\code{plot}", see
#'   \code{\link{plot.discreteQ}}) and summarizing ("\code{summary}", see
#'   \code{\link{summary.discreteQ}}) "\code{discreteQ}" objects. We recommend
#'   using them to analyze the results.
#'
#'   The components contained by an object of class "\code{discreteQ}" depend on the way the function has been called. We can distinguish the same three cases as above:
#'
#'  \enumerate{
#'  \item
#'  If no treatment variable is specified in the argument \code{d}, then \code{discreteQ} contains the following components:
#'  \describe{
#'  \item{\code{Q}}{The empirical quantile function of the outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.Q}}{The lower bound of the uniform confidence band for the unconditional quantile function of the outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.Q}}{The upper bound of the uniform confidence band for the unconditional quantile function of the outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{F}}{The empirical distribution function of the outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.F}}{The lower bound of the uniform confidence band for the unconditional distribution function of the outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.F}}{The upper bound of the uniform confidence band for the unconditional distribution function of the outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{q.range}}{Vector of length 2 that contains the lowest and highest quantile indexes. The uniform bands for the quantile function cover the true quantile function in this quantile range. The uniform bands for the distribution function cover the true function in the range of values of the outcome that are between the quantiles corresponding to these indexes.}
#'  \item{\code{ys}}{Vector containing the thresholds at which the cumulative distribution has been estimated.}
#'  \item{\code{bsrep}}{Scalar containing the number of performed bootstrap replications.}
#'  \item{\code{model}}{String scalar that takes the value "univariate" in this case.}
#'  \item{\code{method}}{String scalar that takes the value "empirical" in this case.}
#'  \item{\code{F.b}}{Matrix with \code{length(ys)} rows and \code{bsrep} columns. Each columns contains the estimated distribution function for the corresponding bootstrap replication. This object, which can be voluminous, is returned only if \code{return.boot = TRUE}.}
#'  \item{\code{seeds}}{List of length \code{bsrep} containing the seeds used for L'Ecuyer's RNG in the bootstrap replications. This object is returned only if \code{return.seeds = TRUE}.}
#'  }
#'
#'  \item
#'  If a treatment variable is specified in the argument \code{d} but \code{decomposition=FALSE}, then \code{discreteQ} contains the components below. Note that the uniform bands \strong{jointly} cover the true Q0, Q1, F0, F1 and QTE functions with probability \code{1-alpha}.
#'  \describe{
#'  \item{\code{Q0}}{The estimated unconditional quantile function of the control outcome. This is the quantile function of the outcome that we would observe if all observations had \code{d = 0}. If \code{x = NULL}, this is simply the empirical quantile function of the outcome for the subsample with \code{d = 0}. If regressors have been provided in the argument \code{x}, then the estimated conditional distribution of the outcome in the sample with \code{d = 0} is integrated over the distribution of the covariates in the whole sample. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.Q0}}{The lower bound of the uniform confidence band for the unconditional quantile function of the control outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.Q0}}{The upper bound of the uniform confidence band for the unconditional quantile function of the control outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{Q1}}{The estimated unconditional quantile function of the treated outcome. This is the quantile function of the outcome that we would observe if all observations had \code{d = 1}. If \code{x = NULL}, this is simply the empirical quantile function of the outcome for the subsample with \code{d = 1}. If regressors have been provided in the argument \code{x}, then the estimated conditional distribution of the outcome in the sample with \code{d = 1} is integrated over the distribution of the covariates in the whole sample. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.Q1}}{The lower bound of the uniform confidence band for the unconditional quantile function of the treated outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.Q1}}{The upper bound of the uniform confidence band for the unconditional quantile function of the treated outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{QTE}}{The estimated quantile treatment effect function: \code{QTE = Q1 - Q0}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.QTE}}{The lower bound of the uniform confidence band for the quantile treatment effect.  This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.QTE}}{The upper bound of the uniform confidence band for the quantile treatment effect.  This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{F0}}{The estimated unconditional distribution function of the control outcome. This is the distribution function of the outcome that we would observe if all observations had \code{d = 0}. If \code{x = NULL}, this is simply the empirical distribution function of the outcome for the subsample with \code{d = 0}. If regressors have been provided in the argument \code{x}, then the estimated conditional distribution of the outcome in the sample with \code{d = 0} is integrated over the distribution of the covariates in the whole sample. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.F0}}{The lower bound of the uniform confidence band for the unconditional distribution function of the control outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.F0}}{The upper bound of the uniform confidence band for the unconditional distribution function of the control outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{F1}}{The estimated unconditional distribution function of the treated outcome. This is the distribution function of the outcome that we would observe if all observations had \code{d = 1}. If \code{x = NULL}, this is simply the empirical distribution function of the outcome for the subsample with \code{d = 1}. If regressors have been provided in the argument \code{x}, then the estimated conditional distribution of the outcome in the sample with \code{d = 1} is integrated over the distribution of the covariates in the whole sample. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.F1}}{The lower bound of the uniform confidence band for the unconditional distribution function of the treated outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.F1}}{The upper bound of the uniform confidence band for the unconditional distribution function of the treated outcome. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{q.range}}{Vector of length 2 that contains the lowest and highest quantile indexes. The uniform bands for the quantile functions cover the true quantile function in this quantile range. The uniform bands for the distribution functions cover the true function in the range of values of the outcome that are between the quantiles corresponding to thes indexes.}
#'  \item{\code{ys0}}{Vector containing the thresholds at which the cumulative distribution of the control outcome has been estimated.}
#'  \item{\code{ys1}}{Vector containing the thresholds at which the cumulative distribution of the treated outcome has been estimated.}
#'  \item{\code{bsrep}}{Scalar containing the number of performed bootstrap replications.}
#'  \item{\code{model}}{String scalar that takes the value "qte" in this case.}
#'  \item{\code{method}}{String scalar. Name of the method used to estimate the conditional distribution functions.}
#'  \item{\code{F.b}}{Matrix with \code{length(ys0) + length(ys1)} rows and \code{bsrep} columns. Each columns contains the estimated distribution functions for the corresponding bootstrap replication. The first \code{length(ys0)} rows contains the estimated distribution function for the control outcome F0. The remaining \code{length(ys1)} rows contains the estimated distribution function for the treated outcome F1. This object, which can be voluminous, is returned only if \code{return.boot = TRUE}.}
#'  \item{\code{seeds}}{List of length \code{bsrep} containing the seeds used for L'Ecuyer's RNG in the bootstrap replications. This object is returned only if \code{return.seeds = TRUE}.}
#'  }
#'
#'  \item
#'  If a treatment variable is specified in the argument \code{d} and \code{decomposition=TRUE}, then \code{discreteQ} contains the components below. Note that the uniform bands \strong{jointly} cover the true functions Q0, Q1, Qc, F0, F1, and Fc as well as the difference between any two of these functions with probability \code{1-alpha}.
#'  \describe{
#'  \item{\code{Q0}}{The empirical quantile function of the outcome for the group with \code{d = 0}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.Q0}}{The lower bound of the uniform confidence band for the unconditional quantile function of outcome in the group with \code{d = 0}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.Q0}}{The upper bound of the uniform confidence band for the unconditional quantile function of outcome in the group with \code{d = 0}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{Q1}}{The empirical quantile function of the outcome for the group with \code{d = 1}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.Q1}}{The lower bound of the uniform confidence band for the unconditional quantile function of outcome in the group with \code{d = 1}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.Q1}}{The upper bound of the uniform confidence band for the unconditional quantile function of outcome in the group with \code{d = 1}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{Qc}}{The estimated counterfactual quantile function of the outcome that we would observe if the distribution of the covariates was the same as that of the group with \code{d = 0} and the conditional distribution of the outcome given the covariates was the same as that of the group with \code{d = 1}.  This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.Qc}}{The lower bound of the uniform confidence band for the counterfactual quantile function. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.Qc}}{The upper bound of the uniform confidence band for the counterfactual quantile function. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{F0}}{The empirical distribution function of the outcome for the group with \code{d = 0}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.F0}}{The lower bound of the uniform confidence band for the unconditional distribution function of outcome in the group with \code{d = 0}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.F0}}{The upper bound of the uniform confidence band for the unconditional distribution function of outcome in the group with \code{d = 0}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{F1}}{The empirical distribution function of the outcome for the group with \code{d = 1}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.F1}}{The lower bound of the uniform confidence band for the unconditional distribution function of outcome in the group with \code{d = 1}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.F1}}{The upper bound of the uniform confidence band for the unconditional distribution function of outcome in the group with \code{d = 1}. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{Fc}}{The estimated counterfactual distribution function of the outcome that we would observe if the distribution of the covariates was the same as that of the group with \code{d = 0} and the conditional distribution of the outcome given the covariates was the same as that of the group with \code{d = 1}.  This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.Fc}}{The lower bound of the uniform confidence band for the counterfactual distribution function. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.Fc}}{The upper bound of the uniform confidence band for the counterfactual distribution function. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{observed}}{The difference between the observed quantile function for the group \code{d = 1} and and the observed quantile function for the group  with \code{d = 0}: Q1 - Q0. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.observed}}{The lower bound of the uniform confidence band for Q1 - Q0. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.observed}}{The upper bound of the uniform confidence band for Q1 - Q0. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{composition}}{The difference between the observed quantile function for the group with \code{d = 1} and the counterfactual quantile function: Q1 - Qc. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.composition}}{The lower bound of the uniform confidence band for Q1 - Qc. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.composition}}{The upper bound of the uniform confidence band for Q1 - Qc. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{unexplained}}{The difference between the counterfactual quantile function and the quantile function for the group with \code{d = 0}: Qc - Q0. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{lb.unexplained}}{The lower bound of the uniform confidence band for Qc - Q0. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{ub.unexplained}}{The upper bound of the uniform confidence band for Qc - Q0. This is a function of class \code{\link[stats]{stepfun}}.}
#'  \item{\code{q.range}}{Vector of length 2 that contains the lowest and highest quantile indexes. The uniform bands for the quantile functions cover the true quantile function in this quantile range. The uniform bands for the distribution functions cover the true function in the range of values of the outcome that are between the quantiles corresponding to thes indexes.}
#'  \item{\code{ys0}}{Vector containing the thresholds at which the cumulative distribution of the outcome for the group with \code{d = 0} has been estimated.}
#'  \item{\code{ys1}}{Vector containing the thresholds at which the cumulative distribution of the outcome for the group with \code{d = 1} and the counterfactual distribution function have been estimated.}
#'  \item{\code{bsrep}}{Scalar containing the number of performed bootstrap replications.}
#'  \item{\code{model}}{String scalar that takes the value "distribution" in this case.}
#'  \item{\code{method}}{String scalar. Name of the method used to estimate the conditional distribution functions.}
#'  \item{\code{F.b}}{Matrix with \code{length(ys0) + 2 * length(ys1)} rows and \code{bsrep} columns. Each columns contains the estimated distribution functions for the corresponding bootstrap replication. The first \code{length(ys0)} rows contains the estimated distribution function for group 0, F0. The next \code{length(ys1)} rows contains the estimated distribution function for group 1, F1. The remaining \code{length(ys1)} rows contains the estimated counterfactual distribution, Fc. This object, which can be voluminous, is returned only if \code{return.boot = TRUE}.}
#'  \item{\code{seeds}}{List of length \code{bsrep} containing the seeds used for L'Ecuyer's RNG in the bootstrap replications. This object is returned only if \code{return.seeds = TRUE}.}
#'  }
#'  }
#'
#' @references
#' Chernozhukov, Victor, Iván Fernández-Val, Blaise Melly, and Kaspar Wüthrich. 2019. “Generic Inference on Quantile and Quantile Effect Functions for Discrete Outcomes.” arXiv Preprint, \url{https://arxiv.org/abs/1608.05142}.
#'
#' @seealso
#' For continous outcomes, see also \code{\link[Counterfactual]{counterfactual}} from the package Counterfactual.
#'
#' @examples
#' ##Example 1: univariate quantile function
#' #Generate the data
#' set.seed(1234)
#' outcome <- rpois(100, 3)
#' #Estimate the functions and the confidence bands
#' results1 <- discreteQ(outcome)
#' #Table containing the estimated quantile function with its confidence band
#' summary(results1)
#' #Plot the estimated quantile function with its confidence band
#' plot(results1)
#'
#' ##Example 2: quantile treatment effect function (QTE)
#' #Generate the data
#' set.seed(1234)
#' treatment <- c(rep(0,100), rep(1,100))
#' reg <- rbinom(200, 1, 0.4 + treatment*0.2)
#' outcome <- rpois(200, lambda = 1+reg)
#' #Estimate the functions and the confidence bands (takes about 1 minute)
#' results2 <- discreteQ(outcome, treatment, cbind(1, reg))
#' #Table containing the estimated QTE function with its confidence band
#' summary(results2)
#' #Plot the QTE with its confidence band
#' plot(results2)
#' #Plot the quantile function of the control outcome
#' plot(results2, which="Q0")
#' #Add the quantile function of the treated outcome
#' plot(results2, which="Q1", add=TRUE, shift=0.2, col.l="dark green", col.b="light green")
#'
#' ##Example 3: decomposition
#' #Generate the data
#' set.seed(1234)
#' group <- c(rep(0,100), rep(1,100))
#' reg <- rbinom(200, 1, 0.2 + group*0.6)
#' outcome <- rpois(200, lambda = exp(-3+4*reg))
#' #Estimate the functions and the confidence bands (takes about 30 seconds)
#' results3 <- discreteQ(outcome, group, cbind(1, reg), decomposition=TRUE)
#' #Table containing the unexplained component with its confidence band
#' summary(results3)
#' #Table containing the difference between the observed quantile functions
#' summary(results3, which="observed")
#' #Plot the observed quantile functions and their decomposition
#' plot(results3)
#' #Plot only the composition component
#' plot(results3, which="composition")

#' @importFrom graphics legend par plot segments polygon lines
#' @importFrom stats IQR binomial glm.fit knots lm.wfit quantile rexp stepfun
#'   weighted.mean glm optim poisson ppois
#' @importFrom plyr join
#' @importFrom foreach %dopar%
#' @importFrom rngtools RNGseed
#' @importFrom parallel clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom fastglm fastglm

#' @export
discreteQ <-
  function(y,
           d = NULL,
           x = NULL,
           w = NULL,
           decomposition = FALSE,
           q.range = c(0.05, 0.95),
           method = NULL,
           bsrep = 200,
           alpha = 0.05,
           ys = NULL,
           cl = NULL,
           cluster = NULL,
           old.res = NULL,
           return.boot = FALSE,
           list_of_seeds = NULL,
           return.seeds = FALSE,
           estim.glm = fastglm::fastglm,
           par.estim = NULL) {
    rng_old <- RNGkind()
    on.exit(RNGkind(rng_old[1], rng_old[2]), add = TRUE)
    RNGkind(kind = "L'Ecuyer-CMRG")
    if(!is.null(d)){
      if(length(d)!=length(y)) stop("y and d must have the same length.")
    }
    if(!is.null(x)){
      if(nrow(as.matrix(x))!=length(y)){
        stop("y and x must have the same number of rows.")
      }
    }
    if(!is.null(w)){
      if(length(w)!=length(y)) stop("y and w must have the same length.")
    }
    if(!is.null(cluster)){
      if(length(cluster)!=length(y)) stop("y and cluster must have the same length.")
    }
    if(!is.null(old.res)){
      bsrep <- old.res$bsrep
    }
    if (is.null(d)) {
      if(is.null(method)) method <- "empirical"
      if(method != "empirical")
        stop("The selected method has not yet been implemented.")
      if (!is.null(x))
        stop("The argument x cannot be specified if d = NULL.")
      fit <- dq_univariate(y, q.range, w, bsrep, alpha, ys, old.res, return.boot, list_of_seeds, return.seeds)
      fit$model <- "univariate"
    } else if (!decomposition) {
      if(is.null(method)) method <- "logit"
      if(!(method %in% c("logit", "probit", "cloglog", "poisson", "lpm", "drp", "cauchit", "log")))
        stop("The selected method has not yet been implemented.")
      fit <- dq_qte(y, d, x, w, q.range, method, bsrep, alpha, ys, cl, cluster, old.res, return.boot, list_of_seeds, return.seeds, estim.glm, par.estim)
      fit$model <- "qte"
    } else {
      if (is.null(x))
        stop("The argument x must be specified if decomposition = TRUE.")
      if(is.null(method)) method <- "logit"
      if(!(method %in% c("logit", "probit", "cloglog", "poisson", "lpm", "drp", "cauchit", "log")))
        stop("The selected method has not yet been implemented.")
      fit <-
        dq_decomposition(y, x, d, w, q.range, method, bsrep, alpha, ys, cl, cluster, old.res, return.boot, list_of_seeds, return.seeds, estim.glm, par.estim)
      fit$model <- "decomposition"
    }
    fit$method <- method
    class(fit) <- "discreteQ"
    fit
  }
