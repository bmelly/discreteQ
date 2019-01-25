#' Uniform inference on quantile and quantile effects for discrete outcomes
#'
#' \code{discreteQ} provides uniform confidence bands for the unconditional
#' quantile function, the quantile treatment effect function or the
#' decomposition of the observed difference between the quantile function of an
#' outcome for two groups.
#'
#' \code{discreteQ} can be used in three different ways. First, if no treatment
#' variable in \code{d} and no regressor in \code{x} are provided, then the the
#' command will provide uniform bands for the unconditional quantile function of
#' \code{y}. Second, if a treatment variable is provided in \code{d} but
#' \code{decomposition=FALSE}, then the command will provide uniform bands that
#' cover the quantile functions of the treated and control outcomes and the
#' quantile treatment effect functions. Third, if a treatment variable is
#' provided in \code{d} and \code{decomposition=TRUE}, then the command provides
#' uniforma bands that cover both quantile functions as well as the
#' decomposition of their difference into an explained (by the regressors
#' \code{x}) and unexplained component.
#'
#' The output of the function is a list of step functions. We recommend to use
#' the \code{\link{plot.discreteQ}} and \code{\link{summary.discreteQ}} functions.
#'
#' @param y outcome (vector).
#' @param d treatment or group variable (vector of 0-1 binary values).
#' @param x matrix of regressors (which must include a constant if appropriate).
#' @param w optional sampling weights (vector).
#' @param decomposition logical, indicating if the decomposition of the observed difference between the control and treated quantile functions should be performed. By default, inference about the quantile treatment effect function is performed.
#' @param q.range vector of length 2 that provides the lowest and highest quantile indexes. The uniform bands will cover the whole QF and QE functionin this range. Default: c(0.05,0.95).
#' @param method link function for the distribution regression model. Possible values: "logit" (the default), "probit", "cloglog", "lpm" (linear probability model, i.e. OLS estimation of binary regressions), "cauchit". This argument is relevant only if there are regressors in \code{x}.
#' @param bsrep number of bootstrap replications. Default: 200.
#' @param alpha confidence level. Default: 0.05.
#' @param ys specifies the thresholds at which the cumulative distribution function will be estimated. This argument can be specified either as a scalar that will be interpreted as the number of thresholds or as a vector that will contain the values of the thresholds. By default, the cdf is estimated at all distinct observed values of the outcome in the sample if there are less than 100 unique values and at 99 different values if there are more than 100 distinct values.
#' @param cl a cluster object as returned by makeCluster. Parallel computing is not used if this argument is not specified.
#' @param cluster vector that specifies to which group each observation belongs. The cluster bootstrap is used if this argument is specified. Otherwise, simple random sampling is assumed.
#' @param old.res a discreteQ object (obtained with the argument return.boot set to TRUE). This argument allows for instance to change the size alpha without recomputing the estimates.
#' @param return.boot logical scalar. The results of the bootstrap are return in the matrix F.b when this argument is set to TRUE.
#' @param list_of_seeds list of seeds for L'Ecuyer RNG. The length of this list must be the same as the value of the argument bsrep.
#' @param return.seeds logical scalar. The list of seeds is return by the function if this argument is set to TRUE.
#' @return A list of step functions (of class \code{stepfun}). For each disctribution function, quantile function or quantile effect function of interest (2 in case (i), 5 in case (ii) and 9 in case (iii)), three step functions are returned: one for the point estimates, one for the lower bound of the confidence band and one for the upper bound of the confidence band. There are methods available for plotting ("\code{plot}", see \code{\link{plot.discreteQ}}) and summarizing ("\code{summary}", see \code{\link{summary.discreteQ}}) "\code{discreteQ}" objects. We recommend using them to analyze the results.
#' @examples
#' set.seed(1234)
#' outcome <- rpois(100, 3)
#' results1 <- discreteQ(outcome)
#' summary(results1)
#' plot(results1)
#'
#' set.seed(1234)
#' treatment <- c(rep(0,1000), rep(1,1000))
#' reg <- rbinom(2000, 1, 0.4+treatment*0.2)
#' outcome <- rpois(2000, lambda = 2+4*reg)
#' results2 <- discreteQ(outcome, treatment, cbind(1, reg))
#' summary(results2)
#' plot(results2)
#' plot(results2, which="Q0")
#' plot(results2, which="Q1", add=TRUE, shift=0.2, col.l="dark green", col.b="light green")
#'
#' set.seed(1234)
#' group <- c(rep(0,1000), rep(1,1000))
#' reg <- rbinom(2000, 1, 0.4+group*0.2)
#' outcome <- rpois(2000, lambda = exp(-2+4*reg))
#' results3 <- discreteQ(outcome, group, cbind(1, reg), decomposition=TRUE)
#' summary(results3)
#' summary(results3, which="observed")
#' plot(results3)
#' plot(results3, which="composition")

#' @importFrom graphics legend par plot segments polygon lines
#' @importFrom stats IQR binomial glm.fit knots lm.wfit quantile rexp stepfun
#'   weighted.mean glm optim poisson ppois
#' @importFrom plyr join
#' @importFrom foreach %dopar%
#' @importFrom rngtools RNGseed
#' @importFrom parallel clusterExport
#' @importFrom doParallel registerDoParallel

#' @export
discreteQ <-
  function(y,
           d = NULL,
           x = NULL,
           w = NULL,
           decomposition = FALSE,
           q.range = c(0.05, 0.95),
           method = "logit",
           bsrep = 200,
           alpha = 0.05,
           ys = NULL,
           cl = NULL,
           cluster = NULL,
           old.res = NULL,
           return.boot = FALSE,
           list_of_seeds = NULL,
           return.seeds = FALSE) {
    if (is.null(d)) {
      if (!is.null(x))
        stop("The argument x cannot be specified if d=NULL.")
      fit <- dq_univariate(y, q.range, w, bsrep, alpha, ys, old.res, return.boot, list_of_seeds, return.seeds)
      fit$model <- "univariate"
    } else if (!decomposition) {
      fit <- dq_qte(y, d, x, w, q.range, method, bsrep, alpha, ys, cl, cluster, old.res, return.boot, list_of_seeds, return.seeds)
      fit$model <- "qte"
    } else {
      if (is.null(x))
        stop("The argument x must be specified if decomposition=TRUE.")
      fit <-
        dq_decomposition(y, x, d, w, q.range, method, bsrep, alpha, ys, cl, cluster, old.res, return.boot, list_of_seeds, return.seeds)
      fit$model <- "decomposition"
    }
    class(fit) <- "discreteQ"
    fit
  }
