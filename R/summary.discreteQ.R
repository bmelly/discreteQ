#' Summary method for class "discreteQ"
#'
#' Return a summary table for \code{\link{discreteQ}} object.
#'
#' \code{summary.discreteQ} produces a matrix with the number of rows equal to
#' the number of requested quantile indexes and 4 columns (one for the quantile
#' index, one for the point estimate and two for the bounds of the uniform
#' bands). The function for which the point estimates and the bands are reported
#' depends on the type of the \code{discreteQ} object and on the value of the
#' argument \code{which}. If the unconditional quantile function of the outcome
#' (i.e. no treatment was provided) has been estimated, then the results for
#' this functions are shown. If a treatment variable has been provided by
#' \code{decomposition=FALSE}, then by default the quantile treatment effect
#' function is shown but it is also possible to tabulate the quantile functions
#' of the control (\code{which="Q0"}) and treated outcomes (\code{which="Q1"}).
#' If \code{decomposition=TRUE}, by default the unexplained component is shown
#' but it is also possible to set the argument \code{which} to one of the
#' following values: "Q0", "Q1", "Qc", "observed", "composition", "unexplained".
#'
#' @param object an object produced by \code{discreteQ}.
#' @param taus defines the quantile indexes at which the results will be shown. There are two ways to specify this argument. (i) If taus is a scalar: the results are provided for a sequence of quantiles indexes in the range defined by q.range with an increment defined by taus. (ii) taus can also be a vector containing the requested quantile indexes. Default: taus=0.05.
#' @param which specifies the function for which the results are plotted. Possible values are (depending on the characteristics of the \code{discreteQ} object): "Q0", "Q1", "Qc", "QTE", observed", "composition", "unexplained". See `Details' below.
#' @param ... additional optional arguments.
#' @return A matrix with the same number of rows as specified with the argument taus and 4 columns. The first column contains the quantile indexes, the second column the point estimates, the third and the fourth column the uniform bands evaluated at this quantile index.
#' @examples
#' set.seed(1234)
#' outcome <- rpois(100, 3)
#' results1 <- discreteQ(outcome)
#' summary(results1)
#'
#' set.seed(1234)
#' treatment <- c(rep(0,100), rep(1,100))
#' reg <- rbinom(200, 1, 0.4+treatment*0.2)
#' outcome <- rpois(200, lambda = 2+4*reg)
#' results2 <- discreteQ(outcome, treatment, cbind(1, reg))
#' summary(results2)
#' summary(results2, which="Q0")
#'
#' set.seed(1234)
#' group <- c(rep(0,100), rep(1,100))
#' reg <- rbinom(200, 1, 0.4+group*0.2)
#' outcome <- rpois(200, lambda = exp(-2+4*reg))
#' results3 <- discreteQ(outcome, group, cbind(1, reg), decomposition=TRUE)
#' summary(results3)
#' summary(results3, which="composition")

#' @export
summary.discreteQ <- function(object, which=NULL, taus=0.05, ...){
  if (object$model=="decomposition") {
    if (is.null(which)) which <- "unexplained"
    dq_summary.decomposition(object, taus, which)
  } else if (object$model=="qte"){
    if (is.null(which)) which <- "QTE"
    dq_summary.qte(object, taus, which)
  } else if (object$model=="univariate"){
    dq_summary.univariate(object, taus)
  }
}
