#' Plotting uniform bands for quantile and quantile effect functions for
#' discrete outcomes
#'
#' Plots point estimates and uniform bands for quantile functions and quantile
#' effect functions of possibly discrete outcomes.
#'
#' \code{plot.discreteQ} will produce different types of plots depending on the
#' type of the \code{discreteQ} object and on the value of the argument
#' \code{which}. If the unconditional quantile function of the outcome (i.e. no
#' treatment was provided) has been estimated, then the only possibility
#' consists in plotting the estimated quantile function and its confidence
#' bands. If a treatment variable has been provided by
#' \code{decomposition=FALSE}, then by default the quantile treatment effect
#' function is plotted but it is also possible to plot the quantile functions of
#' the control (\code{which="Q0"}) and treated outcomes (\code{which="Q1"}). If
#' \code{decomposition=TRUE}, by default a matrix of 4 plots is produced: all
#' three quantile functions in the top-left panel, the observed difference in
#' the top-right panel, the composition effect in the bottom-right panel and the
#' unexplained component in the bottom-right panel. It is possible to plot any
#' of these components by setting the argument \code{which} to one of the
#' following values: "Q0", "Q1", "Qc", "observed", "composition", "unexplained".
#'
#' None.
#'
#' @param x an object produced by \code{discreteQ}.
#' @param which specifies the function for which the results are plotted. Possible values are (depending on the characteristics of the \code{discreteQ} object): "Q0", "Q1", "Qc", "QTE", "decomposition", "observed", "composition", "unexplained". See `Details' below.
#' @param xlim the x limits (x1, x2) of the plot. By default, this corresponds to the range of quantiles specified by the argument q.range when calling \code{discreteQ}. It is not allowed to specify wider x limits than q.range because the uniform bands cover the functions only over this range.
#' @param ylim the y limits of the plot.
#' @param main the title of the plot. If \code{which="decomposition"} then main must be a character vector of length 4 for the 4 panels of the figure. Sensible default values are provided.
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param add logical; if TRUE only add to an existing plot.
#' @param col.l color used for the line representing the point estimates. If \code{which="decomposition"} then col.l must be a color vector of length 3 for the three quantile functions.
#' @param lty.l line type.
#' @param lwd.l line width.
#' @param col.b color used for the confidence bands. If \code{which="decomposition"} then col.b must be a color vector of length 3 for the three quantile functions.
#' @param lty.b line types for the uniform bands.
#' @param lwd.b line width for the confidence bands.
#' @param shift scalar shifting the QF and the bands. Useful if add=TRUE and the functions overlap.
#' @param support type of support restriction. If \code{support="empirical"}, it is assumed that the support of the function(s) of interest (Q0, Q1, Qc, QTE, etc.) is equal to the support of the outcome in the sample. Alternatively, no restriction on the support is assumed if \code{support="continuous"}. Finally, \code{support} may directly contain the vector of all values that the relevant function (Q0, Q1 or QTE) can take. By default, \code{support="empirical"} if there are less than 200 points in the empirical support and \code{support="continuous"} otherwise.
#' @param ... other graphical parameters passed to \code{plot}.
#' @return None.
#' @examples
#' set.seed(1234, kind = "L'Ecuyer-CMRG")
#' outcome <- rpois(100, 3)
#' results1 <- discreteQ(outcome)
#' plot(results1)
#'
#' set.seed(1234, kind = "L'Ecuyer-CMRG")
#' treatment <- c(rep(0,100), rep(1,100))
#' reg <- rbinom(200, 1, 0.4+treatment*0.2)
#' outcome <- rpois(200, lambda = 2+4*reg)
#' results2 <- discreteQ(outcome, treatment, cbind(1, reg))
#' plot(results2)
#' plot(results2, which="Q0")
#' plot(results2, which="Q1", add=TRUE, shift=0.2, col.l="dark green", col.b="light green")
#'
#' set.seed(1234, kind = "L'Ecuyer-CMRG")
#' group <- c(rep(0,100), rep(1,100))
#' reg <- rbinom(200, 1, 0.4+group*0.2)
#' outcome <- rpois(200, lambda = exp(-2+4*reg))
#' results3 <- discreteQ(outcome, group, cbind(1, reg), decomposition=TRUE)
#' plot(results3)
#' plot(results3, which="composition")

#' @export
plot.discreteQ <- function(x, ..., which=NULL, xlim=NULL, ylim=NULL, main=NULL, xlab=NULL, ylab=NULL, add=FALSE, col.l="dark blue", col.b="light blue", shift=NULL, lty.l=1, lwd.l=1, lty.b=1, lwd.b=5, support=NULL){
  if (x$model=="decomposition") {
    dq_plot.decomposition(x, which, xlim, ylim, main, xlab, ylab, add, col.l, col.b, shift, lty.l, lwd.l, lty.b, lwd.b, support, ...)
  } else if (x$model=="qte"){
    dq_plot.qte(x, which, xlim, ylim, main, xlab, ylab, add, col.l, col.b, shift, lty.l, lwd.l, lty.b, lwd.b, support, ...)
  } else if (x$model=="univariate"){
    dq_plot.univariate(x, xlim, ylim, main, xlab, ylab, add, col.l, col.b, shift, lty.l, lwd.l, lty.b, lwd.b, support, ...)
  }
  invisible(x)
}
