## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)

## ----example1------------------------------------------------------------
library(discreteQ)
set.seed(1234)
outcome <- rpois(100, 3)
results1 <- discreteQ(outcome)

## ----example2------------------------------------------------------------
results1

## ----example3, fig.width=10, fig.height=10-------------------------------
plot(results1)

## ----example4------------------------------------------------------------
summary(results1)

## ----example5------------------------------------------------------------
set.seed(1234)
treatment <- c(rep(0,1000), rep(1,1000))
reg <- rbinom(2000, 1, 0.4+treatment*0.2)
outcome <- rpois(2000, lambda = 2+4*reg)

## ----example6, fig.width=10, fig.height=10-------------------------------
results2 <- discreteQ(outcome, treatment)
plot(results2, main="Difference between the unconditional quantile functions")

## ----example7------------------------------------------------------------
results3 <- discreteQ(outcome, treatment, cbind(1, reg))
plot(results3)

## ----example8, fig.width=10, fig.height=10-------------------------------
plot(results3, which="Q0")
plot(results3, which="Q1", add=TRUE, shift=0.2, col.l="dark green", col.b="light green")

## ----example9, fig.width=10, fig.height=10-------------------------------
results4 <- discreteQ(outcome, treatment, cbind(1, reg), decomposition=TRUE)
plot(results4)

## ----example10, fig.width=10, fig.height=10------------------------------
set.seed(1234)
outcome <- rnorm(500, 3)
results5 <- discreteQ(outcome, ys = Inf)
plot(results5, support = "continuous")

## ----example11-----------------------------------------------------------
set.seed(1234)
treatment <- c(rep(0,1000), rep(1,1000))
reg <- rbinom(2000, 1, 0.4+treatment*0.2)
outcome <- rpois(2000, lambda = 2+4*reg)
#Without parallel computing
set.seed(42)
system.time(results6 <- discreteQ(outcome, treatment, reg))
my_cl <- parallel::makePSOCKcluster(2)
#With parallel computing
set.seed(42)
system.time(results7 <- discreteQ(outcome, treatment, reg, cl = my_cl ))

## ----example12-----------------------------------------------------------
#Results with and without parallel computing are equal
all.equal(results6, results7)

## ----example13-----------------------------------------------------------
#95% confidence bands (default value)
results8 <- discreteQ(outcome, treatment, reg, return.boot = TRUE)
#90% confidence bands
results9 <- discreteQ(outcome, treatment, reg, old.res = results8, alpha = 0.1)


