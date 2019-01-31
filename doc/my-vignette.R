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

## ----example3------------------------------------------------------------
plot(results1)

## ----example 4-----------------------------------------------------------
summary(results1)

## ----example 5-----------------------------------------------------------
set.seed(1234)
treatment <- c(rep(0,1000), rep(1,1000))
reg <- rbinom(2000, 1, 0.4+treatment*0.2)
outcome <- rpois(2000, lambda = 2+4*reg)

## ----example 6-----------------------------------------------------------
results2 <- discreteQ(outcome, treatment)
plot(results2, main="Difference between the unconditional quantile functions")

## ----example 7-----------------------------------------------------------
results3 <- discreteQ(outcome, treatment, cbind(1, reg))
plot(results3)

## ----example 8-----------------------------------------------------------
plot(results3, which="Q0")
plot(results3, which="Q1", add=TRUE, shift=0.2, col.l="dark green", col.b="light green")

## ----example 9, fig.width=10, fig.height=10------------------------------
results4 <- discreteQ(outcome, treatment, cbind(1, reg), decomposition=TRUE)
plot(results4)

