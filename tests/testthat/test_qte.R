context("QTE")
library(discreteQ)

set.seed(41, kind = "L'Ecuyer-CMRG")
treatment <- c(rep(0,50), rep(1,50))
reg <- rbinom(100, 1, 0.4+treatment*0.2)
outcome <- rpois(100, lambda = 2+4*reg)
weights <- rbinom(100, 1, 0.5)
quants <- runif(10)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results1 <- discreteQ(outcome, treatment, bsrep = 10)

test_that(
  "Check that the quantile and distribution function are correct for qte", {
  expect_equal(results1$F0(0:12), sapply(0:12, function(y) mean(outcome[treatment==0]<=y)), check.attributes = FALSE, check.names = FALSE)
  expect_equal(results1$F1(0:12), sapply(0:12, function(y) mean(outcome[treatment==1]<=y)), check.attributes = FALSE, check.names = FALSE)
  for(q in quants) expect_equal(results1$Q0(q), quantile(outcome[treatment==0], q, type=1), check.attributes = FALSE, check.names = FALSE)
  for(q in quants) expect_equal(results1$Q1(q), quantile(outcome[treatment==1], q, type=1), check.attributes = FALSE, check.names = FALSE)
  })

test_that(
  "Check that the lower bound is lower than upper bound for qte", {
    for(y in (-1):15) expect_lte(results1$lb.F0(y), results1$ub.F0(y))
    for(y in (-1):15) expect_lte(results1$lb.F1(y), results1$ub.F1(y))
    for(q in quants) expect_lte(results1$lb.Q0(q),results1$ub.Q0(q))
    for(q in quants) expect_lte(results1$lb.Q1(q),results1$ub.Q1(q))
    for(q in quants) expect_lte(results1$lb.QTE(q),results1$ub.QTE(q))
  })

test_that(
  "Check the formal output for qte",{
  expect_s3_class(results1, "discreteQ")
  expect_true(results1$model=="qte")
  expect_true(is.list(results1))
  expect_silent(plot(results1))
  expect_true(is.matrix(summary(results1)))
  expect_length(summary(results1), 76)
  }
)

set.seed(1234, kind = "L'Ecuyer-CMRG")
results2 <- discreteQ(outcome, treatment, bsrep = 10, return.boot = TRUE, return.seeds = TRUE)
test_that(
  "Check the bootstrap for qte",{
    expect_false("F.b" %in% names(results1))
    expect_true("F.b" %in% names(results2))
    expect_true(is.matrix(results2$F.b))
    expect_length(results2$F.b, (length(results2$ys0)+length(results2$ys1))*10)
  }
)

set.seed(1234, kind = "L'Ecuyer-CMRG")
list_of_seeds <- rngtools::RNGseq(10, simplify=FALSE)
set.seed(999, kind = "L'Ecuyer-CMRG")
results3 <- discreteQ(outcome, treatment, bsrep=10, list_of_seeds = list_of_seeds)
results4 <- discreteQ(outcome, treatment, bsrep=10, list_of_seeds = results2$seeds)

test_that(
  "Check the replicability of the bootstrap for qte",{
    expect_equal(results1$lb.F0(0:10),results2$lb.F0(0:10))
    expect_equal(results1$lb.F0(0:10),results3$lb.F0(0:10))
    expect_equal(results1$lb.F0(0:10),results4$lb.F0(0:10))
    expect_equal(results1$lb.F1(0:10),results2$lb.F1(0:10))
    expect_equal(results1$lb.F1(0:10),results3$lb.F1(0:10))
    expect_equal(results1$lb.F1(0:10),results4$lb.F1(0:10))
  }
)

ncores <- min(2,parallel::detectCores())
cl <- parallel::makePSOCKcluster(ncores)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results5 <- discreteQ(outcome, treatment, bsrep=10, cl=cl)
results6 <- discreteQ(outcome, treatment, bsrep=10, list_of_seeds = list_of_seeds)

test_that(
  "Check parallel computing for qte",{
    expect_equal(results1$lb.F0(0:10),results5$lb.F0(0:10))
    expect_equal(results1$lb.F0(0:10),results6$lb.F0(0:10))
    expect_equal(results1$lb.F0(0:10),results4$lb.F0(0:10))
    expect_equal(results1$lb.F1(0:10),results5$lb.F1(0:10))
    expect_equal(results1$lb.F1(0:10),results6$lb.F1(0:10))
    expect_equal(results1$lb.F1(0:10),results4$lb.F1(0:10))
  }
)

cluster <- rep(1:20, each=5)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results7 <- discreteQ(outcome, treatment, bsrep=10, cluster = cluster)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results8 <- discreteQ(outcome, treatment, bsrep=10, cluster = cluster, cl=cl)
results9 <- discreteQ(outcome, treatment, bsrep=10, list_of_seeds = list_of_seeds, cluster=cluster, cl=cl)
parallel::stopCluster(cl)
test_that(
  "Check clustered s.e. for qte",{
    expect_equal(results1$F0(0:10), results7$F0(0:10))
    expect_equal(results1$F0(0:10), results8$F0(0:10))
    expect_equal(results1$F1(0:10), results7$F1(0:10))
    expect_equal(results1$F1(0:10), results8$F1(0:10))
    for(q in quants) expect_equal(results1$Q0(q),results7$Q0(q))
    for(q in quants) expect_equal(results1$Q0(q),results8$Q0(q))
    for(q in quants) expect_equal(results1$Q1(q),results7$Q1(q))
    for(q in quants) expect_equal(results1$Q1(q),results8$Q1(q))
    expect_equal(results7$lb.F0(0:10),results8$lb.F0(0:10))
    expect_equal(results7$lb.F0(0:10),results9$lb.F0(0:10))
    expect_equal(results7$lb.F1(0:10),results8$lb.F1(0:10))
    expect_equal(results7$lb.F1(0:10),results9$lb.F1(0:10))
  }
)

set.seed(1234, kind = "L'Ecuyer-CMRG")
results10 <- discreteQ(outcome, treatment, w = weights, bsrep = 10)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results11 <- discreteQ(outcome[weights==1], treatment[weights==1], bsrep=10)
test_that(
  "Weights for qte",{
    expect_equal(results10$F0(0:10), results11$F0(0:10))
    expect_equal(results10$F1(0:10), results11$F1(0:10))
  }
)

set.seed(1234, kind = "L'Ecuyer-CMRG")
results12 <- discreteQ(outcome, treatment, alpha = 0.1, bsrep = 10)
results13 <- discreteQ(outcome, treatment, alpha=0.1, old.res = results2, bsrep=10)
test_that(
  "Check that the size alpha is correct", {
    for(y in (-1):15) expect_lte(results1$lb.F0(y), results12$lb.F0(y))
    for(y in (-1):15) expect_gte(results1$ub.F0(y), results12$ub.F0(y))
    for(y in (-1):15) expect_lte(results1$lb.F1(y), results12$lb.F1(y))
    for(y in (-1):15) expect_gte(results1$ub.F1(y), results12$ub.F1(y))
    for(q in quants) expect_lte(results1$lb.Q0(q),results12$lb.Q0(q))
    for(q in quants) expect_gte(results1$ub.Q0(q),results12$ub.Q0(q))
    for(q in quants) expect_lte(results1$lb.Q1(q),results12$lb.Q1(q))
    for(q in quants) expect_gte(results1$ub.Q1(q),results12$ub.Q1(q))
    expect_true(any(results1$lb.F0(0:15)-results12$lb.F0(0:15)<0))
    expect_true(any(results1$ub.F0(0:15)-results12$ub.F0(0:15)>0))
    expect_true(any(results1$lb.F1(0:15)-results12$lb.F1(0:15)<0))
    expect_true(any(results1$ub.F1(0:15)-results12$ub.F1(0:15)>0))
    expect_equal(results12$lb.F0(0:10),results13$lb.F0(0:10))
    expect_equal(results12$ub.F0(0:10),results13$ub.F0(0:10))
    expect_equal(results12$lb.F1(0:10),results13$lb.F1(0:10))
    expect_equal(results12$ub.F1(0:10),results13$ub.F1(0:10))
  })

set.seed(1234, kind = "L'Ecuyer-CMRG")
results14 <- discreteQ(outcome, treatment, bsrep = 10, q.range = c(0.2,0.8))
test_that(
  "Check that q.range is correct", {
    for(y in (-1):15) expect_lte(results1$lb.F0(y), results14$lb.F0(y))
    for(y in (-1):15) expect_gte(results1$ub.F0(y), results14$ub.F0(y))
    for(y in (-1):15) expect_lte(results1$lb.F1(y), results14$lb.F1(y))
    for(y in (-1):15) expect_gte(results1$ub.F1(y), results14$ub.F1(y))
    for(q in quants) expect_lte(results1$lb.Q0(q), results14$lb.Q0(q))
    for(q in quants) expect_gte(results1$ub.Q0(q), results14$ub.Q0(q))
    for(q in quants) expect_lte(results1$lb.Q1(q), results14$lb.Q1(q))
    for(q in quants) expect_gte(results1$ub.Q1(q), results14$ub.Q1(q))
  })

for(method in c("logit", "probit", "lpm", "cloglog", "poisson", "drp")){
  set.seed(1234, kind = "L'Ecuyer-CMRG")
  results1 <- discreteQ(outcome, treatment, cbind(1,reg), method = method, bsrep = 10)

  test_that(
    paste("Check that the lower bound is lower than upper bound for qte, method: ", method), {
      for(y in (-1):15) expect_lte(results1$lb.F0(y), results1$ub.F0(y))
      for(y in (-1):15) expect_lte(results1$lb.F1(y), results1$ub.F1(y))
      for(q in quants) expect_lte(results1$lb.Q0(q),results1$ub.Q0(q))
      for(q in quants) expect_lte(results1$lb.Q1(q),results1$ub.Q1(q))
      for(q in quants) expect_lte(results1$lb.QTE(q),results1$ub.QTE(q))
    }
  )

  test_that(
    paste("Check the formal output for qte, method:", method),{
      expect_s3_class(results1, "discreteQ")
      expect_true(results1$model=="qte")
      expect_true(is.list(results1))
      expect_silent(plot(results1))
      expect_true(is.matrix(summary(results1)))
      expect_length(summary(results1), 76)
    }
  )

  set.seed(1234, kind = "L'Ecuyer-CMRG")
  results2 <- discreteQ(outcome, treatment, cbind(1,reg), method=method, bsrep = 10, return.boot = TRUE, return.seeds = TRUE)
  test_that(
    paste("Check the bootstrap for qte, method: ", method),{
      expect_false("F.b" %in% names(results1))
      expect_true("F.b" %in% names(results2))
      expect_true(is.matrix(results2$F.b))
      expect_length(results2$F.b, (length(results2$ys0)+length(results2$ys1))*10)
    }
  )

  set.seed(1234, kind = "L'Ecuyer-CMRG")
  list_of_seeds <- rngtools::RNGseq(10, simplify=FALSE)
  set.seed(999, kind = "L'Ecuyer-CMRG")
  results3 <- discreteQ(outcome, treatment, cbind(1,reg), method=method, bsrep=10, list_of_seeds = list_of_seeds)
  results4 <- discreteQ(outcome, treatment, cbind(1,reg), method=method, bsrep=10, list_of_seeds = results2$seeds)

  test_that(
    paste("Check the replicability of the bootstrap for qte, method: ", method),{
      expect_equal(results1$lb.F0(0:10),results2$lb.F0(0:10))
      expect_equal(results1$lb.F0(0:10),results3$lb.F0(0:10))
      expect_equal(results1$lb.F0(0:10),results4$lb.F0(0:10))
      expect_equal(results1$lb.F1(0:10),results2$lb.F1(0:10))
      expect_equal(results1$lb.F1(0:10),results3$lb.F1(0:10))
      expect_equal(results1$lb.F1(0:10),results4$lb.F1(0:10))
    }
  )

  ncores <- min(2,parallel::detectCores())
  cl <- parallel::makePSOCKcluster(ncores)
  set.seed(1234, kind = "L'Ecuyer-CMRG")
  results5 <- discreteQ(outcome, treatment, cbind(1,reg), method=method, bsrep=10, cl=cl)
  results6 <- discreteQ(outcome, treatment, cbind(1,reg), method=method, bsrep=10, list_of_seeds = list_of_seeds)

  test_that(
    paste("Check parallel computing for qte, method: ", method),{
      expect_equal(results1$lb.F0(0:10),results5$lb.F0(0:10))
      expect_equal(results1$lb.F0(0:10),results6$lb.F0(0:10))
      expect_equal(results1$lb.F0(0:10),results4$lb.F0(0:10))
      expect_equal(results1$lb.F1(0:10),results5$lb.F1(0:10))
      expect_equal(results1$lb.F1(0:10),results6$lb.F1(0:10))
      expect_equal(results1$lb.F1(0:10),results4$lb.F1(0:10))
    }
  )

  cluster <- rep(1:20, each=5)
  set.seed(1234, kind = "L'Ecuyer-CMRG")
  results7 <- discreteQ(outcome, treatment, cbind(1,reg), method = method, bsrep=10, cluster = cluster)
  set.seed(1234, kind = "L'Ecuyer-CMRG")
  results8 <- discreteQ(outcome, treatment, cbind(1,reg), method = method, bsrep=10, cluster = cluster, cl=cl)
  results9 <- discreteQ(outcome, treatment, cbind(1,reg), method = method, bsrep=10, list_of_seeds = list_of_seeds, cluster=cluster, cl=cl)
  parallel::stopCluster(cl)
  test_that(
    paste("Check clustered s.e. for qte, method: ", method),{
      expect_equal(results1$F0(0:10), results7$F0(0:10))
      expect_equal(results1$F0(0:10), results8$F0(0:10))
      expect_equal(results1$F1(0:10), results7$F1(0:10))
      expect_equal(results1$F1(0:10), results8$F1(0:10))
      for(q in quants) expect_equal(results1$Q0(q),results7$Q0(q))
      for(q in quants) expect_equal(results1$Q0(q),results8$Q0(q))
      for(q in quants) expect_equal(results1$Q1(q),results7$Q1(q))
      for(q in quants) expect_equal(results1$Q1(q),results8$Q1(q))
      expect_equal(results7$lb.F0(0:10),results8$lb.F0(0:10))
      expect_equal(results7$lb.F0(0:10),results9$lb.F0(0:10))
      expect_equal(results7$lb.F1(0:10),results8$lb.F1(0:10))
      expect_equal(results7$lb.F1(0:10),results9$lb.F1(0:10))
    }
  )

  set.seed(1234, kind = "L'Ecuyer-CMRG")
  results10 <- discreteQ(outcome, treatment, cbind(1,reg), method = method, w = weights, bsrep = 10)
  set.seed(1234, kind = "L'Ecuyer-CMRG")
  results11 <- discreteQ(outcome[weights==1], treatment[weights==1], cbind(1,reg[weights==1]), method = method, bsrep=10)
  test_that(
    paste("Weights for qte, method: ", method),{
      expect_equal(results10$F0(results11$ys0), results11$F0(results11$ys0), tolerance=0.001)
      expect_equal(results10$F1(results11$ys1), results11$F1(results11$ys1), tolerance=0.001)
    }
  )

  set.seed(1234, kind = "L'Ecuyer-CMRG")
  results12 <- discreteQ(outcome, treatment, cbind(1,reg), method = method, alpha = 0.1, bsrep = 10)
  results13 <- discreteQ(outcome, treatment, cbind(1,reg), method = method, alpha=0.1, old.res = results2, bsrep=10)
  test_that(
    paste("Check that the size alpha is correct, method: ", method), {
      for(y in (-1):15) expect_lte(results1$lb.F0(y), results12$lb.F0(y))
      for(y in (-1):15) expect_gte(results1$ub.F0(y), results12$ub.F0(y))
      for(y in (-1):15) expect_lte(results1$lb.F1(y), results12$lb.F1(y))
      for(y in (-1):15) expect_gte(results1$ub.F1(y), results12$ub.F1(y))
      for(q in quants) expect_lte(results1$lb.Q0(q),results12$lb.Q0(q))
      for(q in quants) expect_gte(results1$ub.Q0(q),results12$ub.Q0(q))
      for(q in quants) expect_lte(results1$lb.Q1(q),results12$lb.Q1(q))
      for(q in quants) expect_gte(results1$ub.Q1(q),results12$ub.Q1(q))
      expect_true(any(results1$lb.F0(0:15)-results12$lb.F0(0:15)<0))
      expect_true(any(results1$ub.F0(0:15)-results12$ub.F0(0:15)>0))
      expect_true(any(results1$lb.F1(0:15)-results12$lb.F1(0:15)<0))
      expect_true(any(results1$ub.F1(0:15)-results12$ub.F1(0:15)>0))
      expect_equal(results12$lb.F0(0:10),results13$lb.F0(0:10))
      expect_equal(results12$ub.F0(0:10),results13$ub.F0(0:10))
      expect_equal(results12$lb.F1(0:10),results13$lb.F1(0:10))
      expect_equal(results12$ub.F1(0:10),results13$ub.F1(0:10))
    })

  set.seed(1234, kind = "L'Ecuyer-CMRG")
  results14 <- discreteQ(outcome, treatment, cbind(1,reg), method = method, bsrep = 10, q.range = c(0.2,0.8))
  test_that(
    paste("Check that q.range is correct, method: ", method), {
      for(y in (-1):15) expect_lte(results1$lb.F0(y), results14$lb.F0(y))
      for(y in (-1):15) expect_gte(results1$ub.F0(y), results14$ub.F0(y))
      for(y in (-1):15) expect_lte(results1$lb.F1(y), results14$lb.F1(y))
      for(y in (-1):15) expect_gte(results1$ub.F1(y), results14$ub.F1(y))
      for(q in quants) expect_lte(results1$lb.Q0(q), results14$lb.Q0(q))
      for(q in quants) expect_gte(results1$ub.Q0(q), results14$ub.Q0(q))
      for(q in quants) expect_lte(results1$lb.Q1(q), results14$lb.Q1(q))
      for(q in quants) expect_gte(results1$ub.Q1(q), results14$ub.Q1(q))
    })

}
