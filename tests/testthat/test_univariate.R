context("Univariate")
library(discreteQ)

set.seed(41, kind = "L'Ecuyer-CMRG")
outcome <- rpois(100, 3)
quants <- runif(10)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results1 <- discreteQ(outcome, bsrep = 10)

test_that(
  "Check that the quantile and distribution function are correct for univariate", {
  expect_equal(results1$F((-1):20),sapply((-1):20, function(y) mean(outcome<=y)), check.attributes = FALSE, check.names = FALSE)
  for(q in quants) expect_equal(results1$Q(q), quantile(outcome, q, type=1), check.attributes = FALSE, check.names = FALSE)
})

test_that(
  "Check that the lower bound is lower than upper bound for univariate", {
    for(y in (-1):15) expect_lte(results1$lb.F(y), results1$ub.F(y))
    for(q in quants) expect_lte(results1$lb.Q(q),results1$ub.Q(q))
  })

test_that(
  "Check the formal output for univariate",{
  expect_s3_class(results1, "discreteQ")
  expect_true(results1$model=="univariate")
  expect_true(is.list(results1))
  expect_silent(plot(results1))
  expect_true(is.matrix(summary(results1)))
  expect_length(summary(results1), 76)
  }
)

set.seed(1234, kind = "L'Ecuyer-CMRG")
results2 <- discreteQ(outcome, bsrep = 10, return.boot = TRUE, return.seeds = TRUE)
test_that(
  "Check the bootstrap for univariate",{
    expect_false("F.b" %in% names(results1))
    expect_true("F.b" %in% names(results2))
    expect_true(is.matrix(results2$F.b))
    expect_length(results2$F.b, length(results2$ys)*10)
  }
)

set.seed(1234, kind = "L'Ecuyer-CMRG")
list_of_seeds <- rngtools::RNGseq(10, simplify=FALSE)
set.seed(999, kind = "L'Ecuyer-CMRG")
results3 <- discreteQ(outcome, bsrep=10, list_of_seeds = list_of_seeds)
results4 <- discreteQ(outcome, bsrep=10, list_of_seeds = results2$seeds)

test_that(
  "Check the replicability of the bootstrap",{
    expect_equal(results1$lb.F(0:10),results2$lb.F(0:10))
    expect_equal(results1$lb.F(0:10),results3$lb.F(0:10))
    expect_equal(results1$lb.F(0:10),results4$lb.F(0:10))
  }
)

ncores <- min(2,parallel::detectCores())
cl <- parallel::makePSOCKcluster(ncores)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results5 <- discreteQ(outcome, bsrep=10, cl=cl)
results6 <- discreteQ(outcome, bsrep=10, list_of_seeds = list_of_seeds)

test_that(
  "Check parallel computing for univariate",{
    expect_equal(results1$lb.F(0:10),results5$lb.F(0:10))
    expect_equal(results1$lb.F(0:10),results6$lb.F(0:10))
    expect_equal(results1$lb.F(0:10),results4$lb.F(0:10))
  }
)

cluster <- rep(1:20, each=5)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results7 <- discreteQ(outcome, bsrep=10, cluster = cluster)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results8 <- discreteQ(outcome, bsrep=10, cluster = cluster, cl=cl)
results9 <- discreteQ(outcome, bsrep=10, list_of_seeds = list_of_seeds, cluster=cluster, cl=cl)
parallel::stopCluster(cl)
test_that(
  "Check clustered s.e. for univariate",{
    expect_equal(results1$F(0:10), results7$F(0:10))
    expect_equal(results1$F(0:10), results8$F(0:10))
    for(q in quants) expect_equal(results1$Q(q),results7$Q(q))
    for(q in quants) expect_equal(results1$Q(q),results8$Q(q))
    expect_equal(results7$lb.F(0:10),results8$lb.F(0:10))
    expect_equal(results7$lb.F(0:10),results9$lb.F(0:10))
  }
)

weights <- c(rep(0,50), rep(2,50))
set.seed(1234, kind = "L'Ecuyer-CMRG")
results10 <- discreteQ(outcome, w = weights, bsrep = 10)
set.seed(1234, kind = "L'Ecuyer-CMRG")
results11 <- discreteQ(outcome[51:100], bsrep=10)
test_that(
  "Weights for univariate",{
    expect_equal(results10$F(0:10),results11$F(0:10))
  }
)

set.seed(1234, kind = "L'Ecuyer-CMRG")
results12 <- discreteQ(outcome, alpha = 0.1, bsrep = 10)
results13 <- discreteQ(outcome, alpha=0.1, old.res = results2, bsrep=10)
test_that(
  "Check that the size alpha is correct", {
    for(y in (-1):15) expect_lte(results1$lb.F(y), results12$lb.F(y))
    for(y in (-1):15) expect_gte(results1$ub.F(y), results12$ub.F(y))
    for(q in quants) expect_lte(results1$lb.Q(q),results12$lb.Q(q))
    for(q in quants) expect_gte(results1$ub.Q(q),results12$ub.Q(q))
    expect_true(any(results1$lb.F(0:15)-results12$lb.F(0:15)<0))
    expect_true(any(results1$ub.F(0:15)-results12$ub.F(0:15)>0))
    expect_equal(results12$lb.F(0:10),results13$lb.F(0:10))
    expect_equal(results12$ub.F(0:10),results13$ub.F(0:10))
  })

set.seed(1234, kind = "L'Ecuyer-CMRG")
results14 <- discreteQ(outcome, bsrep = 10, q.range = c(0.2,0.8))
test_that(
  "Check that q.range is correct", {
    for(y in (-1):15) expect_lte(results1$lb.F(y), results14$lb.F(y))
    for(y in (-1):15) expect_gte(results1$ub.F(y), results14$ub.F(y))
    for(q in quants) expect_lte(results1$lb.Q(q), results14$lb.Q(q))
    for(q in quants) expect_gte(results1$ub.Q(q), results14$ub.Q(q))
  })
