###############
# Data
###############

data <- read.dta(paste0(code.dir, "/Data/data_oregon.dta"))
data$treatment <-
  as.numeric(data$treatment) - 1 # Destring treatment
data <-
  na.omit(subset(
    data,
    select = c(
      doc_num_mod_12m,
      weight_12m,
      treatment,
      ddddraw_sur_2,
      ddddraw_sur_3,
      ddddraw_sur_4,
      ddddraw_sur_5,
      ddddraw_sur_6,
      ddddraw_sur_7,
      dddnumhh_li_2,
      dddnumhh_li_3,
      ddddraXnum_2_2,
      ddddraXnum_2_3,
      ddddraXnum_3_2,
      ddddraXnum_3_3,
      ddddraXnum_4_2,
      ddddraXnum_5_2,
      ddddraXnum_6_2,
      ddddraXnum_7_2,
      household_id
    )
  ))

#Dependent variable
y_name <- "doc_num_mod_12m"
#Treatment variable
d_name <- "treatment"
#Control variables
x_name <-
  c(
    "ddddraw_sur_2",
    "ddddraw_sur_3",
    "ddddraw_sur_4",
    "ddddraw_sur_5",
    "ddddraw_sur_6",
    "ddddraw_sur_7",
    "dddnumhh_li_2",
    "dddnumhh_li_3",
    "ddddraXnum_2_2",
    "ddddraXnum_2_3",
    "ddddraXnum_3_2",
    "ddddraXnum_3_3",
    "ddddraXnum_4_2",
    "ddddraXnum_5_2",
    "ddddraXnum_6_2",
    "ddddraXnum_7_2"
  )
#Weights
w_name <- "weight_12m"

###############
# Estimation
###############

#Poisson model
set.seed(42)
figure4_poisson <-
  discreteQ(
    y = data[, y_name],
    d = data[, d_name],
    x = cbind(1, data[, x_name]),
    w = data[, w_name],
    decomposition = FALSE,
    q.range = c(0, 0.97),
    method = "poisson",
    bsrep = 1000,
    alpha = 0.05,
    ys = 0:15,
    cl = cl,
    cluster = data$household_id,
    return.boot = TRUE
  )

#Distribution regression (with Poisson link function)
set.seed(42)
figure4_dr <-
  discreteQ(
    y = data[, y_name],
    d = data[, d_name],
    x = cbind(1, data[, x_name]),
    w = data[, w_name],
    decomposition = FALSE,
    q.range = c(0, 0.97),
    method = "drp",
    bsrep = 1000,
    alpha = 0.05,
    ys = 0:15,
    cl = cl,
    cluster = data$household_id,
    return.boot = TRUE
  )

###############
# Figure 4
###############

pdf(
  paste0(code.dir, "/Results/figure4.pdf"),
  pointsize = 15,
  width = 10,
  height = 16.5
)
par(mfrow = c(1, 1),
    mar = c(5.1, 4.1, 4.1, 2.1))
layout(matrix(c(1, 3, 5, 2, 4, 6), 3, 2), c(1, 1), c(1, 1, 1))

### Unconditional CDFs
##Poisson model (top left panel)
#Treated outcome
plot(
  figure4_poisson,
  which = "F1",
  xlim = c(0, 9),
  ylim = c(0, 1),
  ylab = "Probability",
  xlab = "Number of outpatient visits",
  main = "CDFs - Poisson Model",
  col.b = adjustcolor("light blue", alpha.f = 0.5)
)
#Control outcome
plot(
  figure4_poisson,
  which = "F0",
  col.l = "dark green",
  col.b = adjustcolor("light green", alpha.f = 0.5),
  add = TRUE
)
#Legend
legend(
  -0.2,
  1.05,
  c(" ", " "),
  col = c(
    adjustcolor("light blue", alpha.f = 0.5),
    adjustcolor("light green", alpha.f = 0.5)
  ),
  lwd = c(4, 4, 4),
  horiz = FALSE,
  bty = "n"
)
legend(
  -0.2,
  1.05,
  c("Treatment group", "Control group"),
  col = c("dark blue", "dark green"),
  lwd = c(1, 1, 1),
  horiz = FALSE,
  bty = "n"
)

##Distribution regression (top right panel)
#Treated outcome
plot(
  figure4_dr,
  which = "F1",
  xlim = c(0, 9),
  ylim = c(0, 1),
  ylab = "Probability",
  xlab = "Number of outpatient visits",
  main = "CDFs - Distribution Regression",
  col.b = adjustcolor("light blue", alpha.f = 0.5)
)
#Control outcome
plot(
  figure4_dr,
  which = "F0",
  col.l = "dark green",
  col.b = adjustcolor("light green", alpha.f = 0.5),
  add = TRUE
)
#Legend
legend(
  -0.2,
  1.05,
  c(" ", " "),
  col = c(
    adjustcolor("light blue", alpha.f = 0.5),
    adjustcolor("light green", alpha.f = 0.5)
  ),
  lwd = c(4, 4, 4),
  horiz = FALSE,
  bty = "n"
)
legend(
  -0.2,
  1.05,
  c("Treatment group", "Control group"),
  col = c("dark blue", "dark green"),
  lwd = c(1, 1, 1),
  horiz = FALSE,
  bty = "n"
)

### Unconditional quantile functions
#Poisson model (middle left panel)
#Treated outcome
plot(
  figure4_poisson,
  which = "Q1",
  ylim = c(0, 9),
  col.l = "dark blue",
  ylab = "Number of outpatient visits",
  xlab = "Probability",
  main = "Quantile functions - Poisson Model",
  col.b = adjustcolor("light blue", alpha.f = 0.5),
  shift = -0.1
)
#Control outcome
plot(
  figure4_poisson,
  which = "Q0",
  col.l = "dark green",
  col.b = adjustcolor("light green", alpha.f = 0.5),
  add = TRUE,
  shift = 0.1
)

##Distribution regression (middle right panel)
#Treated outcome
plot(
  figure4_dr,
  which = "Q1",
  ylim = c(0, 9),
  col.l = "dark blue",
  ylab = "Number of outpatient visits",
  xlab = "Probability",
  main = "Quantile functions - Distribution Regression",
  col.b = adjustcolor("light blue", alpha.f = 0.5),
  shift = -0.1
)
#Control outcome
plot(
  figure4_dr,
  which = "Q0",
  col.l = "dark green",
  col.b = adjustcolor("light green", alpha.f = 0.5),
  add = TRUE,
  shift = 0.1
)

### QE functions
##Poisson model (bottom left panel)
plot(
  figure4_poisson,
  ylim = c(-2, 3),
  col.l = "dark blue",
  ylab = "Difference in the number of outpatient visits",
  xlab = "Probability",
  main = "Quantile treatment effects - Poisson Model",
  col.b = "light blue"
)
##Distribution regression (bottom right panel)
plot(
  figure4_dr,
  ylim = c(-2, 3),
  col.l = "dark blue",
  ylab = "Difference in the number of outpatient visits",
  xlab = "Probability",
  main = "Quantile treatment effects - Distribution Regression",
  col.b = "light blue"
)

dev.off()

### Test of equality between Poisson and DR
ys1 <- figure4_poisson$ys1
ys0 <- figure4_poisson$ys0
delta1.dif <-
  figure4_poisson$F.b[(length(ys0) + 1):(length(ys0) + length(ys1)),] - 
  figure4_poisson$F1(ys1) - 
  figure4_dr$F.b[(length(ys0) + 1):(length(ys0) + length(ys1)),] + 
  figure4_dr$F1(ys1)
delta0.dif <-
  figure4_poisson$F.b[(length(ys0) + 1):(length(ys0) + length(ys1)),] - 
  figure4_poisson$F1(ys1) - 
  figure4_dr$F.b[(length(ys0) + 1):(length(ys0) + length(ys1)),] + 
  figure4_dr$F1(ys1)
se.1 <- apply(delta1.dif, 1, function(x) IQR(x) / 1.349)
se.0 <- apply(delta0.dif, 1, function(x) IQR(x) / 1.349)
# test statistic
t1 <- max(abs(figure4_poisson$F1(ys1) - figure4_dr$F1(ys1)) / se.1)
t0 <- max(abs(figure4_poisson$F0(ys0) - figure4_dr$F0(ys0)) / se.0)
# p-values
z1.dif <- apply(abs(delta1.dif) / se.1, 2, max)
mean(z1.dif > t1)
z0.dif <- apply(abs(delta0.dif) / se.0, 2, max)
mean(z0.dif > t0)
