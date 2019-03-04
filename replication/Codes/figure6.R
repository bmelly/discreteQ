#######################################################################################

# This application performs a decomposition of the black-white testscore gap based on
# "Testing for Racial Differences in the Mental Ability of Young Children" 
# by R.G. Fryer and S.D. Levitt

# Authors: V. Chernozhukov, I. Fernandez-Val, B. Melly, K. Wuethrich

# Data source: AER webpage, CPP dataset

#######################################################################################


##############################################################################
# Data
# Prior sample selection: drop all observations with hispanic=1 and other=1
# Code works if black=0 only if white=1
##############################################################################

data <- read.dta(paste0(code.dir, "/Data/cpp_selected.dta"))

#drop if interviewer interviewd only blacks
for(i in setdiff(unique(data$interviewer_7years[data$black==1]), unique(data$interviewer_7years[data$black==0]))){
  data <- data[data$interviewer_7years!=i,]
}

# Factor variables (to shorten formulas)
data$age_7years.f <- factor(data$age_7years) 
data$interviewer_7years.f <- factor(data$interviewer_7years)

#Dependent variable
dep <- data[, "stand_fullscale_iq_7years"]
#Group variable
white <- 1-data[, "black"]
#Control variables  
x_name <- "~age_7years.f+female+dad_hs_dropout+dad_hs_grad+dad_some_college+dad_college_plus+dad_no_occupation+dad_professional+dad_non_professional+mom_hs_dropout+mom_hs_grad+mom_some_college+mom_college_plus+mom_no_occupation+mom_professional+mom_non_professional+inc_less_500+inc_500_1000+inc_1000_1500+inc_1500_2000+inc_2000_2500+inc_2500_plus+siblings_0+siblings_1+siblings_2+siblings_3+siblings_4+siblings_5+siblings_6_plus+both_bio_parents+age_mom+age_mom_2 +age_mom_3+age_mom_4+age_mom_5+miss_age_mom+mother_indifferent+mother_accepting+mother_attentive+mother_over_caring+mother_other+miss_parental_score+w_less_1500+w_1500_2500+w_2500_3500+w_3500_more+weeks_premature_0+weeks_premature_1+weeks_premature_2+weeks_premature_3+weeks_premature_4+weeks_premature_5+ weeks_premature_6+weeks_premature_7+weeks_premature_8+weeks_premature_9+weeks_premature_10+weeks_premature_11 +singleton+twin+high_order_multiple+interviewer_7years.f"
x <- model.matrix(formula(x_name), data = data)
rm(data)
ys <- unique(dep)
ys <- ys[ys>=quantile(dep, 0.01) & ys<=quantile(dep, 0.99)]

### Point estimates
set.seed(1234)
list_of_seeds <- rngtools::RNGseq(1000, simplify=FALSE)
figure6 <- 
  discreteQ(
    y = dep,
    d = white,
    x = x,
    ys = ys,
    decomposition = TRUE,
    q.range = c(0.03, 0.97),
    method = "logit",
    bsrep = 1000,
    alpha = 0.05,
    cl = cl,
    return.boot = TRUE,
    list_of_seeds = list_of_seeds
  )


### Graphs

pdf(paste0(code.dir, "/Results/figure6.pdf"), pointsize=15,width=10.0,height=10.5)
par(mfrow=c(1,1), lend="butt", mar=c(2.1,4.1,4.1,2.1))
layout(matrix(c(1,3,2,4), 2),c(1,1),c(0.9,1))

plot(figure6, which="Qc", col.l="black", col.b=adjustcolor("gray", alpha.f = 0.5), ylab="Test scores", xlab="", ylim= c(-2.5,2.39), main="Quantile functions")
plot(figure6, which="Q0", col.l="dark blue", col.b=adjustcolor("light blue", alpha.f = 0.5), add=TRUE)
plot(figure6, which="Q1", col.l="dark green", col.b=adjustcolor("light green", alpha.f = 0.5), add=TRUE)
legend(0.03, max(ys), c(' ',' ',' '), col = c(adjustcolor("light blue", alpha.f = 0.5),adjustcolor("light green", alpha.f = 0.5),adjustcolor("gray", alpha.f = 0.5)), lwd = c(4,4,4), horiz = FALSE, bty = 'n', cex=0.75)
legend(0.03, max(ys), c('Observed black quantiles', 'Observed white quantiles','Counterfact. quantiles'), col = c('dark blue','dark green','black'), lwd = c(1,1,1), horiz = FALSE, bty = 'n', cex=0.75)

plot(figure6, which="observed", xlab="Quantile index", ylab="Difference in test scores", ylim=c(-0.4,1.2), main="Black white observed gap")

par(mar=c(5.1,4.1,4.1,2.1))
plot(figure6, which="composition", xlab="Probability", ylab="Difference in test scores", ylim=c(-0.4,1.2), main="Composition effect")

plot(figure6, which="unexplained", xlab="Probability", ylab="Difference in test scores", ylim=c(-0.4,1.2), main="Unexplained difference")

dev.off()
