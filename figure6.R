#######################################################################################;

# This application performs a decomposition of the black-white testscore gap based on
# "Testing for Racial Differences in the Mental Ability of Young Children" 
# by R.G. Fryer and S.D. Levitt

# Authors: V. Chernozhukov, I. Fernandez-Val, B. Melly, K. Wuethrich

# Data source: AER webpage, CPP dataset

#######################################################################################;


library(boot)
library(foreign)
library(gdata)
library(MASS)
library(pscl)
library(matrixStats)
library(snow)

#Folder containing the dataset. CHANGE THE PATH!
setwd("D:/Dropbox/papers/DR")
rm(list = ls());


###############
# Functions
###############

cdf<-  function(ys, Fs){;
  ys<- sort(ys);
  Fs<- sort(Fs);
  F<- stepfun(ys, c(0,Fs));
  return(F);
};

left.inv<- function(ys, Fs) {;
  ys<- sort(ys);
  Fs<- sort(Fs);  
  iF<- stepfun(Fs, c(ys,max(ys)), right=TRUE);
  return(iF);
};

# Fuction to estimate empirical distribution;

ecdf.7y <- function(y, data, group)
{;
  Ff    <- weighted.mean((data$stand_fullscale_iq_7years[data$black==group] <= y),data$weight[data$black==group]); 
  return(Ff);
};

# Fuctions to estimate counterfactual distribution using distribution regression;

counter.7y.lpm <- function(y, data)
{;
  form_basic <- "age_7years.f+female"
  form_ses <- "(dad_hs_dropout+dad_hs_grad+dad_some_college+dad_college_plus+dad_no_occupation+dad_professional+dad_non_professional+mom_hs_dropout+mom_hs_grad+mom_some_college+mom_college_plus+mom_no_occupation+mom_professional+mom_non_professional+inc_less_500+inc_500_1000+inc_1000_1500+inc_1500_2000+inc_2000_2500+inc_2500_plus)"
  form_home_environment <- "siblings_0+siblings_1+siblings_2+siblings_3+siblings_4+siblings_5+siblings_6_plus+both_bio_parents+age_mom+age_mom_2 +age_mom_3+age_mom_4+age_mom_5+miss_age_mom+mother_indifferent+mother_accepting+mother_attentive+mother_over_caring+mother_other+miss_parental_score"
  form_prenatal <- "w_less_1500+w_1500_2500+w_2500_3500+w_3500_more+weeks_premature_0+weeks_premature_1+weeks_premature_2+weeks_premature_3+weeks_premature_4+weeks_premature_5+ weeks_premature_6+weeks_premature_7+weeks_premature_8+weeks_premature_9+weeks_premature_10+weeks_premature_11 +singleton+twin+high_order_multiple"
  form <-  as.formula(paste("I(stand_fullscale_iq_7years<= y)", "~", form_basic,"+",form_ses,"+",form_home_environment,"+",form_prenatal,"+interviewer_7years.f"));  
  fit   <- lm(form, data = data, subset=data$black==0,weights=data$weights);
  Fc  <- weighted.mean(predict(fit, newdata = data[data$black==1, ]),data[data$black==1,"weight"]); 
  return(Fc);
};
counter.7y <- function(y, data)
{;
  form_basic <- "age_7years.f+female"
  form_ses <- "(dad_hs_dropout+dad_hs_grad+dad_some_college+dad_college_plus+dad_no_occupation+dad_professional+dad_non_professional+mom_hs_dropout+mom_hs_grad+mom_some_college+mom_college_plus+mom_no_occupation+mom_professional+mom_non_professional+inc_less_500+inc_500_1000+inc_1000_1500+inc_1500_2000+inc_2000_2500+inc_2500_plus)"
  form_home_environment <- "siblings_0+siblings_1+siblings_2+siblings_3+siblings_4+siblings_5+siblings_6_plus+both_bio_parents+age_mom+age_mom_2 +age_mom_3+age_mom_4+age_mom_5+miss_age_mom+mother_indifferent+mother_accepting+mother_attentive+mother_over_caring+mother_other+miss_parental_score"
  form_prenatal <- "w_less_1500+w_1500_2500+w_2500_3500+w_3500_more+weeks_premature_0+weeks_premature_1+weeks_premature_2+weeks_premature_3+weeks_premature_4+weeks_premature_5+ weeks_premature_6+weeks_premature_7+weeks_premature_8+weeks_premature_9+weeks_premature_10+weeks_premature_11 +singleton+twin+high_order_multiple"
  form <-  as.formula(paste("I(stand_fullscale_iq_7years<= y)", "~", form_basic,"+",form_ses,"+",form_home_environment,"+",form_prenatal,"+interviewer_7years.f"));  
  fit   <- glm(form, data = data, subset=data$black==0,weights=data$weights, family = binomial(link = "logit"));
  Fc  <- weighted.mean(predict(fit,  type = "response", newdata = data[data$black==1, ]),data[data$black==1,"weight"]); 
  return(Fc);
};

##############################################################################
# Data
# Prior sample selection: drop all observations with hispanic=1 and other=1
# Code works if black=0 only if white=1
##############################################################################

############################
# Specification Analysis
############################

set.seed(123)
ql <- 0.01
qu <- 0.99
alpha <- 0.05;
nrep <- 1000

############################
### Analysis 7 years
############################

data  <- read.dta("Data/cpp_selected.dta");
#drop if interviewer interviewd only blacks
for(i in setdiff(unique(data$interviewer_7years[data$black==1]), unique(data$interviewer_7years[data$black==0]))){
  data <- data[data$interviewer_7years!=i,]
}
n     <- nrow(data); 

# Factor variables (to shorten formulas)
data$age_7years.f <- factor(data$age_7years) 

data$interviewer_7years.f <- factor(data$interviewer_7years)

data$weight <- 1

### Point estimates
# Multicore
cl <- makeCluster(7)

ys.7y <- sort(unique(data$stand_fullscale_iq_7years))
ys.7y <- ys.7y[ys.7y>=quantile(data$stand_fullscale_iq_7years,ql) & ys.7y<=quantile(data$stand_fullscale_iq_7years,qu)]

Fb.7y <- sapply(ys.7y, ecdf.7y, data=data, group=1)
Fw.7y <- sapply(ys.7y, ecdf.7y, data=data, group=0)
Fc.7y <- sort(unlist(clusterApply(cl,ys.7y,counter.7y,data=data)))

### Bootstrapped CI

Fb.7y.b <- Fw.7y.b <- Fc.7y.b <- matrix(NA,length(ys.7y),nrep)
for(r in 1:nrep){
  print(r)
  data$weight <- rexp(n)
  Fb.7y.b[,r] <- sapply(ys.7y, ecdf.7y, data=data, group=1)
  Fw.7y.b[,r] <- sapply(ys.7y, ecdf.7y, data=data, group=0)
  Fc.7y.b[,r] <- sort(unlist(clusterApply(cl,ys.7y,counter.7y,data=data)))
}

delta.b.7y    <- Fb.7y.b - Fb.7y;
variance.b.7y <- apply(Fb.7y.b,1,FUN=function(x) IQR(x)/1.349)^2
delta.w.7y    <- Fw.7y.b - Fw.7y;
variance.w.7y <- apply(Fw.7y.b,1,FUN=function(x) IQR(x)/1.349)^2
delta.c.7y    <- Fc.7y.b-Fc.7y;
variance.c.7y <- apply(Fc.7y.b,1,FUN=function(x) IQR(x)/1.349)^2

zs.j.7y<- apply(rbind(abs(delta.b.7y) /sqrt(variance.b.7y),abs(delta.w.7y) /sqrt(variance.w.7y),abs(delta.c.7y) /sqrt(variance.c.7y)),  2, max, na.rm = TRUE)  #max abs t-stat
crt.j.7y<- quantile(zs.j.7y, 1-alpha)  #critical value
zs.b.7y<- apply(abs(delta.b.7y) /sqrt(variance.b.7y),  2, max, na.rm = TRUE)  #max abs t-stat
crt.b.7y<- quantile(zs.b.7y, 1-alpha)  #critical value
zs.w.7y<- apply(abs(delta.w.7y) /sqrt(variance.w.7y),  2, max, na.rm = TRUE)  #max abs t-stat
crt.w.7y<- quantile(zs.w.7y, 1-alpha)  #critical value
zs.c.7y<- apply(abs(delta.c.7y) /sqrt(variance.c.7y),  2, max, na.rm = TRUE)  #max abs t-stat
crt.c.7y<- quantile(zs.c.7y, 1-alpha)  #critical value

# Algorithm 1
ub.Fb.7y<-  sort(Fb.7y + crt.b.7y*sqrt(variance.b.7y))  
lb.Fb.7y<-  sort(Fb.7y - crt.b.7y*sqrt(variance.b.7y))
ub.Fb.7y<-  ifelse(ub.Fb.7y <= 1, ub.Fb.7y, 1);  #imposing support restriction
lb.Fb.7y<-  ifelse(lb.Fb.7y >= 0, lb.Fb.7y, 0);  #imposing support restriction

ub.Fw.7y<-  sort(Fw.7y + crt.w.7y*sqrt(variance.w.7y))
lb.Fw.7y<-  sort(Fw.7y - crt.w.7y*sqrt(variance.w.7y))
ub.Fw.7y<-  ifelse(ub.Fw.7y <= 1, ub.Fw.7y, 1);  #imposing support restriction
lb.Fw.7y<-  ifelse(lb.Fw.7y >= 0, lb.Fw.7y, 0);  #imposing support restriction

ub.Fc.7y<-  sort(Fc.7y + crt.c.7y*sqrt(variance.c.7y))
lb.Fc.7y<-  sort(Fc.7y - crt.c.7y*sqrt(variance.c.7y))
ub.Fc.7y<-  ifelse(ub.Fc.7y <= 1, ub.Fc.7y, 1);  #imposing support restriction
lb.Fc.7y<-  ifelse(lb.Fc.7y >= 0, lb.Fc.7y, 0);  #imposing support restriction


# Algorithm 2
ub.Fbj.7y<-  sort(Fb.7y + crt.j.7y*sqrt(variance.b.7y))  
lb.Fbj.7y<-  sort(Fb.7y - crt.j.7y*sqrt(variance.b.7y))
ub.Fbj.7y<-  ifelse(ub.Fb.7y <= 1, ub.Fb.7y, 1);  #imposing support restriction
lb.Fbj.7y<-  ifelse(lb.Fb.7y >= 0, lb.Fb.7y, 0);  #imposing support restriction

ub.Fwj.7y<-  sort(Fw.7y + crt.j.7y*sqrt(variance.w.7y))
lb.Fwj.7y<-  sort(Fw.7y - crt.j.7y*sqrt(variance.w.7y))
ub.Fwj.7y<-  ifelse(ub.Fw.7y <= 1, ub.Fw.7y, 1);  #imposing support restriction
lb.Fwj.7y<-  ifelse(lb.Fw.7y >= 0, lb.Fw.7y, 0);  #imposing support restriction

ub.Fcj.7y<-  sort(Fc.7y + crt.j.7y*sqrt(variance.c.7y))
lb.Fcj.7y<-  sort(Fc.7y - crt.j.7y*sqrt(variance.c.7y))
ub.Fcj.7y<-  ifelse(ub.Fc.7y <= 1, ub.Fc.7y, 1);  #imposing support restriction
lb.Fcj.7y<-  ifelse(lb.Fc.7y >= 0, lb.Fc.7y, 0);  #imposing support restriction

### Create step functions for distributions and quantiles;

# Algorithm 1
Fb.7y.func <- cdf(ys.7y,Fb.7y)
ub.Fb.7y.func <- cdf(ys.7y,lb.Fb.7y)
lb.Fb.7y.func <- cdf(ys.7y,ub.Fb.7y)

Fw.7y.func <- cdf(ys.7y,Fw.7y)
ub.Fw.7y.func <- cdf(ys.7y,lb.Fw.7y)
lb.Fw.7y.func <- cdf(ys.7y,ub.Fw.7y)

Fc.7y.func <- cdf(ys.7y,Fc.7y)
ub.Fc.7y.func <- cdf(ys.7y,lb.Fc.7y)
lb.Fc.7y.func <- cdf(ys.7y,ub.Fc.7y)

Qb.7y.func   <- left.inv(ys.7y, Fb.7y);
ub.Qb.7y.func  <- left.inv(ys.7y, lb.Fb.7y);
lb.Qb.7y.func  <- left.inv(ys.7y, ub.Fb.7y);

Qw.7y.func   <- left.inv(ys.7y, Fw.7y);
ub.Qw.7y.func  <- left.inv(ys.7y, lb.Fw.7y);
lb.Qw.7y.func  <- left.inv(ys.7y, ub.Fw.7y);

Qc.7y.func   <- left.inv(ys.7y, Fc.7y);
ub.Qc.7y.func  <- left.inv(ys.7y, lb.Fc.7y);
lb.Qc.7y.func  <- left.inv(ys.7y, ub.Fc.7y);

# Algorithm 2

ub.Fbj.7y.func <- cdf(ys.7y,ub.Fbj.7y)
lb.Fbj.7y.func <- cdf(ys.7y,lb.Fbj.7y)
ub.Fwj.7y.func <- cdf(ys.7y,ub.Fwj.7y)
lb.Fwj.7y.func <- cdf(ys.7y,lb.Fwj.7y)
ub.Fcj.7y.func <- cdf(ys.7y,ub.Fcj.7y)
lb.Fcj.7y.func <- cdf(ys.7y,lb.Fcj.7y)

ub.Qbj.7y.func  <- left.inv(ys.7y, lb.Fbj.7y);
lb.Qbj.7y.func  <- left.inv(ys.7y, ub.Fbj.7y);
ub.Qwj.7y.func  <- left.inv(ys.7y, lb.Fwj.7y);
lb.Qwj.7y.func  <- left.inv(ys.7y, ub.Fwj.7y);
ub.Qcj.7y.func  <- left.inv(ys.7y, lb.Fcj.7y);
lb.Qcj.7y.func  <- left.inv(ys.7y, ub.Fcj.7y);



### Graphs

qlg <- max(min(ub.Fb.7y),min(ub.Fw.7y),min(ub.Fc.7y),0.03)
qug <- min(max(lb.Fb.7y),max(lb.Fw.7y),max(lb.Fc.7y),0.97)

# QFs
pdf("Results/Test-Scores/figure6.pdf", pointsize=15,width=10.0,height=10.5);
par(mfrow=c(1,1), lend="butt", mar=c(2.1,4.1,4.1,2.1))
layout(matrix(c(1,3,2,4), 2),c(1,1),c(0.9,1))

#redefine the quantile functions to avoid seeing the tails.
Qb.7y.func  <- left.inv(c(-999,ys.7y,999), c(min(Fb.7y)-.Machine$double.eps,Fb.7y,max(Fb.7y)+.Machine$double.eps))
Qw.7y.func  <- left.inv(c(-999,ys.7y,999), c(min(Fw.7y)-.Machine$double.eps,Fw.7y,max(Fw.7y)+.Machine$double.eps))
Qc.7y.func  <- left.inv(c(-999,ys.7y,999), c(min(Fc.7y)-.Machine$double.eps,Fc.7y,max(Fc.7y)+.Machine$double.eps))

plot(Qb.7y.func, xval=c(qlg,Fb.7y[Fb.7y>=qlg & Fb.7y<=qug],qug), xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", ylab="Test scores", xlab="Quantile index", 
     ylim= c(-2.5,2.39), main="Quantile functions",
     sub=" ");
for(i in 2:length(ys.7y)) if(ub.Fb.7y[i]>qlg & lb.Fb.7y[i-1]<qug) segments(max(qlg,lb.Fb.7y[i-1]),ys.7y[i],min(qug,ub.Fb.7y[i]),ys.7y[i],col="light blue", lty = 1,lwd=5) 
for(i in 2:length(ys.7y)) if(ub.Fw.7y[i]>qlg & lb.Fw.7y[i-1]<qug) segments(max(qlg,lb.Fw.7y[i-1]),ys.7y[i],min(qug,ub.Fw.7y[i]),ys.7y[i],col="light green", lty = 1,lwd=5) 
for(i in 2:length(ys.7y)) if(ub.Fc.7y[i]>qlg & lb.Fc.7y[i-1]<qug) segments(max(qlg,lb.Fc.7y[i-1]),ys.7y[i],min(qug,ub.Fc.7y[i]),ys.7y[i],col="grey", lty = 1,lwd=5) 
lines(Qb.7y.func, xval=c(qlg,Fb.7y[Fb.7y>=qlg & Fb.7y<=qug],qug), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1);
lines(Qw.7y.func, xval=c(qlg,Fw.7y[Fw.7y>=qlg & Fw.7y<=qug],qug), verticals=FALSE, do.points=FALSE, col="dark green", lty = 1);
lines(Qc.7y.func, xval=c(qlg,Fc.7y[Fc.7y>=qlg & Fc.7y<=qug],qug), verticals=FALSE, do.points=FALSE, col="black", lty = 1);

legend(qlg, max(ys.7y), c(' ',' ',' '), col = c('light blue','light green','grey'), lwd = c(4,4,4), horiz = FALSE, bty = 'n', cex=0.75);
legend(qlg, max(ys.7y), c('Observed black quantiles', 'Observed white quantiles','Counterfact. quantiles'), col = c('dark blue','dark green','black'), lwd = c(1,1,1), horiz = FALSE, bty = 'n', cex=0.75);

# Quantile effects

Qb.7y.func  <- left.inv(ys.7y, Fb.7y)
Qw.7y.func  <- left.inv(ys.7y, Fw.7y)
Qc.7y.func  <- left.inv(ys.7y, Fc.7y)

quant <- unique(c(Fb.7y,Fw.7y))
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
QTE <- Qw.7y.func(quant)-Qb.7y.func(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-10,QTE,10),right=TRUE)
plot(QTE.func,xval=c(qlg,quant[quant>qlg & quant<qug],qug),xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Quantile index", ylab="Difference in test scores", 
     ylim=c(-0.4,1.2), main="Black white observed gap",
     sub=" ");
quant <- c(lb.Fbj.7y,lb.Fwj.7y,ub.Fbj.7y,ub.Fwj.7y)
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
quant <- c(qlg,quant[quant>=qlg & quant<=qug],qug)
all.possible <- c()
for(i in 1:length(ys.7y)){
  for(j in 1:length(ys.7y)){
    all.possible <- c(all.possible,ys.7y[i]-ys.7y[j])
  }
}
all.possible <- sort(unique(all.possible))
for(v in all.possible){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=ub.Qwj.7y.func(quant[x]) - lb.Qbj.7y.func(quant[x]) & v >=lb.Qwj.7y.func(quant[x]) - ub.Qbj.7y.func(quant[x])
  }
  if(max(accepted)==1){
    temp_segments <- 1-accepted
    temp_segments[2:length(temp_segments)][temp_segments[1:(length(temp_segments)-1)]==1] <- 0
    temp_segments <- accepted*(cumsum(temp_segments)+accepted[1])
    for(s in 1:max(temp_segments)) segments(quant[temp_segments==s][1],v,tail(quant[temp_segments==s],1),v,col="light blue", lty = 1,lwd=5)
  }
}
lines(QTE.func, xval=c(qlg,knots(QTE.func)[knots(QTE.func)>=qlg & knots(QTE.func)<=qug],qug), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1, lwd=1);

par(mar=c(5.1,4.1,4.1,2.1))
quant <- unique(c(Fc.7y,Fw.7y))
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
QTE <- Qw.7y.func(quant)-Qc.7y.func(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-10,QTE,10),right=TRUE)
plot(QTE.func,xval=c(qlg,quant[quant>=qlg & quant<=qug],qug), xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Probability", ylab="Difference in test scores", 
     ylim=c(-0.4,1.2), main="Composition effect",
     sub=" ");
quant <- c(lb.Fcj.7y,lb.Fwj.7y,ub.Fcj.7y,ub.Fwj.7y)
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
quant <- c(qlg,quant[quant>=qlg & quant<=qug],qug)
for(v in all.possible){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=ub.Qwj.7y.func(quant[x]) - lb.Qcj.7y.func(quant[x]) & v >=lb.Qwj.7y.func(quant[x]) - ub.Qcj.7y.func(quant[x])
  }
  if(max(accepted)==1){
    temp_segments <- 1-accepted
    temp_segments[2:length(temp_segments)][temp_segments[1:(length(temp_segments)-1)]==1] <- 0
    temp_segments <- accepted*(cumsum(temp_segments)+accepted[1])
    for(s in 1:max(temp_segments)) segments(quant[temp_segments==s][1],v,tail(quant[temp_segments==s],1),v,col="light blue", lty = 1,lwd=5)
  }
}
lines(QTE.func, xval=c(qlg,knots(QTE.func)[knots(QTE.func)>=qlg & knots(QTE.func)<=qug],qug), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1, lwd=1);

quant <- unique(c(Fc.7y,Fb.7y))
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
QTE <- Qc.7y.func(quant)-Qb.7y.func(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-10,QTE,10),right=TRUE)
plot(QTE.func,xval=c(qlg,quant[quant>=qlg & quant<=qug],qug), xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Probability", ylab="Difference in test scores", 
     ylim=c(-0.4,1.2), main="Unexplained difference",
     sub=" ");
quant <- c(lb.Fcj.7y,lb.Fbj.7y,ub.Fcj.7y,ub.Fbj.7y)
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
quant <- c(qlg,quant[quant>=qlg & quant<=qug],qug)
for(v in all.possible){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=ub.Qcj.7y.func(quant[x]) - lb.Qbj.7y.func(quant[x]) & v >=lb.Qcj.7y.func(quant[x]) - ub.Qbj.7y.func(quant[x])
  }
  if(max(accepted)==1){
    temp_segments <- 1-accepted
    temp_segments[2:length(temp_segments)][temp_segments[1:(length(temp_segments)-1)]==1] <- 0
    temp_segments <- accepted*(cumsum(temp_segments)+accepted[1])
    for(s in 1:max(temp_segments)) segments(quant[temp_segments==s][1],v,tail(quant[temp_segments==s],1),v,col="light blue", lty = 1,lwd=5)
  }
}
lines(QTE.func, xval=c(qlg,knots(QTE.func)[knots(QTE.func)>=qlg & knots(QTE.func)<=qug],qug), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1, lwd=1);

dev.off();
detach(data)

