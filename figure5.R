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

ecdf.8m <- function(y, data, group)
{;
  Ff    <- weighted.mean((data$stand_mental_score_8months[data$black==group] <= y),data$weight[data$black==group]); 
  return(Ff);
};

# Fuctions to estimate counterfactual distribution using distribution regression;
counter.8m.lpm <- function(y, data)
{;
  form_basic <- "age_8months.f+female"
  form_ses <- "(dad_hs_dropout+dad_hs_grad+dad_some_college+dad_college_plus+dad_no_occupation+dad_professional+dad_non_professional+mom_hs_dropout+mom_hs_grad+mom_some_college+mom_college_plus+mom_no_occupation+mom_professional+mom_non_professional+inc_less_500+inc_500_1000+inc_1000_1500+inc_1500_2000+inc_2000_2500+inc_2500_plus)"
  form_home_environment <- "siblings_0+siblings_1+siblings_2+siblings_3+siblings_4+siblings_5+siblings_6_plus+both_bio_parents+age_mom+age_mom_2 +age_mom_3+age_mom_4+age_mom_5+miss_age_mom+mother_indifferent+mother_accepting+mother_attentive+mother_over_caring+mother_other+miss_parental_score"
  form_prenatal <- "w_less_1500+w_1500_2500+w_2500_3500+w_3500_more+weeks_premature_0+weeks_premature_1+weeks_premature_2+weeks_premature_3+weeks_premature_4+weeks_premature_5+ weeks_premature_6+weeks_premature_7+weeks_premature_8+weeks_premature_9+weeks_premature_10+weeks_premature_11 +singleton+twin+high_order_multiple"
  form <-  as.formula(paste("I(stand_mental_score_8months<= y)", "~", form_basic,"+",form_ses,"+",form_home_environment,"+",form_prenatal,"+interviewer_8months.f"));  
  fit   <- lm(form, data = data, subset=data$black==0,weights=data$weights);
  Fc  <- weighted.mean(predict(fit, newdata = data[data$black==1, ]),data[data$black==1,"weight"]); 
  return(Fc);
};

# Fuctions to estimate counterfactual distribution using distribution regression;
counter.8m <- function(y, data)
{;
  form_basic <- "age_8months.f+female"
  form_ses <- "(dad_hs_dropout+dad_hs_grad+dad_some_college+dad_college_plus+dad_no_occupation+dad_professional+dad_non_professional+mom_hs_dropout+mom_hs_grad+mom_some_college+mom_college_plus+mom_no_occupation+mom_professional+mom_non_professional+inc_less_500+inc_500_1000+inc_1000_1500+inc_1500_2000+inc_2000_2500+inc_2500_plus)"
  form_home_environment <- "siblings_0+siblings_1+siblings_2+siblings_3+siblings_4+siblings_5+siblings_6_plus+both_bio_parents+age_mom+age_mom_2 +age_mom_3+age_mom_4+age_mom_5+miss_age_mom+mother_indifferent+mother_accepting+mother_attentive+mother_over_caring+mother_other+miss_parental_score"
  form_prenatal <- "w_less_1500+w_1500_2500+w_2500_3500+w_3500_more+weeks_premature_0+weeks_premature_1+weeks_premature_2+weeks_premature_3+weeks_premature_4+weeks_premature_5+ weeks_premature_6+weeks_premature_7+weeks_premature_8+weeks_premature_9+weeks_premature_10+weeks_premature_11 +singleton+twin+high_order_multiple"
  form <-  as.formula(paste("I(stand_mental_score_8months<= y)", "~", form_basic,"+",form_ses,"+",form_home_environment,"+",form_prenatal,"+interviewer_8months.f"));  
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
### Analysis 8 months
############################

data  <- read.dta("Data/cpp_selected.dta");
#drop if interviewer interviewd only blacks
for(i in setdiff(unique(data$interviewer_8months[data$black==1]), unique(data$interviewer_8months[data$black==0]))){
  data <- data[data$interviewer_8months!=i,]
}
n     <- nrow(data); 

# Factor variables (to shorten formulas)
data$age_8months.f <- factor(data$age_8months) 

data$interviewer_8months.f <- factor(data$interviewer_8months)

data$weight <- 1

### Point estimates
# Multicore
cl <- makeCluster(7)

ys.8m <- sort(unique(data$stand_mental_score_8months))
ys.8m <- ys.8m[ys.8m>=quantile(data$stand_mental_score_8months,ql) & ys.8m<=quantile(data$stand_mental_score_8months,qu)]

Fb.8m <- sapply(ys.8m, ecdf.8m, data=data, group=1)
Fw.8m <- sapply(ys.8m, ecdf.8m, data=data, group=0)
Fc.8m <- sort(unlist(clusterApply(cl,ys.8m,counter.8m,data=data)))

### Bootstrapped CI

Fb.8m.b <- Fw.8m.b <- Fc.8m.b <- matrix(NA,length(ys.8m),nrep)
for(r in 1:nrep){
  print(r)
  data$weight <- rexp(n)
  Fb.8m.b[,r] <- sapply(ys.8m, ecdf.8m, data=data, group=1)
  Fw.8m.b[,r] <- sapply(ys.8m, ecdf.8m, data=data, group=0)
  Fc.8m.b[,r] <- sort(unlist(clusterApply(cl,ys.8m,counter.8m,data=data)))
}

delta.b.8m    <- Fb.8m.b - Fb.8m;
variance.b.8m <- apply(Fb.8m.b,1,FUN=function(x) IQR(x)/1.349)^2
delta.w.8m    <- Fw.8m.b - Fw.8m;
variance.w.8m <- apply(Fw.8m.b,1,FUN=function(x) IQR(x)/1.349)^2
delta.c.8m    <- Fc.8m.b-Fc.8m;
variance.c.8m <- apply(Fc.8m.b,1,FUN=function(x) IQR(x)/1.349)^2

zs.j.8m<- apply(rbind(abs(delta.b.8m) /sqrt(variance.b.8m),abs(delta.w.8m) /sqrt(variance.w.8m),abs(delta.c.8m) /sqrt(variance.c.8m)),  2, max, na.rm = TRUE)  #max abs t-stat
crt.j.8m<- quantile(zs.j.8m, 1-alpha)  #critical value
zs.b.8m<- apply(abs(delta.b.8m) /sqrt(variance.b.8m),  2, max, na.rm = TRUE)  #max abs t-stat
crt.b.8m<- quantile(zs.b.8m, 1-alpha)  #critical value
zs.w.8m<- apply(abs(delta.w.8m) /sqrt(variance.w.8m),  2, max, na.rm = TRUE)  #max abs t-stat
crt.w.8m<- quantile(zs.w.8m, 1-alpha)  #critical value
zs.c.8m<- apply(abs(delta.c.8m) /sqrt(variance.c.8m),  2, max, na.rm = TRUE)  #max abs t-stat
crt.c.8m<- quantile(zs.c.8m, 1-alpha)  #critical value

# Algorithm 1
ub.Fb.8m<-  sort(Fb.8m + crt.b.8m*sqrt(variance.b.8m))  
lb.Fb.8m<-  sort(Fb.8m - crt.b.8m*sqrt(variance.b.8m))
ub.Fb.8m<-  ifelse(ub.Fb.8m <= 1, ub.Fb.8m, 1);  #imposing support restriction
lb.Fb.8m<-  ifelse(lb.Fb.8m >= 0, lb.Fb.8m, 0);  #imposing support restriction

ub.Fw.8m<-  sort(Fw.8m + crt.w.8m*sqrt(variance.w.8m))
lb.Fw.8m<-  sort(Fw.8m - crt.w.8m*sqrt(variance.w.8m))
ub.Fw.8m<-  ifelse(ub.Fw.8m <= 1, ub.Fw.8m, 1);  #imposing support restriction
lb.Fw.8m<-  ifelse(lb.Fw.8m >= 0, lb.Fw.8m, 0);  #imposing support restriction

ub.Fc.8m<-  sort(Fc.8m + crt.c.8m*sqrt(variance.c.8m))
lb.Fc.8m<-  sort(Fc.8m - crt.c.8m*sqrt(variance.c.8m))
ub.Fc.8m<-  ifelse(ub.Fc.8m <= 1, ub.Fc.8m, 1);  #imposing support restriction
lb.Fc.8m<-  ifelse(lb.Fc.8m >= 0, lb.Fc.8m, 0);  #imposing support restriction


# Algorithm 2
ub.Fbj.8m<-  sort(Fb.8m + crt.j.8m*sqrt(variance.b.8m))  
lb.Fbj.8m<-  sort(Fb.8m - crt.j.8m*sqrt(variance.b.8m))
ub.Fbj.8m<-  ifelse(ub.Fb.8m <= 1, ub.Fb.8m, 1);  #imposing support restriction
lb.Fbj.8m<-  ifelse(lb.Fb.8m >= 0, lb.Fb.8m, 0);  #imposing support restriction

ub.Fwj.8m<-  sort(Fw.8m + crt.j.8m*sqrt(variance.w.8m))
lb.Fwj.8m<-  sort(Fw.8m - crt.j.8m*sqrt(variance.w.8m))
ub.Fwj.8m<-  ifelse(ub.Fw.8m <= 1, ub.Fw.8m, 1);  #imposing support restriction
lb.Fwj.8m<-  ifelse(lb.Fw.8m >= 0, lb.Fw.8m, 0);  #imposing support restriction

ub.Fcj.8m<-  sort(Fc.8m + crt.j.8m*sqrt(variance.c.8m))
lb.Fcj.8m<-  sort(Fc.8m - crt.j.8m*sqrt(variance.c.8m))
ub.Fcj.8m<-  ifelse(ub.Fc.8m <= 1, ub.Fc.8m, 1);  #imposing support restriction
lb.Fcj.8m<-  ifelse(lb.Fc.8m >= 0, lb.Fc.8m, 0);  #imposing support restriction

### Create step functions for distributions and quantiles;

# Algorithm 1
Fb.8m.func <- cdf(ys.8m,Fb.8m)
ub.Fb.8m.func <- cdf(ys.8m,lb.Fb.8m)
lb.Fb.8m.func <- cdf(ys.8m,ub.Fb.8m)

Fw.8m.func <- cdf(ys.8m,Fw.8m)
ub.Fw.8m.func <- cdf(ys.8m,lb.Fw.8m)
lb.Fw.8m.func <- cdf(ys.8m,ub.Fw.8m)

Fc.8m.func <- cdf(ys.8m,Fc.8m)
ub.Fc.8m.func <- cdf(ys.8m,lb.Fc.8m)
lb.Fc.8m.func <- cdf(ys.8m,ub.Fc.8m)

Qb.8m.func   <- left.inv(ys.8m, Fb.8m);
ub.Qb.8m.func  <- left.inv(ys.8m, lb.Fb.8m);
lb.Qb.8m.func  <- left.inv(ys.8m, ub.Fb.8m);

Qw.8m.func   <- left.inv(ys.8m, Fw.8m);
ub.Qw.8m.func  <- left.inv(ys.8m, lb.Fw.8m);
lb.Qw.8m.func  <- left.inv(ys.8m, ub.Fw.8m);

Qc.8m.func   <- left.inv(ys.8m, Fc.8m);
ub.Qc.8m.func  <- left.inv(ys.8m, lb.Fc.8m);
lb.Qc.8m.func  <- left.inv(ys.8m, ub.Fc.8m);

# Algorithm 2

ub.Fbj.8m.func <- cdf(ys.8m,ub.Fbj.8m)
lb.Fbj.8m.func <- cdf(ys.8m,lb.Fbj.8m)
ub.Fwj.8m.func <- cdf(ys.8m,ub.Fwj.8m)
lb.Fwj.8m.func <- cdf(ys.8m,lb.Fwj.8m)
ub.Fcj.8m.func <- cdf(ys.8m,ub.Fcj.8m)
lb.Fcj.8m.func <- cdf(ys.8m,lb.Fcj.8m)

ub.Qbj.8m.func  <- left.inv(ys.8m, lb.Fbj.8m);
lb.Qbj.8m.func  <- left.inv(ys.8m, ub.Fbj.8m);
ub.Qwj.8m.func  <- left.inv(ys.8m, lb.Fwj.8m);
lb.Qwj.8m.func  <- left.inv(ys.8m, ub.Fwj.8m);
ub.Qcj.8m.func  <- left.inv(ys.8m, lb.Fcj.8m);
lb.Qcj.8m.func  <- left.inv(ys.8m, ub.Fcj.8m);



### Graphs

qlg <- max(min(ub.Fb.8m),min(ub.Fw.8m),min(ub.Fc.8m),0.03)
qug <- min(max(lb.Fb.8m),max(lb.Fw.8m),max(lb.Fc.8m),0.97)


# QFs
pdf("Results/Test-Scores/figure5.pdf", pointsize=15,width=10.0,height=10.5);
par(mfrow=c(1,1), lend="butt", mar=c(2.1,4.1,4.1,2.1))
layout(matrix(c(1,3,2,4), 2),c(1,1),c(0.9,1))
#redefine the quantile functions to avoid seeing the tails.
Qb.8m.func  <- left.inv(c(-999,ys.8m-0.05,999), c(min(Fb.8m)-.Machine$double.eps,Fb.8m,max(Fb.8m)+.Machine$double.eps))
Qw.8m.func  <- left.inv(c(-999,ys.8m+0.05,999), c(min(Fw.8m)-.Machine$double.eps,Fw.8m,max(Fw.8m)+.Machine$double.eps))
Qc.8m.func  <- left.inv(c(-999,ys.8m,999), c(min(Fc.8m)-.Machine$double.eps,Fc.8m,max(Fc.8m)+.Machine$double.eps))

plot(Qb.8m.func, xval=c(qlg,Fb.8m[Fb.8m>=qlg & Fb.8m<=qug],qug), xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", ylab="Test scores", xlab="", 
     ylim= c(-2.5,2.39), main="Quantile functions",
     sub=" ")
for(i in 2:length(ys.8m)) if(ub.Fb.8m[i]>qlg & lb.Fb.8m[i-1]<qug) segments(max(qlg,lb.Fb.8m[i-1]),ys.8m[i]-0.05,min(qug,ub.Fb.8m[i]),ys.8m[i]-0.05,col="light blue", lty = 1,lwd=5) 
for(i in 2:length(ys.8m)) if(ub.Fw.8m[i]>qlg & lb.Fw.8m[i-1]<qug) segments(max(qlg,lb.Fw.8m[i-1]),ys.8m[i]+0.05,min(qug,ub.Fw.8m[i]),ys.8m[i]+0.05,col="light green", lty = 1,lwd=5) 
for(i in 2:length(ys.8m)) if(ub.Fc.8m[i]>qlg & lb.Fc.8m[i-1]<qug) segments(max(qlg,lb.Fc.8m[i-1]),ys.8m[i],min(qug,ub.Fc.8m[i]),ys.8m[i],col="grey", lty = 1,lwd=5) 
lines(Qb.8m.func, xval=c(qlg,Fb.8m[Fb.8m>=qlg & Fb.8m<=qug],qug), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1);
lines(Qw.8m.func, xval=c(qlg,Fw.8m[Fw.8m>=qlg & Fw.8m<=qug],qug), verticals=FALSE, do.points=FALSE, col="dark green", lty = 1);
lines(Qc.8m.func, xval=c(qlg,Fc.8m[Fc.8m>=qlg & Fc.8m<=qug],qug), verticals=FALSE, do.points=FALSE, col="black", lty = 1);
legend(qlg, max(ys.8m), c(' ',' ',' '), col = c('light blue','light green','grey'), lwd = c(4,4,4), horiz = FALSE, bty = 'n', cex=0.75);
legend(qlg, max(ys.8m), c('Observed black quantiles', 'Observed white quantiles','Counterfact. quantiles'), col = c('dark blue','dark green','black'), lwd = c(1,1,1), horiz = FALSE, bty = 'n', cex=0.75);
box()

# Quantile effects
Qb.8m.func  <- left.inv(ys.8m, Fb.8m)
Qw.8m.func  <- left.inv(ys.8m, Fw.8m)
Qc.8m.func  <- left.inv(ys.8m, Fc.8m)

quant <- unique(c(Fb.8m,Fw.8m))
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
QTE <- Qw.8m.func(quant)-Qb.8m.func(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-10,QTE,10),right=TRUE)
plot(QTE.func, xval=c(qlg,quant[quant>=qlg & quant<=qug],qug), xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Quantile index", ylab="Difference in test scores", 
     ylim=c(-0.4,1.2), main="Black-white observed gap",
     sub=" ");
quant <- c(lb.Fbj.8m,lb.Fwj.8m,ub.Fbj.8m,ub.Fwj.8m)
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
quant <- c(qlg,quant[quant>=qlg & quant<=qug],qug)
all.possible <- c()
for(i in 1:length(ys.8m)){
  for(j in 1:length(ys.8m)){
    all.possible <- c(all.possible,ys.8m[i]-ys.8m[j])
  }
}
all.possible <- sort(unique(all.possible))
for(v in all.possible){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=ub.Qwj.8m.func(quant[x]) - lb.Qbj.8m.func(quant[x]) & v >=lb.Qwj.8m.func(quant[x]) - ub.Qbj.8m.func(quant[x])
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
quant <- unique(c(Fc.8m,Fw.8m))
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
QTE <- Qw.8m.func(quant)-Qc.8m.func(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-10,QTE,10),right=TRUE)
plot(QTE.func, xval=c(qlg,quant[quant>=qlg & quant<=qug],qug), xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Probability", ylab="Difference in test scores", 
     ylim=c(-0.4,1.2), main="Composition effect",
     sub=" ");
quant <- c(lb.Fcj.8m,lb.Fwj.8m,ub.Fcj.8m,ub.Fwj.8m)
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
quant <- c(qlg,quant[quant>=qlg & quant<=qug],qug)
for(v in all.possible){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=ub.Qwj.8m.func(quant[x]) - lb.Qcj.8m.func(quant[x]) & v >=lb.Qwj.8m.func(quant[x]) - ub.Qcj.8m.func(quant[x])
  }
  if(max(accepted)==1){
    temp_segments <- 1-accepted
    temp_segments[2:length(temp_segments)][temp_segments[1:(length(temp_segments)-1)]==1] <- 0
    temp_segments <- accepted*(cumsum(temp_segments)+accepted[1])
    for(s in 1:max(temp_segments)) segments(quant[temp_segments==s][1],v,tail(quant[temp_segments==s],1),v,col="light blue", lty = 1,lwd=5)
  }
}
lines(QTE.func, xval=c(qlg,knots(QTE.func)[knots(QTE.func)>=qlg & knots(QTE.func)<=qug],qug), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1, lwd=1);

quant <- unique(c(Fc.8m,Fb.8m))
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
QTE <- Qc.8m.func(quant)-Qb.8m.func(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-10,QTE,10),right=TRUE)
plot(QTE.func, xval=c(qlg,quant[quant>=qlg & quant<=qug],qug), xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Probability", ylab="Difference in test scores", 
     ylim=c(-0.4,1.2), main="Unexplained difference",
     sub=" ");
quant <- c(lb.Fcj.8m,lb.Fbj.8m,ub.Fcj.8m,ub.Fbj.8m)
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps,1))
quant <- c(qlg,quant[quant>=qlg & quant<=qug],qug)
for(v in all.possible){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=ub.Qcj.8m.func(quant[x]) - lb.Qbj.8m.func(quant[x]) & v >=lb.Qcj.8m.func(quant[x]) - ub.Qbj.8m.func(quant[x])
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
