
# This empirical example uses data from the Oregon health experiment

# Authors: V. Chernozhukov, I. Fernandez-Val, B. Melly, K. Wuethrich

# Data source: Oregon Study Group 


library(boot)
library(foreign)
library(MASS)
library(pscl)
library(matrixStats)
library(mgcv)
library(snow)

#Folder containing the dataset. CHANGE THE PATH!
setwd("C:/Users/Melly/Dropbox/papers/DR")
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


# Functions Distribution Regression

uncond.cdfs.dr <- function(y, data)
{;
  form <- I(doc_num_mod_12m<=y) ~ ddddraw_sur_2+ddddraw_sur_3+ddddraw_sur_4+ddddraw_sur_5+ddddraw_sur_6+ddddraw_sur_7+dddnumhh_li_2+dddnumhh_li_3+ddddraXnum_2_2+ddddraXnum_2_3+ddddraXnum_3_2+ddddraXnum_3_3+ddddraXnum_4_2+ddddraXnum_5_2+ddddraXnum_6_2+ddddraXnum_7_2
  fit1  <- glm(form, family = binomial(link = "logit"), data = data[data$treatment==1,], weight=data$weights[data$treatment==1]);
  fit0  <- glm(form, family = binomial(link = "logit"), data = data[data$treatment==0,], weight=data$weights[data$treatment==0]);
  F1  <- weighted.mean(predict(fit1, type = "response",newdata = data), data$weights, na.rm=TRUE);
  F0  <- weighted.mean(predict(fit0, type = "response", newdata = data), data$weights, na.rm=TRUE);
  return(c(F0,F1));
};

# Functions Poisson Regression

uncond.cdfs.po <- function(ys,data){
  form <- doc_num_mod_12m~ddddraw_sur_2+ddddraw_sur_3+ddddraw_sur_4+ddddraw_sur_5+ddddraw_sur_6+ddddraw_sur_7+dddnumhh_li_2+dddnumhh_li_3+ddddraXnum_2_2+ddddraXnum_2_3+ddddraXnum_3_2+ddddraXnum_3_3+ddddraXnum_4_2+ddddraXnum_5_2+ddddraXnum_6_2+ddddraXnum_7_2
  fit0 <- glm(form, family=poisson, data=data, weight=data$weights,subset=(treatment==0))
  fit1 <- glm(form, family=poisson, data=data, weight=data$weights,subset=(treatment==1))
  lambda1 <- predict(fit1, type = "response", newdata = data)
  lambda0 <- predict(fit0, type = "response", newdata = data)
  cond.cdf1 <- sapply(cbind(ys), ppois, lambda=lambda1)
  cond.cdf0 <- sapply(cbind(ys), ppois, lambda=lambda0)
  F1  <- colWeightedMeans(cond.cdf1, w=data$weights)
  F0  <- colWeightedMeans(cond.cdf0, w=data$weights)
  return(cbind(F0,F1))
}


boot.uncond.cdfs.po <- function(data,ys,indices){
  data.b <- data[indices,]
  F <- uncond.cdfs.po(ys,data.b)
  return(F)
}

#Distribution with Poisson link
#Maximum likelihood function
objective <- function(beta,y,binary,x,w=1){         
  lambda <- exp(x%*%beta)
  prob <- pmin(pmax(ppois(y,lambda),10^-15),1-10^-15)
  -sum(w*(binary*log(prob)+(1-binary)*log(1-prob)))  
}

#estimate the unconditional cdf
uncond.cdfs.dp <- function(ys, data, cl){
  dep0 <- data[data$treatment==0,"doc_num_mod_12m"]
  dep1 <- data[data$treatment==1,"doc_num_mod_12m"]
  reg0 <- as.matrix(cbind(1,data[data$treatment==0,c("ddddraw_sur_2","ddddraw_sur_3","ddddraw_sur_4","ddddraw_sur_5","ddddraw_sur_6","ddddraw_sur_7","dddnumhh_li_2","dddnumhh_li_3","ddddraXnum_2_2","ddddraXnum_2_3","ddddraXnum_3_2","ddddraXnum_3_3","ddddraXnum_4_2","ddddraXnum_5_2","ddddraXnum_6_2","ddddraXnum_7_2")]))
  reg1 <- as.matrix(cbind(1,data[data$treatment==1,c("ddddraw_sur_2","ddddraw_sur_3","ddddraw_sur_4","ddddraw_sur_5","ddddraw_sur_6","ddddraw_sur_7","dddnumhh_li_2","dddnumhh_li_3","ddddraXnum_2_2","ddddraXnum_2_3","ddddraXnum_3_2","ddddraXnum_3_3","ddddraXnum_4_2","ddddraXnum_5_2","ddddraXnum_6_2","ddddraXnum_7_2")]))
  wei0 <- data[data$treatment==0,"weights"]
  wei1 <- data[data$treatment==1,"weights"]
  F0 <- F1 <- rep(0,length(ys))
  start <- glm(dep0~reg0-1, weight=wei0, family=poisson)$coef
  fit <- matrix(unlist(clusterApply(cl, ys, fun = function(y, dep, reg, wei, start) optim(start, objective, y=y, binary=(dep<=y), x=reg, w=wei)$par, dep=dep0, reg=reg0, wei=wei0, start=start)), ncol=length(ys))
  F0 <- sapply(1:length(ys), FUN=function(index) weighted.mean(ppois(ys[index],exp(rbind(reg0,reg1)%*%fit[,index])),c(wei0,wei1)))
  start <- glm(dep1~reg1-1, weight=wei1, family=poisson)$coef
  fit <- matrix(unlist(clusterApply(cl, ys, fun = function(y, dep, reg, wei, start) optim(start, objective, y=y, binary=(dep<=y), x=reg, w=wei)$par, dep=dep1, reg=reg1, wei=wei1, start=start)), ncol=length(ys))
  F1 <- sapply(1:length(ys), FUN=function(index) weighted.mean(ppois(ys[index],exp(rbind(reg0,reg1)%*%fit[,index])),c(wei0,wei1)))
  cbind(F0,F1)
}

###############
# Data
###############

data <- read.dta("Data/data_oregon.dta")
data$treatment <- as.numeric(data$treatment)-1 # Destring treatment
data <- na.omit(subset(data, select = c(doc_num_mod_12m,weight_12m,treatment,ddddraw_sur_2,ddddraw_sur_3,ddddraw_sur_4,ddddraw_sur_5,ddddraw_sur_6,ddddraw_sur_7,dddnumhh_li_2,dddnumhh_li_3,ddddraXnum_2_2,ddddraXnum_2_3,ddddraXnum_3_2,ddddraXnum_3_3,ddddraXnum_4_2,ddddraXnum_5_2,ddddraXnum_6_2,ddddraXnum_7_2,household_id)))

###############
# Analysis
###############

### Define forms and threhold values: here you can adjust the outcome
ys <- 0:15
qlg <- 0
qug <- 0.97

### OLS estimates as in Table 5 of Finkelstein et al. 2012
form <- doc_num_mod_12m ~ treatment+ddddraw_sur_2+ddddraw_sur_3+ddddraw_sur_4+ddddraw_sur_5+ddddraw_sur_6+ddddraw_sur_7+dddnumhh_li_2+dddnumhh_li_3+ddddraXnum_2_2+ddddraXnum_2_3+ddddraXnum_3_2+ddddraXnum_3_3+ddddraXnum_4_2+ddddraXnum_5_2+ddddraXnum_6_2+ddddraXnum_7_2
summary(lm(form, data=data, weight=weight_12m))

### Poisson benchmark
form <- doc_num_mod_12m ~ treatment+ddddraw_sur_2+ddddraw_sur_3+ddddraw_sur_4+ddddraw_sur_5+ddddraw_sur_6+ddddraw_sur_7+dddnumhh_li_2+dddnumhh_li_3+ddddraXnum_2_2+ddddraXnum_2_3+ddddraXnum_3_2+ddddraXnum_3_3+ddddraXnum_4_2+ddddraXnum_5_2+ddddraXnum_6_2+ddddraXnum_7_2
glm(form,family=poisson,data=data,weight=weight_12m)

### Point estimates
# Multicore
cl <- makeCluster(7)
clusterExport(cl, "objective")
data$weights <- data$weight_12m
F.po <- uncond.cdfs.po(ys, data)
F1.po <- F.po[,2]
F0.po <- F.po[,1]
F.dr.logit    <- sapply(ys, uncond.cdfs.dr, data = data);
F1.dr.logit <- F.dr.logit[2,]
F0.dr.logit <- F.dr.logit[1,]
F.dr <- uncond.cdfs.dp(ys, data, cl)
F1.dr <- F.dr[,2]
F0.dr <- F.dr[,1]

### Bootstrap confidence intervals;
# Clustered bootstrap (with multicore processing, clusterApply)

set.seed(12345)
nreps <- 1000;
alpha <- 0.05;

F.b1.po <- F.b0.po <- F.b1.dr <- F.b0.dr <- matrix(NA,length(ys),nreps)

household_id <- unique(data$household_id)
nc <- length(household_id)

a <- Sys.time()
for(i in 1:nreps){
  print(i)
  #weighted bootstrap with clustering
  boot.weights <- rexp(nc)
  data <- merge(data,data.frame(cbind(household_id,boot.weights)))
  data$weights <- data$weight_12m*data$boot.weight
  # apply functions
  F.po <- uncond.cdfs.po(ys, data)
  F.b1.po[,i] <- F.po[,2]
  F.b0.po[,i] <- F.po[,1]
  F.dr    <- uncond.cdfs.dp(ys, data, cl)
  F.b1.dr[,i] <- F.dr[,2]
  F.b0.dr[,i] <- F.dr[,1]
  data <- data[,names(data)!="boot.weights"]
}
stopCluster(cl)
Sys.time()-a

# centered/scaled draws and critical values, algorithms 1 and 2
delta1.po  <- F.b1.po - F1.po
delta0.po  <- F.b0.po - F0.po
variance1.po <- apply(F.b1.po,1,FUN=function(x) IQR(x)/1.349)^2
variance0.po <- apply(F.b0.po,1,FUN=function(x) IQR(x)/1.349)^2

select.1 <- (F.b1.po>=qlg)*(rbind(0,F.b1.po[1:(nrow(F.b1.po)-1),])<qug)
select.0 <- (F.b0.po>=qlg)*(rbind(0,F.b0.po[1:(nrow(F.b0.po)-1),])<qug)

#crt for quantile range
zsj.po<- apply(rbind(abs(delta1.po*select.1)/sqrt(variance1.po), abs(delta0.po*select.0) /sqrt(variance0.po)),  2, max, na.rm = TRUE)  #max abs t-stat
crtj.po<- quantile(zsj.po, 1-alpha)  #critical value

delta1.dr  <- F.b1.dr - F1.dr
delta0.dr  <- F.b0.dr - F0.dr
variance1.dr <- apply(F.b1.dr,1,FUN=function(x) IQR(x)/1.349)^2
variance0.dr <- apply(F.b0.dr,1,FUN=function(x) IQR(x)/1.349)^2

select.1 <- (F.b1.dr>=qlg)*(rbind(0,F.b1.dr[1:(nrow(F.b1.dr)-1),])<qug)
select.0 <- (F.b0.dr>=qlg)*(rbind(0,F.b0.dr[1:(nrow(F.b0.dr)-1),])<qug)

zsj.dr<- apply(rbind(abs(delta1.dr*select.1)/sqrt(variance1.dr), abs(delta0.dr*select.0) /sqrt(variance0.dr)),  2, max, na.rm = TRUE)  #max abs t-stat
crtj.dr<- quantile(zsj.dr, 1-alpha)  #critical value

### Confidence bands cdfs

# CDFs and QFs and QTE
ubound.F1j.po <- sort(F1.po + crtj.po*sqrt(variance1.po))
lbound.F1j.po <- sort(F1.po - crtj.po*sqrt(variance1.po))
ubound.F1j.po <-  ifelse(ubound.F1j.po <= 1, ubound.F1j.po, 1)  
lbound.F1j.po <-  ifelse(lbound.F1j.po >= 0, lbound.F1j.po, 0) 

ubound.F0j.po <- sort(F0.po + crtj.po*sqrt(variance0.po))
lbound.F0j.po <- sort(F0.po - crtj.po*sqrt(variance0.po))
ubound.F0j.po <-  ifelse(ubound.F0j.po <= 1, ubound.F0j.po, 1) 
lbound.F0j.po <-  ifelse(lbound.F0j.po >= 0, lbound.F0j.po, 0) 

ubound.F1j.dr <- sort(F1.dr + crtj.dr*sqrt(variance1.dr))
lbound.F1j.dr <- sort(F1.dr - crtj.dr*sqrt(variance1.dr))
ubound.F1j.dr <-  ifelse(ubound.F1j.dr <= 1, ubound.F1j.dr, 1)  
lbound.F1j.dr <-  ifelse(lbound.F1j.dr >= 0, lbound.F1j.dr, 0) 

ubound.F0j.dr <- sort(F0.dr + crtj.dr*sqrt(variance0.dr))
lbound.F0j.dr <- sort(F0.dr - crtj.dr*sqrt(variance0.dr))
ubound.F0j.dr <-  ifelse(ubound.F0j.dr <= 1, ubound.F0j.dr, 1) 
lbound.F0j.dr <-  ifelse(lbound.F0j.dr >= 0, lbound.F0j.dr, 0)

### Create step functions for distributions and quantiles;

# CDF and QFs
F1j.func.po  <- cdf(ys, F1.po);
uF1j.func.po <- cdf(ys, ubound.F1j.po);
lF1j.func.po <- cdf(ys, lbound.F1j.po);

Q1j.func.po  <- left.inv(ys, F1.po)
uQ1j.func.po  <- left.inv(ys, lbound.F1j.po)
lQ1j.func.po  <- left.inv(ys, ubound.F1j.po)

F0j.func.po  <- cdf(ys, F0.po);
uF0j.func.po <- cdf(ys, ubound.F0j.po);
lF0j.func.po <- cdf(ys, lbound.F0j.po);

Q0j.func.po   <- left.inv(ys, F0.po)
uQ0j.func.po  <- left.inv(ys, lbound.F0j.po)
lQ0j.func.po  <- left.inv(ys, ubound.F0j.po)

F1j.func.dr  <- cdf(ys, F1.dr);
uF1j.func.dr <- cdf(ys, ubound.F1j.dr);
lF1j.func.dr <- cdf(ys, lbound.F1j.dr);

Q1j.func.dr  <- left.inv(ys, F1.dr)
uQ1j.func.dr  <- left.inv(ys, lbound.F1j.dr)
lQ1j.func.dr  <- left.inv(ys, ubound.F1j.dr)

F0j.func.dr  <- cdf(ys, F0.dr);
uF0j.func.dr <- cdf(ys, ubound.F0j.dr);
lF0j.func.dr <- cdf(ys, lbound.F0j.dr);

Q0j.func.dr   <- left.inv(ys, F0.dr)
uQ0j.func.dr  <- left.inv(ys, lbound.F0j.dr)
lQ0j.func.dr  <- left.inv(ys, ubound.F0j.dr)

# QTE
uQ1j.func.po  <- left.inv(ys, lbound.F1j.po)
lQ1j.func.po  <- left.inv(ys, ubound.F1j.po)
uQ0j.func.po  <- left.inv(ys, lbound.F0j.po)
lQ0j.func.po  <- left.inv(ys, ubound.F0j.po)

uQ1j.func.dr  <- left.inv(ys, lbound.F1j.dr)
lQ1j.func.dr  <- left.inv(ys, ubound.F1j.dr)
uQ0j.func.dr  <- left.inv(ys, lbound.F0j.dr)
lQ0j.func.dr  <- left.inv(ys, ubound.F0j.dr)

### Graphs

# Histogram

pdf("Results/Oregon/Oregon-hist-doc.pdf", pointsize=15,width=6.0,height=6.0);

par(mfrow=c(1,1));

barplot(prop.table(table(doc_num_mod_12m)),xlab="Outpatient visits (last 6 months)",ylab="Fraction",ylim=c(0,0.4))

dev.off();


# Unconditional CDF

pdf("Results/Oregon/figure4.pdf", pointsize=15,width=10,height=16.5);
par(mfrow=c(1,1), lend="butt", mar=c(5.1,4.1,4.1,2.1))
layout(matrix(c(1,3,5,2,4,6), 3,2),c(1,1),c(1,1,1))

F1.func.po <- cdf(c(-999,ys,max(ys)+.Machine$double.eps),c(0,F1.po,1))
F0.func.po <- cdf(c(-999,ys,max(ys)+.Machine$double.eps),c(0,F0.po,1))

plot(F1.func.po, xlim=c(0,9),verticals=FALSE, do.points=FALSE,col="dark blue", ylab="Probability", xlab="Number of outpatient visits", 
     ylim= c(0,1), main="CDFs - Poisson Model",
     sub=" ");
for(v in 1:length(ys)) polygon(c(ys[v],ys[v+1],ys[v+1],ys[v]),c(lbound.F1j.po[v],lbound.F1j.po[v],ubound.F1j.po[v],ubound.F1j.po[v]),col=adjustcolor("light blue",alpha.f=0.5),border=NA)
for(v in 1:length(ys)) polygon(c(ys[v],ys[v+1],ys[v+1],ys[v]),c(lbound.F0j.po[v],lbound.F0j.po[v],ubound.F0j.po[v],ubound.F0j.po[v]),col=adjustcolor("light green",alpha.f=0.5),border=NA)
segments(-10,0,0,0,col="light green", lty = 1, lwd=5, lend=1)
lines(F1.func.po, ys, verticals=FALSE, do.points=FALSE,col="dark blue", lty = 1);
lines(F0.func.po, ys, verticals=FALSE, do.points=FALSE,col="dark green", lty = 1);
box()
legend(-0.2, 1.05, c(' ', ' '), col = c(adjustcolor("light blue",alpha.f=0.5),adjustcolor("light green",alpha.f=0.5)), lwd = c(4,4,4), horiz = FALSE, bty = 'n');
legend(-0.2, 1.05, c('Treatment group', 'Control group'), col = c('dark blue','dark green'), lwd = c(1,1,1), horiz = FALSE, bty = 'n');


F1.func.dr <- cdf(c(-999,ys,max(ys)+.Machine$double.eps),c(0,F1.dr,11))
F0.func.dr <- cdf(c(-999,ys,max(ys)+.Machine$double.eps),c(0,F0.dr,11))

plot(F1.func.dr, xlim=c(0,9), verticals=FALSE, do.points=FALSE, col="dark blue", ylab="Probability", xlab="Number of outpatient visits", 
      ylim= c(0,1), main="CDFs - Distribution Regression",
      sub=" ");
for(v in 1:(length(ys)-1)) polygon(c(ys[v],ys[v+1],ys[v+1],ys[v]),c(lbound.F1j.dr[v],lbound.F1j.dr[v],ubound.F1j.dr[v],ubound.F1j.dr[v]),col=adjustcolor("light blue",alpha.f=0.5),border=NA)
for(v in 1:(length(ys)-1)) polygon(c(ys[v],ys[v+1],ys[v+1],ys[v]),c(lbound.F0j.dr[v],lbound.F0j.dr[v],ubound.F0j.dr[v],ubound.F0j.dr[v]),col=adjustcolor("light green",alpha.f=0.5),border=NA)
segments(-1,0,0,0,col="light green", lty = 1, lwd=5, lend=1)
lines(F1.func.dr, ys, verticals=FALSE, do.points=FALSE,col="dark blue", lty = 1);
lines(F0.func.dr, ys, verticals=FALSE, do.points=FALSE,col="dark green", lty = 1);
box()
legend(-0.2, 1.05, c(' ', ' '), col = c(adjustcolor("light blue",alpha.f=0.5),adjustcolor("light green",alpha.f=0.5)), lwd = c(4,4,4), horiz = FALSE, bty = 'n');
legend(-0.2, 1.05, c('Treatment group', 'Control group'), col = c('dark blue','dark green'), lwd = c(1,1,1), horiz = FALSE, bty = 'n');


# Unconditional quantile functions
Q1.func.po  <- left.inv(c(-999,ys-0.1,9999), c(0,F1.po,max(F1.po)+.Machine$double.eps))
Q0.func.po  <- left.inv(c(-999,ys+0.1,9999), c(0,F0.po,max(F0.po)+.Machine$double.eps))

plot(Q1.func.po, xval=c(0,F1.po[F1.po<qug],qug), xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", ylab="Number of outpatient visits", xlab="Probability", 
      ylim= c(0,9), main="Quantile functions - Poisson Model",
      sub=" ");
segments(0,-0.1,ubound.F1j.po[1],-0.1,col="light blue", lty = 1,lwd=5, lend=1)
for(i in 2:length(ys)) if(ubound.F1j.po[i]>qlg & lbound.F1j.po[i-1]<qug) segments(lbound.F1j.po[i-1],ys[i]-0.1,min(qug,ubound.F1j.po[i]),ys[i]-0.1,col="light blue", lty = 1,lwd=5, lend=1) 
lines(Q1.func.po, xval=c(0,F1.po[F1.po<qug],qug), verticals=FALSE, do.points=FALSE, col="dark green", lty = 1);
segments(0,0.1,ubound.F0j.po[1],0.1,col="light green", lty = 1,lwd=5)
for(i in 2:length(ys)) if(ubound.F0j.po[i]>qlg & lbound.F0j.po[i-1]<qug) segments(lbound.F0j.po[i-1],ys[i]+0.1,min(qug,ubound.F0j.po[i]),ys[i]+0.1,col="light green", lty = 1,lwd=5, lend=1) 
lines(Q0.func.po, xval=c(0,F0.po[F0.po<qug],qug), verticals=FALSE, do.points=FALSE, col="dark green", lty = 1);

Q1.func.dr  <- left.inv(c(-999,ys-0.1,9999), c(0,F1.dr,max(F1.dr)+.Machine$double.eps))
Q0.func.dr  <- left.inv(c(-999,ys+0.1,9999), c(0,F0.dr,max(F0.dr)+.Machine$double.eps))

plot(Q1.func.dr, xval=c(0,F1.dr[F1.dr<qug],qug), xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", ylab="Number of outpatient visits", xlab="Probablity", 
     ylim= c(0,9), main="Quantile functions - Distribution Regression",
     sub=" ");
segments(0,-0.1,ubound.F1j.dr[1],-0.1,col="light blue", lty = 1,lwd=5, lend=1)
for(i in 2:length(ys)) if(ubound.F1j.dr[i]>qlg & lbound.F1j.dr[i-1]<qug) segments(lbound.F1j.dr[i-1],ys[i]-0.1,min(qug,ubound.F1j.dr[i]),ys[i]-0.1,col="light blue", lty = 1,lwd=5, lend=1) 
lines(Q1.func.dr, xval=c(0,F1.dr[F1.dr<qug],qug), verticals=FALSE, do.points=FALSE, col="dark green", lty = 1);
segments(0,0.1,ubound.F0j.dr[1],0.1,col="light green", lty = 1,lwd=5)
for(i in 2:length(ys)) if(ubound.F0j.dr[i]>qlg & lbound.F0j.dr[i-1]<qug) segments(lbound.F0j.dr[i-1],ys[i]+0.1,min(qug,ubound.F0j.dr[i]),ys[i]+0.1,col="light green", lty = 1,lwd=5, lend=1) 
lines(Q0.func.dr, xval=c(0,F0.dr[F0.dr<qug],qug), verticals=FALSE, do.points=FALSE, col="dark green", lty = 1);

# QTE 

Q1.func.po  <- left.inv(ys, F1.po)
Q0.func.po  <- left.inv(ys, F0.po)

quant <- unique(c(F1.po,F0.po))
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps))
QTE <- Q1.func.po(quant)-Q0.func.po(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-10,QTE,10),right=TRUE)
plot(QTE.func,xval=c(qlg,quant[quant>=qlg & quant<=qug],qug),xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Probability", ylab="Difference in the number of outpatient visits", 
     ylim= c(-2,3), main="Quantile treatment effects - Poisson Model",
     sub=" ");
quant <- c(lbound.F1j.po,lbound.F0j.po,ubound.F1j.po,ubound.F0j.po)
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps))
quant <- c(qlg,quant[quant>=qlg & quant<=qug],qug)
for(v in (min(ys)-max(ys)):(max(ys)-min(ys))){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=uQ1j.func.po(quant[x]) - lQ0j.func.po(quant[x]) & v >=lQ1j.func.po(quant[x]) - uQ0j.func.po(quant[x])
  }
  if(max(accepted)==1){
    temp_segments <- 1-accepted
    temp_segments[2:length(temp_segments)][temp_segments[1:(length(temp_segments)-1)]==1] <- 0
    temp_segments <- accepted*(cumsum(temp_segments)+accepted[1])
    for(s in 1:max(temp_segments)) segments(quant[temp_segments==s][1],v,tail(quant[temp_segments==s],1),v,col="light blue", lty = 1,lwd=10, lend=1)
  }
}
lines(QTE.func, xval=c(qlg,knots(QTE.func)[knots(QTE.func)>qlg & knots(QTE.func)<qug],qug), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1, lwd=1);


Q1.func.dr  <- left.inv(ys, F1.dr)
Q0.func.dr  <- left.inv(ys, F0.dr)

quant <- unique(c(F1.dr,F0.dr))
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps))
QTE <- Q1.func.dr(quant)-Q0.func.dr(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-10,QTE,10),right=TRUE)
plot(QTE.func,xval=c(qlg,quant[quant>=qlg & quant<=qug],qug),xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Probability", ylab="Difference in the number of outpatient visits", 
     ylim= c(-2,3), main="Quantile treatment effects - Distribution Regression",
     sub=" ");
quant <- c(lbound.F1j.dr,lbound.F0j.dr,ubound.F1j.dr,ubound.F0j.dr)
quant <- sort(c(0,quant-.Machine$double.eps,quant+.Machine$double.eps))
quant <- c(qlg,quant[quant>=qlg & quant<=qug],qug)
for(v in (min(ys)-max(ys)):(max(ys)-min(ys))){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=uQ1j.func.dr(quant[x]) - lQ0j.func.dr(quant[x]) & v >=lQ1j.func.dr(quant[x]) - uQ0j.func.dr(quant[x])
  }
  if(max(accepted)==1){
    temp_segments <- 1-accepted
    temp_segments[2:length(temp_segments)][temp_segments[1:(length(temp_segments)-1)]==1] <- 0
    temp_segments <- accepted*(cumsum(temp_segments)+accepted[1])
    for(s in 1:max(temp_segments))  segments(quant[temp_segments==s][1],v,tail(quant[temp_segments==s],1),v,col="light blue", lty = 1,lwd=10, lend=1)
  }
}
lines(QTE.func, xval=c(qlg,knots(QTE.func)[knots(QTE.func)>qlg & knots(QTE.func)<qug],qug), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1, lwd=1);

dev.off()

#Test of equality between Poisson and DR
delta1.dif  <- F.b1.po - F.b1.dr - F1.po + F1.dr
delta0.dif  <- F.b0.po - F.b0.dr - F0.po + F0.dr
variance1.dif <- apply(delta1.dif*delta1.dif,1,mean)
variance0.dif <- apply(delta0.dif*delta0.dif,1,mean)
#test statistic
t1 <- max(abs(F1.po-F1.dr)/sqrt(variance1.dif))
t0 <- max(abs(F0.po-F0.dr)/sqrt(variance0.dif))
#p-values
z1.dif<- apply(abs(delta1.dif) /sqrt(variance1.dif),  2, max, na.rm = TRUE)  #max abs t-stat
mean(z1.dif>t1)
z0.dif<- apply(abs(delta0.dif) /sqrt(variance0.dif),  2, max, na.rm = TRUE)  #max abs t-stat
mean(z0.dif>t0)
