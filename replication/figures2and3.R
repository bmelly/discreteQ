#Folder containing the dataset. CHANGE THE PATH!
setwd("D:/Dropbox/papers/DR")

set.seed(1)
ys<- c(1:10)
Fs<- sort(c(runif(9),1))
F1s<- sort(c (rep(0,5), 1, runif(4)))
bound <- (Fs*(1-Fs))^0.5*0.15
Us<-  sort(Fs + bound)
Ls<-  sort(Fs - bound)
bound <- (F1s*(1-F1s))^0.5*0.15
U1s<- sort(F1s + bound)
L1s<- sort(F1s - bound)

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

pdf("Results/figure2.pdf", height=5, width=14)
par(mfrow=c(1,3),lend="butt")
F<- cdf(c(-999,ys,max(ys)+.Machine$double.eps), c(0,Fs,1))
plot(F,xval=-1:10,xlim=c(0,10),verticals=FALSE, do.points=FALSE,col="dark blue", ylab="Probability", xlab="y", 
     ylim= c(0,1), main="Cumulative distribution function with uniform band",
     sub=" ");
for(v in 1:(length(ys)-1)) polygon(c(ys[v],ys[v+1]+0.01,ys[v+1]+0.01,ys[v]),c(Ls[v],Ls[v],min(Us[v],1),min(Us[v],1)),col="light blue",border=NA)
segments(-1,0,1,0,col="light blue",lwd=5)
segments(10,1,11,1,col="light blue",lwd=5)
lines(F, ys, verticals=FALSE, do.points=FALSE,col="dark blue", lty = 1);
box()

Q  <- left.inv(c(-999,ys,9999), c(0,Fs,max(Fs)+.Machine$double.eps))
plot(Q, xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", ylab="y", xlab="Probability", 
     ylim= c(0,10), main="Quantile function with uniform band (no support restrictions)",
     sub=" ");
segments(0,1,Us[1],1,col="light blue",lwd=5)
for(i in 2:length(ys)) segments(Ls[i-1],ys[i],min(1,Us[i]),ys[i],col="light blue", lty = 1,lwd=5) 

for(v in 1:length(ys)) polygon(c(Ls[v],Ls[v],min(Us[v],1),min(Us[v],1)),c(ys[v],ys[v+1],ys[v+1],ys[v]),col="light blue",border=NA)
lines(Q, verticals=FALSE, do.points=FALSE,col="dark blue", lty = 1);

plot(Q, xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", ylab="y", xlab="Probability", 
     ylim= c(0,10), main="Quantile function with uniform band (support restrictions)",
     sub=" ");
segments(0,1,Us[1],1,col="light blue",lwd=5)
for(i in 2:length(ys)) segments(Ls[i-1],ys[i],min(Us[i],1),ys[i],col="light blue", lty = 1,lwd=5) 
lines(Q, verticals=FALSE, do.points=FALSE,col="dark blue", lty = 1);

dev.off()

pdf("Results/figure3.pdf", height=5, width=14)
par(mfrow=c(1,3),lend="butt")

Q0  <- left.inv(c(-999,ys,9999), c(-2*.Machine$double.eps,Fs,max(Fs)+.Machine$double.eps))
Q1  <- left.inv(c(-999,ys[F1s>0],9999), c(-2*.Machine$double.eps,F1s[F1s>0],max(F1s)+.Machine$double.eps))

plot(Q1, xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", ylab="y", xlab="Probability", 
     ylim= c(0,10.5), main="Quantile functions with uniform bands (no support restrictions)",
     sub=" ");
for(i in 2:length(ys)) if(max(0,L1s[i-1])<min(1,U1s[i])) segments(max(0,L1s[i-1]),ys[i],min(1,U1s[i]),ys[i],col="light blue", lty = 1,lwd=5) 
for(v in 1:length(ys)) if(max(0,L1s[v])<min(U1s[v],1)) polygon(c(max(0,L1s[v]),max(0,L1s[v]),min(U1s[v],1),min(U1s[v],1)),c(ys[v],ys[v+1],ys[v+1],ys[v]),col="light blue",border=NA)
segments(0,1,Us[1],1,col="light green",lwd=5)
for(i in 2:length(ys)) segments(max(0,Ls[i-1]),ys[i],min(1,Us[i]),ys[i],col="light green", lty = 1,lwd=5) 
for(v in 1:length(ys)) polygon(c(max(0,Ls[v]),max(0,Ls[v]),min(Us[v],1),min(Us[v],1)),c(ys[v],ys[v+1],ys[v+1],ys[v]),col="light green",border=NA)
plot(Q1, xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", add=TRUE);
plot(Q, xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark green", add=TRUE);
legend(0, 12, c('Treatment group', 'Control group'), col = c('light blue','light green'), lwd = c(4,4,4), horiz = T, bty = 'n');
legend(0, 12, c('Treatment group', 'Control group'), col = c('dark blue','dark green'), lwd = c(1,1,1), horiz = T, bty = 'n');

quant <- unique(c(F1s,Fs))
quant <- sort(c(quant-.Machine$double.eps,quant+.Machine$double.eps))
QTE <- Q1(quant)-Q0(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-100,QTE,100),right=TRUE)
plot(QTE.func,xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Probability", ylab="y", 
     ylim= c(-1,9), main="QE function with uniform band (no support restrictions)",
     sub=" ");
quant <- unique(c(0,1,L1s,Ls,U1s,Us))
quant <- sort(c(quant-.Machine$double.eps,quant+.Machine$double.eps))
quant <- c(quant[quant>0 & quant<1],1)
uQ1  <- left.inv(ys, L1s)
lQ1  <- left.inv(ys, U1s)
uQ  <- left.inv(ys, Ls)
lQ  <- left.inv(ys, Us)
for(v in (min(ys)-max(ys)):(max(ys)-min(ys))){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=uQ1(quant[x]) - lQ(quant[x]) & v >=lQ1(quant[x]) - uQ(quant[x])
  }
  if(max(accepted)==1){
    temp_segments <- 1-accepted
    temp_segments[2:length(temp_segments)][temp_segments[1:(length(temp_segments)-1)]==1] <- 0
    temp_segments <- accepted*(cumsum(temp_segments)+accepted[1])
    for(s in 1:max(temp_segments)) segments(quant[temp_segments==s][1],v,tail(quant[temp_segments==s],1)+0.001,v,col="light blue", lty = 1,lwd=5)
  }
}
for(i in 2:length(quant)) polygon(c(quant[i-1],quant[i-1],quant[i]+0.001,quant[i]+0.001),c(uQ1(quant[i]) - lQ(quant[i]),lQ1(quant[i]) - uQ(quant[i]),lQ1(quant[i]) - uQ(quant[i]),uQ1(quant[i]) - lQ(quant[i])),col="light blue",border=NA)
lines(QTE.func, xval=c(0,quant[quant>=0 & quant<=1],1), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1, lwd=1);

quant <- unique(c(F1s,Fs))
quant <- sort(c(quant-.Machine$double.eps,quant+.Machine$double.eps))
QTE <- Q1(quant)-Q0(quant)
QTE.func <- stepfun(c(quant,max(quant)+.Machine$double.eps),c(-100,QTE,100),right=TRUE)
plot(QTE.func,xlim=c(0,1), verticals=FALSE, do.points=FALSE, col="dark blue", cex=0, xlab="Probability", ylab="y", 
     ylim= c(-1,9), main="QE function with uniform band (support restrictions)",
     sub=" ");
quant <- unique(c(0,1,L1s,Ls,U1s,Us))
quant <- sort(c(quant-.Machine$double.eps,quant+.Machine$double.eps))
quant <- c(quant[quant>0 & quant<1],1)
uQ1  <- left.inv(ys, L1s)
lQ1  <- left.inv(ys, U1s)
uQ  <- left.inv(ys, Ls)
lQ  <- left.inv(ys, Us)
for(v in (min(ys)-max(ys)):(max(ys)-min(ys))){
  accepted <- rep(0,length(quant))
  for(x in 1:length(quant)){
    accepted[x] <-  v<=uQ1(quant[x]) - lQ(quant[x]) & v >=lQ1(quant[x]) - uQ(quant[x])
  }
  if(max(accepted)==1){
    temp_segments <- 1-accepted
    temp_segments[2:length(temp_segments)][temp_segments[1:(length(temp_segments)-1)]==1] <- 0
    temp_segments <- accepted*(cumsum(temp_segments)+accepted[1])
    for(s in 1:max(temp_segments)) segments(quant[temp_segments==s][1],v,tail(quant[temp_segments==s],1),v,col="light blue", lty = 1,lwd=5)
  }
}
lines(QTE.func, xval=c(0,quant[quant>=0 & quant<=1],1), verticals=FALSE, do.points=FALSE, col="dark blue", lty = 1, lwd=1);

dev.off()
