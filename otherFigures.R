
####################################################################
### approximate the discrete Gaussian
####################################################################
T = 1000 #### range the Discrete Gaussian on -T,T
m=0
mu=1
#mu=1.5
#A = (pnorm(mu*(3/2))-pnorm(mu*(1/2)))/(pnorm(mu*(1/2))-pnorm(mu*(-1/2)))
#s= sqrt((1/2)*(-log(A))^(-1)) ###scale parameter. usually using large values for privacy. 
s=1
dDG = function(t){
  return(exp(-(t-m)^2/(2*s^2)))#/(s*sqrt(2*pi)))
}

val = seq(-T,T)
#C = sum(dDG(val))
#C

library(elliptic)
C=theta3(0,q=exp(-1/(2*s^2)))


#k = seq(1,1000)
#U = 1+2*sum(exp(-2*pi*s^2*k^2))
#U
#C=U

c = (1-(1/C))/2#(1-1/(C*s*sqrt(2*pi)))/2
mu = -2*qnorm(c)

sum(dDG(val))/C

v= sum(val^2*dDG(val))/C
v
#csqrt(2*pi*s^2)
### seems to be EQUAL for s>=1 (!)
### for s<=.181, seems that C=1 (!?)
pDG = function(t){
  if(t<0){
    val = seq(t+1,-t-1)
    return(1/2-1/2*sum(dDG(val)/C))
  }
  if(t==0){
    return(pDG(-1)+dDG(0)/C)
  }
  if(t>0){
    return(1-pDG(-t-1))
  }
    
  #val = seq(-T,t)
  #return(sum(dDG(val))/C)
}

t = seq(-100,100)


pdf("discreteGaussian2.pdf",width=7,height=7)
L=3
plot(sapply(t,pDG),sapply((t-1),pDG),type="l",lwd=L,col="red",lty=2,xlab="1-type I", ylab="type II")
abline(a=0,b=1,lty=3,lwd=L)
abline(h=0,lty=3,lwd=L)
abline(v=1,lty=3,lwd=L)
lines(alpha,pnorm(qnorm(alpha)-mu),col="black",lwd=L)
#lines(alpha,pnorm(qnorm(alpha)-.88),col="black",lwd=L)

#lines(alpha,exp(-18)*alpha,col="purple")
lines(pnorm((t)+1/2,s=1/mu),pnorm((t)-1/2,s=1/mu),col="blue",lty=3,lwd=L)
legend("topleft",c("1-GDP","Discrete CND", "Discrete Gaussian"),lty=c(1,3,2),col=c("black","blue","red"),lwd=L)
dev.off()

t = seq(-5,5)
t2 = seq(-5,5,by=.01)
plot(t-1/2,sapply(t,pDG))
points(t+1/2,sapply(t,pDG),col="blue")
lines(t2,pnorm(t2,s=1/mu),col="red")
lines(t2-1,pnorm(t2,s=1/mu),col="red")

t = seq(-10,10)

plot(t,log(dDG(t)))
points(t,log(pnorm(t+1/2,s=1/mu)-pnorm(t-1/2,s=1/mu)),col="blue")

for(i in -1:-10){
  print(i)
print(1-2*pnorm(i+1/2,s=1/mu))
print(1-2*pDG(i))
}

t3 = -1:-10
plot(t3,-log(-log(1-2*pnorm(t3+1/2,s=1/mu))))
points(t3,-log(-log(1-2*sapply(t3,pDG))),col="blue")

### we see that rounded Normal is more concentrated than discrete Gaussian at the same TV.

pnorm(-1+1/2,s=1/mu)
pnorm(-1-1/2,s=1/mu)

pDG(-1)
pDG(-2)## higher (more private).


pdf("discreteGaussian1.pdf",width=7,height=7)
mu=2.29
L=3
plot(sapply(t,pDG),sapply((t-1),pDG),type="l",lwd=L,col="red",lty=2,xlab="1-type I", ylab="type II")
abline(a=0,b=1,lty=3,lwd=L)
abline(h=0,lty=3,lwd=L)
abline(v=1,lty=3,lwd=L)
lines(alpha,pnorm(qnorm(alpha)-mu),col="black",lwd=L)
lines(pnorm((t)+1/2,s=1/mu),pnorm((t)-1/2,s=1/mu),col="blue",lty=3,lwd=L)
legend("topleft",c("1-GDP","discrete CND", "Discrete Gaussian"),lty=c(1,3,2),col=c("black","blue","red"),lwd=L)
dev.off()





###################################################################
###   Staircase DP
ptulap <- function (t, median = 0, lambda = 0, cut=0) {
  lcut=cut/2
  rcut=cut/2
  # Normalize
  t <- t - median
  
  # Split the positive and negsative t calculations, and factor out stuff
  r <- round(t)
  g <- -log(lambda)
  l <- log(1 + lambda)
  k <- 1 - lambda
  negs <- exp((r * g) - l + log(lambda + ((t - r + (1/2)) * k)))
  poss <- 1 - exp((r * (-g)) - l + log(lambda + ((r - t + (1/2)) * k)))
  
  # Check for infinities
  negs[is.infinite(negs)] <- 0
  poss[is.infinite(poss)] <- 0
  
  # Truncate wrt the indicator on t's positivity
  is.leq0 <- t <= 0
  trunc <- (is.leq0 * negs) + ((1 - is.leq0) * poss)
  
  # Handle the cut adjustment and scaling
  cut <- lcut + rcut
  is.mid <- (lcut <= trunc) & (trunc <= (1 - rcut))
  is.rhs <- (1 - rcut) < trunc
  return (((trunc - lcut) / (1 - cut)) * is.mid + is.rhs)
}
Delta = 6
ep=1
b=exp(-ep)
de=.05
pDCND = function(t){
  return(ptulap((t+1/2)/Delta,lambda=b,cut = 2*de*b/(1-b+2*de*b)))
}
dDCND = function(t){
  return(pDCND(t)-pDCND(t-1))
}

thick=3
pdf("staircase.pdf",width=7,height=7)
val = seq(-18,19,by=1)
plot(val-1/2,dDCND(val),type="s",yaxs="i",ylim=c(0,.1),xlab="t",ylab="P(N=t)",lwd=thick,xaxt="n")
points(val-1/2,dDCND(val),type="h",lwd=thick)
abline(h=0,lwd=thick)
axis(1,at=seq(-18,18,by=3),label=seq(-18,18,by=3))
dev.off()



#############################################################
#Plot GDP and discrete CND f-DP
pdf("discreteFDP.pdf",width=7,height=7)
L=3
alpha=seq(0,1,by=.001)
plot(alpha,pnorm(qnorm(alpha)-1),type="l",lwd=L,ylab="type II",xlab="1-type I")
t = seq(-100,100)
lines(pnorm(t+1/2),pnorm(t-1/2),col="red",lty=2,lwd=L)
abline(a=0,b=1,lty=3,lwd=L)
abline(h=0,lty=3,lwd=L)
abline(v=1,lty=3,lwd=L)
legend("topleft",c("1-GDP","discrete CND for 1-GDP"),col=c("black","red"),lwd=L,lty=c(1,2))
dev.off()
###############################################################


######################################################################
###  non-integer centered noise. (uniform)
######################################################################
T = seq(-1000,1000)

delta=.4
cdf = function(t){
  return(punif(t,min=-1/(2*delta),max=1/(2*delta)))
}

quantile = function(t){
  return(qunif(t,min=-1/(2*delta),max=1/(2*delta)))
}

sum(T^2*(cdf(T+1/2)-cdf(T-1/2)))
sum((T+1/2)^2*(cdf(T+1)-cdf(T)))

L=3
alpha=seq(0,1,by=.001)
plot(c(alpha,1),c(cdf(quantile(alpha)-1),1),type="l",lwd=L,ylab="type II",xlab="1-type I",xlim=c(0,1),ylim=c(0,1))
t = seq(-100,100)
lines(cdf(t+1/2),cdf(t-1/2),col="red",lty=2,lwd=L)
lines(cdf(t),cdf(t-1),col="blue",lty=3,lwd=L)
abline(a=0,b=1,lty=3,lwd=L)
abline(h=0,lty=3,lwd=L)
abline(v=1,lty=3,lwd=L)
legend("topleft",c("1-GDP","discrete CND for 1-GDP","floor(N(0,1))"),col=c("black","red","blue"),lwd=L,lty=c(1,2,3))




############################################################
###   Tulap versus Laplace 
############################################################

epsilon = seq(0,5,by=.01)

pdf("tulapLaplacehalf.pdf",width=7,height=7)
plot(epsilon,(exp(epsilon)-1)/(exp(epsilon)+1),type="l",ylim=c(0,1),lwd=3,ylab="P(|.|<=1/2)")
lines(epsilon,1-exp(-epsilon/2),col="red",lty=2,lwd=3)
legend("topleft",c("Tulap","Laplace"),col=c("black","red"),lty=c(1,2),lwd=3)
dev.off()


pdf("tulapLaplacefourth.pdf",width=7,height=7)
plot(epsilon,(exp(epsilon)-1)/(2*(exp(epsilon)+1)),type="l",ylim=c(0,1),lwd=3,ylab="P(|.|<=1/4)")
lines(epsilon, 1-exp(-epsilon/4),col="red",lty=2,lwd=3)
legend("topleft",c("Tulap","Laplace"),col=c("black","red"),lty=c(1,2),lwd=3)
dev.off()


###########################################################
###   Anti-concentration lemma plot
library("sn")
val = seq(-2,4,by=.01)
xi=-.75
t=2
th = 3

pdf("antiLemma.pdf",width=8,height=4)
plot(0,0,type="l",xaxt="n",xlab="",ylab="",yaxt="n",lwd=th,xlim=c(-2,4),ylim=c(0,.7),xaxs="i", yaxs="i")
s=1
d = 25
axis(1,at=c(-1,0,1,2,3),labels=c("a-t/2","a","a+t/2","a+t","a+3t/2"))
v1 = seq(1,4,by=.01)
a1 = dsn(v1,xi=xi,omega=1,alpha=3,tau=0)
polygon(x=c(v1,1),y=c(a1,0),col="red",density=d,angle=90)
v2 = seq(3,6,by=.01)
polygon(x=c(v2,3),y=c(a1,0),col="red",density=d,angle=90)
v3 = seq(-3,-1,by=.01)
a3 = dsn(v3,xi=xi,omega=1,alpha=3,tau=0)
polygon(x=c(v3,-1),y=c(a3,0),col="blue",density=d,angle=0)
v4 = seq(-1,1,by=.01)
polygon(x=c(v4,1),y=c(a3,0),col="blue",density=d,angle=0)
for(i in -2:2){
  abline(v=(s+i),lty=3)
}
lines(val,dsn(val,xi=xi,omega=1,alpha=3,tau=0),lwd=th)
lines(val,dsn(val,xi=xi+t,omega=1,alpha=3,tau=0),col="black",lty=2,lwd=th)
legend("topleft",c("type II error","type I error"),fill=c("blue","red"),density=d,angle=c(0,90),bg="white")
legend("topright",c("density of N","density of N+t"),lty=c(1,2),lwd=th,bg="white")
dev.off()

