# The code in this file replicates the experiments in Example 2 of 
# Optimizing Noise for f-Differential Privacy
# via Anti-Concentration and Stochastic Dominance
# by Jordan Awan and Aishwarya Ramasethu

# Note that the tradeoff function f is defined as the type II error as a 
# function of the 1 - type I error


# Load in basic tools for working with CNDs
source("CNDtools.R")

vals = seq(-10,10,by=.01)
# visualize the log-likelihood function
plot(vals,log(dcauchy(vals-1)/dcauchy(vals)),type="l")
# illustrate the symmetry with respect to the mapping ell(x)=-ell(1-x)
lines(vals,-log(dcauchy(1-vals-1)/dcauchy(1-vals)),col="red")
abline(h=0)
abline(v=1/2)

# given threshold t, determine the x-values of the rejection region 
x_t = function(t){
  return(c(t+1-sqrt(-t^2+t+1),t+1+sqrt(-t^2+t+1))/t)
}       

# given type I error 1-alpha, determine the threshold t (RHS)
t_alpha = function(alpha){
  obj = function(t){
    x=x_t(t)
    return(((1-alpha)-(pcauchy(x[2])-pcauchy(x[1])))^2)
  }
  sol=optim(par=1,fn=obj,method="Brent",lower=0,upper=(1/2+sqrt(5)/2))
  return(sol$par)
}


# for alpha<1-c i.e. type I error >c
t_alpha_left = function(alpha){
  obj = function(t){
    x=1-x_t(t)
    x_l=x[2]
    x_r=x[1]
    return(((1-alpha)-(pcauchy(x_l)+1-pcauchy(x_r)))^2)
  }
  sol=optim(par=1,fn=obj,method="Brent",lower=0,upper=(1/2+sqrt(5)/2))
  return(sol$par)
}


# evaluate the tradeoff function for C_1 using rejection threshold t (RHS)
f_t = function(t){
  x = x_t(t)
  return(1-(pcauchy(x[2]-1)-pcauchy(x[1]-1)))
}
# evaluate the tradeoff function for C_1 using rejection threshold t (LHS)
f_t_left = function(t){
  x=1-x_t(t)
  x_l=x[2]
  x_r=x[1]
  return(pcauchy(x_r-1)-pcauchy(x_l-1))
}

# tradeoff function for C_1
f_alpha = function(alpha){
  # formula for the "c"
  c=1/2-(1/pi)*atan(1/2)
  
  if(alpha>=1-c){
  t=t_alpha(alpha)
  return(f_t(t))
  }
  else{
    t=t_alpha_left(alpha)
    return(f_t_left(t))
  }
}

#vectorized tradeoff function for C_1
f_vect = function(alpha){
  return(sapply(X=alpha,FUN=f_alpha))
}

# Verify the tradeoff function agrees with c
alpha = seq(0,1,by=.001)
plot(alpha,f_vect(alpha),type="l",xlim=c(0,1),ylim=c(0,1))
# formula for the "c"
c=1/2-(1/pi)*atan(1/2)
points(1-c,c)

val = seq(-5,5,by=.01)
plot(val,pCND(x=val,f=f_vect,c=c),type="l",ylim=c(0,1))
lines(val,pcauchy(val),col="red")

# generate right plot of Figure 3 in the paper
pdf("cauchyCND.pdf",width=7,height=7)
L=3
pos = seq(0,10,by=.01)
plot(pos,2*pCND(x=pos,f=f_vect,c=c)-1,type="l",xlab="t",ylab="P(|N|<=t)",lwd=L)
lines(pos,2*pcauchy(pos)-1,col="red",lty=2,lwd=L)
abline(h=1,lty=3,lwd=L)
abline(v=.5,lty=3,lwd=L)
axis(side=1,at=1/2,labels=.5)
legend("bottomright",c(expression("CND for C"[1]),"Cauchy(0,1)"),col=c("black","red"),lty=c(1,2),lwd=L)
dev.off()

# Cauchy has a greater P(|x|<=t) for 0<t<1/2
# maximum gap favoring cauchy is .00425
# P(|x|<=.28)

# generate 100 samples from CND and Cauchy distributions
set.seed(100)
reps = 100
cnds = rCND(n=reps,f=f_vect,c=c)
cauchys = rcauchy(n=reps)
hist(cnds)
hist(cauchys)

# values that appear in Example 2
max(abs(cnds))
max(abs(cauchys))

approxDP = function(alpha,ep,de){
  return(pmax(pmax(0,1-de-exp(ep)+exp(ep)*alpha),exp(-ep)*(alpha-de)))
}

# Left plot of Figure 3
pdf("CauchyDP.pdf",width=7,height=7)
alpha = seq(0,1,by=.001)
m=1
plot(alpha,f_vect(alpha),type="l",lwd=L,xlab="1-type I", ylab="type II")
abline(a=0,b=1,lty=3,lwd=L)
abline(h=0,lty=3,lwd=L)
abline(v=1,lty=3,lwd=L)
ep = log((4+(m+sqrt(m^2+4))^2)/(4+(m-sqrt(m^2+4))^2))
TV=(2/pi)*atan(m/2)
c=(1-TV)/2
ep2=log(1/c-1)
alpha=seq(0,1,by=.001)
lines(alpha,approxDP(alpha,ep,0),col="red",lty=2,lwd=L)
lines(alpha,approxDP(alpha,ep2,0),col="blue",lty=4,lwd=L)
legend("topleft",c(expression("f"[epsilon[U]][",0"]),expression("C"[1]),expression("f"[epsilon[L]][",0"])),
       lty=c(4,1,2),col=c("blue","black","red"),lwd=L)
dev.off()
