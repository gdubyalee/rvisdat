#!/usr/bin/env Rscript

#This isn't correct as such, uses equations that assume 1d drift.

#

library(dplyr)
library(tidyr)
library(rmutil)

#Evolution matrix describing equations
dayEvolution<-function(lambda,Ns){
  #Generate matrix
  A<-matrix(0,nrow=Ns+1,ncol=Ns+1)
  for(i in 1:Ns){
    A[i,i]=-2
    A[i,i+1]=1
    A[i+1,i]=1
  }
  A[,Ns+1]=A[,1]=0
  return(mexp(A,lambda))
}

solveIVPDays<-function(Ns,init,lambda,days){
  A<-dayEvolution(lambda,Ns)
  for(i in 1:days){
    init<-A%*%init
  }
  return(init)
}

solveIVP<-function(init,lambda,Ns,t){
  #Generate matrix
  A<-matrix(0,nrow=Ns+1,ncol=Ns+1)
  for(i in 1:Ns){
    A[i,i]=-2
    A[i,i+1]=1
    A[i+1,i]=1
  }
  A[,Ns+1]=A[,1]=0
  return(mexp(A,lambda*t)%*%init)
}

fracExpressing<-.1
lambda<-.2
Ns<-7
ic<-dbinom(0:Ns,Ns,fracExpressing)
print(ic)
sln<-solveIVP(ic,lambda,Ns,30)
print(sln)

probs<-seq(0,1,0.01)
fracPar<-c()
fracFull<-c()
for(pExpressing in probs){
  ic<-dbinom(0:Ns,Ns,pExpressing)
  sln<-solveIVP(ic,lambda,Ns,30)
  fracPar[length(fracPar)+1]<-sum(sln[2:Ns])
  fracFull[length(fracFull)+1]<-sln[Ns+1]
}

library(ggplot2)
plt<-ggplot(data.frame(probs,fracPar,fracFull))+
  geom_line(aes(x=probs,y=fracFull))+
  geom_line(aes(x=probs,y=fracPar))
print(plt)

