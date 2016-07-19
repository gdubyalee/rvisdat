library(DriftR)

lsfParamsNeutral<-function(dataset){
  mins<-list()
  for(N in 3:20){
    lsCostDataset<-function(lambdatau)lsCost(N,lambdatau[1],lambdatau[2],dataset)
    #min<-optim(c(8,0.1,1),lsCostDataset)
    mins[[N]]<-nlm(lsCostDataset,c(0.1,1),iterlim=1000)
  }
  mins
}

lsCost<-function(N,lambda,tau,dataset){
  #if(tau<0){return(1000000)}
  times<-strtoi(substring(names(dataset),2))
  #Need to factor to the right number of bins - last argument for Crypt_drift_c
  analyticSln<-Crypt_drift_c(lambda,lambda,N,times-tau,1,nrow(dataset))
  #print(analyticSln)
  #print(dataset)
  sum((analyticSln-dataset)^2)
}
