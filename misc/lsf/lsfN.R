library(DriftR)

lsfParamsNeutral<-function(dataset){
  lsCostDataset<-function(Nlambdatau)lsCost(Nlambdatau[1],Nlambdatau[2],Nlambdatau[3],dataset)
  mins<-list()
  for(N in 3:20)mins[[N-2]]<-optim(c(N,0.1,1),lsCostDataset)
  #mins<-nlm(lsCostDataset,c(8,1,0.1),iterlim=1000)
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
