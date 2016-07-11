TIME_INTERVAL<-1


#Assume that the names are of the form 'X<day number>'
getNeutralDriftParams<-function(neutralDriftData){
  if(is.null(neutralDriftData)){print('Got a null...');return(NULL)}
  return(getNeutralDirftParams(fitNeutralDrift(neutralDriftData,strtoi(substring(names(neutralDriftData),2)))))
}

handleUpload<-function(uploadedDataset){
  uploadedObj<-read.csv(uploadedDataset$datapath,sep=' ')
  #Assume filename is [NAME].data for now
  saveName<-paste0(substr(uploadedDataset$name,1,nchar(uploadedDataset$name)-4),'rds')
  saveRDS(uploadedObj,paste0('data/',saveName))
  saveRDS(getNeutralDriftParams(uploadedObj),paste0('cache/mcmc_',saveName))
}

processDataForPlots<-function(selectedDatasets,mouseLife,N,lambda,tau){
  analyticPlotTimes<-seq(TIME_INTERVAL,mouseLife,TIME_INTERVAL)
  #This feels rather hacky...
  if(length(selectedDatasets)){
    rawIn<-readRDS(paste0('data/',selectedDatasets[1]))
    rawIn<-rawIn/colSums(rawIn)
    rawData<-cbind(
      rawIn,
      experiment=selectedDatasets[1],
      proportion=1:nBins
    )

    analyticIn<-readRDS(paste0('cache/mcmc_',selectedDatasets[1]))
    analyticIn<-data.frame(Crypt_drift_c(
      analyticIn$lambda,
      analyticIn$lambda,
      analyticIn$N,
      analyticPlotTimes,
      1,
      nBins
    ))
    names(analyticIn)<-analyticPlotTimes
    analyticData<-cbind(
      analyticIn,
      experiment=selectedDatasets[1],
      proportion=1:nBins
      #time<-analyticPlotTimes
    )
  }
  if(length(selectedDatasets)>1){
    for(i in 2:length(selectedDatasets)){
      #Need rbind.fill - should perhaps actually gather data earlier...
      rawIn<-readRDS(paste0('data/',selectedDatasets[i]))
      rawIn<-rawIn/colSums(rawIn)
      rawData<-rbind.fill(
        rawData,
        cbind(
          rawIn,
          experiment=selectedDatasets[i],
          proportion=1:nBins
        )
      )
      #Analytic stuff
      analyticIn<-readRDS(paste0('cache/mcmc_',selectedDatasets[i]))
      analyticIn<-data.frame(Crypt_drift_c(
        analyticIn$lambda,
        analyticIn$lambda,
        analyticIn$N,
        #Probabilities fixed 
        pmax(analyticPlotTimes-analyticIn$tau,0),
        1,
        nBins
      ))
      names(analyticIn)<-analyticPlotTimes
      analyticIn<-cbind(
        analyticIn,
        experiment=selectedDatasets[i],
        proportion=1:nBins
        #time<-analyticPlotTimes
      )
      analyticData<-rbind(
        analyticData,
        analyticIn
      )
    }
  }

  if(length(selectedDatasets)){
    #Now make the frames into collections of time points rather than having huge numbers of colu,mns in general
    rawData<-gather(rawData,'time','n',1:(ncol(rawData)-2))
    analyticData<-gather(analyticData,'time','p',1:(ncol(analyticData)-2))
    #coerce time points in raw data into numers
    rawData[['time']]<-strtoi(substring(rawData[['time']],2))
    analyticData[['time']]<-as.numeric(analyticData[['time']])
    return(list(rawData,analyticData))
  }
  list(NULL,NULL)
}
