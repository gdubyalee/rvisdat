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

processDataForPlots<-function(selectedDatasets){
  #This feels rather hacky...
  if(length(selectedDatasets)){
    rawData<-cbind(
      readRDS(paste0('data/',selectedDatasets[1])),
      experiment=selectedDatasets[1],
      proportion=1:nBins
    )

    d<-readRDS(paste0('cache/mcmc_',selectedDatasets[1]))
    d<-data.frame(Crypt_drift_c(
      d$lambda,
      d$lambda,
      d$N,
      analyticPlotTimes,
      1,
      nBins
    ))
    names(d)<-analyticPlotTimes
    analyticData<-cbind(
      d,
      experiment=selectedDatasets[1],
      proportion=1:nBins
      #time=analyticPlotTimes
    )
  }
  if(length(selectedDatasets)>1){
    for(i in 2:length(selectedDatasets)){
      #Need rbind.fill - should perhaps actually gather data earlier...
      rawData<-rbind.fill(
        rawData,
        cbind(
          readRDS(paste0('data/',selectedDatasets[i])),
          experiment=selectedDatasets[i],
          proportion=1:nBins
        )
      )
      #Analytic stuff
      d<-readRDS(paste0('cache/mcmc_',selectedDatasets[i]))
      d<-data.frame(Crypt_drift_c(
        d$lambda,
        d$lambda,
        d$N,
        analyticPlotTimes,
        1,
        nBins
      ))
      names(d)<-analyticPlotTimes
      d<-cbind(
        d,
        experiment=selectedDatasets[i],
        proportion=1:nBins
        #time<-analyticPlotTimes
      )
      analyticData<-rbind(
        analyticData,
        d
      )
    }
  }

  if(length(selectedDatasets)){
    #Now make the frames into collections of time points rather than having huge numbers of colu,mns in general
    rawData<-gather(rawData,'time','n',1:(ncol(rawData)-2))
    analyticData<-gather(analyticData,'time','p',1:(ncol(analyticData)-2))
    #coerce time points in raw data into numers
    rawData[['time']]<-strtoi(substring(rawData[['time']],2))
    return(list(rawData,analyticData))
  }
  list(NULL,NULL)
}
