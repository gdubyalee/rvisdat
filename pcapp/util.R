TIME_INTERVAL<-1


availableDatasetList<-function(){
  ret<-dir('data')
  ret<-ret[ret!='user_defined']
  l<-list()
  if(length(ret)){
    for(i in 1:length(ret)){
      r<-readRDS(paste0('cache/mcmc_',ret[i]))
      l[i]<-paste0(
        substr(ret[i],1,nchar(ret[i])-4),
        ' (N=',
        r$N,
        ', lambda=',
        format(r$lambda,digits=3),
        ', tau=',
        format(r$tau,digits=3),
        ')'
      )
    }
  }
  names(ret)=l
  c(ret,'user_defined')
}

#Assume that the names are of the form 'X<day number>'
getNeutralDriftParams<-function(neutralDriftData,name){
  if(is.null(neutralDriftData)){print('Got a null...');return(NULL)}
  mcmc<-fitNeutralDrift(neutralDriftData,strtoi(substring(names(neutralDriftData),2)))
  #HACK
  saveRDS(mcmc,paste0('raw/raw_',name))
  return(getNeutralDirftParams(mcmc))
}

handleUpload<-function(uploadedDataset){
  uploadedObj<-read.csv(uploadedDataset$datapath,sep=' ')
  #Assume filename is [NAME].data for now
  saveName<-paste0(substr(uploadedDataset$name,1,nchar(uploadedDataset$name)-4),'rds')
  if(saveName=='user_defined')return('Please choose another name!')
  saveRDS(uploadedObj,paste0('data/',saveName))
  saveRDS(getNeutralDriftParams(uploadedObj,saveName),paste0('cache/mcmc_',saveName))
}

processDataForPlots<-function(selectedDatasets,mouseLife,N,lambda,tau,errorBars=F){
  analyticPlotTimes<-seq(TIME_INTERVAL,mouseLife,TIME_INTERVAL)
  rawData<-NULL
  analyticData<-NULL
  errData<-NULL
  for(i in 1:length(selectedDatasets)){
    #Need rbind.fill - should perhaps actually gather data earlier...
    rawIn<-readRDS(paste0('data/',selectedDatasets[i]))
    analyticIn<-readRDS(paste0('cache/mcmc_',selectedDatasets[i]))
    errIn<-NULL

    #Get error bars from Ed's fn
    if(errorBars&&(selectedDatasets[i]!='user_defined')){
      errIn<-format_exp_data(
        format_tbl2list(
          rawIn,
          strtoi(substring(names(rawIn),2))
        ),
        clone_fractions=nrow(rawIn)
      )

      #?!?!?!?!?!?!?!?!?!!
      prop<-if(typeof(errIn$Nths[1])=='character') errIn$Nths else errIn$Nths[[1]]
      errIn<-data.frame(
        proportion=prop,
        lo=if(typeof(errIn$low_lim)=='character')errIn$low_lim else errIn$low_lim[[1]],
        hi=if(typeof(errIn$hi_lim)=='character')errIn$hi_lim else errIn$hi_lim[[1]],
        time=if(typeof(errIn$Day)=='character')errIn$Day else errIn$Day[[1]],
        experiment=selectedDatasets[i]
      )
    }

    if(selectedDatasets[i]!='user_defined'){
      rawIn<-cbind(
        #Normalise to proportions at each time rather than raw counts
        rawIn/rep(colSums(rawIn),each=nrow(rawIn)),
        experiment=selectedDatasets[i],
        proportion=1:nBins
      )

      #if(errorBars){
      #  rawIn<-full_join(errIn,rawIn,by=c(experiment,time,frac))
      #}

      errData<-if(i==1){
        errIn
      }else{
        rbind.fill(
          errData,
          errIn
        )
      }
      print(errData)

      rawData<-if(i==1){
        rawIn
      }else{
        rbind.fill(
          rawData,
          rawIn
        )
      }
    }
    #Analytic stuff
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
    )

    analyticData<-if(i==1){
      analyticIn
    }else{
      rbind(
        analyticData,
        analyticIn
      )
    }
  }


  if(length(selectedDatasets)){
    #Now make the frames into collections of time points rather than having huge numbers of colu,mns in general
    analyticData<-gather(analyticData,'time','p',1:(ncol(analyticData)-2))
    analyticData[['time']]<-as.numeric(analyticData[['time']])
    if(exists('rawData')&&length(rawData)){
      rawData<-gather(rawData,'time','n',1:(ncol(rawData)-2))
      #coerce time points in raw data into numers
      rawData[['time']]<-strtoi(substring(rawData[['time']],2))
      if(errorBars){
        errData$time<-as.numeric(errData$time)
        rawData$time<-as.numeric(rawData$time)
        errData$proportion<-as.numeric(errData$proportion)
        rawData$proportion<-as.numeric(rawData$proportion)
        print(rawData)
        rawData<-full_join(rawData,errData,by=c('proportion','experiment','time'))
      }
      return(list(rawData,analyticData))
    }else return(list(NULL,analyticData))
  }
  list(NULL,NULL)
}


genExpPlots<-function(input){
  saveRDS(list(lambda=input$lambda,tau=input$tau,N=input$N),'cache/mcmc_user_defined')
  if(length(input$datasets)){
    renderedData<-processDataForPlots(input$datasets,input$T)
    renderedData[[2]]<-ddply(renderedData[[2]],.(experiment,time),summarize,expectation=sum(proportion*p)/nBins)
    if(!is.null(renderedData[[1]])){
      renderedData[[1]]<-ddply(renderedData[[1]],.(experiment,time),summarize,expectation=sum(proportion*n)/nBins)
      ggplot()+
        geom_point(data=renderedData[[1]],mapping=aes(x=time,y=expectation,col=experiment))+
        geom_line(data=renderedData[[2]],mapping=aes(x=time,y=expectation,col=experiment))+
        ggtitle('Expected proportion of crypt occupied by clones where clones found')
    }else{
      ggplot()+
        geom_line(data=renderedData[[2]],mapping=aes(x=time,y=expectation,col=experiment))+
        ggtitle('Expected proportion of crypt occupied by clones where clones found')
    }
  }
}

genClonPlots<-function(input){
  saveRDS(list(lambda=input$lambda,tau=input$tau,N=input$N),'cache/mcmc_user_defined')
  if(length(input$datasets)){
    renderedData<-processDataForPlots(input$datasets,input$T,input$N,input$lambda,input$tau,T)
    if(!is.null(renderedData[[1]])){
      #print(head(renderedData[[1]])
      ggplot()+
        geom_point(data=renderedData[[1]],mapping=aes(x=time,y=n,col=experiment))+
        geom_line(data=renderedData[[2]],mapping=aes(x=time,y=p,group=experiment,col=experiment))+
        #geom_errorbar(data=renderedData[[1]],mapping=aes(x=time,ymax=hi,ymin=lo))+
        ylab('p')+
        facet_grid(~proportion)+
        ggtitle('Clonal drift profiles')
    }else{
      ggplot()+
        geom_line(data=renderedData[[2]],mapping=aes(x=time,y=p,group=experiment,col=experiment))+
        ylab('p')+
        facet_grid(~proportion)+
        ggtitle('Clonal drift profiles')
    }
  }
}

genPosteriorPlots<-function(input){
  p<-list()
  ds<-input$datasets[input$datasets!='user_defined']
  for(i in 1:length(ds)){
    p[[i]]<-plotPosterior_Neutral(readRDS(paste0('raw/raw_',ds[i])))
  }
  if(i>1)p[[1]]<-do.call(grid.arrange,c(p,ncol=1))
  p[[1]]
}
