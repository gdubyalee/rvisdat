library(dplyr)
library(plyr)
library(InferCryptDrift)
library(DriftR)
#Constants for sliders
LAMBDA_MIN<-0
LAMBDA_MAX<-1
TAU_MIN<-0
TAU_MAX<-2
N_MIN<-2
N_MAX<-10
MOUSE_TIME<-200
TIME_INTERVAL<-1
#For now, anyway...
nBins<-8

populateAvailableDatasets<-function(){
  return(dir('data'))
}

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
  print('Now to work out some params using the mcmc algorithm')
  print(uploadedObj)
  saveRDS(getNeutralDriftParams(uploadedObj),paste0('cache/mcmc_',saveName))
}

processDataForPlots<-function(selectedDatasets){
  if(is.null(selectedDatasets))return(NULL)
  #Tend to get a nicer number of elements like this
  timesForAnalyticFit<-seq(0,MOUSE_TIME-TIME_INTERVAL,TIME_INTERVAL)
  #nBins<-16
  mcmcParams<-list()
  rawData<-list()
  for(i in 1:length(selectedDatasets)){
    mcmcParams[[i]]<-readRDS(paste0('cache/mcmc_',selectedDatasets[i]))
    rawData[[i]]<-readRDS(paste0('data/',selectedDatasets[i]))
    #nBins<-min(nBins,ncol(rawData[[i]]))
  }

  analyticList<-lapply(mcmcParams,function(obj)Crypt_drift_c(obj$lambda,obj$lambda,obj$N,timesForAnalyticFit,1,nBins))
  #Add in a row for times
  analyticList<-lapply(1:length(selectedDatasets),function(i){
    rownames(analyticList[[i]])<-as.character(1:nBins)
    ret<-rbind.data.frame(
      times=timesForAnalyticFit,
      analyticList[[i]]
    )
    ret<-rbind.data.frame(ret,experiment=selectedDatasets[i])
    print(ret)
    #names(ret)<-c('times',as.character(1:nBins))
    ret
  })

  #names(analyticList)<-selectedDatasets
  print(analyticList)
  #frameToPlotAnalytic<-ldply(analyticList)
  #frameToPlotAnalytic<-do.call('rbind',analyticList)
  #print(frameToPlotAnalytic)
  frameToPlotAnalytic<-analyticList[[1]]
  if(length(selectedDatasets)>1){
    for(i in 2:length(selectedDatasets)){
      frameToPlotAnalytic<-cbind.data.frame(frameToPlotAnalytic,analyticList[[i]])
    }
  }
  #print(frameToPlotAnalytic)
  frameToPlotRaw<-NULL
  
  ##Now squish so all have the right number of bins and put in a list
  #frameToPlotRaw<-NULL
  #for(i in 1:length(selectedDatasets)){
  #  dat<-NULL
  #  #Don't care about this yet
  #  if(FALSE&&ncol(rawData[[i]])!=nBins){
  #    tryCatch(someVariableWhichDoesNotExist,error=function(e)print('Haven\'t implemented comparison of multiple bin sized data yet!'))
  #    #squish<-ncol(rawData[i])/nBins
  #    #for(j in 1:nBins){
  #    #  #Need to check indices are arranged the right way here
  #    #  dat[j]<-colSums(rawData[j:(j+squish-1),])
  #    #}
  #  }else{
  #    dat<-rawData[[i]]
  #  }
  #  dat[['experiment']]<-selectedDatasets[i]
  #  dat[['times']]<-strtoi(substring(names(rawData[i]),2))
  #  frameToPlotRaw<-rbind(frameToPlotRaw,dat)
  #}

  c(frameToPlotRaw,frameToPlotAnalytic)
}

serve<-function(input,output){

  

  output$availableDatasets<-renderUI({
    checkboxGroupInput('datasets','Visualise datasets:',populateAvailableDatasets())
  })
  #Data needs to be shared between bits of the form in a sensible way.  Idealy stuff shouldn't get called multiple times...
  #Should try and work out how to do this properly later.
  output$driftPlots<-renderPlot({
  })
  output$expectation<-renderPlot({
    renderedData<-processDataForPlots(input$datasets)
    print(renderedData)
    plot(c(1,1,2,3,5,8,13,21))
  })
  observe({
    input$newDataset
    if(length(input$newDataset)){
      handleUpload(input$newDataset)
      output$availableDatasets<-renderUI({
        checkboxGroupInput('datasets','Visualise datasets:',populateAvailableDatasets())
      })
    }
  })
}

pcApp<-shinyUI(fluidPage(
  #Display available datasets and expected number of clones in first row
  sidebarLayout(
    #List available data in side panel
    sidebarPanel(uiOutput('availableDatasets')),
    mainPanel(
      plotOutput('expectation')
    )
  ),
  fluidRow(
    plotOutput('driftPlots')
  ),
  flowLayout(
    fileInput('newDataset','Upload a new dataset'),
    sliderInput('lambda',HTML('&lambda; (replacement rate)'),LAMBDA_MIN,LAMBDA_MAX,.5*(LAMBDA_MIN+LAMBDA_MAX)),
    sliderInput('tau',HTML('&tau; (initial delay)'),TAU_MIN,TAU_MAX,.5*(TAU_MIN+TAU_MAX),step=.01),
    sliderInput('N','N (#active stem cells/crypt)',N_MIN,N_MAX,.5*(N_MIN+N_MAX))
  )
))
