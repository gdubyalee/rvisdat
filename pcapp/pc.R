library(shiny)
library(plyr)
library(dplyr)
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
analyticPlotTimes<-seq(TIME_INTERVAL,MOUSE_TIME,TIME_INTERVAL)

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
      rawData<-rbind(
        rawData,
        cbind(
          readRDS(paste0('data/',selectedDatasets[i])),
          experiment=selectedDatasets[i],
          proportion=1:nBins
        )
      )
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
