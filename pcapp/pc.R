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

populateAvailableDatasets<-function(){
  return(dir('data'))
}

#Assume that the names are of the form 'X<day number>'
getNeutralDriftParams<-function(neutralDriftData){
  return(getNeutralDirftParams(fitNeutralDrift(neutralDriftData,strtoi(substring(names(neutralDriftData),2)))))
}

handleUpload<-function(uploadedDataset){
  uploadedObj<-read.csv(uploadedDataset$datapath,sep=' ')
  #Assume filename is [NAME].data for now
  saveName<-paste0(substr(uploadedDataset$name,1,nchar(uploadedDataset$name)-4),'rds')
  saveRDS(uploadedObj,paste0('data/',saveName))
  saveRDS(getNeutralDriftParams(uploadedObj),paste0('cache/mcmc_',saveName))
}

expectationPlot<-function(selectedDatasets){
  if(is.null(selectedDatasets))return(NULL)
  mcmcParams<-lapply(paste0('cache/mcmc_',selectedDatasets),readRDS)
  rawData<-lapply(paste0('data/',selectedDatasets),readRDS)
  timesForAnalyticFit<-seq(0,MOUSE_TIME,.1)
  analyticPlot<-lapply(mcmcParams,function(obj)Crypt_drift_c(obj$lambda,obj$lambda,obj$N,timesForAnalyticFit,1))
  analyticPlotExpectation<-lapply(analyticPlot,function(data)t(data)%*%1:nrow(data))
  names(analyticPlotExpectation)<-selectedDatasets
  names(timesForAnalyticFit)<-'times'
  analyticPlotExpectationTimes<-lapply(analyticPlotExpectation,function(data)join(data,timesForAnalyticFit))
  #toPlot<-do.call(rbind.data.frame,analyticPlotExpectation)
  #print(toPlot)
  #ggplot(analyticPlotExpectation,aes(
  frameToPlot<-Reduce(function(...)merge(...,all=T),analyticPlotExpectationTimes)
  #!!!!!!!!! :(
  #ggplot(frameToPlot,aes(x=times,y=selectedDatasets)+geom_line(mapping=aes(
}

serve<-function(input,output){
  output$availableDatasets<-renderUI({
    checkboxGroupInput('datasets','Visualise datasets:',populateAvailableDatasets())
  })
  output$driftPlots<-renderPlot({
  })
  output$expectation<-renderPlot({
    plot(expectationPlot(input$datasets))
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
