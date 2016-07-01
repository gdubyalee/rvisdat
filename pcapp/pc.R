library(InferCryptDrift)
library(DriftR)
#Constants for sliders
LAMBDA_MIN<-0
LAMBDA_MAX<-1
TAU_MIN<-0
TAU_MAX<-2
N_MIN<-2
N_MAX<-10

serve=function(input,output){
  output$availableDatasets<-renderUI({
    checkboxGroupInput('datasets','Visualise datasets:',populateAvailableDatasets())
  })
  output$driftPlots<-renderPlot({
  })
  output$expectedNumberClonesPlot<-renderPlot({
  })
}

populateAvailableDatasets=function(){
  return(dir('data'))
}

#Assume that the names are of the form 'X<day number>'
getNeutralDriftParams=function(neutralDriftData){
  return(getNeutralDirftParams(fitNeutralDrift(neutralDriftData,strtoi(substring(names(neutralDriftData),2)))))
}

#Updates the cache to contain our processed data
#updateAnalysis<-function(updateAll){
#  dataFiles<-dir('data')
#  toUpdate<-ifelse(
#    updateAll,
#    dataFiles,
#    dataFiles[!(dataFiles %in% substring(dir('cache'),5))]
#  )
#  mcmc<-lapply(lapply(paste0('data/',dataFiles),readRDS),fitNeutralDriftWithTimes)
#}

pcApp=shinyUI(fluidPage(
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
    sliderInput('N','N (#stem cells in crypt)',N_MIN,N_MAX,.5*(N_MIN+N_MAX))
  )
))
