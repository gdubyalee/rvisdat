library(shiny)
library(plyr)
library(dplyr)
library(InferCryptDrift)
library(DriftR)
source('util.R')

#Constants for sliders
LAMBDA_MIN<-0
LAMBDA_MAX<-1
TAU_MIN<-0
TAU_MAX<-2
N_MIN<-2
N_MAX<-10
#Paramaters for analytic solution plots
MOUSE_TIME<-200
TIME_INTERVAL<-1
#For now, anyway...
nBins<-8
analyticPlotTimes<-seq(TIME_INTERVAL,MOUSE_TIME,TIME_INTERVAL)

serve<-function(input,output){
  output$availableDatasets<-renderUI({
    checkboxGroupInput('datasets','Visualise datasets:',dir('data'))
  })
  #Data needs to be shared between bits of the form in a sensible way.  Idealy stuff shouldn't get called multiple times...
  #Should try and work out how to do this properly later.
  output$driftPlots<-renderPlot({
  })
  output$expectation<-renderPlot({
    renderedData<-processDataForPlots(input$datasets)
    print(renderedData[[1]])
    ggplot()+
      geom_point(data=renderedData[[1]],mapping=aes(x=time,y=n,col=experiment,size=proportion))#+
      #geom_line(data=renderedData[[2]],mapping=aes(x=time,y=p,col=experiment))
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
