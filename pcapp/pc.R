library(shiny)
library(plyr)
library(dplyr)
library(InferCryptDrift)
library(DriftR)
source('util.R')

#Constants for sliders
LAMBDA_MIN<-0
LAMBDA_MAX<-2
TAU_MIN<-0
TAU_MAX<-10
N_MIN<-2
N_MAX<-20
#Paramaters for analytic solution plots
MOUSE_TIME_MIN<-20
MOUSE_TIME_MAX<-200
#For now, anyway...
nBins<-8
renderedData<-NULL

serve<-function(input,output){
  output$availableDatasets<-renderUI({
    checkboxGroupInput('datasets','Visualise datasets:',dir('data'))
  })
  #Data needs to be shared between bits of the form in a sensible way.  Idealy stuff shouldn't get called multiple times...
  #Should try and work out how to do this properly later.
  output$driftPlots<-renderPlot({
    saveRDS(list(lambda=input$lambda,tau=input$tau,N=input$N),'cache/mcmc_user_defined')
    if(length(input$datasets)){
      renderedData<-processDataForPlots(input$datasets,input$T,input$N,input$lambda,input$tau)
      if(!is.null(renderedData[[1]])){
        ggplot()+
          geom_point(data=renderedData[[1]],mapping=aes(x=time,y=n,col=experiment))+
          geom_line(data=renderedData[[2]],mapping=aes(x=time,y=p,group=experiment,col=experiment))+
          facet_grid(~proportion)+
          ggtitle('Clonal drift profiles')
      }else{
        ggplot()+
          geom_line(data=renderedData[[2]],mapping=aes(x=time,y=p,group=experiment,col=experiment))+
          facet_grid(~proportion)+
          ggtitle('Clonal drift profiles')
      }
    }
    #ggplot()+
    #  geom_point(data=renderedData[[1]],mapping=aes(x=time,y=n,col=experiment,size=proportion))+
    #  geom_line(data=renderedData[[2]],mapping=aes(x=time,y=p,col=experiment,group=interaction(proportion,experiment)))
    #  #geom_line(data=renderedData[[2]],mapping=aes(x=time,y=p,col=experiment))
  })
  output$expectation<-renderPlot({

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
  })
  observe({
    input$newDataset
    if(length(input$newDataset)){
      handleUpload(input$newDataset)
      output$availableDatasets<-renderUI({
        checkboxGroupInput('datasets','Available datasets',dir('data'))
      })
    }
  })
  
}

pcApp<-shinyUI(fluidPage(
  titlePanel('Clonal pulse-chase visualisation'),
  #Display available datasets and expected number of clones in first row
  sidebarLayout(
    #List available data in side panel
    sidebarPanel(
      HTML('<h4>Add new datasets using the file upload input below.
        the \'user\' line can be adjusted, when selected, using the sliders below.
        All curves shown are of the same analytic solution of the uynderlying dynamical equations, but with different choices of parameters (N,&lambda;,&tau;) to fit observed datasets.
        Note then all plots make assumptions of neutral drift for now.</h4>'), uiOutput('availableDatasets')
    ),
    mainPanel(
      plotOutput('expectation')
    )
  ),
  fluidRow(
    plotOutput('driftPlots')
  ),
  flowLayout(
    fileInput('newDataset','Upload a new dataset'),
    sliderInput('lambda',HTML('&lambda; (replacement rate)'),LAMBDA_MIN,LAMBDA_MAX,.5*(LAMBDA_MIN+LAMBDA_MAX),step=.01),
    sliderInput('tau',HTML('&tau; (initial delay)'),TAU_MIN,TAU_MAX,.5*(TAU_MIN+TAU_MAX),step=.01),
    sliderInput('N','N (#active stem cells/crypt)',N_MIN,N_MAX,.5*(N_MIN+N_MAX)),
    sliderInput('T','Maximum time in plot',MOUSE_TIME_MIN,MOUSE_TIME_MAX,.5*(MOUSE_TIME_MIN+MOUSE_TIME_MAX))
  )
))
