library(ggplot2)
library(stringr)
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
  output$posteriorPdf<-downloadHandler(
    filename='posterioranalysis.pdf',
    content=function(file){
      ggsave(
        file,
        plot=genPosteriorPlots(input),
        device='pdf'
      )
    }
  )
  output$genPdfCld<-downloadHandler(
    filename='clonaldriftprofiles.pdf',
    content=function(file){
      ggsave(
        file,
        width=12,
        height=6,
        plot=genClonPlots(input),
        device='pdf'
      )
    }
  )
  output$genPdfExp<-downloadHandler(
    filename='clonaldriftexpectations.pdf',
    content=function(file){
      ggsave(
        file,
        plot=genExpPlots(input),
        device='pdf'
      )
    }
  )
  output$availableDatasets<-renderUI({
    checkboxGroupInput('datasets','Visualise datasets:',availableDatasetList())
  })
  #Data needs to be shared between bits of the form in a sensible way.  Idealy stuff shouldn't get called multiple times...
  #Should try and work out how to do this properly later.
  output$driftPlots<-renderPlot({
    genClonPlots(input)
  })
  output$expectation<-renderPlot({
    genExpPlots(input)
  })
  observe({
    input$newDataset
    if(length(input$newDataset)){
      handleUpload(input$newDataset)
      output$availableDatasets<-renderUI({
        checkboxGroupInput('datasets','Available datasets',availableDatasetList())
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
      HTML(
        '<h4>Add new datasets using the file upload input below.
        the \'user defined\' line can be adjusted, when selected, using the sliders below.
        All curves shown are of the same analytic solution of the uynderlying dynamical equations, but with different choices of parameters (N,&lambda;,&tau;) to fit observed datasets.
        Note all plots make assumptions of neutral drift (mutants have no competitive bias) for now.</h4>'
      ),
      checkboxInput('errbars','Display error bars'),
      uiOutput('availableDatasets')
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
    sliderInput('T','Maximum time in plot',MOUSE_TIME_MIN,MOUSE_TIME_MAX,.5*(MOUSE_TIME_MIN+MOUSE_TIME_MAX)),
    downloadButton('genPdfExp','Download pdf of average clone size profile'),
    downloadButton('genPdfCld','Download pdf of clonal drift profiles'),
    downloadButton('posteriorPdf','Download pdf of fit')
  )
))
