library(shiny)
library(plyr)
library(dplyr)
library(InferCryptDrift)
library(DriftR)

serve<-function(input,output){
  output$controlSets<-renderUI({
    radioButtons('wtToCmp','Select WT set to compare against:',dir('wt'))
  })
  observe({
    input$controlUpload
    if(length(input$controlUpload)){
      wtIn<-read.csv(input$controlUpload$datapath)
      colnames(wtIn)<-paste0('D',substr(colnames(wtIn),2,length(colnames(wtIn))))
      if(nrow(wtIn)!=8)return('Unexpected number of bins.')
    }
  })
}

cmpApp<-shinyUI(fluidPage(
  titlePanel('Neutral/Biased Drift Comparison'),
  sidebarLayout(
    sidebarPanel(
      fileInput('controlUpload','New WT set'),
      uiOutput('controlSets'),
      fileInput('biasedUpload','New set to analyse bias')
    ),
    mainPanel(
    )
  )
))
