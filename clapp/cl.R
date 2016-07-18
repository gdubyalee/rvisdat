library(shiny)
library(plyr)
library(dplyr)
library(InferCryptDrift)
library(DriftR)
MU_MIN<-0
MU_MAX<-10
LAMBDA_MIN<-0
LAMBDA_MAX<-2
PR_MIN<-0
PR_MAX<-1
NS_MIN<-2
NS_MAX<-20
MAXTIME_MIN<-30
MAXTIME_MAX<-200
NUMRUNS_MIN<-200
NUMRUNS_MAX<-20000
NUM_TIME_POINTS<-1000

handleUpload<-function(uploadedDataset){
  uploadedObj<-read.csv(uploadedDataset$datapath)#,sep=' ')
  #Assume filename is [NAME].rds for now
  saveName<-paste0(substr(uploadedDataset$name,1,nchar(uploadedDataset$name)-3),'rds')
  saveRDS(uploadedObj,paste0('data/',saveName))
}

#Constants for sliders

serve<-function(input,output){
  output$availableDatasets<-renderUI({
    radioButtons('dataset','Visualise datasets:',dir('data'))#c('randomly generated simulation',dir('data')))
  })
  output$clPlots<-renderPlot({
    if(length(input$dataset)){
      #Set data to be visualised
      d<-if(input$dataset=='randomly generated simulation'){
        as.data.frame(sim_contLabelling(
          input$mu,
          input$lambda,
          input$Ns,
          seq(0,input$maxTime,length.out=NUM_TIME_POINTS),
          input$numRuns,
          input$Pr
        ))
      }else{
        as.data.frame(readRDS(paste0('data/',input$dataset)))
      }
       ggplot()+
         geom_point(data=d,aes(x=Age,y=Clone.number,color=Location))+
         geom_smooth(data=d,method=input$smoothMethod,aes(x=Age,y=Clone.number,color=Location))+
         facet_wrap(~Type.Clones)         
    }
  })
  observe({
    input$newDataset
    if(length(input$newDataset)){
      handleUpload(input$newDataset)
      output$availableDatasets<-renderUI({
        radioButtons('dataset','Visualise datasets:',dir('data'))
      })
    }
  })
}

clApp<-shinyUI(fluidPage(
  sidebarLayout(
    sidebarPanel(uiOutput('availableDatasets')),
    mainPanel(
      plotOutput('clPlots')
    )
  ),
  flowLayout(
    fileInput('newDataset','Upload a new dataset'),
    radioButtons('smoothMethod','Fit method:',list('lm','glm','gam','loess','none'))
  )#,
  #flowLayout(
  #  checkboxInput('showAnalyticProfile','Display analytic solutions'),
  #  sliderInput('alpha',HTML('&alpha;'),min=0,max=1,step=.0000001),
  #  sliderInput('lambda',HTML('&lambda;'),min=0,max=1,step=.01),
  #  sliderInput('Ns',HTML('N_s'),min=3,max=20),
  #  sliderInput('maxTime','Maximum time to plot analytic solution to',min=200,max=1000)
  #)#,
  #flowLayout(
  #  h2('Simulate data'),
  #  sliderInput('mu',HTML('&mu;'),MU_MIN,MU_MAX,.5*(MU_MIN+MU_MAX),step=.01),
  #  sliderInput('lambda',HTML('&lambda; (cell replacement rate)'),LAMBDA_MIN,LAMBDA_MAX,.5*(LAMBDA_MIN+LAMBDA_MAX),step=.01),
  #  sliderInput('Pr','Pr',PR_MIN,PR_MAX,.5*(PR_MIN+PR_MAX),step=.01),
  #  sliderInput('Ns','N (#stem cells per crypt)',NS_MIN,NS_MAX,.5*(NS_MIN+NS_MAX)),
  #  sliderInput('maxTime','Max time to plot',MAXTIME_MIN,MAXTIME_MAX,.5*(MAXTIME_MIN+MAXTIME_MAX),step=.01),
  #  sliderInput('numRuns','Number of runs to simulate',NUMRUNS_MIN,NUMRUNS_MAX,.5*(NUMRUNS_MIN+NUMRUNS_MAX))
  #)
))
