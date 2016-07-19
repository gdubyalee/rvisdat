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
MAXTIME_MAX<-1000
MAXTIME_DISPLAYTIME<-1000
NUMRUNS_MIN<-10000
NUMRUNS_MAX<-1000000
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
      d<-as.data.frame(readRDS(paste0('data/',input$dataset)))
      simLine<-NULL
      if(input$showSim){
        timePts<-seq(0,MAXTIME_DISPLAYTIME,length.out=NUM_TIME_POINTS)
        sim<-sim_contLabelling(
          input$mu,
          input$lambda,
          input$N,
          timePts,
          input$numSim,
          .5,
          T
        )
        simFrame<-data.frame(t(sim))
        simFrame[['time']]<-timePts
        simFrame[['partial']]<-rowSums(simFrame[,2:input$N])
        simFrame[['whole']]<-simFrame[[input$N+1]]
        print(head(simFrame))
        dSim<-rbind.data.frame(
          data.frame(time=simFrame[['time']],value=simFrame[['partial']],type='partial'),
          data.frame(time=simFrame[['time']],value=simFrame[['whole']],type='whole')
        )
        print(head(dSim))
        simLine<-geom_line(data=dSim,aes(x=time,y=value,color=type))
      }
      p1<-ggplot()+
        geom_point(data=d,aes(x=Age,y=Clone.number,color=Type.Clones))+
        geom_smooth(data=d,method=input$smoothMethod,aes(x=Age,y=Clone.number,color=Type.Clones))+
        geom_vline(xintercept=input$displayTime)+
        facet_wrap(~Location)+
        simLine
      #p2<-NULL
      #grid.arrange(p1,p2,nRow=2)
      p1
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
  titlePanel('Continuous Clonal Labelling Visualisation'),
  sidebarLayout(
    sidebarPanel(uiOutput('availableDatasets')),
    mainPanel(
      plotOutput('clPlots')
    )
  ),
  flowLayout(
    fileInput('newDataset','Upload a new dataset'),
    radioButtons('smoothMethod','Fit method:',list('lm','glm','gam','loess','none'))
  ),
  flowLayout(
    sliderInput('displayTime','Distribution display time',0,MAXTIME_DISPLAYTIME,.5*MAXTIME_DISPLAYTIME),
    sliderInput('mu',HTML('&mu;'),0,.0001,.00005,step=.000001),
    sliderInput('lambda',HTML('&lambda;'),0,.5,.01,step=.001),
    sliderInput('N','N',3,20,10),
    sliderInput('numSim','Number of simulations to run',NUMRUNS_MIN,NUMRUNS_MAX,.5*(NUMRUNS_MIN+NUMRUNS_MAX)),
    checkboxInput('showSim','Show simulated curve')
  #  checkboxInput('showAnalyticProfile','Display analytic solutions'),
  #  sliderInput('alpha',HTML('&alpha;'),min=0,max=1,step=.0000001),
  #  sliderInput('lambda',HTML('&lambda;'),min=0,max=1,step=.01),
  #  sliderInput('Ns',HTML('N_s'),min=3,max=20),
  #  sliderInput('maxTime','Maximum time to plot analytic solution to',min=200,max=1000)
  )#,
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
