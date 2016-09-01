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

generatePlots<-function(input){
  if(length(input$dataset)){
    #Set data to be visualised
    d<-as.data.frame(readRDS(paste0('data/',input$dataset)))
    simLine<-NULL
    p2<-NULL
    if(input$showSim){
      timePts<-seq(0,MAXTIME_DISPLAYTIME,length.out=NUM_TIME_POINTS)
      sim<-sim_contLabelling(
        input$mu,
        input$lambda,
        input$N,
        timePts,
        input$numSim,
        input$P,
        T
      )
      simFrame<-data.frame(t(sim))
      simFrame[['time']]<-timePts
      #Correct for fixed tissue sample size
      simFrame[['partial']]<-rowSums(simFrame[,2:input$N])*input$numCrypt/input$numSim
      simFrame[['whole']]<-simFrame[[input$N+1]]*input$numCrypt/input$numSim
      dSim<-rbind.data.frame(
        data.frame(time=simFrame[['time']],value=simFrame[['partial']],type='Simulated Partial Clones'),
        data.frame(time=simFrame[['time']],value=simFrame[['whole']],type='Simulated Whole Clones')
      )
      simLine<-geom_line(data=dSim,aes(x=time,y=value,color=type))
      
      #Plot of distribution of clone sizes at time line's marked at
      distPlot<-as.numeric(as.vector(simFrame[which.min(abs(timePts-input$displayTime)),2:input$N]))
      p2<-ggplot()+
        geom_line(data=data.frame(cells=distPlot,n=1:(input$N-1)),aes(x=n,y=cells))+
        labs(
          x='Proportion of crypt',
          y='count',
          title='distribution of partially mutant populated crypt population'
        )
    }
    p1<-ggplot()+
      geom_point(data=d,aes(x=Age,y=Clone.number,color=Type.Clones))+
      geom_smooth(data=d,method=input$smoothMethod,aes(x=Age,y=Clone.number,color=Type.Clones))+
      geom_vline(xintercept=input$displayTime)+
      facet_wrap(~Location)+
      simLine
    if(input$showSim)p1<-grid.arrange(p1,p2)
    p1
  }
}

handleUpload<-function(uploadedDataset){
  uploadedObj<-read.csv(uploadedDataset$datapath,stringsAsFactors=F)
  uploadedObj$Location=factor(uploadedObj$location,levels=unique(uploadedObj$Location))
  #Assume filename is [NAME].rds for now
  saveName<-paste0(substr(uploadedDataset$name,1,nchar(uploadedDataset$name)-3),'rds')
  saveRDS(uploadedObj,paste0('data/',saveName))
}

#Constants for sliders

serve<-function(input,output){
  output$availableDatasets<-renderUI({
    radioButtons('dataset','Visualise datasets:',dir('data'))
  })
  output$clPlots<-renderPlot({
    generatePlots(input)
  })
  output$genPdf<-downloadHandler(
    filename='clplots.pdf',
    content=function(file){
      ggsave(
        file,
        width=12,
        height=8,
        plot=generatePlots(input),
        device='pdf'
      )
    }
  )
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

addResourcePath('js','js')
clApp<-shinyUI(fluidPage(
  tags$script(src='js/widthhack.js'),
  #Set the side panel width using js?
  titlePanel('Continuous Clonal Labelling Visualisation'),
  sidebarLayout(
    sidebarPanel(uiOutput('availableDatasets')),
    mainPanel(
      plotOutput('clPlots')
    )
  ),
  flowLayout(
    fileInput('newDataset','Upload a new dataset'),
    radioButtons(
      'smoothMethod',
      'Fit method:',
      list(
        'lm (linear model fit)'='lm',
        'loess (smooth interpolation)'='loess',
        'none'
      )
    )
  ),
  flowLayout(
    sliderInput('displayTime','Distribution display time',0,MAXTIME_DISPLAYTIME,.5*MAXTIME_DISPLAYTIME),
    sliderInput('mu',HTML('&mu; (mutation rate)'),0,.0002,.00005,step=.000001),
    sliderInput('lambda',HTML('&lambda; (replacement rate)'),0,.5,.01,step=.001),
    sliderInput('N','N (#stem cells/crypt)',3,20,10),
    sliderInput('numCrypt','Number of crypts in tissue',10000,200000,100000),
    sliderInput('P','Bias (.5 for neutral)',0,1,.5,step=.01),
    sliderInput('numSim','Number of simulations to run',NUMRUNS_MIN,NUMRUNS_MAX,.5*(NUMRUNS_MIN+NUMRUNS_MAX)),
    checkboxInput('showSim','Show simulated curve'),
    downloadButton('genPdf','Download pdf of plots')
  )
))
