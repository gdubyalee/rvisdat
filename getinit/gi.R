library(ggplot2)
library(stringr)
library(shiny)
library(shinyBS)
library(plyr)
library(dplyr)
library(InferCryptDrift)
library(DriftR)
library(readxl)

serve<-function(input,output){
  output$availableDatasets<-renderUI(
    #If dir('dat') is an empty string then Shiny fails in an opaque manner
    radioButtons('dataset','Visualise datasets:',dir('dat'))
  )
  output$giPlots<-renderPlot({
    generatePlots(input)
  })
  #observe({
  #  input$newDataset
  #  if(length(input$newDataset)){
  #    handleUpload(input$newDataset)
  #    output$availableDatasets<-renderUI(
  #      radioButtons('dataset','Visualise datasets:',dir('dat'))
  #    )
  #  }
  #})
  #observe({
  #  input$dataset
  #  output$giPlots<-renderPlot({
  #    loadPlots(input)
  #  })
  #})
}

giApp<-shinyUI(fluidPage(
  titlePanel('Saturation Labelling Visualisation'),
  br(),
  actionButton('xlsInfo','About'),
  sidebarLayout(
    sidebarPanel(
      fileInput('newDataset','Upload a new dataset'),
      uiOutput('availableDatasets'),
      downloadButton('genPdf','Download pdf of plots')
    ),
    mainPanel(
      plotOutput('giPlots')
    )
  ),
  flowLayout(
    sliderInput('mu',HTML('&alpha; (mutation rate)'),0,.0005,.0002,step=.000001),
    sliderInput('lambda',HTML('&lambda; (replacement rate)'),0,.5,.01,step=.001),
    sliderInput('N','N (#stem cells/crypt)',3,20,10),
    radioButtons('plotToView','View Plot:',c(1,2,3,4))#,
    #sliderInput('numCrypt','Number of crypts in tissue',10000,200000,100000),
    #sliderInput('P','Bias (.5 for neutral)',0,1,.5,step=.01),
    #sliderInput('numSim','Number of simulations to run',NUMRUNS_MIN,NUMRUNS_MAX,.5*(NUMRUNS_MIN+NUMRUNS_MAX))
  )
))

setAvDat<-function(){
  renderUI(
    radioButtons('dataset','Visualise datasets:',dir('dat'))
  )
}
handleUpload<-function(uploadedDataset){
  dat<-data.frame()
  #The readxl package sulks if asked to read a file that doesn't end in xls or xlsx
  #Which is pretty unhelpful behaviour, better hack around it
  file.rename(uploadedDataset$datapath,paste0(uploadedDataset$datapath,'.xls'))
  uploadedDataset$datapath=paste0(uploadedDataset$datapath,'.xls')
  #Read data from all sheets
  for(sheet in excel_sheets(uploadedDataset$datapath)){
    dati<-read_excel(uploadedDataset$datapath,sheet)%>%
      gather(mouse,count,starts_with('Mouse_'))%>%
      spread(Type,count)%>%
      gather(Clone,Counts,Partial,Whole)%>%
      filter(!is.na(Counts))
    dat<-rbind(dat,dati)
  }

  #Check field consistency
  allProbs<-mutate(
      dat,
      p_min=qbeta(0.025,Counts,Total-Counts),
      p_median=qbeta(0.5,Counts,Total-Counts),
      p_max=qbeta(0.975,Counts,Total-Counts)
    )
  ppFieldSI<- allProbs%>%
    filter(Location=="SI")%>%
    mutate(Field=factor(Field))%>%
    ggplot(aes(y = p_median, x = mouse, col = Field, fill = Field))+
    geom_bar(position = position_dodge(0.9), stat = "identity")+
    geom_errorbar(aes(ymin = p_min, ymax = p_max), width = 0.4, position = position_dodge(0.9))+ 
    facet_grid(Clone ~ Genotype + Condition, scales = "free") + ggtitle("Counts - SI")
  ppFieldColon<- allProbs%>%
    filter(Location=="Colon")%>%
    mutate(Field=factor(Field))%>%
    ggplot(aes(y = p_median, x = mouse, col = Field, fill = Field))+
    geom_bar(position = position_dodge(0.9), stat = "identity")+
    geom_errorbar(aes(ymin = p_min, ymax = p_max), width = 0.4, position = position_dodge(0.9))+ 
    facet_grid(Clone ~ Genotype + Condition, scales = "free") + ggtitle("Counts - Colon")
  dataMouse<-dat%>%
    group_by(Genotype,Condition,Location,mouse,Clone)%>%
    summarize(Counts=sum(Counts),Total=sum(Total))%>%
    ungroup()
  allProbsMouse=mutate(
    dataMouse,
    p_min=qbeta(0.025,Counts,Total-Counts),
    p_median=qbeta(0.5,Counts,Total-Counts),
    p_max=qbeta(0.975,Counts,Total-Counts)
  )
  ppMouseColon<-allProbsMouse%>%
    filter(Location=='Colon')%>%
    ggplot(aes(y = p_median, x = mouse, col = mouse, fill = mouse)) + 
    geom_point(position = position_dodge(0.9), stat = "identity", size = 2) +
    geom_errorbar(aes(ymin = p_min, ymax = p_max), width = 0.4, position = position_dodge(0.9)) + 
    scale_x_discrete(breaks=NULL) +
    ylab("Fraction of crypts") +
    xlab("") +
    facet_grid(Clone~Genotype+Condition, scales = "free_x") +
    ggtitle("Counts - Colon")
  ppMouseSI<-allProbsMouse%>%
    filter(Location=='SI')%>%
    ggplot(aes(y = p_median, x = mouse, col = mouse, fill = mouse)) + 
    geom_point(position = position_dodge(0.9), stat = "identity", size = 2) +
    geom_errorbar(aes(ymin = p_min, ymax = p_max), width = 0.4, position = position_dodge(0.9)) + 
    scale_x_discrete(breaks=NULL) +
    ylab("Fraction of crypts") +
    xlab("") +
    facet_grid(Clone~Genotype+Condition, scales = "free_x") +
    ggtitle("Counts - SI")

  #Pulse chase simulation
  saveRDS(
    list(
      ppFieldColon,
      ppFieldSI,
      ppMouseColon,
      ppMouseSI
    ),
    paste0('dat/',uploadedDataset$name,'.rds')
  )
}

generatePlots<-function(input){
  if(!is.null(input$dataset)){
    d<-readRDS(paste0('dat/',input$dataset))
    return(d[as.numeric(input$plotToView)])
  }
}
