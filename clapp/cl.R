library(shiny)
library(plyr)
library(dplyr)
library(InferCryptDrift)
library(DriftR)

handleUpload<-function(uploadedDataset){
  uploadedObj<-read.csv(uploadedDataset$datapath)#,sep=' ')
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
    if(length(input$dataset)){
      d<-as.data.frame(readRDS(paste0('data/',input$dataset)))
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
    radioButtons('smoothMethod','Fit method:',list('loess','lm','glm','gam','rlm'))
  )
))
