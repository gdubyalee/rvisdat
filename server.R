library(dplyr)
library(tidyr)
library(DriftR)
library(stringr)
library(InferCryptDrift)
#source('pchase/pchaseutils.R')
#source('ctslab/ctslabutils.R')

#Start up our server...
library(shiny)
shinyServer(function(input,output){

  #Display data in table
  #output$dataset<-renderDataTable({if(is.character(input$searchStr)&&nchar(input$searchStr)>0){
  #  readRDS(paste0('data/',input$searchStr))
  #}})


  rp<-renderPlot({
    if(input$doMcmc){
      #Cache analysis once it's been done on a dataset
      if(
        !file.exists(paste0('cache/mcmc_',input$searchStr))
        ||
        input$forceRedo
      ){
        cryptData<-readRDS(paste0('data/',input$searchStr))
        times<-strtoi(substring(names(cryptData),2))
        neutralDrift<-fitNeutralDrift(cryptData,times)
        params<-getNeutralDirftParams(neutralDrift)
        neutralDriftFitted<-plotsNeutralDrift_Fit(
          cryptData,
          times,
          params,
          max_x=100,
          input$searchStr
        )
        saveRDS(neutralDriftFitted,paste0('cache/mcmc_',input$searchStr))
      }else{
        neutralDriftFitted<-readRDS(paste0('cache/mcmc_',input$searchStr))
      }
      plot(neutralDriftFitted)
    }
  })
  output$mcmc<-rp
  #output$mcmc2<-rp
})
