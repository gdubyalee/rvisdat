library(DriftR)
library(stringr)
library(InferCryptDrift)

#Start up our server...
library(shiny)
shinyServer(function(input,output){

  #Display data in table
  #output$dataset<-renderDataTable({if(is.character(input$searchStr)&&nchar(input$searchStr)>0){
  #  readRDS(paste0('data/',input$searchStr))
  #}})

  output$convergence<-renderPlot({
    if(input$plotConvergence){
      cryptData<-readRDS(paste0('data/',input$searchStr))
      times<-strtoi(substring(names(cryptData),2))
      if(!file.exists(paste0('data/neutdrift_',input$searchStr))||input$forceRedo){
        neutralDrift<-fitNeutralDrift(cryptData,times)
        saveRDS(neutralDrift,paste0('data/neutdrift_',input$searchStr))
      }else{
        neutralDrift<-readRDS(paste0('data/neutdrift_',input$searchStr))
      }
      plotsConvergence_Neutral(neutralDrift)
    }
  })

  rp<-renderPlot({
    if(input$doMcmc){
      #Cache analysis once it's been done on a dataset
      if(!file.exists(paste0('data/mcmc_',input$searchStr))||input$forceRedo){
        cryptData<-readRDS(paste0('data/',input$searchStr))
        times<-strtoi(substring(names(cryptData),2))
        neutralDrift<-fitNeutralDrift(cryptData,times)
        params<-getNeutralDirftParams(neutralDrift)
        neutralDriftFitted<-plotsNeutralDrift_Fit(
          cryptData,
          times,
          params,
          max_x=100,
          paste0('Pulse chase data for ',input$searchStr)
        )
        saveRDS(neutralDriftFitted,paste0('data/mcmc_',input$searchStr))
      }else{
        neutralDriftFitted<-readRDS(paste0('data/mcmc_',input$searchStr))
      }
      plot(neutralDriftFitted)
    }
  })
  output$mcmc<-rp
  #output$mcmc2<-rp
})
