library(DriftR)
library(InferCryptDrift)

#Start up our server...
library(shiny)
shinyServer(function(input,output){
  output$dataset<-readRDS(paste0('data/',input$searchStr))
})
