library(shiny)
shinyUI(fluidPage(
  sidebarPanel(
    h3('Data visualisation ui'),
    'Pull down data:',
    textInput('searchStr','Dataset name'),
    submitButton('Pull dataset'),
  ),
  uiOutput('dataset')
))
