library(shiny)
library(shinyBS)
shinyUI(fluidPage(
  tags$script(src='togglefs.js'),
  sidebarPanel(
    h3('Pulse Chase Data'),
    'Pull down data:',
    textInput('searchStr','Dataset name'),
    checkboxInput('doMcmc','Perform best fit analysis?'),
    checkboxInput('plotConvergence','Display convergence plot?'),
    checkboxInput('forceRedo','Force fitting and convergence tests to be redone?'),
    submitButton('Go'),
    #dataTableOutput('dataset'),
    plotOutput('mcmc'),
    #actionButton('viewMcmcFS','View MCMC fitting plot'),
    #bsModal('mcmcModal','Mcmc plot','viewMcmcFS',plotOutput('mcmc')),
    plotOutput('convergence'),
    width=6
  ),
  sidebarPanel(
    h3('Continuous labelling data'),
    width=6
  )
))
