library(ggplot2)
library(stringr)
library(shiny)
library(shinyBS)
library(plyr)
library(dplyr)
library(InferCryptDrift)
library(DriftR)
library(readxl)
library(rmutil)
#source('sim.R')

serve<-function(input,output){
  output$availableDatasets<-renderUI(
    #If dir('dat') is an empty string then Shiny fails in an opaque manner
    radioButtons('dataset','Visualise datasets:',dir('dat'))
  )
  output$giPlots<-renderPlot({
    generatePlots(input)
  })
  observe({
    input$newDataset
    if(length(input$newDataset)){
      handleUpload(input$newDataset,input)
      output$availableDatasets<-renderUI(
        #If dir('dat') is an empty string then Shiny fails in an opaque manner
        radioButtons('dataset','Visualise datasets:',dir('dat'))
      )
    }
  })
  observe({
    input$dataset
    output$giPlots<-renderPlot({
      generatePlots(input)
    })
  })
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
    #sliderInput('mu',HTML('&alpha; (mutation rate)'),0,.0005,.0002,step=.000001),
    sliderInput('lambda',HTML('&lambda; (replacement rate)'),0,2,.2,step=.001),
    sliderInput('N','N (#stem cells/crypt)',3,20,5),
    sliderInput('T','Time since treatment',5,100,30),
    radioButtons('plotToView','View Plot:',c('ppFieldSI'=1,'ppFieldColon'=2,'ppMouseColon'=3,'ppMouseSI'=4,'Initial relationship'=5))#,
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
handleUpload<-function(uploadedDataset,input){
  t<-input$T
  Ns<-input$N
  lambda<-input$lambda
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
    group_by(Genotype,Condition,Clone)%>%
    mutate(meanMedian=mean(p_median))%>%
    ggplot(aes(y = p_median, x = mouse, col = mouse, fill = mouse)) + 
    geom_point(position = position_dodge(0.9), stat = "identity", size = 2) +
    geom_errorbar(aes(ymin = p_min, ymax = p_max), width = 0.4, position = position_dodge(0.9))+
    scale_x_discrete(breaks=NULL) +
    ylab("Fraction of crypts") +
    xlab("") +
    facet_grid(Clone~Genotype+Condition, scales = "free_x") +
    ggtitle("Counts - Colon")+
    geom_hline(aes(yintercept=meanMedian),size=.7,lty=2)
  
  ppMouseSI<-allProbsMouse%>%
    filter(Location=='SI')%>%
    group_by(Genotype,Condition,Clone)%>%
    mutate(meanMedian=mean(p_median))%>%
    ggplot(aes(y = p_median, x = mouse, col = mouse, fill = mouse)) + 
    geom_point(position = position_dodge(0.9), stat = "identity", size = 2) +
    geom_errorbar(aes(ymin = p_min, ymax = p_max), width = 0.4, position = position_dodge(0.9))+
    scale_x_discrete(breaks=NULL) +
    ylab("Fraction of crypts") +
    xlab("") +
    facet_grid(Clone~Genotype+Condition, scales = "free_x") +
    ggtitle("Counts - SI")+
    geom_hline(aes(yintercept=meanMedian),size=.7,lty=2)

  initialRelation<-initialRelationPlot(lambda,Ns,t,max(allProbsMouse$p_median))

  saveRDS(
    list(
      ppFieldColon,
      ppFieldSI,
      ppMouseColon,
      ppMouseSI,
      initialRelation
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

solveIVP<-function(init,lambda,Ns,t){
  #Generate matrix
  A<-matrix(0,nrow=Ns+1,ncol=Ns+1)
  for(i in 1:Ns){
    A[i,i]=-2
    A[i,i+1]=1
    A[i+1,i]=1
  }
  A[,Ns+1]=A[,1]=0
  return(mexp(A,lambda*t)%*%init)
}

initialRelationPlot<-function(lambda,Ns,t,probMax){
  probMax<-min(probMax,1)
  probs<-seq(0,probMax,length.out=50)
  fracPar<-c()
  fracFull<-c()
  for(pExpressing in probs){
    ic<-dbinom(0:Ns,Ns,pExpressing)
    sln<-solveIVP(ic,lambda,Ns,30)
    fracPar[length(fracPar)+1]<-sum(sln[2:Ns])
    fracFull[length(fracFull)+1]<-sln[Ns+1]
  }
  bind_rows(
    data.frame('fraction_expressing_initially'=probs,frac=fracPar,type='Fraction Partial'),
    data.frame('fraction_expressing_initially'=probs,frac=fracFull,type='Fraction Clonal')
  )%>%
  ggplot()+
    geom_line(aes(x=fraction_expressing_initially,y=frac,col=type,xlabel='Initial proportion of cells expressing',ylabel='Crypt frequency'))+
    aes(ylabel='Fraction expressing initially')+
    ggtitle(paste0('Proportion of partial and monoclonal crypts for Ns=',Ns,', T=',t,', \\lambda=',lambda))
}
