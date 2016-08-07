library(shiny)
require(igraph)
require(dplyr)
require(reshape)
require(networkD3)
require(d3Network)
require(ggplot2)
require(huge)
require(glasso)
require(scales)
library(plsgenomics)
library(flare)
library(ggm)

require(d3heatmap)
data(leukemia)
leukaemia <- leukemia
data(Colon)
data(SRBCT)
source("help.R")

shinyServer(function(input, output) {
  #______________________________________________________
  output$groups <- reactiveUI(function() {
    sliderInput("groups", 
                "Number of groups:", 
                min = 0,  
                max = input$nodes, 
                value = 5)
  })
  #______________________________________________________
  output$lambda1 <- reactiveUI(function() {
    sliderInput("lambda1", 
                "Lambda to calculate graph:", 
                min = 0,  
                max = input$lambda, 
                value = 4,
                step=1)
  })
  #______________________________________________________
  dataInput <- reactive({
    huge.generator(d = input$nodes ,
                   graph = input$type, 
                   g = input$groups)
  })
  #______________________________________________________
  dataOutput <- reactive({
    huge(dataInput()$data ,
         method = input$method, 
         nlambda = input$lambda)
  })
  #______________________________________________________
  output$distPlot <- renderPlot({
    plot(dataInput())
  })
  #______________________________________________________
  output$heatmap <- renderD3heatmap(
    plot_heatmap(dataInput(),input$nodes)
  )
  #______________________________________________________
  output$network <- renderForceNetwork(
    plot_network(dataInput(),
                 input$nodes, 
                 g = input$groups)
  )
  #______________________________________________________
  output$network1 <- renderForceNetwork(
    plot_network(dataInput(),
                 input$nodes, 
                 g = input$groups)
  )
  #______________________________________________________
  output$network2 <- renderForceNetwork(
    plot_network_emp(dataOutput(), 
                     p = input$nodes, 
                     lambda_n = input$lambda1, 
                     g = input$groups)
  )
  
  # 3 zakladka___________________________________________________________________________________
  output$do <- reactiveUI(function(){
    textInput(inputId = "do",label = "To",value = input$od+input$dlugosc)
    # numericInput("do",label = "to",min = 1,max=input$od+input$dlugosc,value = 125,step=25)
  })
  output$lambdan1 <- reactiveUI(function() {
    sliderInput("lambdan1", 
                "Use lambda:", 
                min = 0,  
                max = input$nlambda1, 
                value = 5,
                step = 1)
  })
  
  
  dataInput1 <- reactive({
    dobor_modelu(input$data,
                 od=input$od,
                 do=input$do,
                 method = input$algorithm1,
                 nlambda = input$nlambda1)
    
  })
  
  
  
  output$heatmap1 <- renderD3heatmap({
    
    plot_heatmap1(dataInput1(),lambda_n = input$lambdan1,names = nazwa(data=input$data,od = input$od,do=input$do,czynazwy = input$nazwy))
  })
  
  output$networkk1 <- renderForceNetwork({
    plot_network1(dataInput1(),lambda_n = input$lambdan1,names = nazwa(data=input$data,od = input$od,do=input$do,czynazwy = input$nazwy),k = input$k )
  })
  
  
  output$lambdan2 <- reactiveUI(function() {
    sliderInput("lambdan2", 
                "Use lambda", 
                min = 0,  
                max = input$nlambda2, 
                value = 5,
                step=1)
  })
  
  dataInput2 <- reactive({
    dobor_modelu(input$data,
                 od=input$od,do=input$do,method = input$algorithm2,nlambda = input$nlambda2)
    
  })
  
  
  
  output$heatmap2 <- renderD3heatmap({
    plot_heatmap1(dataInput2(),lambda_n = input$lambdan2,names = nazwa(data=input$data,od = input$od,do=input$do,czynazwy = input$nazwy))
  })
  
  output$networkk2 <- renderForceNetwork({
    plot_network1(dataInput2(),lambda_n = input$lambdan2,names = nazwa(data=input$data,od = input$od,do=input$do,czynazwy = input$nazwy),k = input$k)
  })
  
  
  
})
