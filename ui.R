library(shiny)
require(d3heatmap)
require(igraph)
require(dplyr)
require(reshape)
require(networkD3)
require(d3Network)
require(ggplot2)
require(huge)
require(glasso)
require(scales)
require(plsgenomics)
library(flare)
library(ggm)
data(leukemia)
leukaemia <- leukemia
data(Colon)
data(SRBCT)
source("help.R")

shinyUI(fluidPage(
  
  titlePanel("Generation and visualization data"),
  
  shinyUI(navbarPage("Gene Application",
  #1st tab______________________________________________________
                     tabPanel("Comparison of gene net types",sidebarPanel(
                       selectInput("type",label = ("Select type of graphs:"),
                                   choices = c("random","hub","band","cluster", "scale-free")
                                   ,selected = "hub"),
                       sliderInput("nodes", 
                                   "Number of genes:", 
                                   min = 0, 
                                   max = 200, 
                                   value = 50, 
                                   step=1),
                       uiOutput("groups")
                       
                      
                     ),
  #1st main Panel______________________________________________________                     
                     mainPanel(
                       #column(width = 15, plotOutput("distPlot")),
                       column(width = 12, h3("Heat map"), d3heatmapOutput("heatmap")),
                       column(width = 12, h3("Network graphs"), forceNetworkOutput("network"))
                       
                     )),
  #2nd tab______________________________________________________
                     tabPanel("Estimation",sidebarPanel(
                       selectInput("method",label = ("Select method of estimation:"),
                                   choices = c("glasso","mb")
                                   ,selected = "glasso"),
                       sliderInput("lambda",
                                   "Number of lambdas:",
                                   min = 1,
                                   max = 10,
                                   value = 5,
                                   step = 1),
                       uiOutput("lambda1")
                     ),
  #2nd main Panel______________________________________________________
                     mainPanel(
                       #column(width = 15, plotOutput("distPlot")),
                       #column(width = 12, d3heatmapOutput("heatmap")),
                       column(width = 12, h3("Network graphs"), forceNetworkOutput("network1")),
                       column(width = 12, h3("Estimated network graphs"), forceNetworkOutput("network2"))
                     ))
  ,
                     tabPanel("The actual data analysis",
                              
                              
                              column(2, 
                                     selectInput("data",label = "Select data: ", choices = list("leukaemia","Colon","SRBCT")),
                                     numericInput("dlugosc",
                                                  label = "Length of data",
                                                  min=50,
                                                  max=150,
                                                  value=25,
                                                  step=25),
                                     strong("Gene numbers:"),
                                     numericInput("od",
                                                  label = "From",
                                                  min = 1,max=3051,
                                                  value = 100,
                                                  step=25),
                                     uiOutput("do"),
                                     sliderInput("k",label = "Confidence ",min = 1,max = 10,value = 3),
                                     strong("Names of genes"),
                                     checkboxInput("nazwy", label = "show", value = TRUE)
                                     ),
                              column(5,selectInput("algorithm1", label = "Select first algorithm:", 
                                                   choices= list("glasso","mb","tiger","clime")),
                                     sliderInput("nlambda1",label = "Number of lambda for first model",
                                                 min = 1,
                                                 max=30,
                                                 value=10),
                                     uiOutput("lambdan1"),
                                     h4("Heat map 1"),
                                     d3heatmapOutput("heatmap1"),
                                     h4("Network graph 1"),
                                     forceNetworkOutput("networkk1")),
                              column(5,      selectInput("algorithm2", label = "Select second algorithm:", 
                                                         choices= list("glasso","mb","tiger","clime")),
                                     sliderInput("nlambda2",label = "Number of lambda for second model",
                                                 min = 1,
                                                 max=30,
                                                 value=10),
                                     uiOutput("lambdan2"),
                                     h4("Heat map 2"),
                                     d3heatmapOutput("heatmap2"),
                                     h4("Network graph 2"),
                                     forceNetworkOutput("networkk2"))
                              
                              
                              
                              
                            )
  ))
    
)
)
