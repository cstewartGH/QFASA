## ui.R ##
library(htmltools)
library(shiny)
library(rCharts)

########################################
## UI
########################################
ui <- shinyUI(
    fluidPage(fluidRow(column(offset=0, width=4, selectInput("prey",
                                                             label = h3("Prey"), 
                                                             choices = getprey(),
                                                             selected = 1,
                                                             multiple=TRUE)),
                       column(offset=0, width=8, selectInput("fattyacid",
                                                             label = h3("Fatty Acids"), 
                                                             choices = getfattyacids(),
                                                             selected = 1,
                                                             multiple=TRUE))),
              actionButton("goButton", "Go"),
              hr(),
              
              # Plots
              fluidRow(chartOutput('dietestbydist', lib='nvd3')),
              fluidRow(chartOutput('dietestbyprey', lib='nvd3')),
              fluidRow(chartOutput('addmeasbydist', lib='nvd3')),
              fluidRow(chartOutput('addmeasbyprey', lib='nvd3'))
              )
)


