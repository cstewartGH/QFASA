## ui.R ##
library(htmltools)
library(shinydashboard)
library(rCharts)

########################################
## UI
########################################
ui <- dashboardPage(
    dashboardHeader(title = "QFASA"),
    dashboardSidebar(disable=TRUE),
    dashboardBody(
        fluidRow(column(width=3,
                        box(width=NULL,
                            selectInput("prey",
                                        label = h3("Prey"), 
                                        choices = getprey(),
                                        selected = 1,
                                        multiple=TRUE),
                            selectInput("fattyacid",
                                        label = h3("Fatty Acids"), 
                                        choices = getfattyacids(),
                                        selected = 1,
                                        multiple=TRUE),
                            actionButton("goButton", "Go")
                            )),
              
                 ## Plots                        
                 column(width=9,
                        box(width=NULL, title='Diet Estimates by Prey Type',
                            chartOutput('dietestbyprey', lib='nvd3')
                           ),
                        box(width=NULL, title='Proportional Contribution by Fatty Acid',
                            chartOutput('addmeasbyprey', lib='nvd3')
                            )
                        )
                 )
    )
)





