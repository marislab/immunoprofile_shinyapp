library(shiny)
library(shinydashboard)
library(shinyIncubator)
library(shinyWidgets)
library(shinyBS)
library(DT)
library(dplyr)
library(reshape2)
library(ggthemes)
library(ggplot2)
library(plotly)
library(pheatmap)
library(ggpubr)
library(grid)
library(heatmaply)
library(sunburstR)
library(d3r)
options(gsubfn.engine = "R")
options(shiny.sanitize.errors = TRUE)

dashboardPage(
  
  # dashboardHeader begins
  dashboardHeader(title = 'Immunoprofile Portal', titleWidth = 300, dropdownMenuOutput("messageMenu")), # dashboardHeader ends
  
  # dashboardSidebar begins
  dashboardSidebar(width = 300,
                   
                   tags$head(
                     tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
                   
                   # enable vertical scrolling
                   div(style="overflow-y: scroll"),
                   
                   # sidebarMenu begin
                   sidebarMenu(
                     menuItem("Dashboard", icon = icon("dashboard"), tabName = "dashboard"),
                     menuItem("Gene Aliases", icon = icon("question-circle"), tabName = "genehelp"),
                     menuItem("Expression Plots", icon = icon("bar-chart"), tabName = "box"),
                     menuItem("GSVA Plots", icon = icon("bar-chart"), tabName = "gsva"),
                     menuItem("Heatmaps", icon = icon("th"), tabName = "heatmap")
                   ) # sidebarMenu
  ), # dashboardSidebar
  
  dashboardBody(
    
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    
    div(style="overflow-x: scroll"),
    
    # tabItems begins
    tabItems(
      
      # dashboard content
      tabItem(tabName = "dashboard",
              fluidRow(
                box(title = "Overview", status = "danger", width = 12, 
                    collapsed = F, collapsible = F, solidHeader = T, 
                    tags$p(style = "font-size: 16px; font-weight: bold;", 
                          HTML("Immunoprofiling Project <br/> 67 Histologies across 4 Datasets (TCGA, TARGET, PNOC and CBTTC) <br/> 31 Gene signatures consisting of 281 unique genes"))),
                box(title = "Data Summary", status = "warning", width = 6, 
                    collapsed = F, collapsible = F, solidHeader = T, 
                    sunburstOutput(outputId = "datsum1")),
                box(title = "Breakdown by Dataset", status = "warning", width = 6,
                    collapsed = F, collapsible = F, solidHeader = T,
                    plotlyOutput(outputId = "datsum2", height = 400))
              ),
              fluidRow(
                box(title = "Gene Signatures", status = "warning", width = 12, height = "500px",
                    collapsed = F, collapsible = F, solidHeader = T,
                    div(style="max-height:400px; overflow-y: scroll; position: relative", plotlyOutput(outputId = "datsum3", height = 600)))
              )
              ),
      tabItem(tabName = "genehelp",
              div(DT::dataTableOutput("genehelptable"), style = "font-size: 12px")
              ),
      tabItem(tabName = "gsva",
              fluidRow(
                box(background = "navy", width = 12,
                  column(2, pickerInput(inputId = "gsvaselectInput1", label = "Dataset", choices = "none", options = list(`actions-box` = TRUE), multiple = TRUE, choicesOpt = list(style = "font-size: 12px;"))),
                  column(3, pickerInput(inputId = "gsvaselectInput2", label = "Histology", choices = "none", options = list(`actions-box` = TRUE), multiple = TRUE)),
                  column(3, selectInput(inputId = "gsvaselectInput3", label = "Signature", choices = "none")),
                  column(2, selectInput(inputId = "gsvaselectInput4", label = "Type", choices = c("log2counts", "z-scores"))), br(),
                  column(2, actionButton(inputId = "gsvasubmit1", label = "Get Boxplot", icon("paper-plane"), style = "font-size: 14px; margin-top: 6px; padding:8px;color: #fff; background-color: #337ab7; border-color: #2e6da4"))
                  )
              ),
              tabsetPanel(type = "tabs",
                          tabPanel("Histology View", div(style="overflow-x: scroll; overflow-y: scroll", plotlyOutput("gsvaplot1", width = "100%", height = 600))),
                          tabPanel("Raw Data", div(DT::dataTableOutput("gsvatable"), style = "font-size: 12px")))
              ),
      tabItem(tabName = "heatmap",
              fluidRow(
                box(background = "navy", width = 12,
                    column(2, pickerInput(inputId = "heatmapselectInput1", label = "Dataset", choices = "none", options = list(`actions-box` = TRUE), multiple = TRUE)),
                    column(3, pickerInput(inputId = "heatmapselectInput2", label = "Histology", choices = "none", options = list(`actions-box` = TRUE), multiple = TRUE)),
                    column(3, selectInput(inputId = "heatmapselectInput3", label = "Signature", choices = "none")),
                    column(2, selectInput(inputId = "heatmapselectInput4", label = "Type", choices = c("Mean FPKM" = "mean.fpkm", "Mean z-score"="mean.zscore", "Median z-score"="median.zscore", "25th Perc (z-score)" = "zscore.25", "75th Perc (z-score)" = "zscore.75", "95th Perc (z-score)" = "zscore.95"))), br(),
                    column(2, actionButton(inputId = "heatmapsubmit1", label = "Get Heatmap", icon("paper-plane"), style = "font-size: 14px; margin-top: 6px; padding:8px;color: #fff; background-color: #337ab7; border-color: #2e6da4"))
                    )
              ),
              tabsetPanel(type = "tabs",
                          tabPanel("Heatmap", div(style="overflow-x: scroll; overflow-y: scroll", plotlyOutput(outputId = "heatmap1", width = "100%", height = 1200))),
                          tabPanel("Raw Data", div(DT::dataTableOutput("heatmaptable"), style = "font-size: 12px"))
                          )
              ),
      tabItem(tabName = "box",
              fluidRow(
                box(background = "navy", width = 12,
                    column(2, pickerInput(inputId = "boxselectInput1", label = "Dataset", choices = "none", options = list(`actions-box` = TRUE), multiple = TRUE)),
                    column(3, pickerInput(inputId = "boxselectInput2", label = "Histology", choices = "none", options = list(`actions-box` = TRUE), multiple = TRUE)),
                    column(2, selectInput(inputId = "boxselectInput3", label = "Gene", choices = "none")), 
                    column(2, selectInput(inputId = "boxselectInput4", label = "Type", choices = c("logFPKM", "FPKM", "z-score"))), br(),
                    column(2, actionButton(inputId = 'boxsubmit1', label = "Get Boxplot", icon("paper-plane"), style = "font-size: 14px; margin-top: 6px; padding:8px;color: #fff; background-color: #337ab7; border-color: #2e6da4")))
              ),
              tabsetPanel(type = "tabs",
                          tabPanel("Histology View", div(plotlyOutput("boxplot1", width = "100%", height = 600), style="overflow-x: scroll")),
                          tabPanel("Raw Data", div(DT::dataTableOutput("boxtable"), style = "font-size: 12px"))
              ))
    ) # tabItems ends
  ) # dashboardBody ends
) # dashboardPage ends