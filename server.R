# load datasets
load('data/filtered_expr_with_meta.RData')
load('data/histology_labs.RData')
gsva.zscore <- readRDS('data/Immunoprofiling_signatures.RDS')
gsva.log2counts <- readRDS('data/Immunoprofiling_signatures_log2counts.RDS')

# source R code
source('R/themes.R')
source('R/create.box.disease.R')
source('R/create.hm.R')
source('R/viewDataTable.R')
source('R/create.box.gsva.R')

shinyServer(function(input, output, session){
  
  # gene names and signatures
  dataset <- read.delim('data/ImmuneSig_GeneList_OBD_final_2_19_2019.txt', stringsAsFactors = F)
  dataset$Signature <- trimws(dataset$Signature)
  genelist.ct <- plyr::count(dataset$Signature)
  genes <- unique(dataset$newName)
  sigs <- unique(dataset$Signature)
  
  # histology breakdown by study
  load('data/histology_labs.RData')
  studies <- unique(lb$study_id)
  
  # update all select boxes with datasets, gene names and signatures
  observe({
    
    # dataset
    updatePickerInput(session = session, inputId = "boxselectInput1", choices = studies, choicesOpt = list(style = rep_len("font-size: 12px;", length(studies))))
    updatePickerInput(session = session, inputId = "heatmapselectInput1", choices = studies, choicesOpt = list(style = rep_len("font-size: 12px;", length(studies))))
    updatePickerInput(session = session, inputId = "gsvaselectInput1", choices = studies, choicesOpt = list(style = rep_len("font-size: 12px;", length(studies))))
    
    # signatures
    updateSelectizeInput(session = session, inputId = "gsvaselectInput3", choices = sigs, server = TRUE)
    updateSelectizeInput(session = session, inputId = "heatmapselectInput3", choices = sigs, server = TRUE)
    
    # genes
    updateSelectizeInput(session = session, inputId = "boxselectInput3", choices = genes, server = TRUE)
    
  })
  
  # expression boxplot
  observe({
    selected.data <- input$boxselectInput1
    hist <- lb[which(lb$study_id %in% selected.data),"disease"]
    updatePickerInput(session = session, inputId = "boxselectInput2", choices = hist, choicesOpt = list(style = rep_len("font-size: 12px;", length(hist))))
  })

  # heatmap  
  observe({
    selected.data <- input$heatmapselectInput1
    hist <- lb[which(lb$study_id %in% selected.data),"disease"]
    updatePickerInput(session = session, inputId = "heatmapselectInput2", choices = hist, choicesOpt = list(style = rep_len("font-size: 12px;", length(hist))))
  })

  # gsva boxplot  
  observe({
    selected.data <- input$gsvaselectInput1
    hist <- lb[which(lb$study_id %in% selected.data),"disease"]
    updatePickerInput(session = session, inputId = "gsvaselectInput2", choices = hist, choicesOpt = list(style = rep_len("font-size: 12px;", length(hist))))
  })
  
  # dashboard data summary
  output$datsum1 <- renderSunburst({
    sb.dat <- plyr::count(x[,c('study_id','disease_name','definition')])
    colnames(sb.dat) <- c("level1","level2","level3","size")
    sn <- sunburst(
      d3_nest(sb.dat, value_cols = "size"),
      count = T, percent = F, 
      legend = F, withD3 = T)
    sn
  })
  
  output$datsum2 <- renderPlotly({
    plot_ly() %>%
      add_pie(data = count(x[x$study_id == "CBTTC",], study_id, disease_name), 
              labels = ~disease_name, values = ~n,
              name = "CBTTC", domain = list(row = 0, column = 0),
              textposition = 'inside') %>%
      add_pie(data = count(x[x$study_id == "TCGA",], study_id, disease_name), 
              labels = ~disease_name, values = ~n,
              name = "TCGA", domain = list(row = 0, column = 1),
              textposition = "inside") %>%
      add_pie(data = count(x[x$study_id == "TARGET",], study_id, disease_name), 
              labels = ~disease_name, values = ~n,
              name = "TARGET", domain = list(row = 1, column = 0),
              textposition = "inside") %>%
      add_pie(data = count(x[x$study_id == "PNOC",], study_id, disease_name), 
              labels = ~disease_name, values = ~n,
              name = "PNOC", domain = list(row = 1, column = 1),
              textposition = "inside") %>%
      layout(showlegend = F, grid=list(rows=2, columns=2),
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             annotations = list(list(x = 0.2, y = 1.05, text = "CBTTC", showarrow = F, xref='paper', yref='paper'),
                                list(x = 0.8, y = 1.05, text = "TCGA", showarrow = F, xref='paper', yref='paper'),
                                list(x = 0.2, y = -0.05, text = "TARGET", showarrow = F, xref = 'paper', yref = 'paper'),
                                list(x = 0.8, y = -0.05, text = "PNOC", showarrow = F, xref = 'paper', yref = 'paper')))
  })
  
  output$datsum3 <- renderPlotly({
    sig.plot <- dataset %>% 
      group_by(Signature) %>% 
      summarise(count = n(), newName = toString(newName)) %>%
      as.data.frame()
    sig.plot$Signature <- reorder(sig.plot$Signature, sig.plot$count)
    p <- plot_ly(data = sig.plot, 
                 x = ~Signature,
                 y = ~count, 
                 text = ~newName,
                 type = "bar") %>%
      layout(xaxis = list(title = "", tickangle = -90, tickfont = list(size = 14)),
             yaxis = list(title = "Number of Genes", tickfont = list(size = 14)),
             title = "Gene Signatures")
    p
  })
  
  # gene alias table
  output$genehelptable <- DT::renderDataTable({
    viewDataTable(dat = dataset, pageLength = 5)
  })
  
  # create expression plot
  output$boxplot1 <- renderPlotly({
    if(input$boxsubmit1 == 0){
      return()
    }
    withProgress(session = session, message = "Plotting Data...", detail = "Takes a while...", min = 1, value = 10, max = 10,{
      isolate({
        hist = input$boxselectInput2
        gene = input$boxselectInput3
        type = input$boxselectInput4
        bx <<- create.box.disease(total.sub, gene = gene, type = type, hist = hist)
        bx[[1]]
      })
    })
  })
  
  # get raw data expression data
  output$boxtable <- renderDataTable({
    if(input$boxsubmit1 == 0){
      return()
    }
    isolate({
      dat <- bx[[2]]
      viewDataTable(dat = dat, pageLength = 5)
    })
  })
  
  # create heatmap
  output$heatmap1 <- renderPlotly({
    if(input$heatmapsubmit1 == 0){
      return()
    }
    withProgress(session = session, message = "Plotting Data...", detail = "Takes a while...", min = 1, value = 10, max = 10,{
      isolate({
        hist = input$heatmapselectInput2
        sig = input$heatmapselectInput3
        type = input$heatmapselectInput4
        hm <<- create.hm(dataset, total.sub = total.sub, genelist.ct = genelist.ct, sig = sig, type = type, hist = hist)
        hm[[1]]
      })
    })
  })

  # get raw data expression data
  output$heatmaptable <- renderDataTable({
    if(input$heatmapsubmit1 == 0){
      return()
    }
    isolate({
      dat <- hm[[2]]
      viewDataTable(dat = dat, pageLength = 5)
    })
  })
  
  # gsva scores boxplot
  output$gsvaplot1 <- renderPlotly({
    if(input$gsvasubmit1 == 0){
      return()
    }
    isolate({
      hist <- input$gsvaselectInput2
      sig <- input$gsvaselectInput3
      type <- input$gsvaselectInput4
      if(type == "z-scores"){
        dat <- gsva.zscore
      } else {
        dat <- gsva.log2counts
      }
      gsva.obj <<- create.box.gsva(dat = dat, meta = x, genelist.ct = genelist.ct, hist = hist, sig = sig, type = type)
      gsva.obj[[1]]
    })
  })
  
  # get raw gsva data
  output$gsvatable <- renderDataTable({
    if(input$gsvasubmit1 == 0){
      return()
    }
    isolate({
      dat <- gsva.obj[[2]]
      viewDataTable(dat = dat, pageLength = 5)
    })
  })
    
})