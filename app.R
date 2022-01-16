#shiny MGP with custom list support
library(shiny)
library(shinydashboard)
library(ggplot2)
library(Morpho)
library(rgl)
library(geomorph)
library(plotly)
library(rmarkdown)
library(shapes)
library(ddsPLS)
library(Jovid)
library(shinycssloaders)
library(shinyjs)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(dplyr)
library(dbplyr)
library(future)
library(promises)
library(data.table)
future::plan("multisession")

#save(combined.markers, DO.go, giga.pca, mutant.db, mutant.lms, shape.mean, Y, DO.probs, file = "/data/MGP_data/offline_data.Rdata")
#save(combined.markers, giga.pca, mutant.db, mutant.lms, shape.mean, Y, file = "~/shiny/shinyapps/MGP/shiny_data2.Rdata")

#local dirs
mmusculusEnsembl <- loadDb(file="~/shiny/shinyapps/MGP/ensemble.sqlite")
load("~/shiny/shinyapps/MGP/shiny_data.Rdata")
# load("~/shiny/shinyapps/MGP/cached.results.Rdata")
# DO_probs_DB <- src_sqlite("~/shiny/shinyapps/MGP/MGP_genotypes.sqlite")
# # DO_probs_DB <- s3read_using(FUN = src_sqlite, object = "s3://mgpgenotypes/MGP_genotypes.sqlite") #src_sqlite("mgpgenotypes.s3.ca-central-1.amazonaws.com/MGP_genotypes.sqlite")

#docker deployment dirs
# mmusculusEnsembl <- loadDb(file="/srv/shiny-server/MGP/ensemble.sqlite")
# # # load("/srv/shiny-server/offline_data.Rdata")
# # # load("/srv/shiny-server/cached.results.Rdata")
# # # mmusculusEnsembl <- loadDb(file="/data/MGP_data/ensemble.sqlite")
# load("/srv/shiny-server/MGP/shiny_data.Rdata")
# load("/srv/shiny-server/MGP/cached.results.Rdata")
# DO_probs_DB <- src_sqlite("/srv/shiny-server/MGP/MGP_genotypes.sqlite")

#genopheno deployment dirs
# mmusculusEnsembl <- loadDb(file="/data/MGP_data/ensemble.sqlite")
# load("/data/MGP_data/shiny_data.Rdata")
# load("/data/MGP_data/cached.results.Rdata")
# DO_probs_DB <- src_sqlite("/data/MGP_data/MGP_genotypes.sqlite")


#get gene numbers for process field
  with.counts <- DO.go[,3:4]
  with.counts[,1] <- as.character(with.counts[,1])
  with.counts[,2] <- as.character(with.counts[,2])
  with.counts.char <- as.character(apply(with.counts, 1, FUN = function(x){ paste0(x[1], " (", x[2], " genes)")}))


#ui####
body <- dashboardBody(useShinyjs(),
                      tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
                      tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "mgp_style.css")),
                      tags$script(HTML('
                           $(document).ready(function() {
                           $("header").find("nav").append(\'<div id="pageHeader" class="titleStyle"></div>\');
                           })
                           ')),
        fluidRow(
           column(width = 12,
                  conditionalPanel(condition = "input.tabs1 == 'About this app'",
                                   tags$b("Pick a tab to see some options"))
           ),
           conditionalPanel(condition = "input.tabs1 == 'Process MGP'",
                            column(width = 4,
                                   # textInput("process", label = "Process", value = c("chondrocyte", "BMP", "fibroblast", "cohesin", "apoptosis")[sample(1:5,1)]),
                                   selectInput("variables2", label = "Process", multiple = T, choices = with.counts.char),
                                   sliderInput("lambda", "Sparsity parameter", min = 0, max = .15, value = .06, step = .01)
                                   # uiOutput('variables')
                            ),
                            column(width = 4,
                                   selectInput("facet", "Type of plot", choices = c("Simple","Messy, but informative", "Just the allele ranges", "Facet by founders")),
                                   # numericInput("lambda", "Sparsity parameter", value = .06, min = 0, max = 1)
                                   sliderInput("pls_axis", "PLS axis", min = 1, max = 4, value = 1, step = 1 )
                            ),
                            column(width = 4,
                                   numericInput("mag", "Magnification", value = 4, min = 1, max = 10),
                                   selectInput("mutant", "Make a comparison?", choices = c(" ","Whole genome", as.character(unique(mutant.db$Genotype))))
                            ),
                            column(width = 4,
                                   disabled(downloadButton("report", "Send me the results!")),
                                   br(),
                                   br(),
                                   actionButton("update_process", "Update process!")
                            ),
           ),
           conditionalPanel(condition = "input.tabs1 == 'Custom MGP'",
                              column(width = 4,
                                   textInput("custom_process", label = "Comma-separated genes", value = NULL, placeholder = c("Eg. Bmp7, Bmp2, Bmp4, Ankrd11")),
                                   sliderInput("pls_axis2", "PLS axis", min = 1, max = 4, value = 1, step = 1 )
                              ),
                              column(width = 4,
                                   selectInput("facet2", "Type of plot", choices = c("Simple","Messy, but informative", "Just the allele ranges", "Facet by founders")),
                                   # numericInput("lambda2", "Sparsity parameter", value = .06, min = 0, max = 1)
                                   sliderInput("lambda2", "Sparsity parameter", min = 0, max = .15, value = .06, step = .01)
                              ),
                              column(width = 4,
                                   numericInput("mag2", "Magnification", value = 4, min = 1, max = 10),
                                   selectInput("mutant2", "Make a comparison?", choices = c(" ","Whole genome", as.character(unique(mutant.db$Genotype))))
                              ),
                              column(width = 4,
                                   disabled(downloadButton("report2", "Send me the results!"))
                              ),
                              column(width = 4,
                                   actionButton("update_process2", "Update process!")
                              )
           ),
           conditionalPanel(condition = "input.tabs1 == 'Recent searches'"
           )
        ),
           column(width = 12, align = "center",
                    br(),
                    tabsetPanel(id = "tabs1",
                              tabPanel("Process MGP",
                                       br(),
                                       verbatimTextOutput("MGP_mutant_cor"),
                                       box(title = tags$b("Effect size plot"),
                                           status = "warning",
                                           solidHeader = F,
                                           width = 12,
                                           withSpinner(plotlyOutput("process_effect_size", width = "95%"), type = 6, color = "#a6192e")),
                                       box(title = tags$b("Phenotype plot"),
                                           status = "warning",
                                           solidHeader = F,
                                           width = 12,
                                           withSpinner(rglwidgetOutput("process_heatmap", width = "95%"), type = 6, color = "#a6192e"))
                              ),
                              tabPanel("Custom MGP",
                                       br(),
                                       verbatimTextOutput("value"),
                                       box(title = tags$b("Effect size plot"),
                                           status = "warning",
                                           solidHeader = F,
                                           width = 12,
                                           withSpinner(plotlyOutput("process_effect_size2", width = "95%"), type = 6, color = "#a6192e")),
                                       box(title = tags$b("Phenotype plot"),
                                           status = "warning",
                                           solidHeader = F,
                                           width = 12,
                                           withSpinner(rglwidgetOutput("process_heatmap2", width = "95%"), type = 6, color = "#a6192e"))
                              ),
                              tabPanel("Recent searches", br(),
                                       box(title = tags$b("Previous process correlations"),
                                           status = "warning",
                                           solidHeader = F,
                                           width = 12,
                                           actionButton("redraw_cors", "Draw correlation heatmap"),
                                           br(),
                                           br(),
                                           withSpinner(plotlyOutput("process_correlations", width = "95%"), type = 6, color = "#a6192e")),
                                       br(),
                                       tableOutput("recents")
                              ),
                              tabPanel("Mutant DB",
                                       br(),
                                       p("Here we provide our mutant database in full for you to explore before selecting options to make comparisons with."),
                                       br(),
                                       box(title = tags$b("Mutant DB"),
                                           status = "warning",
                                           solidHeader = F,
                                           width = 12,
                                           collapsed = T,
                                           collapsible = T,
                                           dataTableOutput("mutantdb"))
                              ),
                              tabPanel("About this app",
                                       br(),
                                       p("This app is designed for a developmentally-focused genomic analysis. Instead of thinking about how all the variants in the genome contribute to the shape of the face, we focus on a more direct question; what do variants associated with a process I’m interested in do to the face? The figure below shows what’s going on under the hood to make that happen. When you select a process, this program looks up which genes are known to be associated with that process. We then look at the genomic markers in our sample to find the nearest ones. Those close markers will represent that gene’s effect. We put those variants together with 54 3-dimensional measurements of craniofacial shape and estimate how they covary. "),
                                       img(src='process_mgp_diagram 2.png', height = 610 , width = 247, ALIGN = "center", VSPACE = 30),
                                       p("In the end you’ll see a plot that describes the relative strength of allelic effects on face shape in our sample. You can go on to compare mouse mutants from our collection of studies and collaborations over the years. The mutant database is always growing! We hope this program helps us better understand how genetic variation creates craniofacial variation. There are lots of contributing loci, and they often overlap in their effects. The dataset that’s driving all of this is a collection of 1145 mice from the Diversity outbred project at Jackson laboratories. Each of these mice have been genotyped at tens of thousands of loci, enabling us to do the marker selection with pretty good fidelity. If you’re interested in more details about this work, have a look at the ",  a("paper.", href="https://www.biorxiv.org/content/10.1101/2020.11.12.378513v2"), "If you'd like to see how this app was put together, here's a link to the", a(" Github repository.", href="https://github.com/J0vid/MGP_shiny"))
                              )
                  )
           )
  )

dbHeader <- dashboardHeader()
 dbHeader$children[[2]]$children <-  tags$a(href='../', tags$img(src="uc_logo.png",height='30',width='116'))

ui <- dashboardPage(title = "Process MGP",
                    dbHeader,
                    dashboardSidebar(disable = T),
                    body
      )

#server####
server <- function(input, output){
  
  #stringify custom process list####
  output$value <- renderText({str(input$custom_process)})
  
  #process filtering reactive####
  # outVar <- reactive({
  #   if(nchar(input$process) > 2){
  #   with.counts <- DO.go[grep(DO.go[,3], pattern = tolower(input$process)),3:4]
  #   with.counts[,1] <- as.character(with.counts[,1])
  #   with.counts[,2] <- as.character(with.counts[,2])
  #   process.ano <- with.counts[,1]
  #   
  #   with.counts.char <- as.character(apply(with.counts, 1, FUN = function(x){ paste0(x[1], " (", x[2], " genes)")}))
  #   
  #   return(list(process.ano, with.counts.char))
  # }
  # })
  
  #slow down reactivity of text input
  # process.list <- debounce(outVar, 300)
  
  # output$variables <- renderUI({
  #   selectInput('variables2', 'Process filter', process.list()[[2]], multiple = T)
  # })
  
  #main process reactive####
  process.svd <- eventReactive(input$update_process, {
    print(input$variables2)
    
    selection.vector <- paste(with.counts[with.counts.char %in% input$variables2, 1], collapse = ",") #process.list()[[1]][process.list()[[2]] %in% input$variables2]
    print(selection.vector)
    print(length(selection.vector))
    tmp.lambda <- input$lambda
    #debug: process.svd <- reactive({
    future::future({
      raw_api_res <- httr::GET(url = paste0("http://127.0.0.1:3636/MGP", "/mgp"),
                             query = list(GO.term = selection.vector, lambda = tmp.lambda),
                             encode = "json")
      print(jsonlite::fromJSON(httr::content(raw_api_res, "text")))
      jsonlite::fromJSON(httr::content(raw_api_res, "text"))
    }) 
  })
  
  old_processes <- reactiveValues("Previous analyses" = list())
  
  observeEvent(input$update_process,
               {
                 old_processes$x <- c(old_processes$x, input$variables2)
                 op.df <- as.data.frame(old_processes$x)
                 names(op.df) <- "Previous analyses"
                 old_processes$x <- op.df
                 enable("report")
               })
  
  output$recents <- renderTable({
    old_processes$x
  })
  
  #effect size plot Process MGP####
  output$process_effect_size <- renderPlotly({
    
    do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
    do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
  
    process.svd() %...>% `[[`("loadings") %...>% 
      {

    if(input$facet == "Simple"){ ggplot() +
      geom_bar(data = .,
               aes(x = gnames, y = gloadings),
               stat = "identity",
               width = .75,
               position=position_dodge()) +
      geom_point(data = .,
                 aes(x = gnames, y = gloadings, color = founders),
                 shape = "-",
                 size = 1) +
      scale_color_manual(values=do.colors,
                         guide = guide_legend(title = "Founder\nGenotype", override.aes = list(shape = rep(19, 8), size = 1))) +
      xlab("Gene") +
      ylab("Genetic marker loading") +
      theme(text = element_text(size=6),
            axis.text.x = element_text(angle = 75, hjust = 1),
            axis.title.x = element_text(margin = margin(t = 20)),
            axis.text = element_text(angle = 55, hjust = 1, size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8, face = "bold", hjust = .5),
            plot.background = element_rect(fill = rgb(245/255, 245/255, 245/255, .9), colour = rgb(245/255, 245/255, 245/255, .9)),
            legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9)))
    } else if(input$facet == "Messy, but informative"){
      ggplot(data = ., aes(x = gnames, y = gloadings, fill = founders)) +
        geom_bar(stat = "identity", width = .75, position=position_dodge()) +
        theme(text = element_text(size=6),
              axis.text.x = element_text(angle = 70, hjust = 1),
              axis.text.y = element_text(size = .2),
              axis.title.x = element_text(margin = margin(t = 20))) +
        scale_fill_manual(values=do.colors) +
        xlab("Gene") +
        ylab("Genetic marker loading")
    } else if(input$facet == "Just the allele ranges"){
      ggplot(data = ., aes(x = gnames, y = gloadings)) +
        geom_bar(stat = "identity", width = .75, position=position_dodge()) +
        theme(text = element_text(size=4),
              axis.text.x = element_text(angle = 75, hjust = 1, size = .5),
              axis.text.y = element_text(size = .4),
              axis.title.x = element_text(margin = margin(t = 20))) +
        scale_fill_manual(values=do.colors) +
        xlab("Gene") +
        ylab("Genetic marker loading")
    } else if(input$facet == "Facet by founders"){
      ggplot(data = ., aes(x = gnames, y = gloadings, fill = founders)) +
        geom_bar(stat = "identity", width = .75, position=position_dodge()) +
        theme(text = element_text(size=3),
              axis.text.x = element_text(angle = 70, hjust = 1),
              axis.title.x = element_text(margin = margin(t = 20))) +
        scale_fill_manual(values=do.colors) +
        xlab("") +
        ylab("Genetic coefficient") +
        facet_wrap( ~ founders, nrow = 3)
    } %...>% ggplotly(. +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                     axis.text.y = element_text(size = 8),
                     axis.title = element_text(size = 12, face = "bold"),
                     legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9)))) %...>% layout(margin = list(b = 100, l = 50)
                     )
    }
  })
  
  #process MGP phenotypic projection####
  #mutant registration reactive####
  mutant.comparison <- reactive({
    tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == input$mutant,])), 54, 3))$coords
    #debug: tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
    do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
    mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])
    return(list(do.mean, mutant.mean, tmp.mutant.registration))
  })

  output$process_heatmap <- renderRglwidget({
   
    do.mean <- matrix(colMeans(Y), nrow = 54, ncol = 3, byrow = T)

    par3d(zoom = .65)
    aspect3d("iso")

    #vectors from DO mean to mutant
    shape.warp <-  plot3d(shape.mean$mesh, col = adjustcolor("lightgrey", .3), alpha = .2, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
    bg3d(rgb(245/255,245/255,245/255, .9))


    if(input$mutant != " "){

      if(input$mutant != "Whole genome"){
        
        do.mean <- mutant.comparison()[[1]]
        mutant.mean <- mutant.comparison()[[2]]
        spheres3d(do.mean, radius = .002, color = 2)
        segments3d(x = (rbind(as.vector(do.mean[,1]), as.vector((mutant.mean + (mutant.mean - do.mean) * (input$mag - 1))[,1]))),
                   y = (rbind(as.vector(do.mean[,2]), as.vector((mutant.mean + (mutant.mean - do.mean) * (input$mag - 1))[,2]))),
                   z = (rbind(as.vector(do.mean[,3]), as.vector((mutant.mean + (mutant.mean - do.mean) * (input$mag - 1))[,3]))),
                   lwd = 3,
                   col = 2)
      }

      if(input$mutant == "Whole genome"){
      }
    }

    then(process.svd(),
         function(phenos) {
           spheres3d(phenos$pheno1, radius = .003, color = 1)
               segments3d(x = rbind(phenos$pheno1[,1], phenos$pheno2[,1]  + ((phenos$pheno2[,1] - phenos$pheno1[,1]) * (input$mag - 1))),
                        y = rbind(phenos$pheno1[,2], phenos$pheno2[,2]  + ((phenos$pheno2[,2] - phenos$pheno1[,2]) * (input$mag - 1))),
                        z = rbind(phenos$pheno1[,3], phenos$pheno2[,3]  + ((phenos$pheno2[,3] - phenos$pheno1[,3]) * (input$mag - 1))),
                        lwd = 3)
               
               rglwidget()
         })
  })
  
  #mutant db
  output$mutantdb <- renderDataTable({
    mutant.db[,2:10]
  }, options = list(
    autoWidth = FALSE, scrollX = TRUE #prevents overflow of datatable out of box
  ))

  
  # #custom MGP gene list code####
  custom.process.svd <- eventReactive(input$update_process2, {
    # selection.vector <- c("Bmp7, Bmp2, Bmp4, Ankrd11")
    tmp.list <- input$custom_process
    tmp.lambda <- input$lambda2
    future::future({
      raw_api_res <- httr::GET(url = paste0("http://127.0.0.1:3636/MGP", "/custom_mgp"),
                               query = list(genelist = tmp.list, lambda = tmp.lambda),
                               encode = "json")
      jsonlite::fromJSON(httr::content(raw_api_res, "text"))
    }) 

  })
  
  #Custom process effect size plot####
  output$process_effect_size2 <- renderPlotly({
    
    do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
    do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
    
    custom.process.svd() %...>% `[[`("loadings") %...>% 
      {
        
        if(input$facet == "Simple"){ ggplot() +
            geom_bar(data = .,
                     aes(x = gnames, y = gloadings),
                     stat = "identity",
                     width = .75,
                     position=position_dodge()) +
            geom_point(data = .,
                       aes(x = gnames, y = gloadings, color = founders),
                       shape = "-",
                       size = 1) +
            scale_color_manual(values=do.colors,
                               guide = guide_legend(title = "Founder\nGenotype", override.aes = list(shape = rep(19, 8), size = 1))) +
            xlab("Gene") +
            ylab("Genetic marker loading") +
            theme(text = element_text(size=6),
                  axis.text.x = element_text(angle = 75, hjust = 1),
                  axis.title.x = element_text(margin = margin(t = 20)),
                  axis.text = element_text(angle = 55, hjust = 1, size = 12),
                  axis.title = element_text(size = 12, face = "bold"),
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 8, face = "bold", hjust = .5),
                  plot.background = element_rect(fill = rgb(245/255, 245/255, 245/255, .9), colour = rgb(245/255, 245/255, 245/255, .9)),
                  legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9)))
        } else if(input$facet == "Messy, but informative"){
          ggplot(data = ., aes(x = gnames, y = gloadings, fill = founders)) +
            geom_bar(stat = "identity", width = .75, position=position_dodge()) +
            theme(text = element_text(size=6),
                  axis.text.x = element_text(angle = 70, hjust = 1),
                  axis.text.y = element_text(size = .2),
                  axis.title.x = element_text(margin = margin(t = 20))) +
            scale_fill_manual(values=do.colors) +
            xlab("Gene") +
            ylab("Genetic marker loading")
        } else if(input$facet == "Just the allele ranges"){
          ggplot(data = ., aes(x = gnames, y = gloadings)) +
            geom_bar(stat = "identity", width = .75, position=position_dodge()) +
            theme(text = element_text(size=4),
                  axis.text.x = element_text(angle = 75, hjust = 1, size = .5),
                  axis.text.y = element_text(size = .4),
                  axis.title.x = element_text(margin = margin(t = 20))) +
            scale_fill_manual(values=do.colors) +
            xlab("Gene") +
            ylab("Genetic marker loading")
        } else if(input$facet == "Facet by founders"){
          ggplot(data = ., aes(x = gnames, y = gloadings, fill = founders)) +
            geom_bar(stat = "identity", width = .75, position=position_dodge()) +
            theme(text = element_text(size=3),
                  axis.text.x = element_text(angle = 70, hjust = 1),
                  axis.title.x = element_text(margin = margin(t = 20))) +
            scale_fill_manual(values=do.colors) +
            xlab("") +
            ylab("Genetic coefficient") +
            facet_wrap( ~ founders, nrow = 3)
        } %...>% ggplotly(. +
                            theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                                  axis.text.y = element_text(size = 8),
                                  axis.title = element_text(size = 12, face = "bold"),
                                  legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9)))) %...>% layout(margin = list(b = 100, l = 50)
                                  )
      }
  })

  #custom mutant registration reactive####
  custom.mutant.comparison <- reactive({

    tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == input$mutant2,])), 54, 3))$coords
    #debug: tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
    do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
    mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])
    return(list(do.mean, mutant.mean, tmp.mutant.registration))
  })

  #render custom process phenotype####
  output$process_heatmap2 <- renderRglwidget({
    
    do.mean <- matrix(colMeans(Y), nrow = 54, ncol = 3, byrow = T)
    
    par3d(zoom = .65)
    aspect3d("iso")
    
    #vectors from DO mean to mutant
    shape.warp <-  plot3d(shape.mean$mesh, col = adjustcolor("lightgrey", .3), alpha = .2, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
    bg3d(rgb(245/255,245/255,245/255, .9))
    
    if(input$mutant2 != " "){
      
      if(input$mutant2 != "Whole genome"){
        do.mean <- custom.mutant.comparison()[[1]]
        mutant.mean <- custom.mutant.comparison()[[2]]
        spheres3d(do.mean, radius = .002, color = 2)
        segments3d(x = (rbind(as.vector(do.mean[,1]), as.vector((mutant.mean + (mutant.mean - do.mean) * (input$mag2 - 1))[,1]))),
                   y = (rbind(as.vector(do.mean[,2]), as.vector((mutant.mean + (mutant.mean - do.mean) * (input$mag2 - 1))[,2]))),
                   z = (rbind(as.vector(do.mean[,3]), as.vector((mutant.mean + (mutant.mean - do.mean) * (input$mag2 - 1))[,3]))),
                   lwd = 3,
                   col = 2)
      }
      
      if(input$mutant2 == "Whole genome"){
      }
    }
    
    then(custom.process.svd(),
         function(phenos) {
           spheres3d(phenos$pheno1, radius = .003, color = 1)
           segments3d(x = rbind(phenos$pheno1[,1], phenos$pheno2[,1]  + ((phenos$pheno2[,1] - phenos$pheno1[,1]) * (input$mag2 - 1))),
                      y = rbind(phenos$pheno1[,2], phenos$pheno2[,2]  + ((phenos$pheno2[,2] - phenos$pheno1[,2]) * (input$mag2 - 1))),
                      z = rbind(phenos$pheno1[,3], phenos$pheno2[,3]  + ((phenos$pheno2[,3] - phenos$pheno1[,3]) * (input$mag2 - 1))),
                      lwd = 3)
           
           rglwidget()
         })

  })
  
  #dynamic title needs to have a reactive expression that looks for changes in process and mutant!
  title.change <- reactive({
    paste(input$update_process , input$mutant)
  })
  
  #dynamic title with correlation b/t mutant effect and current MGP vector####
  title.react <- eventReactive(title.change(), {
    first.title <- "sup"#paste0(paste(DO.go[grep(DO.go[,3], pattern = tolower(input$variables2)), 3], collapse = ", ", sep = ", "), " MGP")#paste0(paste(process.list()[[1]][process.list()[[2]] %in% input$variables2], collapse = ", ", sep = ", "), " MGP")
    tmp.mutant <- input$mutant
    
    
    if(tmp.mutant == " "){my.title <- first.title}
    
    if(tmp.mutant != " "){
      #change title name to selected process
      if(tmp.mutant != "Whole genome"){
        do.mean <- mutant.comparison()[[1]]
        
        my.cor <- then(process.svd(),
                               function(phenos) {
                                 
                                 mutant.cor <- cor(as.numeric(phenos$pheno_loadings), as.numeric(manova(two.d.array(mutant.comparison()[[3]]) ~ c(rep(0, nrow(Y)), rep(1, sum(mutant.db$Genotype == tmp.mutant))))$coef[2,]))
                                 # cor(as.numeric(phenos$pheno1 - phenos$pheno2), manova(two.d.array(mutant.comparison()[[3]]) ~ c(rep(0, nrow(Y)), rep(1, sum(mutant.db$Genotype == tmp.mutant))))$coef[2,])
                                 print(mutant.cor)
                                 
                                 return(paste0("Correlation between ", paste(first.title, collapse = ", ", sep = ", "), " and ", tmp.mutant, " mutant: ", round(mutant.cor, digits = 3)))
                               }
                              )
        print(my.cor)
        my.title <- my.cor
       
        # MGP.mutant.cor <- cor(process.svd()[[1]]$mod$v[,1], prcomp(rbind(as.numeric(do.mean), as.numeric(mutant.comparison()[[2]])))$rotation[,1])
        
      }
      if(tmp.mutant == "Whole genome"){
        MGP.mutant.cor <- cor(process.svd()[[1]]$mod$V_super[,1], wgbetas)
        my.title <- paste0("Correlation between ", paste(first.title, collapse = ", ", sep = ", "), " and whole genome MGP: ", round(MGP.mutant.cor, digits = 3))
      }
    }

  })
  
  output$title <- renderText({
    print(title.react())
  })
  
  #meta-analysis of user selections####
  #convert cached.params from GO to process names
  correlation_draw <- eventReactive(input$redraw_cors, {  
    raw_api_res <- httr::POST(url = paste0("http://127.0.0.1:3636/MGP", "/mgp_cor"), encode = "json")
    cormat <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))
    cormat$cormat <- abs(cormat$cormat)
    colnames(cormat$cormat) <- cormat$process.names
    rownames(cormat$cormat) <- cormat$process.names
    heatmaply::heatmaply(cormat$cormat, symm = T, dendrogram = "none") 
    })
  
  output$process_correlations <- renderPlotly({
    correlation_draw()
  })

  #report generator####
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = paste0("Process_MGP_report.html"),
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("/srv/shiny-server/MGP/report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      parameters <- list(variables2 = input$variables2, lambda = input$lambda, mutant = input$mutant, mag = input$mag, pls_axis = input$pls_axis)
      # parameters <- list(variables2 = "chondrocyte differentiation", lambda = .06, mutant = "Alk6", mag = 4)
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = parameters#,
                        # envir = new.env(parent = globalenv())
      )
    }
  )
  
  #thanks to https://stackoverflow.com/questions/45176030/add-text-on-right-of-shinydashboard-header
  observe({
  #place title to the right of the logo
  shinyjs::html("pageHeader", title.react())
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)

