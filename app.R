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


#save(combined.markers, DO.go, giga.pca, mutant.db, mutant.lms, shape.mean, Y, DO.probs, file = "/data/MGP_data/offline_data.Rdata")
#save(combined.markers, giga.pca, mutant.db, mutant.lms, shape.mean, Y, file = "~/shiny/shinyapps/MGP/shiny_data2.Rdata")

# #local dirs
# mmusculusEnsembl <- loadDb(file="~/shiny/shinyapps/MGP/ensemble.sqlite")
# load("~/shiny/shinyapps/MGP/shiny_data.Rdata")
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
# # # mmusculusEnsembl <- loadDb(file="/data/MGP_data/ensemble.sqlite")
# load("/data/MGP_data/shiny_data.Rdata")
# load("/data/MGP_data/cached.results.Rdata")
# DO_probs_DB <- src_sqlite("/data/MGP_data/MGP_genotypes.sqlite")


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
           conditionalPanel(condition = "input.tabs1 == 'Process MGP'",
                              column(width = 4,
                                   textInput("process", label = "Process", value = c("chondrocyte", "BMP", "fibroblast", "cohesin", "apoptosis")[sample(1:5,1)]),
                                   uiOutput('variables')
                              ),
                              column(width = 4,
                                   selectInput("facet", "Type of plot", choices = c("Simple","Messy, but informative", "Just the allele ranges", "Facet by founders")),
                                   # numericInput("lambda", "Sparsity parameter", value = .06, min = 0, max = 1)
                                   sliderInput("lambda", "Sparsity parameter", min = 0, max = .15, value = .06, step = .01),
                                   sliderInput("pls_axis", "PLS axis", min = 1, max = 4, value = 1, step = 1 )
                              ),
                              column(width = 4,
                                   numericInput("mag", "Magnification", value = 4, min = 1, max = 10),
                                   selectInput("mutant", "Make a comparison?", choices = c(" ","Whole genome", as.character(unique(mutant.db$Genotype))))
                              ),
                             column(width = 4,
                                   disabled(downloadButton("report", "Send me the results!"))
                              ),
                              column(width = 4,
                                   br(),
                                   actionButton("update_process", "Update process!")
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
                                           withSpinner(plotlyOutput("process_correlations", width = "95%"), type = 6, color = "#a6192e")),
                                       br(),
                                       tableOutput("recents")
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


server <- function(input, output){
  
  #stringify custom process list####
  output$value <- renderText({str(input$custom_process)})
  
  #process filtering reactive####
  outVar <- reactive({
    with.counts <- DO.go[grep(DO.go[,3], pattern = tolower(input$process)),3:4]
    with.counts[,1] <- as.character(with.counts[,1])
    with.counts[,2] <- as.character(with.counts[,2])
    process.ano <- with.counts[,1]
    
    with.counts.char <- as.character(apply(with.counts, 1, FUN = function(x){ paste0(x[1], " (", x[2], " genes)")}))
    
    return(list(process.ano, with.counts.char))
  })
  
  #slow down reactivity of text input
  process.list <- debounce(outVar, 1500)
  
  output$variables <- renderUI({
    selectInput('variables2', 'Process filter', process.list()[[2]], multiple = T)
  })
  
  #main process reactive####
  process.svd <- eventReactive(input$update_process, {
    #debug: process.svd <- reactive({
    #selection.vector <- c('chondrocyte differentiation')
    
    #instead of process.ano doing pattern matching, we need to use it to match GO terms precisely
    selection.vector <- process.list()[[1]][process.list()[[2]] %in% input$variables2]
    
    process.ano <- NULL
    for(i in 1: length(selection.vector)) process.ano <- c(process.ano, as.character(DO.go[DO.go[,3] == selection.vector[i], 2]))
    
    #check the cache to see if the analysis has been called before. If so, just give the cached result
    if(sum(cached.params[,1] == paste(process.ano, collapse = '|') & cached.params[,2] == input$lambda & input$pls_axis == 1) == 1){
      results <- cached.results[[which(cached.params[,1] == paste(process.ano, collapse = '|') & cached.params[,2] == input$lambda)]][[1]]
    } else {
      
      #offline method for getting gene metadata
      #pull gene names from process.ano (go terms)
      coi <- c("ENSEMBL", "SYMBOL")
      go2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = process.ano, columns = coi, keytype = "GO")[,-2:-3]))
      
      coi2 <- c("TXCHROM", "TXSTART", "TXEND")
      symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = go2symbol[,2], columns = coi2, keytype="GENEID")
      transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
      
      #symbol, chr, start, end
      chr_name <- rep(NA,  nrow(go2symbol))
      gene.start <- rep(NA,  nrow(go2symbol))
      gene.end <- rep(NA,  nrow(go2symbol))
      
      for(i in 1:length(unique(symbol2info$GENEID))){
        
          tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
        
          chr_name[i] <- tmp.transcript$TXCHROM
          gene.start[i] <- tmp.transcript$TXSTART
          gene.end[i] <- tmp.transcript$TXEND
        
      }
      
      seq.info <- data.frame(mgi_symbol = go2symbol$SYMBOL, chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
      seq.info[,2] <- as.character(seq.info[,2])
      seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6  
      
      #get rid of weird chromosome names
      if(length(grep(seq.info$chromosome_name, pattern = "CHR")) > 0) seq.info <- seq.info[-grep(seq.info$chromosome_name, pattern = "CHR"),]
      
      gene.names <- seq.info[,1]
      seq.indexes <- matrix(NA, ncol = 3, nrow = dim(seq.info)[1])
      
      for(j in 1 : dim(seq.info)[1]){
          tmp.indexes <-  combined.markers[which(combined.markers$chr == seq.info[j,2] & combined.markers$Mbp_mm10 > mean(as.numeric(seq.info[j,3:4])) - 2 & combined.markers$Mbp_mm10 < mean(as.numeric(seq.info[j,3:4])) + 2), c(1,3)]
          #for each gene, select the marker closest to the middle of the gene
          seq.indexes[j,] <- as.matrix(cbind(seq.info[j,1],tmp.indexes[which.min(abs(tmp.indexes[,2] - mean(as.numeric(seq.info[j,3:4])))),]))
      }
      
      #put together selected genotype data
      probs.rows <- matrix(NA, nrow = nrow(Y), ncol = nrow(seq.indexes) * 8)
      probrowseq <- seq(1, ncol(probs.rows) + 8, by = 8)
      for(i in 1:nrow(seq.indexes)) probs.rows[, probrowseq[i]:(probrowseq[i+1] - 1) ] <- as.matrix(collect(tbl(DO_probs_DB, seq.indexes[i,2])))
      
      #fit lasso pls2B
      pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = input$pls_axis, lambda = input$lambda)
      
      do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
      pathway.loadings <- data.frame(gloadings = pls.svd$mod$u[[1]][, input$pls_axis], gnames = as.character(rep(seq.info[,1], each = 8)), founders = rep(do.names, nrow(seq.info)))
      
      #cache to pls list for new analyses if they've asked for the first axis only
      if(input$pls_axis == 1){
        results <- list(pls.svd, gene.names, seq.info, probs.rows, pathway.loadings)
        cached.results[[length(cached.results) + 1]] <<- list(results)
        cached.params <<- rbind(cached.params, c(paste(process.ano, collapse = '|'), input$lambda))
        save(cached.results, cached.params, file = "~/shiny/shinyapps/MGP/cached.results.Rdata")
      
      } else { results <- list(pls.svd, gene.names, seq.info, probs.rows, pathway.loadings)}
    }#end caching if
  return(results)
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
    tmp.reactive <- process.svd()
    reactive.svd <- tmp.reactive[[1]]$mod$u[[1]]
    gene.names <- tmp.reactive[[2]]
    seq.info <- tmp.reactive[[3]]
    
    #now we should be able to take pls.svd directly and maybe label them by founder in a new column, then barplot by family, by gene
    do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
    do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
    print(dim(reactive.svd))
    pathway.loadings <- data.frame(gloadings = reactive.svd[,input$pls_axis], gnames = as.character(rep(seq.info[,1], each = 8)), founders = rep(do.names, nrow(seq.info)))

    p <-  ggplot() +
      geom_bar(data = pathway.loadings, 
               aes(x = gnames, y = gloadings), 
               stat = "identity", 
               width = .75, 
               position=position_dodge()) +
      geom_point(data = pathway.loadings,
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
    
    
    if(input$facet == "Messy, but informative"){
      p <- ggplot(data = pathway.loadings, aes(x = gnames, y = gloadings, fill = founders)) +
        geom_bar(stat = "identity", width = .75, position=position_dodge()) +
        theme(text = element_text(size=6),
              axis.text.x = element_text(angle = 70, hjust = 1),
              axis.text.y = element_text(size = .2),
              axis.title.x = element_text(margin = margin(t = 20))) +
        scale_fill_manual(values=do.colors) +
        xlab("Gene") +
        ylab("Genetic marker loading")
    }
    
    if(input$facet == "Just the allele ranges"){
      p <- ggplot(data = pathway.loadings, aes(x = gnames, y = gloadings)) +
        geom_bar(stat = "identity", width = .75, position=position_dodge()) +
        theme(text = element_text(size=4),
              axis.text.x = element_text(angle = 75, hjust = 1, size = .5),
              axis.text.y = element_text(size = .4),
              axis.title.x = element_text(margin = margin(t = 20))) +
        scale_fill_manual(values=do.colors) +
        xlab("Gene") +
        ylab("Genetic marker loading")
    }
    
    
    if(input$facet == "Facet by founders"){
      
      p <- ggplot(data = pathway.loadings, aes(x = gnames, y = gloadings, fill = founders)) +
        geom_bar(stat = "identity", width = .75, position=position_dodge()) +
        theme(text = element_text(size=3),
              axis.text.x = element_text(angle = 70, hjust = 1),
              axis.title.x = element_text(margin = margin(t = 20))) +
        scale_fill_manual(values=do.colors) +
        xlab("") +
        ylab("Genetic coefficient") +
        facet_wrap( ~ founders, nrow = 3)
    }
    
    ggplotly(p +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                     axis.text.y = element_text(size = 8),
                     axis.title = element_text(size = 12, face = "bold")),
             legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9))) %>% layout(
                       margin = list(b = 100, l = 50) # to fully display the x and y axis labels
                     )
  })
  
  #process MGP phenotypic projection####
  pheno.proj.process <- reactive({
    tmp.reactive <- process.svd()
    probs.rows <- tmp.reactive[[4]]
    
    snp.dim <- input$pls_axis
    
    #calculate projection
    proj.coords.a1 <- row2array3d(predict(tmp.reactive[[1]], probs.rows[c(which.min(tmp.reactive[[1]]$mod$ts[[snp.dim]][, 1]), which.max(tmp.reactive[[1]]$mod$ts[[snp.dim]][,1])),])$y, Nlandmarks = 54)
    proj.coords.a2 <- proj.coords.a1[,,2]
    proj.coords.a1 <- proj.coords.a1[,,1]
    
    return(list(proj.coords.a1, proj.coords.a2))
  })
  
  #mutant registration reactive####
  mutant.comparison <- reactive({
    
    tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == input$mutant,])), 54, 3))$coords
    #debug: tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
    do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
    mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])
    return(list(do.mean, mutant.mean, tmp.mutant.registration))
  })
  
  output$process_heatmap <- renderRglwidget({
    print(input$mutant)
    do.mean <- matrix(colMeans(Y), nrow = 54, ncol = 3, byrow = T)

    par3d(zoom = .65)
    aspect3d("iso")
    
    #vectors from DO mean to mutant
    shape.warp <-  plot3d(shape.mean$mesh, col = adjustcolor("lightgrey", .3), alpha = .2, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
    spheres3d(pheno.proj.process()[[1]], radius = .003, color = 1)
    bg3d(rgb(245/255,245/255,245/255, .9))
    
    
    if(input$mutant != " "){
      
      if(input$mutant != "Whole genome"){
        print(mutant.comparison())
        do.mean <- mutant.comparison()[[1]]
        mutant.mean <- mutant.comparison()[[2]]
        spheres3d(do.mean, radius = .002, color = 2)
        segments3d(x = (rbind(as.vector(do.mean[,1]), as.vector((mutant.mean + (mutant.mean - do.mean) * (input$mag - 1))[,1]))),
                   y = (rbind(as.vector(do.mean[,2]), as.vector((mutant.mean + (mutant.mean - do.mean) * (input$mag - 1))[,2]))),
                   z = (rbind(as.vector(do.mean[,3]), as.vector((mutant.mean + (mutant.mean - do.mean) * (input$mag - 1))[,3]))),
                   lwd = 3,
                   col = 2)
        # for(i in 1:54) arrow3d(do.mean[i,] - (mutant.mean[i,] - do.mean[i,]), mutant.mean[i,], type = "lines", col = "red", barblen = 0.005, lwd = 2)
      }
      
      if(input$mutant == "Whole genome"){
        #for(i in 1:54) arrow3d(proj.pca1[i,], proj.pca2[i,], type = "lines", col = "red", barblen = 0, lwd = 2)
      }
    }
    
    # for(i in 1:54) arrow3d(pheno.proj.process()[[2]][i,], pheno.proj.process()[[1]][i,] + (pheno.proj.process()[[1]][i,] - pheno.proj.process()[[2]][i,]) * (input$mag - 1), type = "lines", col = "black", barblen = 0.005, lwd = 3)
    segments3d(x = (rbind(as.vector(pheno.proj.process()[[1]][,1]), as.vector((pheno.proj.process()[[2]] + (pheno.proj.process()[[2]] - pheno.proj.process()[[1]]) * (input$mag - 1))[,1]))),
               y = (rbind(as.vector(pheno.proj.process()[[1]][,2]), as.vector((pheno.proj.process()[[2]] + (pheno.proj.process()[[2]] - pheno.proj.process()[[1]]) * (input$mag - 1))[,2]))),
               z = (rbind(as.vector(pheno.proj.process()[[1]][,3]), as.vector((pheno.proj.process()[[2]] + (pheno.proj.process()[[2]] - pheno.proj.process()[[1]]) * (input$mag - 1))[,3]))),
               lwd = 3)
    
    rglwidget()
    
  })
  
  #custom MGP gene list code####
  custom.process.svd <- eventReactive(input$update_process2, {
    # process.svd <- reactive({
    #instead of process.ano doing pattern matching, we need to use it to match GO terms precisely
    # selection.vector <- c("Bmp7, Bmp2, Bmp4, Ankrd11")
    selection.vector <- as.character(strsplit(input$custom_process, ", ")[[1]])
    print(str(selection.vector))

    coi <- c("ENSEMBL", "SYMBOL")
    gene2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = selection.vector, columns = coi, keytype = "SYMBOL")))
      
    coi2 <- c("TXCHROM", "TXSTART", "TXEND")
    symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = gene2symbol[,2], columns = coi2, keytype="GENEID")
    transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
      
    chr_name <- rep(NA,  length(unique(symbol2info$GENEID)))
    gene.start <- rep(NA,  length(unique(symbol2info$GENEID)))
    gene.end <- rep(NA,  length(unique(symbol2info$GENEID)))
      
    for(i in 1:length(unique(symbol2info$GENEID))){
        tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
        chr_name[i] <- tmp.transcript$TXCHROM
        gene.start[i] <- tmp.transcript$TXSTART
        gene.end[i] <- tmp.transcript$TXEND
      }
      
    seq.info <- data.frame(mgi_symbol = gene2symbol[match(unique(symbol2info$GENEID), gene2symbol$ENSEMBL),1], chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
    seq.info[,2] <- as.character(seq.info[,2])
    seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6  
    gene.names <- seq.info[,1]
      
    if(length(grep(seq.info$chromosome_name, pattern = "CHR")) > 0) seq.info <- seq.info[-grep(seq.info$chromosome_name, pattern = "CHR"),]
      
      seq.indexes <- matrix(NA, ncol = 3, nrow = dim(seq.info)[1])
      #we have seq.info which gives us a gene name and its location on the chromosome
      
      for(j in 1 : dim(seq.info)[1]){
        tmp.indexes <-  combined.markers[which(combined.markers$chr == seq.info[j,2] & combined.markers$Mbp_mm10 > mean(as.numeric(seq.info[j,3:4])) - 2 & combined.markers$Mbp_mm10 < mean(as.numeric(seq.info[j,3:4])) + 2), c(1,3)]
        #for each gene, select the marker closest to the middle of the gene
        seq.indexes[j,] <- as.matrix(cbind(seq.info[j,1],tmp.indexes[which.min(abs(tmp.indexes[,2] - mean(as.numeric(seq.info[j,3:4])))),]))
    }

    #put together selected genotype data
    probs.rows <- matrix(NA, nrow = nrow(Y), ncol = nrow(seq.indexes) * 8)
    probrowseq <- seq(1, ncol(probs.rows) + 8, by = 8)
    for(i in 1:nrow(seq.indexes)) probs.rows[, probrowseq[i]:(probrowseq[i+1] - 1) ] <- as.matrix(collect(tbl(DO_probs_DB, seq.indexes[i,2])))
    pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = input$pls_axis2, lambda = input$lambda2)
    results <- list(pls.svd, gene.names, seq.info, probs.rows)
    print(seq.info)
    return(results)
  })
  
  #Custom process effect size plot####
  output$process_effect_size2 <- renderPlotly({
    tmp.reactive <- custom.process.svd()
    reactive.svd <- tmp.reactive[[1]]$mod$u[[1]]
    gene.names <- tmp.reactive[[2]]
    seq.info <- tmp.reactive[[3]]
    
    #now we should be able to take pls.svd directly and maybe label them by founder in a new column, then barplot by family, by gene
    do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
    do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
    
    pathway.loadings <- data.frame(gloadings = reactive.svd[,input$pls_axis2], gnames = as.character(rep(seq.info[,1], each = 8)), founders = rep(do.names, nrow(seq.info)))

    p <-  ggplot() +
      geom_bar(data = pathway.loadings, 
               aes(x = gnames, y = gloadings), 
               stat = "identity", 
               width = .75, 
               position=position_dodge()) +
      geom_point(data = pathway.loadings,
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
    
    
    if(input$facet2 == "Messy, but informative"){
      p <- ggplot(data = pathway.loadings, aes(x = gnames, y = gloadings, fill = founders)) +
        geom_bar(stat = "identity", width = .75, position=position_dodge()) +
        theme(text = element_text(size=6),
              axis.text.x = element_text(angle = 70, hjust = 1),
              axis.text.y = element_text(size = .2),
              axis.title.x = element_text(margin = margin(t = 20))) +
        scale_fill_manual(values=do.colors) +
        xlab("Gene") +
        ylab("Genetic marker loading")
    }
    
    if(input$facet2 == "Just the allele ranges"){
      p <- ggplot(data = pathway.loadings, aes(x = gnames, y = gloadings)) +
        geom_bar(stat = "identity", width = .75, position=position_dodge()) +
        theme(text = element_text(size=4),
              axis.text.x = element_text(angle = 75, hjust = 1, size = .5),
              axis.text.y = element_text(size = .4),
              axis.title.x = element_text(margin = margin(t = 20))) +
        scale_fill_manual(values=do.colors) +
        xlab("Gene") +
        ylab("Genetic marker loading")
    }
    
    
    if(input$facet2 == "Facet by founders"){
      
      p <- ggplot(data = pathway.loadings, aes(x = gnames, y = gloadings, fill = founders)) +
        geom_bar(stat = "identity", width = .75, position=position_dodge()) +
        theme(text = element_text(size=3),
              axis.text.x = element_text(angle = 70, hjust = 1),
              axis.title.x = element_text(margin = margin(t = 20))) +
        scale_fill_manual(values=do.colors) +
        xlab("") +
        ylab("Genetic coefficient") +
        facet_wrap( ~ founders, nrow = 3)
    }
    
    ggplotly(p +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                     axis.text.y = element_text(size = 8),
                     axis.title = element_text(size = 12, face = "bold")),
             legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9))) %>% layout(
               margin = list(b = 100, l = 50) # to fully display the x and y axis labels
             )
  })
  
  #custom process phenotypic projection####
  custom.pheno.proj.process <- reactive({
    tmp.reactive <- custom.process.svd()
    probs.rows <- tmp.reactive[[4]]
    snp.dim <- input$pls_axis2
    #calculate projection
    proj.coords.a1 <- row2array3d(predict(tmp.reactive[[1]], probs.rows[c(which.min(tmp.reactive[[1]]$mod$ts[[snp.dim]][,1]), which.max(tmp.reactive[[1]]$mod$ts[[snp.dim]][,1])),])$y, Nlandmarks = 54)
    proj.coords.a2 <- proj.coords.a1[,,2]
    proj.coords.a1 <- proj.coords.a1[,,1]
    
    return(list(proj.coords.a1, proj.coords.a2))
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
    spheres3d(custom.pheno.proj.process()[[1]], radius = .003, color = 1)
    bg3d(rgb(245/255,245/255,245/255, .9))
    
    # shape.warp <- plot3d(do.mean, typ = "s", radius = .001, col = adjustcolor("red", .3), aspect = "iso")
    # shade3d(shape.mean$mesh)
    # rglwidget()
    
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
        #for(i in 1:54) arrow3d(proj.pca1[i,], proj.pca2[i,], type = "lines", col = "red", barblen = 0, lwd = 2)
      }
    }
    
    # for(i in 1:54) arrow3d(custom.custom.pheno.proj.process()[[2]][i,], custom.custom.pheno.proj.process()[[1]][i,] + (custom.custom.pheno.proj.process()[[1]][i,] - custom.custom.pheno.proj.process()[[2]][i,]) * (input$mag - 1), type = "lines", col = "black", barblen = 0.005, lwd = 3)
    segments3d(x = (rbind(as.vector(custom.pheno.proj.process()[[1]][,1]), as.vector((custom.pheno.proj.process()[[2]] + (custom.pheno.proj.process()[[2]] - custom.pheno.proj.process()[[1]]) * (input$mag2 - 1))[,1]))),
               y = (rbind(as.vector(custom.pheno.proj.process()[[1]][,2]), as.vector((custom.pheno.proj.process()[[2]] + (custom.pheno.proj.process()[[2]] - custom.pheno.proj.process()[[1]]) * (input$mag2 - 1))[,2]))),
               z = (rbind(as.vector(custom.pheno.proj.process()[[1]][,3]), as.vector((custom.pheno.proj.process()[[2]] + (custom.pheno.proj.process()[[2]] - custom.pheno.proj.process()[[1]]) * (input$mag2 - 1))[,3]))),
               lwd = 3)
    
    rglwidget()
    
  })
  
  #dynamic title needs to have a reactive expression that looks for changes in process and mutant!
  title.change <- reactive({
    paste(input$update_process , input$mutant)
  })
  
  #dynamic title with correlation b/t mutant effect and current MGP vector####
  title.react <- eventReactive(title.change(), {
    first.title <- paste0(paste(process.list()[[1]][process.list()[[2]] %in% input$variables2], collapse = ", ", sep = ", "), " MGP")
    
    if(input$mutant == " "){my.title <- first.title}
    
    if(input$mutant != " "){
      #change title name to selected process
      if(input$mutant != "Whole genome"){
        do.mean <- mutant.comparison()[[1]]
        MGP.mutant.cor <- cor(process.svd()[[1]]$mod$V_super[,input$pls_axis], manova(two.d.array(mutant.comparison()[[3]]) ~ c(rep(0, nrow(Y)), rep(1, sum(mutant.db$Genotype == input$mutant))))$coef[2,])
        # MGP.mutant.cor <- cor(process.svd()[[1]]$mod$v[,1], prcomp(rbind(as.numeric(do.mean), as.numeric(mutant.comparison()[[2]])))$rotation[,1])
        my.title <- paste0("Correlation between ", paste(first.title, collapse = ", ", sep = ", "), " and ", input$mutant, " mutant: ", round(MGP.mutant.cor, digits = 3))
      }
      if(input$mutant == "Whole genome"){
        MGP.mutant.cor <- cor(process.svd()[[1]]$mod$V_super[,1], wgbetas)
        my.title <- paste0("Correlation between ", paste(first.title, collapse = ", ", sep = ", "), " and whole genome MGP: ", round(MGP.mutant.cor, digits = 3))
      }
    }
    
    return(my.title)
  })
  
  output$title <- renderText({
    print(title.react())
    
  })
  
  
  #meta-analysis of user selections####
  #convert cached.params from GO to process names
  output$process_correlations <- renderPlotly({
    chosen.processes <- sample(1:nrow(cached.params), 10)
    tmp.proc.names <- rep(NA, 10)
    tmp.cor <- matrix(NA, nrow = 162, ncol = 10)
    for(i in 1:10){
      tmp.cor[,i] <- cached.results[[chosen.processes[i]]][[1]][[1]]$mod$V_super[,1]
      tmp.split <- strsplit(cached.params[chosen.processes[i],1], split = "\\|")[[1]]
      name.buffer <- NULL
      for(j in 1 : length(tmp.split)) {
        name.buffer <- c(name.buffer, as.character(DO.go[DO.go[,2] == tmp.split[j],3]))
      }
      tmp.proc.names[i] <- paste(name.buffer, collapse = " + ") 
    }
    heatmaply::heatmaply(cor(tmp.cor), labRow = tmp.proc.names, symm = T)
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

