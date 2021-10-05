#
# This is a Plumber API. You can run the API by clicking
# the 'Run API' button above.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#

library(plumber)
future::plan("multicore")

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

#local dirs
mmusculusEnsembl <- loadDb(file="~/shiny/shinyapps/MGP/ensemble.sqlite")
load("~/shiny/shinyapps/MGP/shiny_data.Rdata")
load("~/shiny/shinyapps/MGP/cached.results.Rdata")
DO_probs_DB <- src_sqlite("~/shiny/shinyapps/MGP/MGP_genotypes.sqlite")

#mgp function####

mgp <- function(GO.term = "chondrocyte differentiation", mutant = NULL, Y, cv = F, lambda = .06){
  
  selection.vector <- c(GO.term)
  # selection.vector <- process.list()[[1]][process.list()[[2]] %in% input$variables2]
  
  process.ano <- NULL
  for(i in 1: length(selection.vector)) process.ano <- c(process.ano, as.character(DO.go[DO.go[,3] == selection.vector[i], 2]))
  
  coi <- c("ENSEMBL", "SYMBOL")
  go2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = process.ano, columns = coi, keytype = "GO")[,-2:-3]))
  coi2 <- c("TXCHROM", "TXSTART", "TXEND")
  symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = go2symbol[,2], columns = coi2, keytype="GENEID")
  
  transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
  
  #symbol, chr, start, end
  chr_name <- rep(NA,  length(unique(symbol2info$GENEID)))
  gene.start <- rep(NA,  length(unique(symbol2info$GENEID)))
  gene.end <- rep(NA,  length(unique(symbol2info$GENEID)))
  
  for(i in 1:length(unique(symbol2info$GENEID))){
    
    tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
    
    chr_name[i] <- tmp.transcript$TXCHROM
    gene.start[i] <- tmp.transcript$TXSTART
    gene.end[i] <- tmp.transcript$TXEND
    
  }
  
  seq.info <- data.frame(mgi_symbol = go2symbol$SYMBOL, chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
  seq.info[,2] <- as.character(seq.info[,2])
  seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6  
  
  #biomart method for getting gene metadata
  # seq.info <- getBM(attributes = c("mgi_symbol", "chromosome_name", "start_position", "end_position") , filters = "go" , values = process.ano ,mart = mouse)
  # seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6
  #get rid of weird chromosome names
  if(length(grep(seq.info$chromosome_name, pattern = "CHR")) > 0) seq.info <- seq.info[-grep(seq.info$chromosome_name, pattern = "CHR"),]
  
  seq.indexes <- matrix(NA, ncol = 3, nrow = dim(seq.info)[1])
  #we have seq.info which gives us a gene name and its location on the chromosome
  
  for(j in 1 : dim(seq.info)[1]){
    #seq.indexes <- rbind(seq.indexes, cbind(seq.info[j,1],MM_snps[which(MM_snps$chr == seq.info[j,2] & MM_snps$pos > mean(as.numeric(seq.info[j,3:4])) - .07 & MM_snps$pos < mean(as.numeric(seq.info[j,3:4])) + .07), c(1,3)]))
    tmp.indexes <-  combined.markers[which(combined.markers$chr == seq.info[j,2] & combined.markers$Mbp_mm10 > mean(as.numeric(seq.info[j,3:4])) - 2 & combined.markers$Mbp_mm10 < mean(as.numeric(seq.info[j,3:4])) + 2), c(1,3)]
    #for each gene, select the marker closest to the middle of the gene
    seq.indexes[j,] <- as.matrix(cbind(seq.info[j,1],tmp.indexes[which.min(abs(tmp.indexes[,2] - mean(as.numeric(seq.info[j,3:4])))),]))
  }
  
  probs.rows <- matrix(NA, nrow = nrow(Y), ncol = nrow(seq.indexes) * 8)
  probrowseq <- seq(1, ncol(probs.rows) + 8, by = 8)
  
  gene.names <- seq.info[,1]
  
  for(i in 1:nrow(seq.indexes)) probs.rows[, probrowseq[i]:(probrowseq[i+1] - 1) ] <- as.matrix(collect(tbl(DO_probs_DB, seq.indexes[i,2])))
  #fit pls2B, need duv, gene names, seq.indexes
  if(cv){
    pls.svd.cv <- perf_mddsPLS(Xs = probs.rows, Y = Y, lambda_min = .03, lambda_max = .15, n_lambda = 4, kfolds = 10, R = 1, mode = "reg", NCORES = 11)
    pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = 1, lambda = pls.svd.cv$Optim$optim_para_one[1])
  } else {pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = 1, lambda = lambda)
  }
  #cache to pls list for new analyses
  results <- list(pls.svd, gene.names, seq.info, probs.rows)
  
  tmp.reactive <- results
  reactive.svd <- tmp.reactive[[1]]$mod$u[[1]]
  gene.names <- tmp.reactive[[2]]
  seq.info <- tmp.reactive[[3]]
  
  #now we should be able to take pls.svd directly and maybe label them by founder in a new column, then barplot by family, by gene
  do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
  do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
  
  pathway.loadings <- data.frame(gloadings = reactive.svd[,1], gnames = as.character(rep(seq.info[,1], each = 8)), founders = rep(do.names, nrow(seq.info)))
  
  
  tmp.reactive <- results
  probs.rows <- tmp.reactive[[4]]
  
  snp.dim = 1#1
  
  #calculate projection
  proj.coords.a1 <- row2array3d(predict(tmp.reactive[[1]], probs.rows[c(which.min(tmp.reactive[[1]]$mod$ts[[1]][,1]), which.max(tmp.reactive[[1]]$mod$ts[[1]][,1])),])$y, Nlandmarks = ncol(Y)/3)
  proj.coords.a2 <- proj.coords.a1[,,2]
  proj.coords.a1 <- proj.coords.a1[,,1]
  
  
  return(list(loadings = pathway.loadings, pheno1 = proj.coords.a1, pheno2 = proj.coords.a2))
}




#* @apiTitle MGP API

#* Run an MGP model
#* @param GO.term GO term to run
#* @param mutant Mutant for comparison
#* @param lambda Regularization strength
#* @get /mgp_loadings

function(GO.term = "chondrocyte differentiation", lambda = .06) {
  future::future({
    mgp(GO.term = GO.term, lambda = as.numeric(lambda), Y = Y)
  })
}

# mgp <- function(GO.term = "GO:0002062", mutant = "Alk6", cv = F, lambda = .06){
# 
#   process.ano <- as.character(GO.term)
#   coi <- c("ENSEMBL", "SYMBOL")
#   go2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = process.ano, columns = coi, keytype = "GO")[,-2:-3]))
#   coi2 <- c("TXCHROM", "TXSTART", "TXEND")
#   symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = go2symbol[,2], columns = coi2, keytype="GENEID")
# 
# 
#   transcipt.size <- abs(symbol2info[,3] - symbol2info[,4])
# 
#   #symbol, chr, start, end
#   chr_name <- rep(NA,  length(unique(symbol2info$GENEID)))
#   gene.start <- rep(NA,  length(unique(symbol2info$GENEID)))
#   gene.end <- rep(NA,  length(unique(symbol2info$GENEID)))
# 
#   for(i in 1:length(unique(symbol2info$GENEID))){
# 
#     tmp.transcript <- symbol2info[symbol2info[,1] == unique(symbol2info$GENEID)[i],][which.max(transcipt.size[symbol2info[,1] == unique(symbol2info$GENEID)[i]]),]
# 
#     chr_name[i] <- tmp.transcript$TXCHROM
#     gene.start[i] <- tmp.transcript$TXSTART
#     gene.end[i] <- tmp.transcript$TXEND
# 
#   }
# 
#   seq.info <- data.frame(mgi_symbol = go2symbol$SYMBOL, chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
#   seq.info[,2] <- as.character(seq.info[,2])
#   seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6
# 
#   #biomart method for getting gene metadata
#   # seq.info <- getBM(attributes = c("mgi_symbol", "chromosome_name", "start_position", "end_position") , filters = "go" , values = process.ano ,mart = mouse)
#   # seq.info[,3:4] <- as.matrix(seq.info[,3:4])/1e6
#   #get rid of weird chromosome names
#   if(length(grep(seq.info$chromosome_name, pattern = "CHR")) > 0) seq.info <- seq.info[-grep(seq.info$chromosome_name, pattern = "CHR"),]
# 
#   seq.indexes <- matrix(NA, ncol = 3, nrow = dim(seq.info)[1])
#   #we have seq.info which gives us a gene name and its location on the chromosome
# 
#   for(j in 1 : dim(seq.info)[1]){
#     #seq.indexes <- rbind(seq.indexes, cbind(seq.info[j,1],MM_snps[which(MM_snps$chr == seq.info[j,2] & MM_snps$pos > mean(as.numeric(seq.info[j,3:4])) - .07 & MM_snps$pos < mean(as.numeric(seq.info[j,3:4])) + .07), c(1,3)]))
#     tmp.indexes <-  combined.markers[which(combined.markers$chr == seq.info[j,2] & combined.markers$Mbp_mm10 > mean(as.numeric(seq.info[j,3:4])) - 2 & combined.markers$Mbp_mm10 < mean(as.numeric(seq.info[j,3:4])) + 2), c(1,3)]
#     #for each gene, select the marker closest to the middle of the gene
#     seq.indexes[j,] <- as.matrix(cbind(seq.info[j,1],tmp.indexes[which.min(abs(tmp.indexes[,2] - mean(as.numeric(seq.info[j,3:4])))),]))
#   }
# 
#   probs.rows <- matrix(NA, nrow = nrow(Y), ncol = nrow(seq.indexes) * 8)
#   probrowseq <- seq(1, ncol(probs.rows) + 8, by = 8)
# 
#   gene.names <- seq.info[,1]
# 
#   for(i in 1:nrow(seq.indexes)) probs.rows[, probrowseq[i]:(probrowseq[i+1] - 1) ] <- as.matrix(collect(tbl(DO_probs_DB, seq.indexes[i,2])))
#   #fit pls2B, need duv, gene names, seq.indexes
#   if(cv){
#     pls.svd.cv <- perf_mddsPLS(Xs = probs.rows, Y = Y, lambda_min = .03, lambda_max = .15, n_lambda = 4, kfolds = 10, R = 1, mode = "reg", NCORES = 11)
#     pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = 1, lambda = pls.svd.cv$Optim$optim_para_one[1])
#   } else {pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = 1, lambda = as.numeric(lambda))
#   }
#   #cache to pls list for new analyses
#   results <- list(pls.svd, gene.names, seq.info, probs.rows)
#    # return(jsonlite::toJSON(results, force = TRUE))
# 
#   return(gene.names)
# 
#  #  tmp.reactive <- results
#  #  reactive.svd <- tmp.reactive[[1]]$mod$u[[1]]
#  #  gene.names <- tmp.reactive[[2]]
#  #  seq.info <- tmp.reactive[[3]]
#  # 
#  #  #now we should be able to take pls.svd directly and maybe label them by founder in a new column, then barplot by family, by gene
#  #  do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
#  #  do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
#  # 
#  #  pathway.loadings <- data.frame(gloadings = reactive.svd[,1], gnames = as.character(rep(seq.info[,1], each = 8)), founders = rep(do.names, nrow(seq.info)))
#  # 
#  # 
#  # p <- ggplot() +
#  #    geom_bar(data = pathway.loadings,
#  #             aes(x = gnames, y = gloadings),
#  #             stat = "identity",
#  #             width = .75,
#  #             position=position_dodge()) +
#  #    geom_point(data = pathway.loadings,
#  #               aes(x = gnames, y = gloadings, color = founders),
#  #               shape = "-",
#  #               size = 1) +
#  #    scale_color_manual(values=do.colors,
#  #                       guide = guide_legend(title = "Founder\nGenotype", override.aes = list(shape = rep(19, 8), size = 1))) +
#  #    xlab("Gene") +
#  #    ylab("Genetic marker loading") +
#  #    theme(text = element_text(size=6),
#  #          axis.text.x = element_text(angle = 75, hjust = 1),
#  #          axis.title.x = element_text(margin = margin(t = 20)),
#  #          axis.text = element_text(angle = 55, hjust = 1, size = 12),
#  #          axis.title = element_text(size = 12, face = "bold"),
#  #          legend.text = element_text(size = 8),
#  #          legend.title = element_text(size = 8, face = "bold", hjust = .5),
#  #          plot.background = element_rect(fill = rgb(245/255, 245/255, 245/255, .9), colour = rgb(245/255, 245/255, 245/255, .9)),
#  #          legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9)))
#  #  
#  # 
#  # 
#  #  ggplotly(p +
#  #                   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
#  #                         axis.text.y = element_text(size = 8),
#  #                         axis.title = element_text(size = 12, face = "bold")),
#  #                 legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9))) %>% layout(
#  #                   margin = list(b = 100, l = 50) # to fully display the x and y axis labels
#  #                 )
# 
#   # 
#   # tmp.reactive <- results
#   # probs.rows <- tmp.reactive[[4]]
#   # 
#   # snp.dim = 1#1
#   # 
#   # #calculate projection
#   # if(is.null(mutant) == F){
#   #   proj.coords.a1 <- row2array3d(predict(tmp.reactive[[1]], probs.rows[c(which.min(tmp.reactive[[1]]$mod$ts[[1]][,1]), which.max(tmp.reactive[[1]]$mod$ts[[1]][,1])),])$y, Nlandmarks = ncol(Y)/3)
#   #   proj.coords.a2 <- proj.coords.a1[,,2]
#   #   proj.coords.a1 <- proj.coords.a1[,,1]
#   #   
#   #   tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(array2row3d(row2array3d(mutant.lms[mutant.db$Genotype == mutant,])[,,]))), ncol(Y)/3, 3))$coords
#   #   # tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
#   #   do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
#   #   mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])
#   #   
#   #   do.mean <- matrix(colMeans(Y), nrow = ncol(Y)/3, ncol = 3, byrow = T)
#   #   
#   #   # shape.mean <- rotmesh.onto(mouse.ply, refmat = as.matrix(consensus.skull), tarmat = do.mean, scale = T, reflection = T)
#   #   
#   #   par3d(zoom = .65)
#   #   aspect3d("iso")
#   #   
#   #   #vectors from DO mean to mutant
#   #   plot3d(shape.mean$mesh, col = "lightgrey", alpha = .3, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
#   #   spheres3d(do.mean, radius = .001, color = adjustcolor("red", .3))
#   #   
#   #   
#   #   for(i in 1:nrow(do.mean)) arrow3d(do.mean[i,] - (mutant.mean[i,] - do.mean[i,]), mutant.mean[i,], type = "lines", col = "red", barblen = 0.005, lwd = 4.5)
#   #   for(i in 1:nrow(do.mean)) arrow3d(proj.coords.a1[i,], proj.coords.a2[i,] + (proj.coords.a2[i,] - proj.coords.a1[i,]) * (2 - 1), type = "lines", col = "black", barblen = 0.005, lwd = 5.5)
#   #   
#   #   rglwidget()
#   #   
#   #   MGP.mutant.cor <- cor(tmp.reactive[[1]]$mod$V_super[,1], manova(two.d.array(tmp.mutant.registration) ~ c(rep(0, nrow(Y)), rep(1, sum(mutant.db$Genotype == mutant))))$coef[2,])
#   # } else {MGP.mutant.cor <- NULL}
#   # return(list(mgp = pls.svd, gene.names = gene.names, seq.info = seq.info, probs.rows = probs.rows, cor = MGP.mutant.cor, gene.plot = p))
# 
# }



#* Return the sum of two numbers
#* @param a The first number to add
#* @param b The second number to add
#* @post /sum
function(a, b) {
  future::future({
    as.numeric(a) + as.numeric(b)
  })
}
