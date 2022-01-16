library(plumber)
library(Morpho)
library(ddsPLS)
library(Jovid)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(dplyr)
library(dbplyr)
library(future)
library(promises)
future::plan("multicore")


#save(combined.markers, DO.go, giga.pca, mutant.db, mutant.lms, shape.mean, Y, DO.probs, file = "/data/MGP_data/offline_data.Rdata")
#save(combined.markers, giga.pca, mutant.db, mutant.lms, shape.mean, Y, file = "~/shiny/shinyapps/MGP/shiny_data2.Rdata")

#local dirs
mmusculusEnsembl <- loadDb(file="~/shiny/shinyapps/MGP/ensemble.sqlite")
load("~/shiny/shinyapps/MGP/shiny_data.Rdata")
load("~/shiny/shinyapps/MGP/cached.results.Rdata")
DO_probs_DB <- src_sqlite("~/shiny/shinyapps/MGP/MGP_genotypes.sqlite")

#docker dirs
# setwd("/srv/shiny-server/")
# mmusculusEnsembl <- loadDb(file="ensemble.sqlite")
# load("shiny_data.Rdata")
# load("cached.results.Rdata")
# DO_probs_DB <- src_sqlite("MGP_genotypes.sqlite")

#genopheno deployment dirs
# mmusculusEnsembl <- loadDb(file="/data/MGP_data/ensemble.sqlite")
# load("/data/MGP_data/shiny_data.Rdata")
# load("/data/MGP_data/cached.results.Rdata")
# DO_probs_DB <- src_sqlite("/data/MGP_data/MGP_genotypes.sqlite")

#core mgp function####
mgp <- function(GO.term = "chondrocyte differentiation", Y, cv = F, lambda = .06, pls_axis = 1, permutation = F){
  
  selection.vector <- GO.term
  # selection.vector <- process.list()[[1]][process.list()[[2]] %in% input$variables2]
  print(length(selection.vector))
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
    print(i)
    chr_name[i] <- tmp.transcript$TXCHROM
    gene.start[i] <- tmp.transcript$TXSTART
    gene.end[i] <- tmp.transcript$TXEND
    
  }
  
  seq.info <- data.frame(mgi_symbol = unique(go2symbol$SYMBOL), chromosome_name = chr_name, start_position = gene.start, end_position = gene.end)
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
    pls.svd.cv <- perf_mddsPLS(Xs = probs.rows, Y = Y, lambda_min = .03, lambda_max = .15, n_lambda = 4, kfolds = 10, R = pls_axis, mode = "reg", NCORES = 11)
    pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = pls_axis, lambda = pls.svd.cv$Optim$optim_para_one[1])
  } else {pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = pls_axis, lambda = lambda)
  }
  
  approximate.p <- NULL
  if(permutation > 0){
    perm.result <- rep(NA, permutation)
    for(i in 1:length(perm.result)){
    process.svd <- mddsPLS(Xs = probs.rows, Y = as.matrix(Y[sample(1:nrow(Y), size = nrow(Y)),]), R = pls_axis, lambda = lambda)
    
    full.pred <- predict(process.svd, probs.rows)$y
    ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y))^2))
    rss <- sum(apply(Y, 1, function(x) (x - colMeans(full.pred))^2))
    perm.result[i] <- ess/(rss + ess)
    }
    
    full.pred <- predict(pls.svd, probs.rows)$y
    ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y))^2))
    rss <- sum(apply(Y, 1, function(x) (x - colMeans(full.pred))^2))
    approximate.p <- sum(ess/(rss + ess) < perm.result)/length(perm.result)
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
  
  pathway.loadings <- data.frame(gloadings = reactive.svd[,pls_axis], gnames = as.character(rep(seq.info[,1], each = 8)), founders = rep(do.names, nrow(seq.info)))
  
  bar_order <- pathway.loadings %>% 
    group_by(gnames) %>%
    summarise(test = diff(range(gloadings))) %>%
    arrange(-test) 
  
  pathway.loadings$gnames <- factor(pathway.loadings$gnames, levels = lapply(bar_order, as.character)$gnames)

  tmp.reactive <- results
  probs.rows <- tmp.reactive[[4]]
  
  snp.dim = pls_axis
  
  #calculate projection
  proj.coords.a1 <- row2array3d(predict(tmp.reactive[[1]], probs.rows[c(which.min(tmp.reactive[[1]]$mod$ts[[snp.dim]][,1]), which.max(tmp.reactive[[1]]$mod$ts[[snp.dim]][,1])),])$y, Nlandmarks = ncol(Y)/3)
  proj.coords.a2 <- proj.coords.a1[,,2]
  proj.coords.a1 <- proj.coords.a1[,,1]

  return(list(loadings = pathway.loadings, pheno1 = proj.coords.a1, pheno2 = proj.coords.a2, pheno_loadings = pls.svd$mod$V_super[,pls_axis], p_value = approximate.p))
}

#custom MGP function####
custom.mgp <- function(genelist = c("Bmp7, Bmp2, Bmp4, Ankrd11"), Y, cv = F, lambda = .06, pls_axis = 1, permutation = 0){
  
  selection.vector <- as.character(strsplit(genelist, ", ")[[1]])
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
    pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = pls_axis, lambda = lambda)
    
    approximate.p <- NULL
    if(permutation > 0){
      perm.result <- rep(NA, permutation)
      for(i in 1:length(perm.result)){
        process.svd <- mddsPLS(Xs = probs.rows, Y = as.matrix(Y[sample(1:nrow(Y), size = nrow(Y)),]), R = pls_axis, lambda = lambda)
        
        full.pred <- predict(process.svd, probs.rows)$y
        ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y))^2))
        rss <- sum(apply(Y, 1, function(x) (x - colMeans(full.pred))^2))
        perm.result[i] <- ess/(rss + ess)
      }
      
      full.pred <- predict(pls.svd, probs.rows)$y
      ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y))^2))
      rss <- sum(apply(Y, 1, function(x) (x - colMeans(full.pred))^2))
      approximate.p <- sum(ess/(rss + ess) < perm.result)/length(perm.result)
    }
    
    
    results <- list(pls.svd, gene.names, seq.info, probs.rows)
    
    tmp.reactive <- results
    reactive.svd <- tmp.reactive[[1]]$mod$u[[1]]
    gene.names <- tmp.reactive[[2]]
    seq.info <- tmp.reactive[[3]]
    
    #now we should be able to take pls.svd directly and maybe label them by founder in a new column, then barplot by family, by gene
    do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
    do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
    
    pathway.loadings <- data.frame(gloadings = reactive.svd[,pls_axis], gnames = as.character(rep(seq.info[,1], each = 8)), founders = rep(do.names, nrow(seq.info)))
    
    bar_order <- pathway.loadings %>% 
      group_by(gnames) %>%
      summarise(test = diff(range(gloadings))) %>%
      arrange(-test) 
    
    pathway.loadings$gnames <- factor(pathway.loadings$gnames, levels = lapply(bar_order, as.character)$gnames)
    
    tmp.reactive <- results
    probs.rows <- tmp.reactive[[4]]
    
    snp.dim = pls_axis
    
    #calculate projection
    proj.coords.a1 <- row2array3d(predict(tmp.reactive[[1]], probs.rows[c(which.min(tmp.reactive[[1]]$mod$ts[[snp.dim]][,1]), which.max(tmp.reactive[[1]]$mod$ts[[snp.dim]][,1])),])$y, Nlandmarks = ncol(Y)/3)
    proj.coords.a2 <- proj.coords.a1[,,2]
    proj.coords.a1 <- proj.coords.a1[,,1]
    
    return(list(loadings = pathway.loadings, pheno1 = proj.coords.a1, pheno2 = proj.coords.a2, pheno_loadings = pls.svd$mod$V_super[,pls_axis], p_value = approximate.p))
}

#set CORS parameters####
#* @filter cors
cors <- function(res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "Origin, X-Requested-With, Content-Type, Accept, Authorization")
  res$setHeader("Access-Control-Allow-Credentials", "true")
  plumber::forward()
}

#API definition####
#* @apiTitle MGP API

#* Run an MGP model
#* @param GO.term GO term to run
#* @param lambda Regularization strength
#* @param pls_axis how many axes to return
#* @param pheno what phenotype to use
#* @param permutation how many permutations should we use for testing? Max 200.
#* @get /mgp

function(GO.term = "chondrocyte differentiation", lambda = .06, pls_axis = 1, pheno = "Y", pheno_index = "1:54", permutation = 0) {
  future::future({
    print(GO.term)
    
    GO.term <- strsplit(GO.term, split = ",")[[1]]
    coordinate.table <- matrix(1:(54 * 3), ncol = 3, byrow = T)
    selected.pheno <- eval(parse(text = pheno_index))
    pheno.xyz <- as.numeric(t(coordinate.table[selected.pheno,]))
    
    mgp(GO.term = GO.term, lambda = as.numeric(lambda), Y = subset(get(pheno), select = pheno.xyz), pls_axis = as.numeric(pls_axis), permutation = permutation)
    
  })
}

#* Run a custom MGP model
#* @param genelist comma-separated list of gene names
#* @param lambda Regularization strength
#* @param pls_axis how many axes to return
#* @param pheno what phenotype to use
#* @param permutation should we run a permutation test?
#* @get /custom_mgp


function(genelist = c("Bmp7, Bmp2, Bmp4, Ankrd11"), lambda = .06, pls_axis = 1,  pheno = "Y", pheno_index = "1:54", permutation = 0) {
  future::future({
    
    coordinate.table <- matrix(1:(54 * 3), ncol = 3, byrow = T)
    selected.pheno <- eval(parse(text = pheno_index))
    pheno.xyz <- as.numeric(t(coordinate.table[selected.pheno,]))
    
    custom.mgp(genelist = genelist, lambda = as.numeric(lambda), Y = subset(get(pheno), select = pheno.xyz), pls_axis = as.numeric(pls_axis), permutation = permutation)
  })

}


#* Draw a random correlation matrix for MGP effects in the cache
#* @post /mgp_cor

function(num.processes = 5){
  future::future({
    chosen.processes <- sample(1:nrow(cached.params), as.numeric(num.processes))
    tmp.proc.names <- rep(NA, num.processes)
    tmp.cor <- matrix(NA, nrow = 162, ncol = as.numeric(num.processes))
    for(i in 1:num.processes){
      tmp.cor[,i] <- cached.results[[chosen.processes[i]]][[1]][[1]]$mod$V_super[,1]
      tmp.split <- strsplit(cached.params[chosen.processes[i],1], split = "\\|")[[1]]
      name.buffer <- NULL
      for(j in 1 : length(tmp.split)) {
        name.buffer <- c(name.buffer, as.character(DO.go[DO.go[,2] == tmp.split[j],3]))
      }
      tmp.proc.names[i] <- paste(name.buffer, collapse = " + ") 
    }
    return(list(cormat = cor(tmp.cor), process.names = tmp.proc.names))
  }, seed = NULL)
}

#* Query the Ensembl database
#* @param process
#* @get /mgi

function(process = "chondrocyte differentiation") {
  future::future({
    process.call <- DO.go[DO.go[,3] == process, -1]
    colnames(process.call) <- c("go_id", "go_term", "ngenes")
    
    selection.vector <- c(process)
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
    
    return(list(process_call = process.call, seq_info = seq.info))
    
  })
}


#* Get the entire GO network for a process or gene list
#* @param genelist custom comma-separated list of genes
#* @param process process name
#* @get /GO_network
function(genelist = NULL, process = NULL){
future_promise({
  if(is.null(genelist) == F){
selection.vector <- as.character(strsplit(genelist, ", ")[[1]])
coi <- c("ENSEMBL", "SYMBOL")
gene2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = selection.vector, columns = coi, keytype = "SYMBOL")))

coi2 <- c("GO", "SYMBOL")
ensembl2go <- AnnotationDbi::select(org.Mm.eg.db, keys = gene2symbol[,2], columns = coi2, keytype="ENSEMBL")

GO_ensembl_join <- right_join(DO.go, ensembl2go, by = c("V2" = "GO"))
}

  if(is.null(process) == F){
    
    selection.vector <- c(process)
  
    process.ano <- NULL
    for(i in 1: length(selection.vector)) process.ano <- c(process.ano, as.character(DO.go[DO.go[,3] == selection.vector[i], 2]))
    coi <- c("ENSEMBL", "SYMBOL")
    go2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = process.ano, columns = coi, keytype = "GO")[,-2:-3]))
    coi2 <- c("GO", "SYMBOL")
    ensembl2go <- AnnotationDbi::select(org.Mm.eg.db, keys = go2symbol[,2], columns = coi2, keytype="ENSEMBL")
    GO_ensembl_join <- right_join(DO.go, ensembl2go, by = c("V2" = "GO"))
  }
  
  return(GO_ensembl_join)
  
})
}

#* Get the list of mutants to make comparisons with
#* @get /mutant_list
function(){
  future_promise({
    as.character(unique(mutant.db$Genotype))
  })
}

#* mutant vector correlation to MGP vector
#* @param MGP_pheno_loadings phenotype loadings from MGP model
#* @param mutant selected mutant from /mutant_list
#* @serializer print
#* @get /mutant_comparison
function(MGP_pheno_loadings = NULL, mutant = "Bmp2"){
  future_promise({
    MGP_pheno_loadings <- as.numeric(strsplit(MGP_pheno_loadings, ", ")[[1]])
    tmp.mutant.registration <- geomorph::gpagen(geomorph::arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == mutant,])), 54, 3))$coords
    mutant.loadings <- as.numeric(manova(geomorph::two.d.array(tmp.mutant.registration) ~ c(rep(0, nrow(Y)), rep(1, sum(mutant.db$Genotype == mutant))))$coef[2,])
    
    mutant.cor <- cor.test(as.numeric(MGP_pheno_loadings), mutant.loadings)
    
    return(mutant.cor)
  })
}

