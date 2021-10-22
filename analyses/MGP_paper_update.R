###are fgf10 mutants asymmetric?

library(shiny)
library(shinydashboard)
library(ggplot2)
# library(biomaRt)
library(Morpho)
library(rgl)
library(geomorph)
library(plotly)
library(rmarkdown)
library(shapes)
library(ddsPLS)
library(Jovid)
#paper resubimssion analyses

library(shinycssloaders)
library(shinyjs)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(dplyr)
library(dbplyr)



mgp <- function(GO.term, mutant = NULL, Y, cv = T, lambda = .06){
  
  process.ano <- GO.term
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
  
  p2 <- ggplotly(p +
             theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                   axis.text.y = element_text(size = 8),
                   axis.title = element_text(size = 12, face = "bold")),
           legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9))) %>% layout(
             margin = list(b = 100, l = 50) # to fully display the x and y axis labels
           )
  
  
  tmp.reactive <- results
  probs.rows <- tmp.reactive[[4]]
  
  snp.dim = 1#1
  
  #calculate projection
  if(is.null(mutant) == F){
  proj.coords.a1 <- row2array3d(predict(tmp.reactive[[1]], probs.rows[c(which.min(tmp.reactive[[1]]$mod$ts[[1]][,1]), which.max(tmp.reactive[[1]]$mod$ts[[1]][,1])),])$y, Nlandmarks = ncol(Y)/3)
  proj.coords.a2 <- proj.coords.a1[,,2]
  proj.coords.a1 <- proj.coords.a1[,,1]
  
  tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(array2row3d(row2array3d(mutant.lms[mutant.db$Genotype == mutant,])[,,]))), ncol(Y)/3, 3))$coords
  # tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
  do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
  mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])
  
  do.mean <- matrix(colMeans(Y), nrow = ncol(Y)/3, ncol = 3, byrow = T)
  
  # shape.mean <- rotmesh.onto(mouse.ply, refmat = as.matrix(consensus.skull), tarmat = do.mean, scale = T, reflection = T)
  
  par3d(zoom = .65)
  aspect3d("iso")
  
  #vectors from DO mean to mutant
  plot3d(shape.mean$mesh, col = "lightgrey", alpha = .3, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
  spheres3d(do.mean, radius = .001, color = adjustcolor("red", .3))
  
  
  for(i in 1:nrow(do.mean)) arrow3d(do.mean[i,] - (mutant.mean[i,] - do.mean[i,]), mutant.mean[i,], type = "lines", col = "red", barblen = 0.005, lwd = 4.5)
  for(i in 1:nrow(do.mean)) arrow3d(proj.coords.a1[i,], proj.coords.a2[i,] + (proj.coords.a2[i,] - proj.coords.a1[i,]) * (2 - 1), type = "lines", col = "black", barblen = 0.005, lwd = 5.5)
  
  rglwidget()
  
  MGP.mutant.cor <- cor(tmp.reactive[[1]]$mod$V_super[,1], manova(two.d.array(tmp.mutant.registration) ~ c(rep(0, nrow(Y)), rep(1, sum(mutant.db$Genotype == mutant))))$coef[2,])
  } else {MGP.mutant.cor <- NULL}
  return(list(mgp = pls.svd, gene.names = gene.names, seq.info = seq.info, probs.rows = probs.rows, cor = MGP.mutant.cor, gene.plot = p))
  
}

loading.plot <- function(mgp, point.size = 10, axis.text.size = 10, axis.title.size = 10, legend.size = 12, gene2highlight){
  do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
  do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")
  
  pathway.loadings <- data.frame(gloadings = mgp$mgp$mod$u[[1]][,1], gnames = as.character(rep(mgp$seq.info[,1], each = 8)), founders = rep(do.names, nrow(mgp$seq.info)))
  
  gname.colors <- rep(1, nrow(mgp$seq.info))
  gname.colors[which(levels(mgp$seq.info[,1]) == gene2highlight)] <- 2
  
  ultimate_theme <-   theme(
    text = element_text(size=6), 
    axis.text.x = element_text(angle = 75, hjust = 1, color = gname.colors),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.text = element_text(angle = 55, hjust = 1, size = axis.text.size),
    axis.title = element_text(size = axis.title.size, face = "bold"),
    legend.text = element_text(size = legend.size), 
    legend.title = element_text(size = legend.size, face = "bold", hjust = .5))
  
  p <- ggplot() +
    geom_bar(data = pathway.loadings, 
             aes(x = gnames, y = gloadings), 
             stat = "identity", 
             width = .75, 
             position=position_dodge()) +
    geom_point(data = pathway.loadings,
               aes(x = gnames, y = gloadings, color = founders),
               shape = "-",
               size = point.size) +
    scale_color_manual(values=do.colors, 
                       guide = guide_legend(title = "Founder\nGenotype", override.aes = list(shape = rep(19, 8), size = 2))) +
    xlab("Gene") +
    ylab("Genetic marker loading") + 
    ultimate_theme
  
  # p2 <- ggplotly(p +
  #                  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
  #                        axis.text.y = element_text(size = 8),
  #                        axis.title = element_text(size = 12, face = "bold")),
  #                legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9))) %>% layout(
  #                  margin = list(b = 100, l = 50) # to fully display the x and y axis labels
  #                )
  return(p)
}

mmusculusEnsembl <- loadDb(file="~/shiny/shinyapps/MGP/ensemble.sqlite")
# load("~/shiny/shinyapps/MGP/offline_data.Rdata")
load("~/shiny/shinyapps/MGP/new_mutants.Rdata")

# DO_probs_DB <- s3read_using(FUN = src_sqlite, object = "s3://mgpgenotypes/MGP_genotypes.sqlite") #src_sqlite("mgpgenotypes.s3.ca-central-1.amazonaws.com/MGP_genotypes.sqlite")
# load("~/shiny/shinyapps/MGP/cached.results.Rdata")

# save(combined.markers, DO.go, giga.pca, mutant.db, mutant.lms, shape.mean, Y, file = "~/shiny/shinyapps/MGP/shiny_data.Rdata")
load("~/shiny/shinyapps/MGP/shiny_data.Rdata")
DO_probs_DB <- src_sqlite("~/shiny/shinyapps/MGP/MGP_genotypes.sqlite")

#is there significant asymmetry in FGF10 mutants?####
#define sides/midline
colnames(mutant.db)[12:length(colnames(mutant.db))]
right <- 1:(69/3)
left <- 32:54
pairedLM <- cbind(left,right)

#symmetric procrustes
mutant.sym <- procSym(arrayspecs(mutant.lms, 54,3)[,,which(mutant.db$Genotype == "Fgf10CuER__R(-P+TA)__Fgf10_reton")], pairedLM = pairedLM)

procAOVsym(mutant.sym)

do.mean <- matrix(colMeans(Y), nrow = 54, ncol = 3, byrow = T)

plot3d(shape.mean$mesh, color = "lightgrey", alpha = .1, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
points3d(do.mean[-pairedLM,], color = 2)
text3d(do.mean[-pairedLM,], texts = c(1:54)[-pairedLM])
#good plane: 31, 25, 27
mirror.right <- mirror2plane(do.mean[left,], v1 = do.mean[31,], v2 = do.mean[25,], v3 = do.mean[26,])

plot3d(shape.mean$mesh, color = "lightgrey", alpha = .1, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
points3d(do.mean[left,], color =2)
points3d(mirror.right, color =4)

#now do this for fgf mutant
# fgf.mesh <- rotmesh.onto(shape.mean$mesh, shape.mean$yrot, mutant.sym$mshape, scale = T)
# we should TPS deform the mesh for a better figure
fgf.mesh <- tps3d(shape.mean$mesh, shape.mean$yrot, mutant.sym$mshape, scale = T)
asym.fgf <- Jovid::array.mean(mutant.sym$Sym + mutant.sym$Asym)
fgf.mirror.right <- mirror2plane(asym.fgf[left,], v1 = asym.fgf[31,], v2 = asym.fgf[25,], v3 = asym.fgf[27,])

plot3d(fgf.mesh, color = "lightgrey", alpha = .2, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
mag <- 0
for(i in 1:23) arrow3d(fgf.mirror.right[i,] - (asym.fgf[right,][i,] - fgf.mirror.right[i,]) * mag, asym.fgf[right,][i,], type = "lines", col = "black", barblen = 0, lwd = 3)

#check that asymmetry differs in magnitude at diff points: rowSums(sqrt((asym.fgf[left,] - fgf.mirror.right)^2))

spheres3d(asym.fgf[right,], color =2, radius = .001)
spheres3d(fgf.mirror.right[,] - (asym.fgf[right,][,] - fgf.mirror.right[,]) * mag, color = 4, radius = .001)
# arrow3d(fgf.mirror.right, asym.fgf[right,])



#left/right symmetry
# GO:0007368
lr.symm <- mgp("GO:0007368", "Fgf10CuER__R(-P+TA)__Fgf10_reton", Y, cv = T)

lr.symm$gene.plot

png("~/Desktop/lrsymm_gloadings.png", height = 750, width = 1100)
loading.plot(lr.symm, point.size = 11.25, axis.text.size = 11.5, legend.size = 18.5, gene2highlight = "Fgf10")
dev.off()

#are the DO mice asymmetric?
#load asymmetric dataset
Y_asym <- read.csv("~/Downloads/asymmetic_DO_lms.csv", row.names = 1)
#make sure Y_asym is in the same order as Y
Y_asym <- Y_asym[rownames(Y_asym) %in% rownames(Y),]
Y_asym <- Y_asym[sort(rownames(Y_asym), index.return = T)$ix,]

lr.symm <- mgp("GO:0007368", "Fgf10CuER__R(-P+TA)__Fgf10_reton", as.matrix(Y_asym), cv = T)


DO.sym <- procSym(arrayspecs(Y_asym, 54,3), pairedLM = pairedLM)

lr.symm <- mgp("GO:0007368", "Fgf10CuER__R(-P+TA)__Fgf10_reton", two.d.array(DO.sym$Asym), cv = T)


procAOVsym(DO.sym)


#which DO mice are homozygous for the cast fgf10 allele? FGF10 is number 11 on the list of snps
fgf10.probs <- lr.symm$probs.rows[,((8 * 11)+1):(8*12)]

#cast is column 4
sum(fgf10.probs[,4] > .6)
sum(fgf10.probs[,4] < .6 & fgf10.probs[,4] > .2)
sum(fgf10.probs[,4] < .2)

#procAOV uses grep to find indnames--that throws an error for our naming converntions, so i'm setting the rowname to null
# rownames(Y_asym) <- NULL

cast.homo <- procSym(arrayspecs(Y_asym, 54,3)[,,which(fgf10.probs[,4] > .6)], pairedLM = pairedLM)
cast.het <- procSym(arrayspecs(Y_asym, 54,3)[,,which(fgf10.probs[,4] < .6 & fgf10.probs[,4] > .2)], pairedLM = pairedLM)
cast.nomo <- procSym(arrayspecs(Y_asym, 54,3)[,,which(fgf10.probs[,4] < .2)], pairedLM = pairedLM)

procAOVsym <- function (symproc, indnames = NULL){
  if (class(symproc) != "symproc") {
    stop("input is not of class symproc")
  }
  m <- dim(symproc$rotated)[2]
  k <- dim(symproc$rotated)[1]
  n <- dim(symproc$rotated)[3]
  if (is.null(indnames)) 
    indnames <- as.factor(rownames(symproc$PCscore_sym))
  else indnames <- as.factor(indnames)
  indlev <- levels(indnames)
  nlev <- length(indlev)
  alist <- list()
  for (i in 1:length(indlev)) {
    alist[[i]] <- (indnames)[i]
  }
  checkinds <- lapply(alist, function(x) {
    length(x) == length(alist[[1]])
  })
  if (prod(as.integer(unlist(checkinds))) == 0) 
    stop("same number of digitizations is needed for all specimen")
  r <- length(alist[[1]])
  pl <- dim(symproc$pairedLM)[1]
  sl <- k - 2 * pl
  if (k == 3) {
    df.ind <- (nlev - 1) * (3 * pl + 2 * sl - 4)
    df.side <- 3 * pl + sl - 3
    df.indxside <- (nlev - 1) * (3 * pl + sl - 3)
    df.res <- (r - 1) * nlev * (3 * (2 * pl + sl) - 7)
  }
  else {
    df.ind <- (nlev - 1) * (2 * pl + sl - 2)
    df.side <- 2 * pl + sl - 2
    df.indxside <- (nlev - 1) * (2 * pl + sl - 2)
    df.res <- (r - 1) * nlev * (2 * (2 * pl + sl) - 4)
  }
  side <- sum(symproc$asymmean^2) * n
  asymfun <- function(x) {
    x <- sum(colMeans(symproc$Asymtan[x, ])^2) * length(x)
    return(x)
  }
  symfun <- function(x) {
    x <- sum(colMeans(symproc$Symtan[x, ])^2) * length(x)
    return(x)
  }
  if (length(alist[[1]]) > 1) {
    indxside <- sum(unlist(lapply(alist, asymfun)))
    ind <- sum(unlist(lapply(alist, symfun)))
  }
  else {
    indxside <- sum(symproc$Asymtan^2)
    ind <- sum(symproc$Symtan^2)
  }
  allsq <- sum(symproc$Asymtan[, ]^2 + symproc$Symtan[, ]^2) + 
    side
  res <- allsq - (ind + indxside + side)
  outss <- c(ind, side, indxside, res)
  outdf <- as.integer(c(df.ind, df.side, df.indxside, df.res))
  exVar <- outss/allsq
  outms <- outss/outdf
  out <- data.frame(outss, outms, exVar, outdf)
  F.values <- c(outms[1]/outms[3], outms[2]/outms[3], outms[3]/outms[4], 
                NA)
  sig <- c(pf(F.values[1], outdf[1], outdf[3], lower.tail = F), 
           pf(F.values[2], outdf[2], outdf[3], lower.tail = F), 
           pf(F.values[3], outdf[3], outdf[4], lower.tail = F), 
           NA)
  out <- data.frame(out, F.values, sig)
  rownames(out) <- c("ind", "side", "ind.x.side", "error")
  colnames(out) <- c("SS", "MS", "exVar", "df", "F", "p-value")
  return(out)
}


procAOVsym(cast.homo)
procAOVsym(cast.het)
procAOVsym(cast.nomo)

cast.homo.asym <- rep(NA, dim(cast.homo$rotated)[3])
for(i in 1:dim(cast.homo$rotated)[3]){
cast.homo.asym[i] <- sum(sqrt(((cast.homo$Sym + cast.homo$Asym)[,,i] - cast.homo$symmean)^2))
}

cast.het.asym <- rep(NA, dim(cast.het$rotated)[3])
for(i in 1:dim(cast.het$rotated)[3]){
  cast.het.asym[i] <- sum(sqrt(((cast.het$Sym + cast.het$Asym)[,,i] - cast.het$symmean)^2))
}

cast.nomo.asym <- rep(NA, dim(cast.nomo$rotated)[3])
for(i in 1:dim(cast.nomo$rotated)[3]){
  cast.nomo.asym[i] <- sum(sqrt(((cast.nomo$Sym + cast.nomo$Asym)[,,i] - cast.nomo$symmean)^2))
}

fgf.cast.genos <- c(rep("Homozygous", dim(cast.homo$rotated)[3]), rep("Heterozygous", dim(cast.het$rotated)[3]), rep("Cast null", dim(cast.nomo$rotated)[3]))
fgf.cast.scores <- c(cast.homo.asym, cast.het.asym, cast.nomo.asym)

boxplot(fgf.cast.scores ~ fgf.cast.genos)

#c57 allele asymmetry####

c57.homo <- procSym(arrayspecs(Y_asym, 54,3)[,,which(fgf10.probs[,6] > .6)], pairedLM = pairedLM)
c57.het <- procSym(arrayspecs(Y_asym, 54,3)[,,which(fgf10.probs[,6] < .6 & fgf10.probs[,6] > .2)], pairedLM = pairedLM)
c57.nomo <- procSym(arrayspecs(Y_asym, 54,3)[,,which(fgf10.probs[,6] < .2)], pairedLM = pairedLM)

c57.homo.asym <- rep(NA, dim(c57.homo$rotated)[3])
for(i in 1:dim(c57.homo$rotated)[3]){
  c57.homo.asym[i] <- sum(sqrt(((c57.homo$Sym + c57.homo$Asym)[,,i] - c57.homo$symmean)^2))
}

c57.het.asym <- rep(NA, dim(c57.het$rotated)[3])
for(i in 1:dim(c57.het$rotated)[3]){
  c57.het.asym[i] <- sum(sqrt(((c57.het$Sym + c57.het$Asym)[,,i] - c57.het$symmean)^2))
}

c57.nomo.asym <- rep(NA, dim(c57.nomo$rotated)[3])
for(i in 1:dim(c57.nomo$rotated)[3]){
  c57.nomo.asym[i] <- sum(sqrt(((c57.nomo$Sym + c57.nomo$Asym)[,,i] - c57.nomo$symmean)^2))
}

fgf.c57.genos <- c(rep("Homozygous", dim(c57.homo$rotated)[3]), rep("Heterozygous", dim(c57.het$rotated)[3]), rep("c57 null", dim(c57.nomo$rotated)[3]))
fgf.c57.scores <- c(c57.homo.asym, c57.het.asym, c57.nomo.asym)

boxplot(fgf.c57.scores ~ fgf.c57.genos)

#asymmetry MGP?####
lr.asymm <- mgp("GO:0007368", "Fgf10CuER__R(-P+TA)__Fgf10_reton", as.matrix(Y_asym))
loading.plot(lr.asymm, gene2highlight = "Fgf10")

#lr variance explained and permutation####
lr.perm.r2 <- rep(NA, 1000)
for(i in 1:1000){
  process.svd <- mddsPLS(Xs = as.matrix(lr.asymm$probs.rows), Y = as.matrix(Y_asym[sample(1:nrow(Y), size = nrow(Y)),]), R = 1, lambda = .06)
  
  full.pred <- predict(process.svd, as.matrix(lr.asymm$probs.rows))$y
  ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y_asym))^2))
  rss <- sum(apply(Y_asym, 1, function(x) (x - colMeans(full.pred))^2))
  lr.perm.r2[i] <- ess/(rss + ess)
  print(i)
}

full.pred <- predict(lr.asymm$mgp, as.matrix(lr.asymm$probs.rows))$y
ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y_asym))^2))
rss <- sum(apply(Y_asym, 1, function(x) (x - colMeans(full.pred))^2))
ess/(rss + ess)

hist(lr.perm.r2)

#other permutation tests####
set.seed(38)
chond_perm <- MGP_permuter(38)

set.seed(73)

palate_perm <- MGP_permuter(73)



png("~/Desktop/lrsymm_perm.png", height = 375*2, width = 550*2)
ggplot(data.frame(r2 = lr.perm.r2), mapping = aes(x = r2, fill = "#adcae6")) + 
  geom_histogram(alpha = 1, show.legend = FALSE) + 
  geom_vline(xintercept = ess/(rss + ess)) +
  xlab(expression(paste("81 marker MGP R"^"2"))) +
  ylab("Frequency") +
  scale_fill_manual(values = "#adcae6") +
  theme_bw()+ theme(text = element_text(size = 20))  
dev.off()


png("~/Desktop/chond_perm.png", height = 375*2, width = 550*2)
ggplot(data.frame(r2 = chond_perm), mapping = aes(x = r2, fill = "#adcae6")) + 
  geom_histogram(alpha = 1, show.legend = FALSE) + 
  geom_vline(xintercept = .0215) +
  xlab(expression(paste("38 marker MGP R"^"2"))) +
  ylab("Frequency") +
  scale_fill_manual(values = "#adcae6") +
  theme_bw()+ theme(text = element_text(size = 20))  
dev.off()


png("~/Desktop/palate_perm.png", height = 375*2, width = 550*2)
ggplot(data.frame(r2 = palate_perm), mapping = aes(x = r2, fill = "#adcae6")) + 
  geom_histogram(alpha = 1, show.legend = FALSE) + 
  geom_vline(xintercept = .024) +
  xlab(expression(paste("73 marker MGP R"^"2"))) +
  ylab("Frequency") +
  scale_fill_manual(values = "#adcae6") +
  theme_bw()+ theme(text = element_text(size = 20))  
dev.off()


png("~/Desktop/lrsymm_gloadings.png", height = 750, width = 1100)
loading.plot(lr.asymm, point.size = 11.25, axis.text.size = 11.5, axis.title.size = 16.5, legend.size = 18.5, "Fgf10")
dev.off()

#get projected phenotypes and measure asymmetry####
proj.coords.a1 <- row2array3d(predict(lr.asymm$mgp, lr.asymm$probs.rows[c(which.min(lr.asymm$mgp$mod$ts[[1]][,1]), which.max(lr.asymm$mgp$mod$ts[[1]][,1])),])$y, Nlandmarks = ncol(Y)/3)
proj.coords.a2 <- proj.coords.a1[,,2]
proj.coords.a1 <- proj.coords.a1[,,1]

colnames(Y_asym) <- colnames(as.matrix(mutant.lms[mutant.db$Genotype == "Fgf10CuER__R(-P+TA)__Fgf10_reton",]))
tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y_asym, as.matrix(mutant.lms[mutant.db$Genotype == "Fgf10CuER__R(-P+TA)__Fgf10_reton",])), 54, 3))$coords
# tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])


do.mean.asym <- matrix(colMeans(Y_asym), nrow = ncol(Y)/3, ncol = 3, byrow = T)

asym.shape.mean <- rotmesh.onto(shape.mean$mesh, refmat = as.matrix(shape.mean$yrot), tarmat = do.mean.asym, scale = T, reflection = T)

par3d(zoom = .65, userMatrix = lat.view)
aspect3d("iso")

plot3d(asym.shape.mean$mesh, col = adjustcolor("lightgrey", .3), alpha = .2, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
spheres3d(do.mean.asym, radius = .001, color = adjustcolor("red", .3))

for(i in 1:nrow(do.mean)) arrow3d(proj.coords.a1[i,], proj.coords.a2[i,] + (proj.coords.a2[i,] - proj.coords.a1[i,]) * (4 - 1), type = "lines", col = "black", barblen = 0.005, lwd = 4)
for(i in 1:54) arrow3d(do.mean.asym[i,] - (mutant.mean[i,] - do.mean.asym[i,]), mutant.mean[i,], type = "lines", col = "red", barblen = 0.005, lwd = 3)
rglwidget()




#visualize asymmetry in the MGP response david's terrible code####

mgp.mirror.right <- mirror2plane(proj.coords.a1[left,], v1 = proj.coords.a1[31,], v2 = proj.coords.a1[25,], v3 = proj.coords.a1[26,])

plot3d(asym.shape.mean$mesh, color = "lightgrey", alpha = .1, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
mag <- 0
for(i in 1:23) arrow3d(mgp.mirror.right[i,] - (proj.coords.a2[right,][i,] - mgp.mirror.right[i,]) * mag, proj.coords.a2[right,][i,], type = "lines", col = "black", barblen = 0, lwd = 3)

#check that asymmetry differs in magnitude at diff points: rowSums(sqrt((asym.fgf[left,] - mgp.mirror.right)^2))

spheres3d(proj.coords.a2[right,], color =2, radius = .001)
# spheres3d(mgp.mirror.right[,] - (proj.coords.a1[right,][,] - mgp.mirror.right[,]) * mag, color = 4, radius = .001)
spheres3d(mgp.mirror.right[,], color = 4, radius = .001)
# arrow3d(mgp.mirror.right, asym.fgf[right,])

#marta's code for asymmetry visualization####

#### 0. Load R packages ####
library(rgl)
library(Morpho) 
library(geomorph)
library(devtools)
library(Rvcg)
library(magick)
library(Evomorph)
library(ggplot2)
library(ggtittle)
library(vegan)
# install_github("marta-vidalgarcia/morpho.tools.GM", force = TRUE)
library(morpho.tools.GM)
# install_github("marta-vidalgarcia/symmetry", force = TRUE)
library(symmetry)


#### 1. LOAD DATA ####
load("~/Desktop/MGP_FGF_asymmetry.Rdata")

ls()

fgf.mesh # is a mesh already fit to the fgf mutant data
asym.fgf # is the mean of the non-symmetrized fgf mutants
fgf.mirror.right # is the mirrored lm set for plotting
proj.coords.a1 # is the MGP phenotype that I'm trying to visualize
pairedLM # are landmark indices
asym.shape.mean # is the mesh fit to the MGP projected data
mgp.mirror.right # are the mirrored lms for the MGP phenotype


#### 3. PLOTTING ####
right <- pairedLM[,2]
left <- pairedLM[,1]
midline <- setdiff(c(1:dim(proj.coords.a1)[1]), c(left, right))

detect.symmetry(array(data = c(proj.coords.a1, proj.coords.a1), dim = c(dim(proj.coords.a1), 2)), plot = TRUE)
# Currently it just works with arrays, hence the duplicated matrix, and forcing them to an array
# According to symmetry::detect.symmetry the reflection plane is the xz plane, so y = 0
# Trying to find the most suitable landmarks int he midline (the closest to the xz-plane)
midline_y <- proj.coords.a1[midline,2]
names(midline_y) <- midline
sorted_midline <- as.numeric(names(sort(abs(midline_y), decreasing = FALSE)))
# so here I just sorted them in order from smallest (absolute values to biggest)

# Morpho approach with a twist
fgf.mirror.right <- mirror2plane(asym.fgf[left,], v1 = asym.fgf[sorted_midline[1],], 
                                 v2 = asym.fgf[sorted_midline[2],], 
                                 v3 = asym.fgf[sorted_midline[3],])

# open3d(zoom=0.9, windowRect = c(0,0,1000,700))
plot3d(fgf.mesh, color = "lightgrey", alpha = .2, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
mag <- 3
symmetric.fgf <- array.mean(abind::abind(asym.fgf[right,], fgf.mirror.right, along = 3))
for (i in 1:length(right)){
  arrow3d(fgf.mirror.right[i,] , asym.fgf[right,][i,] - (fgf.mirror.right[i,] - asym.fgf[right,][i,]) * mag, 
          type = "lines", col = "black", barblen = 0.005, lwd = 6)
  
  # arrow3d(symmetric.fgf[i,] , fgf.mirror.right[i,] - (symmetric.fgf[i,] - fgf.mirror.right[i,]) * mag, 
  #         type = "lines", col = "red", barblen = 0.005, lwd = 3)
} 

spheres3d(fgf.mirror.right, color = 2, radius = .001)
rglwidget()
# 
# spheres3d(asym.fgf[right,], color =2, radius = .001)
# spheres3d(fgf.mirror.right[,] - (asym.fgf[right,][,] - fgf.mirror.right[,]) * mag, color = 4, radius = .001)

#mgp viz####
mgp.mirror.right <- mirror2plane(proj.coords.a1[left,], v1 = proj.coords.a1[sorted_midline[1],], 
                                 v2 = proj.coords.a1[sorted_midline[2],], v3 = proj.coords.a1[sorted_midline[3],])

# open3d(zoom=0.9, windowRect = c(0,0,1000,700), userMatrix = cor.view)
plot3d(asym.shape.mean$mesh, color = "lightgrey", alpha = .2, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
#there's a magnification that's commented out, but if the asymmetry effects are too small to see, that may be useful
# Sorry David, not sure what you mean here! I changed it to 1 because it looked tiny
mag <- 3
for (i in 1:length(right)){
  arrow3d(mgp.mirror.right[i,], proj.coords.a1[right,][i,] - (mgp.mirror.right[i,] - proj.coords.a1[right,][i,] ) * mag, 
          type = "lines", col = "black", barblen = 0.005, lwd = 6)
} 

# spheres3d(proj.coords.a1[right,], color =2, radius = .001)
# spheres3d(mgp.mirror.right[,] - (proj.coords.a1[right,][,] - mgp.mirror.right[,]) * mag, color = 4, radius = .001)
# spheres3d(mgp.mirror.right[,] - (proj.coords.a1[right,][,] - mgp.mirror.right[,]) * mag, color = 4, radius = .001)
spheres3d(mgp.mirror.right, color =2, radius = .001)
rglwidget()


#what do alk6 mutants look like####
tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk6",])), 54, 3))$coords
tmp.mutant.registration <- gpagen(arrayspecs(rbind(as.matrix(mutant.lms[mutant.db$Genotype == "Alk6" & mutant.db$Experimental_Group == "Wildtype",]), as.matrix(mutant.lms[mutant.db$Genotype == "Alk6" & mutant.db$Experimental_Group == "Homozygote",])), 54, 3))$coords
# tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])
mutant.mean <- array.mean(tmp.mutant.registration[,,-(1)])
alk6.wt <- tmp.mutant.registration[,,1]

#heatmap for alk6 pheno
#make alk6 mesh
alk6.mesh <- rotmesh.onto(shape.mean$mesh, shape.mean$yrot, mutant.mean)

#make alk6 wt mesh
alk6.mesh.wt <- rotmesh.onto(shape.mean$mesh, shape.mean$yrot, alk6.wt)
#meshdist
meshDist(alk6.mesh.wt$mesh, alk6.mesh$mesh)
plot3d(alk6.mesh$mesh, color = 2, lit = F, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
plot3d(alk6.mesh.wt$mesh, color = 4, lit = F, specular = 1, add = T, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")

do.mean <- matrix(colMeans(Y), nrow = 54, ncol = 3, byrow = T)

# shape.mean <- rotmesh.onto(mouse.ply, refmat = as.matrix(consensus.skull), tarmat = do.mean, scale = T, reflection = T)

par3d(zoom = .65)
aspect3d("iso")

#vectors from DO mean to mutant
shape.warp <-  plot3d(shape.mean$mesh, col = adjustcolor("lightgrey", .3), alpha = .1, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
spheres3d(do.mean, radius = .001, color = adjustcolor("red", .3))
bg3d(rgb(245/255,245/255,245/255, .9))

for(i in 1:54) arrow3d(do.mean[i,] - (mutant.mean[i,] - do.mean[i,]), mutant.mean[i,], type = "lines", col = "red", barblen = 0.005, lwd = 2)



#alk6 cell morpho analysis####

cell.morpho <- read.csv("~/Downloads/ALK6_chond_morpho.csv", header = T)
cell.morpho$geno <- as.character(cell.morpho$geno)
cell.morpho$geno[cell.morpho$geno == "Alk6+/-"] <- "Bmpr1b control"
cell.morpho$geno[cell.morpho$geno == "Alk6-/-"] <- "Bmpr1b -/-"
cell.morpho$geno <- factor(cell.morpho$geno)

cell.morpho$Area <- cell.morpho$Area/2.1




#initial visualizations
#ISS
cell.morpho %>% 
  filter(synchondrosis == "ISS", magnification == "high") %>%
  ggplot(aes(x = Distance, y = Area, color = geno, group = geno)) +
  geom_point(aes(size = Area), alpha = .2) + 
  geom_smooth(aes(fill = geno)) + 
  scale_color_manual(values=c("red", "black")) +
  scale_fill_manual(values=c("red", "black")) +
  scale_size(guide=FALSE)

#SOS

cell.morpho %>% 
  filter(synchondrosis == "SOS", magnification == "high") %>%
  ggplot(aes(x = Distance, y = Area, color = geno, group = geno)) +
  geom_point(aes(size = Area), alpha = .2) + 
  geom_smooth(aes(fill = geno)) + 
  scale_color_manual(values=c("red", "black")) +
  scale_fill_manual(values=c("red", "black")) +
  scale_size(guide=FALSE)


#let's try normalizing distance for high and low and then putting the data together
iss.high <- cell.morpho %>% 
  filter(synchondrosis == "ISS", magnification == "high") 
scale(iss.high$Distance)

iss.low <- cell.morpho %>% 
  filter(synchondrosis == "ISS", magnification == "low") 
hist(scale(iss.low$Distance))

iss.combined <- rbind(iss.high, iss.low)
iss.combined$Distance <- c(as.numeric(scale(iss.high$Distance)), as.numeric(scale(iss.low$Distance)))

iss.combined %>% 
  ggplot(aes(x = Distance, y = Area, color = geno, group = geno)) +
  geom_point(aes(size = Area), alpha = .2) + 
  geom_smooth(aes(fill = geno)) + 
  scale_color_manual(values=c("red", "black")) +
  scale_fill_manual(values=c("red", "black")) +
  scale_size(guide=FALSE)

#ok the y values are also really different. let's just work with the high mag for now and move on.

library(lme4)
m0 <- lmer(Area ~ poly(Distance,2) + (1| ID/slide), data = iss.high)
m1 <- lmer(Area ~ poly(Distance,2) * geno + (1| ID/slide), data = iss.high)

anova(m0, m1)
random.effects(m1)
summary(m1)

sos.high <- cell.morpho %>% 
  filter(synchondrosis == "SOS", magnification == "high") 

#remove outlier
sos.high <- sos.high[sos.high$Area < 11000,]

sos0 <- lmer(Area ~ poly(Distance,2) + (1| ID/slide), data = sos.high)
sos1 <- lmer(Area ~ poly(Distance,2) * geno + (1| ID/slide), data = sos.high)

anova(sos0, sos1)

random.effects(sos1)
summary(sos1)

#models without distance interaction
m0 <- lmer(Area ~ (1| ID/slide), data = iss.high)
m1 <- lmer(Area ~ geno + (1| ID/slide), data = iss.high)

anova(m0, m1)

sos0 <- lmer(Area ~  (1| ID/slide), data = sos.high)
sos1 <- lmer(Area ~  geno + (1| ID/slide), data = sos.high)

anova(sos0, sos1)

#let's try plotting the models####
#ISS
pdf("alk6_iss_morpho.pdf", width = 10, height = 7)
cell.morpho %>% 
  filter(synchondrosis == "ISS", magnification == "high") %>%
  ggplot(aes(x = Distance, y = Area, color = geno, group = geno)) +
  geom_point(aes(size = Area), alpha = .2) + 
  stat_smooth(aes(fill = geno), method = "lm", formula = y ~ poly(x, 2)) + 
  scale_color_manual(values=c("red", "black")) +
  scale_fill_manual(values=c("red", "black")) +
  scale_size(guide=FALSE) + 
  xlab("Relative distance from synchondrosis midline") +
  ylab(expression(paste("Estimated cell area (",mu, "m)")))
dev.off()

#SOS

cell.morpho %>% 
  filter(synchondrosis == "SOS", magnification == "high", Area < 11000) %>%
  ggplot(aes(x = Distance, y = Area, color = geno, group = geno)) +
  geom_point(aes(size = Area), alpha = .2) + 
  stat_smooth(aes(fill = geno), method = "lm", formula = y ~ poly(x, 2)) + 
  scale_color_manual(values=c("red", "black")) +
  scale_fill_manual(values=c("red", "black")) +
  scale_size(guide=FALSE)



#alk6 cell morphology plots####

cell.morpho <- read.csv("~/Downloads/alk6_cell_morpho.csv", header = T)
cell.morpho$geno <- as.character(cell.morpho$geno)
cell.morpho$geno[cell.morpho$geno == "Alk6"] <- "Alk6-/-"
cell.morpho$geno[cell.morpho$geno == "Control"] <- "Alk6+/-"
cell.morpho$geno <- factor(cell.morpho$geno)

plot(cell.morpho[,1:2], col = cell.morpho$geno)

ggplot(data = cell.morpho, aes(x = Distance, y = Area, color = geno, group = geno)) +
  geom_point(aes(size = Area), alpha = .2) + 
  geom_smooth(aes(fill = geno)) + 
  scale_color_manual(values=c("red", "black")) +
  scale_fill_manual(values=c("red", "black")) +
  scale_size(guide=FALSE)

#alk6 suppplemental cell morpho figure####
library(plotrix)
#load up coordinates
mut.cells <- read.table("~/Downloads/1860A-8 ko +sut H&E ISS 1.txt", header = F)
#load up corresponding image
cell.image <- imager::load.image("~/Downloads/1860A-8 ko +sut H&E ISS 1.jpg")
#plot coords on image
plot(cell.image)
# points(mut.cells[, 1], -mut.cells[, 2] + dim(cell.image)[2], pch =19, cex = .5)
 for(i in 1:248){ 
  segments(mut.cells[seq(1,248, 2)[i], 1], mut.cells[seq(1,248, 2)[i], 2], mut.cells[seq(2,248, 2)[i], 1], mut.cells[seq(2,248, 2)[i], 2], col = 1, lwd = 2.5)
 
   # draw.ellipse(t(colMeans(mut.cells[1:2,])), a = sqrt((mut.cells[2,1] - mut.cells[2,2])^2)/2, b = sqrt((mut.cells[1,1] - mut.cells[1,2])^2)/2, angle = -15)
  
 }
segments(mut.cells[249, 1], mut.cells[249, 2], mut.cells[250, 1], mut.cells[250, 2], col = "black", lwd = 1.5, lty = 2)
segments(mut.cells[249,1], mut.cells[249,2] - 200, mut.cells[249,1], mut.cells[249,2] + 200, col = "black", lty = 2)
segments(mut.cells[250,1], mut.cells[250,2] - 200, mut.cells[250,1], mut.cells[250,2] + 200, col = "black", lty = 2)
midline <- colMeans(mut.cells[249:250,])
segments(midline[1], midline[2] - 20, midline[1], midline[2] + 20, col = "black", lty = 2)
#is there a custom MGP list that highlight ankrd11 effects? Bone remodeling?####
# GO:0046849

#palate development: GO:0060021
#face morphogenesis: GO:0060325 #comment out rbind with ankrd11 code


process.ano <- "GO:0060021"
coi <- c("ENSEMBL", "SYMBOL")
go2symbol_GO <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = process.ano, columns = coi, keytype = "GO")[,-2:-3]))
coi2 <- c("TXCHROM", "TXSTART", "TXEND")
symbol2info_GO <- AnnotationDbi::select(mmusculusEnsembl, keys = go2symbol_GO[,2], columns = coi2, keytype="GENEID")


selection.vector <- c("Ankrd11")
print(str(selection.vector))
#offline method for getting gene metadata
#pull gene names from process.ano (go terms)
coi <- c("ENSEMBL", "SYMBOL")

gene2symbol <- unique(na.omit(AnnotationDbi::select(org.Mm.eg.db, keys = selection.vector, columns = coi, keytype = "SYMBOL")))

coi2 <- c("TXCHROM", "TXSTART", "TXEND")

symbol2info <- AnnotationDbi::select(mmusculusEnsembl, keys = gene2symbol[,2], columns = coi2, keytype="GENEID")

# combine GO with custom genes
gene2symbol <- c(GO = process.ano, ENSEMBL = gene2symbol[,2], SYMBOL = gene2symbol[,1])
go2symbol <- rbind(gene2symbol, go2symbol_GO)
symbol2info <- rbind(symbol2info, symbol2info_GO)

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
    pls.svd.cv <- perf_mddsPLS(Xs = probs.rows, Y = Y, lambda_min = 0, lambda_max = .2, n_lambda = 20, kfolds = 10, R = 1, mode = "reg", NCORES = 11)
    pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = 1, lambda = .05)

#palate permutation and var explained####
full.pred <- predict(pls.svd, as.matrix(probs.rows))$y
ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y))^2))
rss <- sum(apply(Y, 1, function(x) (x - colMeans(full.pred))^2))
ess/(rss + ess)

hist(lr.perm.r2)

#cache to pls list for new analyses
results <- list(pls.svd, gene.names, seq.info, probs.rows)

tmp.reactive <- results
reactive.svd <- tmp.reactive[[1]]$mod$u[[1]]
gene.names <- tmp.reactive[[2]]
seq.info <- tmp.reactive[[3]]

#now we should be able to take pls.svd directly and maybe label them by founder in a new column, then barplot by family, by gene
do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")

pathway.loadings <- data.frame(gloadings = reactive.svd[,1], gnames = sort(as.character(rep(seq.info[,1], each = 8))), founders = rep(do.names, nrow(seq.info)))

gene2highlight = "Ankrd11"
gname.colors <- rep(1, nrow(seq.info))
gname.colors[which(levels(seq.info[,1]) == gene2highlight)] <- 2



point.size = 11.55
axis.text.size = 13
axis.title.size = 16.5
legend.size = 16.5

ultimate_theme <-   theme(
  text = element_text(size=6), 
  axis.text.x = element_text(angle = 75, hjust = 1, color = gname.colors),
  axis.title.x = element_text(margin = margin(t = 20)),
  axis.text = element_text(angle = 55, hjust = 1, size = axis.text.size),
  axis.title = element_text(size = axis.title.size, face = "bold"),
  legend.text = element_text(size = legend.size), 
  legend.title = element_text(size = legend.size, face = "bold", hjust = .5))

p <- ggplot() +
  geom_bar(data = pathway.loadings, 
           aes(x = gnames, y = gloadings), 
           stat = "identity", 
           width = .75, 
           position=position_dodge()) +
  geom_point(data = pathway.loadings,
             aes(x = gnames, y = gloadings, color = founders),
             shape = "-",
             size = point.size) +
  scale_color_manual(values=do.colors, 
                     guide = guide_legend(title = "Founder\nGenotype", override.aes = list(shape = rep(19, 8), size = 2))) +
  xlab("Gene") +
  ylab("Genetic marker loading") + 
  ultimate_theme

png("~/Desktop/ankrd11_gloadings.png", height = 750, width = 1100)
p
dev.off()



tmp.reactive <- results
probs.rows <- tmp.reactive[[4]]

snp.dim = 1#1

#calculate projection
proj.coords.a1 <- row2array3d(predict(tmp.reactive[[1]], probs.rows[c(which.min(tmp.reactive[[1]]$mod$ts[[1]][,1]), which.max(tmp.reactive[[1]]$mod$ts[[1]][,1])),])$y, Nlandmarks = 54)
proj.coords.a2 <- proj.coords.a1[,,2]
proj.coords.a1 <- proj.coords.a1[,,1]

tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Ankrd11",])), 54, 3))$coords
# tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])

# do.mean <- matrix(colMeans(Y), nrow = 54, ncol = 3, byrow = T)

# shape.mean <- rotmesh.onto(mouse.ply, refmat = as.matrix(consensus.skull), tarmat = do.mean, scale = T, reflection = T)

par3d(zoom = .65)
aspect3d("iso")

#vectors from DO mean to mutant
shape.warp <-  plot3d(shape.mean$mesh, col = "lightgrey", alpha = .2, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
spheres3d(proj.coords.a1, radius = .001, color = adjustcolor("red", .3))
# bg3d(rgb(245/255,245/255,245/255, .9))

for(i in 1:54) arrow3d(do.mean[i,] - (mutant.mean[i,] - do.mean[i,]), mutant.mean[i,], type = "lines", col = "red", barblen = 0.005, lwd = 8)
for(i in 1:54) arrow3d(proj.coords.a1[i,], proj.coords.a2[i,] + (proj.coords.a2[i,] - proj.coords.a1[i,]) * (4 - 1), type = "lines", col = "black", barblen = 0.005, lwd = 10)
rglwidget()

cor(tmp.reactive[[1]]$mod$V_super[,1], manova(two.d.array(tmp.mutant.registration) ~ c(rep(0, nrow(Y)), rep(1, sum(mutant.db$Genotype == "Ankrd11"))))$coef[2,])



palate.pheno <-two.d.array(arrayspecs((rbind(tmp.reactive[[1]]$mod$V_super[,1], tmp.reactive[[1]]$mod$V_super[,1])), 54, 3)[palate.lms,,])[1,]

cor(palate.pheno, manova(two.d.array(tmp.mutant.registration[palate.lms,,]) ~ c(rep(0, nrow(Y)), rep(1, sum(mutant.db$Genotype == "Ankrd11"))))$coef[2,])

#MGP with ankrd11 already annotated####
#face morphogenesis: GO:0060325 #comment out rbind with ankrd11 code
#head morphogenesis: GO:0060323 #not interesting
#bone development: GO:0060348 #ankrd11 isn't a high loader
#odontogenesis:GO:0042475 #not interesting
#skeletal system morpho: GO:0048705 #this may be the winner
#tissue homeostasis: GO:0001894 #not a particularly large effect
#anatomical structure morpho: GO:0009653 # ankr is not on this list in my db
#in utero: GO:0001701 #long list
#multicellular organism growth: GO:0035264 #too many genes, ankr not a big factor

process.ano <- "GO:0048705"
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

probs.rows <- NULL

gene.names <- seq.info[,1]

# Y <- big.Y()[[1]][1:1140,]
#use list of marker names to call on probs and build probs.rows for the custom set
#reformat correct dims of probs.rows
for(i in 1: dim(seq.indexes)[1]) probs.rows <- cbind(probs.rows, DO.probs[,,dimnames(DO.probs)[[3]] == seq.indexes[i,2]])
#fit pls2B, need duv, gene names, seq.indexes
if(cv){
pls.svd.cv <- perf_mddsPLS(Xs = probs.rows, Y = Y, lambda_min = 0, lambda_max = .2, n_lambda = 20, kfolds = 10, R = 1, mode = "reg", NCORES = 11)
pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = 1, lambda = pls.svd.cv$Optim$optim_para_all$Lambdas)
} else {pls.svd <- mddsPLS(Xs = probs.rows, Y = Y, R = 1, lambda = .06)
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

pathway.loadings <- data.frame(gloadings = reactive.svd[,1], gnames = sort(as.character(rep(seq.info[,1], each = 8))), founders = rep(do.names, nrow(seq.info)))


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
        )

ggplotly(p +
           theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                 axis.text.y = element_text(size = 8),
                 axis.title = element_text(size = 12, face = "bold")),
         legend.key = element_rect(fill = rgb(245/255, 245/255, 245/255, .9))) %>% layout(
           margin = list(b = 100, l = 50) # to fully display the x and y axis labels
         )


tmp.reactive <- results
probs.rows <- tmp.reactive[[4]]

snp.dim = 1#1

#calculate projection
proj.coords.a1 <- row2array3d(predict(tmp.reactive[[1]], probs.rows[c(which.min(tmp.reactive[[1]]$mod$ts[[1]][,1]), which.max(tmp.reactive[[1]]$mod$ts[[1]][,1])),])$y, Nlandmarks = 54)
proj.coords.a2 <- proj.coords.a1[,,2]
proj.coords.a1 <- proj.coords.a1[,,1]

tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Ankrd11",])), 54, 3))$coords
# tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])

do.mean <- matrix(colMeans(Y), nrow = 54, ncol = 3, byrow = T)

# shape.mean <- rotmesh.onto(mouse.ply, refmat = as.matrix(consensus.skull), tarmat = do.mean, scale = T, reflection = T)

par3d(zoom = .65)
aspect3d("iso")

#vectors from DO mean to mutant
plot3d(shape.mean$mesh, col = "black", alpha = .1, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
spheres3d(do.mean, radius = .001, color = adjustcolor("red", .3))
text3d(do.mean, texts = 1:54, color = adjustcolor("red", .3))
# bg3d(rgb(245/255,245/255,245/255, .9))

for(i in 1:54) arrow3d(do.mean[i,] - (mutant.mean[i,] - do.mean[i,]), mutant.mean[i,], type = "lines", col = "red", barblen = 0.005, lwd = 2)
for(i in 1:54) arrow3d(proj.coords.a1[i,], proj.coords.a2[i,] + (proj.coords.a2[i,] - proj.coords.a1[i,]) * (4 - 1), type = "lines", col = "black", barblen = 0.005, lwd = 3)

cor(tmp.reactive[[1]]$mod$V_super[,1], manova(two.d.array(tmp.mutant.registration) ~ c(rep(0, nrow(Y)), rep(1, sum(mutant.db$Genotype == "Ankrd11"))))$coef[2,])


# let's do palate development on just palate lms/no ankrd11 added
#palate.lms <- c(9,40,8,39,21,52,22,53, 23, 54)

Y.palate <- array2row3d(row2array3d(Y)[palate.lms,,])

palate <- mgp("GO:0060021", "Ankrd11", Y)


#updated chond diff fig####
chond <- mgp("GO:0002062", "Alk6", Y, cv = F)

png("~/Downloads/chond_gloadings.png", height = 750, width = 1100)
loading.plot(chond, point.size = 23, axis.text.size = 16.5, legend.size = 18.5, gene2highlight = "Bmpr1b")
dev.off()

proj.coords.a1 <- row2array3d(predict(chond$mgp, chond$probs.rows[c(which.min(chond$mgp$mod$ts[[1]][,1]), which.max(chond$mgp$mod$ts[[1]][,1])),])$y, Nlandmarks = 54)
proj.coords.a2 <- 1 * proj.coords.a1[,,2]
proj.coords.a1 <- 1 *proj.coords.a1[,,1]

tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk6",])), 54, 3))$coords
# tmp.mutant.registration <- gpagen(arrayspecs(rbind(Y, as.matrix(mutant.lms[mutant.db$Genotype == "Alk2",])), 54, 3))$coords
do.mean <- array.mean(tmp.mutant.registration[,,1:nrow(Y)])
mutant.mean <- array.mean(tmp.mutant.registration[,,-(1:nrow(Y))])

do.mean <- matrix(colMeans(Y), nrow = 54, ncol = 3, byrow = T)

# shape.mean <- rotmesh.onto(mouse.ply, refmat = as.matrix(consensus.skull), tarmat = do.mean, scale = T, reflection = T)

par3d(zoom = .65)
aspect3d("iso")

#vectors from DO mean to mutant
plot3d(shape.mean$mesh, col = "lightgrey", alpha = .3, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
spheres3d(do.mean, radius = .001, color = adjustcolor("red", .3))

# bg3d(rgb(245/255,245/255,245/255, .9))

for(i in 1:54) arrow3d(do.mean[i,] - (mutant.mean[i,] - do.mean[i,]), mutant.mean[i,], type = "lines", col = "red", barblen = 0.005, lwd = 2)
flip.vector <- proj.coords.a1[,] - (proj.coords.a2[,] - proj.coords.a1[,])
vectormag <- proj.coords.a2[,] + (proj.coords.a2[,] - proj.coords.a1[,]) * (2 - 1)
for(i in 1:54) arrow3d(proj.coords.a1[i,] + (proj.coords.a1[i,] - flip.vector[i,]) * (2 - 1), flip.vector[i,], type = "lines", col = "black", barblen = 0.005, lwd = 3)

rglwidget()


#new analyses####
chond <- mgp("GO:0002062", "Alk6", Y, cv = T)

chond$gene.plot

#sort gene loadings by size
do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")

seq.names <- factor(chond$seq.info[,1], levels = chond$seq.info[,1])

pathway.loadings <- data.frame(gloadings = chond$mgp$mod$u[[1]][,1], gnames = (rep(seq.names, each = 8)), founders = rep(do.names, nrow(chond$seq.info)))

ggplot() +
  geom_bar(data = pathway.loadings, 
           aes(x = gnames, y = gloadings), #reorder(gnames, -gloadings)
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
  )


#permutation tests for process examples####
#chond

#symmetry
lr.perm.r2 <- rep(NA, 10000)
for(i in 1:10000){
  process.svd <- mddsPLS(Xs = as.matrix(lr.symm$probs.rows), Y = as.matrix(Y[sample(1:nrow(Y), size = nrow(Y)),]), R = 1, lambda = .04, getVariances = T)
  
  full.pred <- predict(process.svd, as.matrix(lr.symm$probs.rows))$y
  ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y))^2))
  rss <- sum(apply(Y, 1, function(x) (x - colMeans(full.pred))^2))
  lr.perm.r2[i] <- ess/(rss + ess)
  print(i)
}

#palate




pal.perm.r2 <- rep(NA, 10000)
for(i in 1:10000){
  process.svd <- mddsPLS(Xs = as.matrix(probs.rows), Y = as.matrix(Y[sample(1:nrow(Y), size = nrow(Y)),]), R = 1, lambda = .05)
  
  full.pred <- predict(process.svd, as.matrix(probs.rows))$y
  ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y))^2))
  rss <- sum(apply(Y, 1, function(x) (x - colMeans(full.pred))^2))
  pal.perm.r2[i] <- ess/(rss + ess)
  print(i)
}


chond.perm.r2 <- rep(NA, 10000)
for(i in 1:10000){
  process.svd <- mddsPLS(Xs = as.matrix(probs.rows), Y = as.matrix(Y[sample(1:nrow(Y), size = nrow(Y)),]), R = 1, lambda = .08)
  
  full.pred <- predict(process.svd, as.matrix(probs.rows))$y
  ess <- sum(apply(full.pred, 1, function(x) (x - colMeans(Y))^2))
  rss <- sum(apply(Y, 1, function(x) (x - colMeans(full.pred))^2))
  chond.perm.r2[i] <- ess/(rss + ess)
  print(i)
}




