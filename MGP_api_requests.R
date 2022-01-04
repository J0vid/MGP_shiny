library(Morpho)
library(rgl)
library(ggplot2)

genelist <- c("Shh, Hhat, Disp1, Disp2, Ptch1, Ptch2, Gas1, Gas2, Gas3, Gli1, Gli2, Lrp2, HhrpGli2,  Gpc3, Gli3, Smo, Cdon")

raw_api_res <- httr::GET(url = paste0("https://genopheno.ucalgary.ca/api", "/custom_mgp"),
                           query = list(genelist = genelist),
                           encode = "json")

MGP_result <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))

#have a look at the structure of the reponse. $loadings helps with plotting genetic effects, $pheno1 and $pheno2 show the extremes of the PLS axis. $pheno_loadings will help if you want to make a vector correlation to a mutant
str(MGP_result)

#plotting marker effects####
do.names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
do.colors <- c("A/J" = "#F0F000","C57BL/6J" = "#808080", "129S1/SvImJ"= "#F08080", "NOD/ShiLtJ" = "#1010F0","NZO/HlLtJ" = "#00A0F0","CAST/EiJ" = "#00A000", "PWK/PhJ" = "#F00000", "WSB/EiJ" = "#9000E0")

ggplot() +
geom_bar(data = MGP_result$loadings,
         aes(x = gnames, y = gloadings),
         stat = "identity",
         width = .75,
         position=position_dodge()) +
  geom_point(data = MGP_result$loadings,
             aes(x = gnames, y = gloadings, color = founders),
             shape = "-",
             size = 15) +
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
        legend.title = element_text(size = 8, face = "bold", hjust = .5))

#plotting phenotypic effects####
pheno_mag <- 4
par3d(zoom = .65)

#vectors from min PLS score to max
plot3d(MGP_result$pheno1, radius = .005, typ = "s", col = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
segments3d(x = rbind(MGP_result$pheno1[,1], MGP_result$pheno2[,1]  + ((MGP_result$pheno2[,1] - MGP_result$pheno1[,1]) * (pheno_mag - 1))),
           y = rbind(MGP_result$pheno1[,2], MGP_result$pheno2[,2]  + ((MGP_result$pheno2[,2] - MGP_result$pheno1[,2]) * (pheno_mag - 1))),
           z = rbind(MGP_result$pheno1[,3], MGP_result$pheno2[,3]  + ((MGP_result$pheno2[,3] - MGP_result$pheno1[,3]) * (pheno_mag - 1))),
           lwd = 3)

#swap in your mesh wirh landmarks scaled to the space
# shape.warp <-  plot3d(mesh, col = adjustcolor("lightgrey", .3), alpha = .2, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")
# spheres3d(MGP_result$pheno1, radius = .003, color = 1)
# segments3d(x = rbind(MGP_result$pheno1[,1], MGP_result$pheno2[,1]  + ((MGP_result$pheno2[,1] - MGP_result$pheno1[,1]) * (pheno_mag - 1))),
#            y = rbind(MGP_result$pheno1[,2], MGP_result$pheno2[,2]  + ((MGP_result$pheno2[,2] - MGP_result$pheno1[,2]) * (pheno_mag - 1))),
#            z = rbind(MGP_result$pheno1[,3], MGP_result$pheno2[,3]  + ((MGP_result$pheno2[,3] - MGP_result$pheno1[,3]) * (pheno_mag - 1))),
#            lwd = 3)

#generate phenotype morphs####
# my_morph <- tps3d(mesh, mesh_lms, MGP_result$pheno1)
# 
# plot3d(my_morph, col = adjustcolor("lightgrey", .3), alpha = .9, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso")

#elife style ordered plots####
pathway.loadings <- MGP_result$loadings
bar_order <- pathway.loadings %>% 
  group_by(gnames) %>%
  summarise(test = diff(range(gloadings))) %>%
  arrange(-test) 

pathway.loadings$gnames <- factor(pathway.loadings$gnames, levels = lapply(bar_order, as.character)$gnames)


ggplot() +
  geom_bar(data = pathway.loadings,
           aes(x = gnames, y = gloadings),
           stat = "identity",
           width = .75,
           position=position_dodge()) +
  geom_point(data = pathway.loadings,
             aes(x = gnames, y = gloadings, color = founders),
             shape = "-",
             size = 15) +
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
        legend.title = element_text(size = 8, face = "bold", hjust = .5))



