#alk6 cell morpho analysis####
library(dplyr)
library(ggplot2)

cell.morpho <- read.csv("~/Downloads/alk6_chondrocyte_morpho.csv", header = T)
cell.morpho$geno <- factor(cell.morpho$geno)
cell.morpho$Area <- cell.morpho$Area/2.1
cell.morpho$Norm_X <- 2*cell.morpho$Norm_X

#try renormalizing to -1,1 comment out if you don't want to change that scale!
cell.morpho$Norm_X <- scale(cell.morpho$Norm_X, center = 1, scale = F)



#initial visualizations fitting a GAM####
#ISS
cell.morpho %>% 
  filter(magnification == "high") %>%
  ggplot(aes(x = Norm_X, y = Area, color = geno, group = geno)) +
  geom_point(aes(size = Area), alpha = .2) + 
  geom_smooth(aes(fill = geno)) + 
  scale_color_manual(values=c("black", "red")) +
  scale_fill_manual(values=c("black", "red")) +
  scale_size(guide=FALSE)


#subset just high mag data
iss.high <- cell.morpho %>% 
  filter(magnification == "high") 

library(lme4)
m0 <- lmer(Area ~ poly(Norm_X,2) + (1| ID/SLIDE), data = iss.high)
m1 <- lmer(Area ~ poly(Norm_X,2) * geno + (1| ID/SLIDE), data = iss.high)

anova(m0, m1)
random.effects(m1)
summary(m1)

#models without distance interaction####
m0 <- lmer(Area ~ Norm_X + (1| ID/SLIDE), data = iss.high)
m1 <- lmer(Area ~ Norm_X + geno + Norm_X:geno + (1| ID/SLIDE), data = iss.high)

anova(m0, m1)

#models without distance####
m0 <- lmer(Area ~ (1| ID/SLIDE), data = iss.high)
m1 <- lmer(Area ~ geno + (1| ID/SLIDE), data = iss.high)

anova(m0, m1)

#let's try plotting the models####
#ISS
pdf("alk6_iss_morpho.pdf", width = 10, height = 7)
cell.morpho %>% 
  filter(magnification == "high") %>%
  ggplot(aes(x = Norm_X, y = Area, color = geno, group = geno)) +
  geom_point(aes(size = Area), alpha = .2) + 
  stat_smooth(aes(fill = geno), method = "lm", formula = y ~ poly(x, 2)) + 
  scale_color_manual(values=c("black", "red")) +
  scale_fill_manual(values=c("black", "red")) +
  scale_size(guide=FALSE) + 
  xlab("Relative distance from synchondrosis midline") +
  ylab(expression(paste("Estimated cell area (",mu, "m)")))
dev.off()
