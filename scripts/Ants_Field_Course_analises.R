library(vegan)
library(dplyr)
library(hillR)
library(GGally)
library(vegan)
library(iNEXT)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(lme4)
library(dplyr)

setwd("C:/Users/boloq/Box/Dissertation/FCAT/FCAT_data")



ants_data <- read.table("Ants_students_data_june2023.txt", h=T, stringsAsFactors = T)
head(ants_data) #Primeiras linhas da tabela
summary(ants_data) #Sum?rio da tabela
str(ants_data) #Estrutura da tabela
names(ants_data) #Nomes das colunas da tabela

#Preparing the data frame 

#change the table form with species in columns
library(reshape)
ants_data_com <- dcast(ants_data, Site+Transect~Species, value.var = "Species")

# We going to  use the bait as a repetition so we going to creat a code of Transect+bait in a new colun

ants_data$Transect_bait <- paste(ants_data$Transect,ants_data$Bait, sep="_")


#change the table form with species in columns
library(reshape)
ants_data_com_Transect_bait <- dcast(ants_data, Site+Transect_bait~Species, value.var = "Species")

#Just for transects 
ants_data_com <- dcast(ants_data, Site+Transect~Species, value.var = "Species")

summary(ants_data_com)

#creatging a object only with the species 
com <- ants_data_com[,-c(1,2)] 


# Amostral coverge ====

# here we want to know how well we sample this comunity so we are more interested in the sites to do that  we creat a data.frame only with sites and with species in the colomns
library(reshape2)
ants_data_com_site <-dcast(ants_data, Site~Species, value.var = "Species")
summary(ants_data_com_site)



#creatging a object only with the species 
com.site <- ants_data_com_site[,-c(1)] 
com.site



# because ants are base on frequency we need to inclued the numer of sample stations we had
# inclucing a colomn with the number of samples
summary(ants_data_com_Transect_bait$Site)

com.site$n_samples<- c(14,16)

#give name to rows
row.names(com.site) <- ants_data_com_site[,1]

# Reorder to have n_samples in the first coloum
ncol(com.site)
com.site <- com.site[, c(41,1:40)]
com.site

#We going to user the inext to do the sample corve and rarefaction
library(iNEXT)
out1 <- iNEXT(t(com.site), q=0, datatype="incidence_freq", endpoint = 40)
out1$DataInfo
out1$DataInfo$SC
#Number of species
(plot_sp <- ggiNEXT(out1, type=1, se=TRUE, grey=FALSE) 
  + theme_classic(base_size = 18) 
  + theme(legend.position="right")
  + labs(y="Number of species", x = "Number of Samples"))+
  scale_fill_discrete(labels = c("Forest", "Reforestation"))+
  scale_shape_discrete(labels = c("Forest", "Reforestation"))+
  scale_color_discrete(labels = c("Forest", "Reforestation"))+ 
  scale_colour_manual(values=c('#1B9E77', '#D95F02')) +
  scale_fill_manual(values=c('#1B9E77', '#D95F02'))


#sample coverge 
(plot_C_amnostral <- ggiNEXT(out1, type=2, se=TRUE, grey=FALSE) 
  + theme_classic(base_size = 18) 
  + theme(legend.position="right")
  + labs(y="Sample coverage", x = "Number of Samples"))+
  scale_fill_discrete(labels = c("Forest", "Reforestation"))+
  scale_shape_discrete(labels = c("Forest", "Reforestation"))+
  scale_color_discrete(labels = c("Forest", "Reforestation"))+ 
  scale_colour_manual(values=c('#1B9E77', '#D95F02')) +
  scale_fill_manual(values=c('#1B9E77', '#D95F02')) 



#_diversity metrics====
#_N?meros de Hill====

# species richness
(Richness <- sp_rich <- specnumber(com))
(tab_result <- data.frame(ants_data_com[,1:2], Richness))
aggregate(tab_result$Richness, list(tab_result$Site), FUN=mean)
#Shannon diversity
(tab_result$Shannon <- diversity(x = com, index = "shannon"))
aggregate(tab_result$Shannon, list(tab_result$Site), FUN=mean)
# Simpson diversity
(tab_result$Simpson <- diversity(x = com, index = "simpson"))
aggregate(tab_result$Simpson, list(tab_result$Site), FUN=mean)

##_Plot====
library(reshape2)
summary(tab_result)
tab_ggplot <- melt(tab_result, id= c("Site","Transect")) # crarte a table for the plot

cor_plot <- c('#1B9E77', '#D95F02') #choose colours

library(ggplot2)
library(ggpubr)
ggplot(tab_ggplot, aes(x=Site, y=value, fill=Site)) +
  geom_boxplot()+
  theme_bw(base_size = 15)+
  scale_fill_manual(values = cor_plot)+ 
  facet_wrap(~variable, scales = "free")+
  theme(legend.position =  "none")+ 
  labs(y="Diversity", x = "")

# models
#Richens
mod1 <- glm(Richness~Site, poisson, data=tab_result)
(s <- summary(mod1))
anova(mod1, test="F")


#Shannon
mod2 <- lm(Shannon~Site, data=tab_result)
(s <- summary(mod2))
anova(mod2, test="F")


#Simpson
mod2 <- lm(Simpson~Site, data=tab_result)
(s <- summary(mod2))
anova(mod2, test="F")


#### BETA DIVERSITY 
library(vegan) 
library(betapart) 

# Turnover and Nestedness

spp.pa <- decostand(com, method="pa")  # mudar para presencia e ausencia 
(jac <- beta.multi(spp.pa, index.family="jaccard")) # index.family="sorensen" - 
round(jac$beta.JTU / jac$beta.JAC, 3) # 96,3 of turnover/replacement of spp
round(jac$beta.JNE / jac$beta.JAC, 3) # 3.7% of nestedness of spp
1 - jac$beta.JAC # 6,7% similarity




#################____________#############___________############
#______________________________NMDS_________________________________#
#___________________________PERMANOVA_______________________________#
################____________#############___________############

ants_data_nmds <- dcast(ants_data, Site+Transect_bait~Species, value.var = "Species")

summary(ants_data_nmds)

#creatging a object only with the species 
com.nmds <- ants_data_nmds[,-c(1,2)]


library(vegan)
nmds <- metaMDS(com.nmds,  distance = "jaccard", k = 2, trymax = 999, trace = F)

nmds

summary(nmds)

par(mar=c(6,6,2,2))

#Pacote de paleta de cores 
library(RColorBrewer)
display.brewer.all()

#
cores <- brewer.pal(7, "Dark2")

plot(nmds, axes=F, cex.lab=3, cex.axis=3, xlab="", ylab="")
axis(side=1, lwd=4,cex.axis=1.3 )
axis(side=2, lwd=4, las=1.5, cex.axis=1.3)
box(lwd=3, bty="l")
mtext("NMDS1", side=1, line=3, cex=1.5)
mtext("NMDS2", side=2, line=3, cex=1.5)
mypch <- c(15,16,17,18,19,20,21)
mycol <- c(cores)
points(nmds, pch=mypch[as.numeric(ants_data_com$Site)], col=mycol[as.numeric(ants_data_com$Site)], cex=2)


ants_data_com 


#Plotando o elipse de tendencia 
summary(ants_data_com$Site)

treat=c(rep("Forest",14),rep("Reforestation",16))
colors<- c(rep(cores[1],14),rep(cores[2],16))
#Plot convex hulls with colors baesd on treatment
for(i in unique(treat)) {
  ordiellipse(nmds$point[grep(i,treat),],draw="polygon",
              groups=treat[treat==i],col=colors[grep(i,treat)],label=F) }

###legand

legend(3, 2 , c("Forest","Reforestation"), pch=c(15,16), col=cores, cex=1.3, bty="n")
#text(4, 2.05,"Stress = 0.18")


#ANOSIM - o anossim faz a mesma coisa mas calcuala a similaridade ao invez de analise de variancia
anosim(com.nmds, ants_data_nmds$Site, distance="bray",permutations = 1000 )
