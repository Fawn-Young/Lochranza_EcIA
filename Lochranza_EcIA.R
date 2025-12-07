rm(list=ls())
library(dplyr)
library(stringr)
library(ggplot2)
library(vegan)
#Read in the combined associated occurrences and terrestrial data from the Scottish Biodiversity List
allassocc<- read.csv('data/All_occurrences.csv', header=T, sep=',')
terrestrialspecies<- read.csv('data/terrestrial species.csv', header=T)

#Attempt to get rid of spacing errors by trimming white space
allassocc$scientificName<- trimws(allassocc$scientificName)

#Using the terrestrial data, filter based on threatened status and cross check it with the associated occurrences
terrestrialspeciesred<- filter(terrestrialspecies, terrestrialspecies$Threatened.species == 'Red')
terrestrialspeciesamber<- filter(terrestrialspecies, terrestrialspecies$Threatened.species == 'Amber')

allassoccfilt<- filter(allassocc, allassocc$scientificName %in% terrestrialspeciesred$Scientific.Name | 
                         allassocc$scientificName %in% terrestrialspeciesamber$Scientific.Name)
allassoccfilt$hill<-ifelse(str_detect(allassoccfilt$eventID, regex('N')), 'North', 'South')

#Plot abundance data of these threatened species
theme_set(theme_grey())
gr<- ggplot(data=allassoccfilt, aes(x=Commonname, y=individualCount, fill=hill))+ geom_col()+
  labs(x='Common Name', y='Abundance')+
  scale_fill_manual(values=c('purple4','yellow4'))
gr

## ANALYSIS OF TAXA GROUPS ##

#load in the required datasets
birds<- read.csv('data/Birds.csv', header = T)
fwinverts<- read.csv('data/freshwaterinverts.csv', header=T)
bats<- read.csv('data/Batdetect.csv', header = T)
terinverts<- read.csv('data/Terrestrialinverts.csv', header=T)
mam<- read.csv('data/Lochranza_mammals.csv', header=T)
moth<- read.csv('data/Mothtraps.csv', header=T)
plants<- read.csv('data/Plants.csv', header=T)

##Richness Based Analysis

#Firstly, due to the data collected, freshwater invertebrates and mammals can be 
#analysed using richness

#Freshwater Invertebrates
#as most individuals were only identified to scientificName, using scientificName we can get the gamma diversity of the area
fwinvertsgamma<- length(unique(fwinverts$Order))
fwinvertsgamma

#from looking at the data there is a spelling mistake, correct it
fwinverts[31,6] = 'Ephemeroptera'
fwinverts[13,16] = '1'
fwinverts[14,16] = '5, 10'
fwinverts[15,16] = '5, 10'
#Separate them into separate hills to make alpha and beta richness easier
fwinverts$hill<- as.factor(ifelse(str_detect(fwinverts$eventID, regex('N')),'North', 'South'))
Nfwinverts<- filter(fwinverts, fwinverts$hill=='North')
Sfwinverts<- filter(fwinverts, fwinverts$hill=='South')

#For alpha diversities
Nfwinvertsalpha<- length(unique(Nfwinverts$Order))
Sfwinvertsalpha<- length(unique(Sfwinverts$Order))

#Beta diversity, using the Jaccard index
fwbeta<- length(intersect(Nfwinverts$scientificName, Sfwinverts$scientificName))/length(unique(c(Nfwinverts$scientificName, Sfwinverts$scientificName)))
fwbeta

ggplot(fwinverts, aes(x=hill, y=Order, color=BMWPscore, size = 4)) +geom_point()+
  labs(x= 'Hill', y='Order', color = 'ASPT') + theme_gray()

#Now terrestrial mammals#Now terOrderrestrial mammals
#Now we see an issue with the data as it counts "no records" as a different species
mam<-filter(mam, mam$scientificName != 'No records')
mamgamma<- length(unique(mam$scientificName))
mamgamma

#Separate them into separate hills to make alpha and beta richness easier
mam$hill<- as.factor(ifelse(str_detect(mam$eventID, regex('N')),'North', 'South'))
Nmam<- filter(mam, mam$hill=='North')
Smam<- filter(mam, mam$hill=='South')

#For alpha diversities
Nmamalpha<- length(unique(Nmam$scientificName))
Smamalpha<- length(unique(Smam$scientificName))

#Beta diversity, using the Jaccard index
mambeta<- length(intersect(Nmam$scientificName, Smam$scientificName))/length(unique(c(Nmam$scientificName, Smam$scientificName)))
mambeta

ggplot(mam, aes(x=hill, y=scientificName, color= Commonname,size = 4)) +geom_point()+
  scale_color_manual(values=c('black', 'blue3', 'orange', 'red3'))+
  labs(x= 'Hill', y='Species', color= 'Common name', title='(Fig.6) Mammal Presence') + theme_gray()

#Abundance based
#In order to show alpha diversity 
moth$hill<- as.factor(ifelse(str_detect(moth$eventID, regex('N')),'North', 'South'))
moth[42,6] = 'Lepidoptera'
moth<-filter(moth, moth$Order=='Lepidoptera ')
moth[18,6] = 'Araneae'
moth[24,6] = 'Diptera'
moth<-filter(moth, moth$Order=='Lepidoptera ')


Nmoth<- filter(moth, moth$hill=='North')
Smoth<- filter(moth, moth$hill=='South')

Nmothdiv<-diversity(Nmoth$organismQuantity, index='shannon')
Smothdiv<-diversity(Smoth$organismQuantity, index='shannon')

par(mfrow=c(2,1))
ggplot(Nmoth, aes(y= scientificName, x=individualCount, fill=Commonname)) +geom_col()+
  labs(x='Abundance', y='Species', fill = 'Common Name')+
  ggtitle('(Fig.8a) North hill')+
  xlim(0,15)
ggplot(Smoth, aes(x= individualCount, y=scientificName, fill=Commonname)) +geom_col()+
  labs(x='Abundance', y='Species', fill='Common Name')+
  ggtitle('(Fig.8b) South hill')+
  xlim(0,15)

#gamma diversity 
mothgamma<- length(unique(moth$scientificName))

#For non moth invertebrates
#reset the moth data
moth<- read.csv('data/Mothtraps.csv', header=T)
moth[42,6] = 'Lepidoptera '
moth[23,6] = 'Araneae'
moth[38,6] = 'Diptera'
moth$hill<- as.factor(ifelse(str_detect(moth$eventID, regex('N')),'North', 'South'))
nonmoth<-filter(moth, moth$Order!='Lepidoptera ')
Nnonmoth<- filter(nonmoth, nonmoth$hill=='North')
Snonmoth<- filter(nonmoth, nonmoth$hill=='South')

Nnonmothdiv<-diversity(Nnonmoth$organismQuantity, index='shannon')
Snonmothdiv<-diversity(Snonmoth$organismQuantity, index='shannon')

par(mfrow=c(1,1))
ggplot(Nnonmoth, aes(y= scientificName, x=individualCount, fill=Commonname)) +geom_col()+
  labs(x='Abundance', y='Species', fill = 'Common Name')+
  ggtitle('(Fig.8c) North hill')+
  xlim(0,15)
ggplot(Snonmoth, aes(x= individualCount, y=scientificName, fill=Commonname)) +geom_col()+
  labs(x='Abundance', y='Species', fill='Common Name')+
  ggtitle('(Fig.8d) South hill')+
  xlim(0,15)

#For Birds
birds$hill<- as.factor(ifelse(str_detect(birds$eventID, regex('N')),'North', 'South'))
Nbirds<- filter(birds, birds$hill=='North')
Sbirds<- filter(birds, birds$hill=='South')

Nbirdsdiv<-diversity(Nbirds$organismQuantity, index='shannon')
Sbirdsdiv<-diversity(Sbirds$organismQuantity, index='shannon')

par(mfrow=c(2,1))
ggplot(Nbirds, aes(y= scientificName, x=individualCount, fill=Commonname)) +geom_col()+
  labs(x='Abundance', y='Species', fill = 'Common Name')+
  ggtitle('(Fig.5a) North hill')+
  xlim(0,15)
ggplot(Sbirds, aes(x= individualCount, y=scientificName, fill=Commonname)) +geom_col()+
  labs(x='Abundance', y='Species', fill='Common Name')+
  ggtitle('(Fig.5b) South hill')+
  xlim(0,15)

#gamma diversity 
birdsgamma<- length(unique(birds$scientificName))

#For terrestrial invertebrates
terinverts$hill<- as.factor(ifelse(str_detect(terinverts$eventID, regex('S')),'South', 'North'))
Nterinverts<- filter(terinverts, terinverts$hill=='North')
Sterinverts<- filter(terinverts, terinverts$hill=='South')

Nterinvertsdiv<-diversity(Nterinverts$organismQuantity, index='shannon')
Sterinvertsdiv<-diversity(Sterinverts$organismQuantity, index='shannon')

par(mfrow=c(2,1))
ggplot(Nterinverts, aes(y= scientificName, x=individualCount, fill=Commonname)) +geom_col()+
  labs(x='Abundance', y='Species', fill = 'Common Name')+
  ggtitle('(Fig.7a) North hill')+
  xlim(0,15)
ggplot(Sterinverts, aes(x= individualCount, y=scientificName, fill=Commonname)) +geom_col()+
  labs(x='Abundance', y='Species', fill='Common Name')+
  ggtitle('(Fig.7b) South hill')+
  xlim(0,40)

#gamma diversity 
terinvertsgamma<- length(unique(terinverts$scientificName))

#For Plants
plants$hill<- as.factor(ifelse(str_detect(plants$eventID, regex('N')),'North', 'South'))
Nplants<- filter(plants, plants$hill=='North')
Splants<- filter(plants, plants$hill=='South')

Nplantsdiv<-diversity(Nplants$organismQuantity, index='shannon')
Splantsdiv<-diversity(Splants$organismQuantity, index='shannon')

par(mfrow=c(2,1))
#As this is plants, we use organism quantity (percentage cover) instead of count 
ggplot(Nplants, aes(y= scientificName, x=organismQuantity, color=Commonname, size = 3)) +geom_point()+
  labs(x='Percentage Cover', y='Species', color = 'Common Name')+
  ggtitle('(Fig.10a) North hill')+
  xlim(0,100)+
  ylim(unique(Nplants$scientificName))
ggplot(Splants, aes(y= scientificName, x=organismQuantity, color=Commonname, size = 3)) +geom_point()+
  labs(x='Percentage Cover', y='Species', color = 'Common Name')+
  ggtitle('(Fig.10b) South hill')+
  xlim(0,100)+
  ylim(unique(Splants$scientificName))

#gamma diversity 
plantsgamma<- length(unique(plants$scientificName))

#For Bats
bats$hill<- as.factor(ifelse(str_detect(bats$eventID, regex('N')),'North', 'South'))
Nbats<- filter(bats, bats$hill=='North')
Sbats<- filter(bats, bats$hill=='South')

Nbatsdiv<-diversity(Nbats$individualCount, index='shannon')
Sbatsdiv<-diversity(Sbats$individualCount, index='shannon')

par(mfrow=c(2,1))

#gamma diversity 
batsgamma<- length(unique(bats$scientificName))

