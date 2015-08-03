# Read data and merge using BiomasaFP and function mergefp
# Note, the new veresion of FP will have a different format for the advanced search downloads

library (BiomasaFP)
Heights <- mergefp('IndvDataV10.csv', 'MetadataR10A.csv','WDHeightsV10.csv')
head (Heights)

# This plot type file probably needs to be incorporated with the metadata for R
# The plot type file includes the 3 ForestStatus categories needed for the height protocol


PlotType <- read.csv('PlotTypeV10A.csv', header = TRUE, stringsAsFactors = FALSE) 

names(PlotType)<- c('PlotID','PlotCode_B','ClusterID','ForestMoistureID','ForestMoistureName',
                      'ForestEdaphicID','ForestEdaphicName', 'EdaphicHeight',	
                      'EdaphicHeightCode','ForestElevationID', 'ForestElevationName',
                      'ElevationHeight','ElevationHeightCode',	'BiogeographicalRegionID',
                       'BiogeographicalRegionName', 'ContinentID','CountryID')

head(PlotType)
nrow(PlotType)

PlotType$ClusterID<- as.numeric(gsub("=","",PlotType$ClusterID))
PlotType$ClusterID<- as.numeric(gsub("NULL","NA",PlotType$ClusterID))

PlotType$ForestMoistureID<- as.numeric(gsub("=","",PlotType$ForestMoistureID))
PlotType$ForestMoistureID<- as.numeric(gsub("NULL","NA",PlotType$ForestMoistureID))

PlotType$ForestEdaphicID<- as.numeric(gsub("=","",PlotType$ForestEdaphicID))
PlotType$ForestEdaphicID<- as.numeric(gsub("NULL","NA",PlotType$ForestEdaphicID))

PlotType$EdaphicHeightCode<- as.numeric(gsub("=","",PlotType$EdaphicHeightCode))
PlotType$EdaphicHeightCode<- as.numeric(gsub("NULL","NA",PlotType$EdaphicHeightCode))

PlotType$ForestElevationID<- as.numeric(gsub("=","",PlotType$ForestElevationID))
PlotType$ForestElevationID<- as.numeric(gsub("NULL","NA",PlotType$ForestElevationID))

PlotType$ElevationHeightCode<- as.numeric(gsub("=","",PlotType$ElevationHeightCode))
PlotType$ElevationHeightCode<- as.numeric(gsub("NULL","NA",PlotType$ElevationHeightCode))

PlotType$ContinentID <- as.numeric(gsub("=","",PlotType$ContinentID))
PlotType$ContinentID <- as.numeric(gsub("NULL","NA",PlotType$ContinentID))

PlotType$CountryID <- as.numeric(gsub("=","",PlotType$CountryID))
PlotType$CountryID <- as.numeric(gsub("NULL","NA",PlotType$CountryID))


#AGB indviduals Chave Moist. Use this output to compare later Feldpaush estimated heights vs new heights protocol.
#HtF is the Height with Weibull model using regional parameters as in Feldpaush.
#In the next version of the R package need to correct the output for snapped trees so Htf is not 0.

Heights1<-AGBChv05MH(Heights)

#nrow(Heights1)
#head(Heights1)
#Example of a snapped tree
#TreeIDk <- Heights1[Heights1$TreeID== 203374,]
#TreeIDk[order(TreeIDk$Census.No),]

## Recategorize F5 Methods. This bit of code is useful when doing the comparisons of laser vs clinometer

Heights1$Method<- ifelse(Heights1$F5==1,1,
                        ifelse(Heights1$F5==6,6,
                               ifelse( Heights1$F5==2 | Heights1$F5==3, 3,4
                               )
                        )
)

## CHeck number of directly measured trees
#aggregate(TreeID ~ Method, data=Heights1, FUN= length)


# add useful columns. These columns should be incorporated in the next version of the Rpackage      

Heights1$Palm<- ifelse ( grepl('Arecaceae',Heights1$Family)==TRUE |
                        grepl('Strelitziaceae',Heights1$Family)==TRUE |
                          grepl('Poaceae',Heights1$Family)==TRUE |
                          grepl('Cyatheaceae',Heights1$Family)==TRUE
                        ,1,0)

#to check number of records in each family
#aggregate(TreeID ~ Family +Palm, data = Heights1, FUN = length )


Heights1$PomChange <- ifelse(grepl('6',Heights1$F4)==TRUE,1,0)


#Merge with PlotType
Heights2<-merge(Heights1, PlotType, by ='PlotID')
#head(Heights2)
#nrow(Heights2)



# Select data with diameters and height 
# Remove trees without heights or with height=0.  Remove PAlms. To generate D:H relationships main plot views should be used.

TreesHt <- Heights2[Heights2$Height>0 & Heights2$Alive==1 & Heights2$DBH1>90 & Heights2$DBH1<5000    & !is.na(Heights2$Height) & Heights2$Height<90 
                   & Heights2$Palm==0 &  !is.na(Heights2$F5), ]

# Exclude F3-3
TreesHt <- TreesHt[ grepl('3',TreesHt$F3)==FALSE,  ]
# Exclude F4 not like '60' or '0'
TreesHt <- TreesHt[ grepl('0',TreesHt$F4)==TRUE|grepl('06',TreesHt$F4)==TRUE ,  ]
# Exclude f1 flag1= b, c, d, f,g,j,k,m,o,p
TreesHt <- TreesHt[ grepl('b',TreesHt$F1)==FALSE & grepl('c',TreesHt$F1)==FALSE &
                            grepl('d',TreesHt$F1)==FALSE & grepl('f',TreesHt$F1)==FALSE&
                            grepl('g',TreesHt$F1)==FALSE & grepl('j',TreesHt$F1)==FALSE &
                            grepl('k',TreesHt$F1)==FALSE & grepl('m',TreesHt$F1)==FALSE &
                            grepl('o',TreesHt$F1)==FALSE & grepl('p',TreesHt$F1)==FALSE
                    ,  ]

nrow(TreesHt)
#Exclude treeswith method 1
TreesHt <- TreesHt[!(TreesHt$Method==1),]

## Aggregate to check rules are applied correctly
aggregate(TreeID ~ F1, data=TreesHt, FUN=length)
aggregate(TreeID ~ F3, data=TreesHt, FUN=length)
aggregate(TreeID ~ F4, data=TreesHt, FUN=length)
aggregate(TreeID ~ Method, data=TreesHt, FUN=length)
aggregate(TreeID ~ Palm, data=TreesHt, FUN=length)




###1,ContinentData## 
HtCont<- TreesHt[,c('ContinentID','PlotID','Height','DBH1')]

library (nlme)
weib1 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|ContinentID,
                  data=HtCont,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
summary(weib1)
ContinentCoef<-data.frame(coef(weib1))
colnames(ContinentCoef) <- c('a_Continent','b_Continent', 'c_Continent')
ContinentCoef$ContinentID = rownames(ContinentCoef)
ContinentCoef
plot(weib1)

#2. ContinentData and Forest Type
HtCont_Type<- (TreesHt[,c('ContinentID','ForestMoistureID', 'EdaphicHeightCode', 'ElevationHeightCode','Height','DBH1')])
good<-complete.cases(HtCont_Type)
HtCont_Typea<-HtCont_Type[good,]
#head(HtCont_Typea)

library (nlme)
weib2 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|ContinentID/ForestMoistureID/EdaphicHeightCode/ElevationHeightCode,
                  data=HtCont_Typea,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
summary(weib2)
ContinentCoef_Type<-data.frame(coef(weib2))
colnames(ContinentCoef_Type) <- c('a_Continent_T','b_Continent_T', 'c_Continent_T')

ContinentCoef_Type$ContType <-rownames(ContinentCoef_Type)
head(ContinentCoef_Type)

ContinentCoef_Type$ContinentID<-unlist(lapply(strsplit(ContinentCoef_Type$ContType, "/"),"[", 1))
ContinentCoef_Type$ForestMoistureID<-unlist(lapply(strsplit(ContinentCoef_Type$ContType, "/"),"[", 2))
ContinentCoef_Type$EdaphicHeightCode<-unlist(lapply(strsplit(ContinentCoef_Type$ContType, "/"),"[", 3))
ContinentCoef_Type$ElevationHeightCode<-unlist(lapply(strsplit(ContinentCoef_Type$ContType, "/"),"[", 4))
ContinentCoef_Type
plot (weib2)

#3. Biogeographic Region

HtBiogeo<- TreesHt[,c('ContinentID','CountryID','BiogeographicalRegionID', 'PlotID','Height','DBH1')]
library (nlme)
weib3 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|'BiogeographicalRegionID',
                  data=HtBiogeo,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
summary(weib3)
BioRCoef<-data.frame(coef(weib3))
colnames(BioRCoef) <- c('a_BioR','b_BioR', 'c_BioR')
BioRCoef$BiogeographicalRegionID = rownames(BioRCoef)
BioRCoef
plot(weib3)
#head(BioRCoef)
##write.csv(BioRCoef,'BioRCoef.csv')


#4. Biogeographic Region/ForestType
HtBiogeoFt<- TreesHt[,c('ContinentID','CountryID','BiogeographicalRegionID', 'PlotID','ForestMoistureID', 'EdaphicHeightCode', 'ElevationHeightCode','Height','DBH1')]
weib4 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|BiogeographicalRegionID/ForestMoistureID/EdaphicHeightCode/ElevationHeightCode,
                  data=HtBiogeoFt,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
summary(weib4)
BioRCoefFt<-data.frame(coef(weib4))
colnames(BioRCoefFt) <- c('a_BioRF','b_BioRF', 'c_BioRF')
BioRCoefFt$BioRCoefFtype = rownames(BioRCoefFt)
head(BioRCoefFt)
BioRCoefFt$BiogeographicalRegionID<-unlist(lapply(strsplit(BioRCoefFt$BioRCoefFtype, "/"),"[", 1))
BioRCoefFt$ForestMoistureID<-unlist(lapply(strsplit(BioRCoefFt$BioRCoefFtype, "/"),"[", 2))
BioRCoefFt$EdaphicHeightCode<-unlist(lapply(strsplit(BioRCoefFt$BioRCoefFtype, "/"),"[", 3))
BioRCoefFt$ElevationHeightCode<-unlist(lapply(strsplit(BioRCoefFt$BioRCoefFtype, "/"),"[", 4))
BioRCoefFt
plot(weib4)


#write.csv(BioRCoefFt,'BioRCoefFt.csv')



##5 Analyze data by country
HtCountry<- TreesHt[,c('ContinentID','CountryID', 'PlotID','Height','DBH1')]
library (nlme)
weib5 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|CountryID,
                  data=HtCountry,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
summary(weib5)
CountryCoef<-data.frame(coef(weib5))
colnames(CountryCoef) <- c('a_Country','b_Country', 'c_Country')
CountryCoef$CountryID = rownames(CountryCoef)
CountryCoef
plot(weib5)
#write.csv(CountryCoef,'CountryCoef.csv')

##6 Analyze data by countryand ForestType
HtCountryF<- TreesHt[,c('ContinentID','CountryID', 'PlotID','ForestMoistureID', 'EdaphicHeightCode', 'ElevationHeightCode','Height','DBH1')]
library (nlme)
weib6 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|CountryID/ForestMoistureID/EdaphicHeightCode/ElevationHeightCode,
                  data=HtCountryF,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
summary(weib6)

HtCountryFt<-data.frame(coef(weib6))
colnames(HtCountryFt) <- c('a_CountryF','b_CountryF', 'c_CountryF')
HtCountryFt$HtCountryFtypw = rownames(HtCountryFt)
head(HtCountryFt)

HtCountryFt$CountryID<-unlist(lapply(strsplit(HtCountryFt$HtCountryFtypw, "/"),"[", 1))
HtCountryFt$ForestMoistureID<-unlist(lapply(strsplit(HtCountryFt$HtCountryFtypw, "/"),"[", 2))
HtCountryFt$EdaphicHeightCode<-unlist(lapply(strsplit(HtCountryFt$HtCountryFtypw, "/"),"[", 3))
HtCountryFt$ElevationHeightCode<-unlist(lapply(strsplit(HtCountryFt$HtCountryFtypw, "/"),"[", 4))
HtCountryFt

plot(weib6)


#write.csv(HtCountryFt,'HtCountryFt.csv')



#7.Analyze data by clusterid
HtCluster<- (TreesHt[,c('ClusterID','Height','DBH1')])
goodCl<-complete.cases(HtCluster)
HtClusterA<-HtCluster[goodCl,]
head(HtClusterA)

library (nlme)
weib7 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|ClusterID,
                  data=HtClusterA,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)

summary(weib7)

ClusterCoef<-data.frame(coef(weib7))
colnames(ClusterCoef) <- c('a_Cluster','b_Cluster', 'c_Cluster')
ClusterCoef$ClusterID = rownames(ClusterCoef)
ClusterCoef
plot(weib7)
#write.csv(ClusterCoef,'ClusterCoef.csv')


#8. Cluster id and ForestType
HtClusterFt <- TreesHt[,c('ClusterID','ForestMoistureID', 'EdaphicHeightCode', 'ElevationHeightCode','Height','DBH1')]
goodClFt<-complete.cases(HtClusterFt)
HtClusterFtA<-HtClusterFt[goodClFt,]
head(HtClusterFtA)
##write.csv (HtClusterFA, 'HtClusterFA.csv')

#aggregate (DBH1 ~ ClusterID +ForestMoistureID + EdaphicHeightCode +ElevationHeightCode, data = HtClusterFA, FUN = length)
##

library (nlme)
weib8 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))| ClusterID/ForestMoistureID/EdaphicHeightCode/ElevationHeightCode,
                  data=HtClusterFtA,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
summary(weib8)
ClusterFCoef<-data.frame(coef(weib8))
colnames(ClusterFCoef) <- c('a_ClusterF','b_ClusterF', 'c_ClusterF')
ClusterFCoef$ClusterF = rownames(ClusterFCoef)
ClusterFCoef

ClusterFCoef$ClusterID<-unlist(lapply(strsplit(ClusterFCoef$ClusterF, "/"),"[", 1))
ClusterFCoef$ForestMoistureID<-unlist(lapply(strsplit(ClusterFCoef$ClusterF, "/"),"[", 2))
ClusterFCoef$EdaphicHeightCode<-unlist(lapply(strsplit(ClusterFCoef$ClusterF, "/"),"[", 3))
ClusterFCoef$ElevationHeightCode<-unlist(lapply(strsplit(ClusterFCoef$ClusterF, "/"),"[", 4))
head(ClusterFCoef)
#write.csv(ClusterFCoef,'ClusterFCoef.csv')
plot (weib8)


#Analyze data by  plot
HtPlot<- TreesHt[,c('PlotID','Height','DBH1')]
library (nlme)
weib9 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))| PlotID,
                  data=HtPlot,
                  na.action=na.omit,
                  start = c(a =25, b= 0.05, c= 0.7),
                  pool=FALSE)
summary(weib9)
PlotCoef<-data.frame(coef(weib9))
colnames(PlotCoef) <- c('a_Plot','b_Plot', 'c_Plot')
PlotCoef$PlotID = rownames(PlotCoef)
head(PlotCoef)
nrow(PlotCoef)

plot (weib9)
#write.csv(PlotCoef,'PlotCoef.csv')



##join all datawithsubset of data where trees are alive
#Height PlotID
head(Heights2)

nrow (Heights2)
HtAl<-Heights2[Heights2$Alive==1, c('PlotID','ContinentID', 'CountryID', 'ContinentName','ClusterID', 'PlotViewID', 'PlotCode', 
                            'ForestMoistureID','ForestEdaphicID','ForestElevationID','EdaphicHeightCode', 'ElevationHeightCode',
                            'BiogeographicalRegionID','AllometricRegionID', 'PlotArea',
                            'TreeID','Family', 'Species','DBH1','DBH2','DBH3', 'DBH4', 'F1','F2','F3', 'F4','F5' ,'Census.No' ,'Census.Mean.Date','Alive','Snapped','Palm', 'CensusStemDied', 
                            'Height','HtF','AGBind','WD')]
head (HtAl)
nrow(HtAl)

## merging plot coefficients with alive records.D4 is used to estimate height. For trees without POM change D1=D4

Ht1 <- merge(HtAl,PlotCoef, by='PlotID', all.x=TRUE)

Ht1$HtPlot<-Ht1$a_Plot*(1-exp(-Ht1$b_Plot *(Ht1$DBH4/10)^Ht1$c_Plot))

nrow (Ht1)
head(Ht1)

#Cluster and Forest type
Ht2 <- merge(Ht1,ClusterFCoef, by=c('ClusterID','ForestMoistureID','EdaphicHeightCode','ElevationHeightCode'), all.x=TRUE)

Ht2$HtClFt<-Ht2$a_ClusterF*(1-exp(-Ht2$b_ClusterF *(Ht2$DBH4/10)^Ht2$c_ClusterF))

nrow (Ht2)
head(Ht2)

# CLUSTER
Ht3<- merge(Ht2, ClusterCoef, by='ClusterID', all.x=TRUE)

Ht3$HtCl<-Ht3$a_Cluster*(1-exp(-Ht3$b_Cluster *(Ht3$DBH4/10)^Ht3$c_Cluster))

nrow(Ht3)
head(Ht3)

#Country/ForestType
Ht4 <- merge(Ht3,HtCountryFt, by=c('CountryID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'), all.x=TRUE)

Ht4$HtCtF<-Ht4$a_CountryF*(1-exp(-Ht4$b_CountryF *(Ht4$DBH4/10)^Ht4$c_CountryF))

nrow(Ht4)
head(Ht4)

#Country
Ht5 <- merge(Ht4,CountryCoef, by='CountryID', all.x=TRUE)

Ht5$HtCt<-Ht5$a_Country*(1-exp(-Ht5$b_Country *(Ht5$DBH4/10)^Ht5$c_Country))

nrow(Ht5)
head(Ht5)

#Biogeographic region and Forest type
Ht6<- merge (Ht5,BioRCoefFt, by=c('BiogeographicalRegionID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'), all.x=TRUE)

Ht6$HtBf<- Ht6$a_BioRF*(1-exp(-Ht6$b_BioRF*(Ht6$DBH4/10)^Ht6$c_BioRF))

nrow(Ht6)
head(Ht6)

#Biogeographic region 
Ht7<- merge (Ht6,BioRCoef, by='BiogeographicalRegionID', all.x=TRUE)


Ht7$HtB<- Ht7$a_BioR*(1-exp(-Ht7$b_BioR*(Ht7$DBH4/10)^Ht7$c_BioR))

nrow(Ht7)
head(Ht7)


#Continent and ForestType

Ht8<- merge (Ht7,ContinentCoef_Type, by=c('ContinentID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'), all.x=TRUE)

Ht8$HtCoF<- Ht8$a_Continent_T*(1-exp(-Ht8$b_Continent_T*(Ht8$DBH4/10)^Ht8$c_Continent_T))

nrow(Ht8)
head(Ht8)

#Continent
Ht9 <- merge(Ht8, ContinentCoef, by='ContinentID', all.x=TRUE)

Ht9$HtCo<-Ht9$a_Continent*(1-exp(-Ht9$b_Continent*(Ht9$DBH4/10)^Ht9$c_Continent))

nrow(Ht9)
head(Ht9)

#Select Reference Height to Use, for trees with multiple heights. Later on the last measured height will be selected for offset.
#For offset only for trees that where selected for the analysis will be included
# for Palms heights is NA


            
Ht9$HtRef<-   ifelse (Ht9$Palm==1, NA,     #Palm NA
               ifelse (!is.na(Ht9$Height), Ht9$Height,   #Measured Height 
                ifelse (!is.na(Ht9$HtPlot),Ht9$HtPlot,   # Plot height     
                  ifelse(!is.na(Ht9$HtClFt),Ht9$HtClFt, #ClusterFt
                          ifelse(!is.na(Ht9$HtCl),  Ht9$HtCl, #Cluster
                                 ifelse(!is.na( Ht9$HtCtF), Ht9$HtCtF,# CountryFt
                                        ifelse(!is.na( Ht9$HtCt),Ht9$HtCt,  #Country
                                               ifelse(!is.na(Ht9$HtBf), Ht9$HtBf, #Biogeographical Region ForestType
                                                      ifelse(!is.na(Ht9$HtB), Ht9$HtB,# Biogeographical Region
                                                             ifelse(!is.na(Ht9$HtCoF),Ht9$HtCoF,#ContinentForestType
                                                                    Ht9$HtCo # Continent
                                                             ))))))))))

Ht9$HtRefType<- ifelse(Ht9$Palm == 1, NA, #Palm NA
                    ifelse(!is.na(Ht9$Height),'MeasuredHt', 
                       ifelse(!is.na(Ht9$HtPlot),'HtPlot' ,
                              ifelse(!is.na(Ht9$HtClFt),'HtClFt', #ClusterFt
                                     ifelse(!is.na(Ht9$HtCl),  'HtCl', #Cluster
                                            ifelse(!is.na( Ht9$HtCtF), 'HtCtF',# CountryFt
                                                   ifelse(!is.na( Ht9$HtCt),'HtCt',  #Country
                                                          ifelse(!is.na(Ht9$HtBf), 'HtBf', #Biogeographical Region ForestType
                                                                 ifelse(!is.na(Ht9$HtB), 'HtB',# Biogeographical Region
                                                                        ifelse(!is.na(Ht9$HtCoF),'HtCoF',#ContinentForestType
                                                                               'HtCo' # Continent
                                                                        ))))))))))


Ht9$HtRefEst<- ifelse (Ht9$Palm == 1, NA, # Palm NA 
                  ifelse(!is.na(Ht9$HtPlot),Ht9$HtPlot, #Plot
                     ifelse(!is.na(Ht9$HtClFt),Ht9$HtClFt,#ClusterFt
                            ifelse(!is.na(Ht9$HtCl),  Ht9$HtCl, #Cluster
                                   ifelse(!is.na( Ht9$HtCtF), Ht9$HtCtF,# CountryFt
                                          ifelse(!is.na( Ht9$HtCt),Ht9$HtCt,  #Country
                                                 ifelse(!is.na(Ht9$HtBf), Ht9$HtBf, #Biogeographical Region ForestType
                                                        ifelse(!is.na(Ht9$HtB), Ht9$HtB,# Biogeographical Region
                                                               ifelse(!is.na(Ht9$HtCoF),Ht9$HtCoF,#ContinentForestType
                                                                      Ht9$HtCo # continent
                                                               )))))))))

Ht9$HtRefEstType<- ifelse (Ht9$Palm == 1, NA, # Palm NA 
                       ifelse(!is.na(Ht9$HtPlot),'Plot', #Plot
                              ifelse(!is.na(Ht9$HtClFt),'ClusterFt',#ClusterFt
                                     ifelse(!is.na(Ht9$HtCl),  'Cluster', #Cluster
                                            ifelse(!is.na( Ht9$HtCtF), 'CountryFt',# CountryFt
                                                   ifelse(!is.na( Ht9$HtCt),'Country',  #Country
                                                          ifelse(!is.na(Ht9$HtBf), 'BioFt', #Biogeographical Region ForestType
                                                                 ifelse(!is.na(Ht9$HtB), 'Bio',# Biogeographical Region
                                                                        ifelse(!is.na(Ht9$HtCoF),'ContF',#ContinentForestType
                                                                               'Cont' # continent
                                                                        )))))))))



head(Ht9)
aggregate(TreeID ~ HtRefType + Palm, data=Ht9, FUN=length)

#to check levels have been allocated properly

head ( Ht9[grepl('Plot', Ht9$HtRefType)== TRUE, ],3)
head ( Ht9[grepl('HtClFt', Ht9$HtRefType)== TRUE, ],3)
head ( Ht9[grepl('HtB', Ht9$HtRefType)== TRUE, ],3)
#head(Heights9[(Heights9$HtRef>90),],3)

#list of trees with multiple heights, not Palms
#find number of heights by treeid

TreeswithHeights <- Ht9[Ht9$Height>0 & Ht9$Alive==1 & Ht9$DBH1>90 & Ht9$DBH1<5000    & !is.na(Ht9$Height) & Ht9$Height<90 
                    & Ht9$Palm==0 &  !is.na(Ht9$F5), ]

# Exclude F3-3
TreeswithHeights<- TreeswithHeights[ grepl('3',TreeswithHeights$F3)==FALSE,  ]
# Exclude F4 not like '60' or '0'
TreeswithHeights <-TreeswithHeights[ grepl('0',TreeswithHeights$F4)==TRUE|grepl('06',TreeswithHeights$F4)==TRUE ,  ]
# Exclude f1 flag1= b, c, d, f,g,j,k,m,o,p
TreeswithHeights <- TreeswithHeights[ grepl('b',TreeswithHeights$F1)==FALSE & grepl('c',TreeswithHeights$F1)==FALSE &
                      grepl('d',TreeswithHeights$F1)==FALSE & grepl('f',TreeswithHeights$F1)==FALSE&
                      grepl('g',TreeswithHeights$F1)==FALSE & grepl('j',TreeswithHeights$F1)==FALSE &
                      grepl('k',TreeswithHeights$F1)==FALSE & grepl('m',TreeswithHeights$F1)==FALSE &
                      grepl('o',TreeswithHeights$F1)==FALSE & grepl('p',TreeswithHeights$F1)==FALSE
                    ,  ]

#remove tree height estimated by eye
TreeswithHeights<- TreeswithHeights[!(TreeswithHeights$F5 ==1),]

nrow(TreeswithHeights)

#TreeswithHeights[TreeswithHeights$TreeID== 55856,] # test tree with two censuses with height measured

#nrow(TreeswithHeights)

HtM1<-aggregate (Height ~ TreeID+PlotID +PlotViewID, data= TreeswithHeights, FUN=length)
colnames(HtM1) <-c('TreeID','PlotID','PlotViewID','Records')
nrow(HtM1)
#Add number of tree height records information to the table
HtM2<-merge (HtM1, TreeswithHeights, by=c('TreeID','PlotID','PlotViewID'))
head(HtM2)
nrow(HtM2)

#HtM2[HtM2$TreeID==55856,]

#maximum censusnumber for all tree_census_heights
mxHt<- aggregate(Census.No ~TreeID+PlotID+PlotViewID, data=HtM2,FUN=max)

#mxHt[mxHt$TreeID == 55856,]

HtCns<-merge(mxHt,TreeswithHeights, by= c('TreeID', 'PlotID','PlotViewID','Census.No') )

nrow(mxHt)
nrow(HtCns)
head(HtCns)

#HtCns[HtCns$TreeID==55856,]

#Select information necessary for Offset

HtCns2<-HtCns[ , c('TreeID','PlotViewID','PlotID','Census.No','Height','HtRefEst','HtRefEstType')]

head(HtCns2)
#HtCns2[HtCns2$TreeID==112142,]

#HtCns2[HtCns2$TreeID==100024,]

#HtCns2[HtCns2$TreeID==55856,]

#Offset table

colnames(HtCns2)<- c('TreeID','PlotViewID','PlotID','CensusNoHtLast','HtLast','HtRefEstLast','HtRefEstTypeLast')
Toffset<- HtCns2
Toffset$OffsetH<-Toffset$HtLast-Toffset$HtRefEstLast
head(Toffset)
##Toffset[Toffset$TreeID== 55856,]


#merge with heights list
#nrow(Heights[Heights$Alive==1,])


##Ht9[Ht9$TreeID==55856,]

##

Heightsb<-merge(Ht9,Toffset,by=c('TreeID','PlotID','PlotViewID'),all=TRUE )
head(Heightsb)
Heightsb[Heightsb$TreeID== 55856,]

Heightsb$Ht1<- ifelse (is.na(Heightsb$HtLast), Heightsb$HtRefEst,
                         ifelse(Heightsb$CensusNoHtLast== Heightsb$Census.No, Heightsb$Height, Heightsb$HtRefEst+Heightsb$OffsetH)
                              )

Heightsb$Ht2 <- ifelse (is.na(Heightsb$HtLast), Heightsb$HtRefEst,
                       ifelse(Heightsb$CensusNoHtLast== Heightsb$Census.No, Heightsb$Height, 
                              Heightsb$HtRefEst+(Heightsb$OffsetH*(Heightsb$HtRefEstLast/Heightsb$HtRefEst)) )
                        )

head(Heightsb)
  


#Test2 <- Heightsb[Heightsb$TreeID==100024,]

#Test3 <- Heightsb[Heightsb$TreeID==55856,]
#Test2
#write.csv (Test3, 'Test3.csv') # to test Offest proposed by Oliver


##AGB individuals using Chave Moist only for plots with Height at plot level
PlotHt<-PlotCoef[!is.na(PlotCoef$a_Plot),]
nrow(PlotHt)
PlotIDtouse<-PlotHt[,c('PlotID')]



HtPlotTest<- Heightsb[Heightsb$PlotID %in% PlotIDtouse, ]

head(HtPlotTest)
nrow(HtPlotTest)


##AGBind ## need to add the Palm AGB equation ,example using dbh4 for AGB 
HtPlotTest$AGB_HFed <- ifelse (HtPlotTest$Palm==1 & HtPlotTest$DBH1>0, (exp(-3.3488 + 2.7483*log(HtPlotTest$DBH4/10))*exp((0.588^2)/2))/1000 ,
                               ifelse (HtPlotTest$Palm==0 & HtPlotTest$DBH1>0
                                       , (0.0509*HtPlotTest$WD * ((HtPlotTest$DBH4/10)^2)* HtPlotTest$HtF)/1000, NA))   
                         
##head(HtPlotTest[HtPlotTest$Palm==1, ],3)

HtPlotTest$AGB_Hest <- ifelse (HtPlotTest$Palm==1 & HtPlotTest$DBH1>0, (exp(-3.3488 + 2.7483*log(HtPlotTest$DBH4/10))*exp((0.588^2)/2))/1000 ,
                               ifelse (HtPlotTest$Palm==0 & HtPlotTest$DBH1>0
                                       , (0.0509*HtPlotTest$WD * ((HtPlotTest$DBH4/10)^2)* HtPlotTest$HtRefEst)/1000, NA))  

        
    

HtPlotTest$AGB_Hest_O1 <- ifelse (HtPlotTest$Palm==1 & HtPlotTest$DBH1>0, (exp(-3.3488 + 2.7483*log(HtPlotTest$DBH4/10))*exp((0.588^2)/2))/1000 ,
                                 ifelse (HtPlotTest$Palm==0 & HtPlotTest$DBH1>0
                                         , (0.0509*HtPlotTest$WD * ((HtPlotTest$DBH4/10)^2)* HtPlotTest$Ht1)/1000, NA)) 
        
        
      


HtPlotTest$AGB_Hest_O2 <- ifelse (HtPlotTest$Palm==1 & HtPlotTest$DBH1>0, (exp(-3.3488 + 2.7483*log(HtPlotTest$DBH4/10))*exp((0.588^2)/2))/1000 ,
                                  ifelse (HtPlotTest$Palm==0 & HtPlotTest$DBH1>0
                                          , (0.0509*HtPlotTest$WD * ((HtPlotTest$DBH4/10)^2)* HtPlotTest$Ht2)/1000, NA))
        
    

#to check calculations
#  HtPlotTest[HtPlotTest$TreeID==55856, c('TreeID','Census.No','Height', 'HtRef' ,'HtRefEst','AGBind', 'AGB_HFed','AGB_Hest','AGB_Hest_O1', 'AGB_Hest_O2')]
head(HtPlotTest)

AGBSummary<- aggregate( cbind(AGBind, AGB_HFed,AGB_Hest,AGB_Hest_O1, AGB_Hest_O2) ~ PlotCode + Census.No + ContinentName + PlotArea, data=HtPlotTest, FUN=sum)
write.csv (AGBSummary,'AGBSummary.csv')


## Need to Improve graphs
library(lattice)
xyplot( c(AGBind, AGB_HFed,AGB_Hest,AGB_Hest_O1, AGB_Hest_O2)~ Census.No| PlotCode, 
        data =AGBSummary,
        #layout= c(2,10),
        xlab = "Census.No",
        ylab = "AGB"
)


##


xyplot(c(Height, Ht_PlotMt1, Ht_PlotMt2)~ DBH1/10| PlotCode,  groups = c(Height, Ht_PlotMt1, Ht_PlotMt2),
       data =TreesHt_use,
       #layout= c(2,10),
       xlab = "Diameter (cm)",
       ylab = "Height (m)"
)






