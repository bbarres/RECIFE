##############################################################################/
##############################################################################/
#Data for the analysis and figure production of the RECIFE project
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(ade4)
library(dplyr)
library(drc)
library(factoextra)
library(gdata)
library(genepop)
library(hierfstat)
library(mapplots)
library(plotrix)
library(poppr)
library(raster)
library(RColorBrewer)
library(sf)
library(tidyr)
library(treemap)
library(vioplot)


##############################################################################/
#loading the geographical data####
##############################################################################/

#load geographical data of departements and regions
load("data/DEP_SHP.RData")
load("data/REG_SHP.RData")

#load the barycentre coordinates of departements and regions
load("data/coorddep.RData")
load("data/coordreg.RData")

DEP_SHP.1<-crop(DEP_SHP,extent(114528.2,1132915.3,6500000,7168463))
REG_SHP.1<-crop(REG_SHP,extent(114528.2,1132915.3,6500000,7168463))

DEP_SHP.2<-crop(DEP_SHP,extent(466064.2,854936.4,6741138,7123635))
REG_SHP.2<-crop(REG_SHP,extent(466064.2,854936.4,6741138,7123635))


##############################################################################/
#loading the bioassay results data####
##############################################################################/

#dose response data set
datamyc2<-read.table("data/CRMYC_270821.txt",header=TRUE,
                     sep=";",stringsAsFactors=TRUE)

#load the resistance results for the 2019-2020 campaign
oldSA<-read.delim(
  "data/data_DC_AZ_FH_Carb_2019_2020.txt",
  header=TRUE,
  sep="\t"
)
#load the resistance results for the 2019-2020 campaign
newSA<-read.delim(
  "data/data_DC_classes_2019_2020.txt",
  header=TRUE,
  sep="\t"
)
AllSamp<-rbind(oldSA[,c("prelvt_id","gps_lat","gps_long")],
               newSA[,c("prelvt_id","gps_lat","gps_long")])
AllSamp<-AllSamp[!is.na(AllSamp$gps_lat),]
#turning this dataframe into a spatial dataframe (wgs84)
AllSamp.wgs<-SpatialPointsDataFrame(coords=AllSamp[,c("gps_long","gps_lat")],
                                    data=AllSamp,
                                proj4string=CRS("+proj=longlat +datum=WGS84")
)
AllSamp<-spTransform(AllSamp.wgs,CRS("+init=epsg:2154"))

#CI50 for the different cyp51 haplotypes
haplo51<-read.table("data/haplo_pheno.txt",sep="\t",header=TRUE,quote="",
                    colClasses=c("character","factor","factor","character",
                                 "numeric","numeric","Date","factor",
                                 "factor","factor","character","character",
                                 "character","character","character",
                                 "character","character","character",
                                 "character","character","character",
                                 "character"))


##############################################################################/
#loading the microsatellite data####
##############################################################################/

sugmic<-read.table("data/sugarMicro2.dat",sep="\t",header=TRUE,
                   stringsAsFactors=TRUE)


##############################################################################/
#Writing info session for reproducibility####
##############################################################################/

sink("session_info.txt")
print(sessioninfo::session_info())
sink()
#inspired by an R gist of François Briatte: 
#https://gist.github.com/briatte/14e47fb0cfb8801f25c889edea3fcd9b


##############################################################################/
#END
##############################################################################/