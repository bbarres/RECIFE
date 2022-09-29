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
library(mapplots)
library(plotrix)
library(poppr)
library(rgdal)
library(rgeos)
library(raster)
library(RColorBrewer)
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


##############################################################################/
#loading the microsatellite data####
##############################################################################/

sugmic<-read.table("data/sugarMicro.dat",sep="\t",header=TRUE,
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