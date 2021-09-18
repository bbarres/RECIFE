##############################################################################/
##############################################################################/
#Code for producing a map of Cercospora beticola resistance
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(rgdal)
library(rgeos)
library(plotrix)
library(mapplots)
library(RColorBrewer)


##############################################################################/
#loading and preparing the data####
##############################################################################/

#load geographical data
load("data/commu.RData")
load("data/arrond.RData")
load("data/departe.RData")
load("data/regionsLight.RData")

#some information on the data structure of the geodata
class(commu)
slotNames(commu)
summary(commu@data)

#isolate the information in the spatial data on the communes
db_commu<-commu@data
summary(db_commu)
db_arrond<-arrond@data
db_arrond$DEPARR<-paste(db_arrond$INSEE_DEP,db_arrond$INSEE_ARR)

#load the resistance results for the 2019 campaign
databruteTOT <- read.delim(
  "data/cerco_germ_dec19.txt",
  header = TRUE,
  sep = "\t",
  colClasses = c("character", "character", "character",
                 "character", "factor")
)

#load the resistance results for the 2020 campaign
databruteTOT <- read.delim(
  "data/data_DC_2020.txt",
  header = TRUE,
  sep = "\t",
  colClasses = c("character", "character", "character",
                 "character", "factor")
)

#load the resistance results for the 2020 campaign
databruteTOT <- read.delim(
  "data/data_DC_classes_2020.txt",
  header = TRUE,
  sep = "\t",
  colClasses = c("character", "character", "character",
                 "character", "factor")
)

# #first we merge the resistance table with the commune info
# databrute<-merge(databrute,db_commu,by.x="parcel_cmne",by.y="NOM_COM_M")
# #in order to acces to the arrondissement ID, we create an individual ID for 
# #each arrondissement combining INSEE_DEP and INSEE_ARR
# Raox_list$DEPARR<-paste(Raox_list$INSEE_DEP,Raox_list$INSEE_ARR)
# #then we merge the resistance table with the arrondissement info
# Raox_list<-merge(Raox_list,db_arrond,by.x="DEPARR",by.y="DEPARR")


##############################################################################/
#FENTINE HYDROXYDE MAP####
##############################################################################/

databrute<-databruteTOT[databruteTOT$pest_sa_id=="FENTINE HYDROXYDE" & 
                          databruteTOT$dose!=0,]

#counting the number of sampling according to their status
datacount <- cbind(
  "dep_ID" = row.names(table(
    databrute$parcel_dpt,
    databrute$rslt_03, exclude =
      ""
  )),
  "Resistant" = table(databrute$parcel_dpt,
                      as.numeric(databrute$rslt_03)>0, exclude = "")[,2],
  "Sensible" = table(databrute$parcel_dpt,
                     as.numeric(databrute$rslt_03)>0, exclude = "")[,1],
  "Total" = rowSums(table(
    databrute$parcel_dpt,
    as.numeric(databrute$rslt_03)>0, exclude = ""
  ))
)
#export the count data by departement
write.table(
  datacount,
  file = "output/datacountFENTHY.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#extract the department coordinates
ind_list <- departe$INSEE_DEP
coorddep <- data.frame(
  "longitude" = departe@polygons[1][[1]]@labpt[1],
  "latitude" = departe@polygons[1][[1]]@labpt[2]
)
for (i in 2:length(ind_list)) {
  coorddep <- rbind(
    coorddep,
    cbind(
      "longitude" = departe@polygons[i][[1]]@labpt[1],
      "latitude" = departe@polygons[i][[1]]@labpt[2]
    )
  )
}
coorddep <- cbind("dep_ID" = ind_list, coorddep)

#combining the information of the geographic
data2map <- merge(datacount, coorddep, by = "dep_ID")
data2map <- data2map[order(as.numeric(as.character(data2map$Total)),
                           decreasing = TRUE), ]

#defining the colors of the pies
colovec <- c(brewer.pal(9, "Reds")[6], brewer.pal(9, "Blues")[7])

#actual plotting
op <- par(mar = c(0, 0, 2, 0))
plot(departe,main="FENTINE HYDROXIDE")
draw.pie(
  x = data2map$longitude,
  y = data2map$latitude,
  z = cbind((as.numeric(
    as.character(data2map$Resistant)
  )),
  (as.numeric(
    as.character(data2map$Sensible)
  ))),
  col = colovec,         #colors of the pie
  lty = 1,               #line type of the pie
  border = "grey60",     #color of the border of the pie
  lwd = 0.1,             #control the width of the border
  radius = 35000, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels = NA,
  scale=FALSE # should the radius be scaled according to sample size
)

#writing the number of samples for each departement
text(
  x = data2map$longitude,
  y = data2map$latitude,
  col = "black",
  font = 2,
  labels = as.character(data2map$Total),
  cex = 1.5              #size of the text
)
par(op)

#export the map to a pdf file 7 x 7 inches (for examples)

#boxplot for the resistant populations
boxplot(as.numeric(databrute$rslt_03[databrute$rslt_03!=0]),
        boxwex=0.5,las=1,ylim=c(0,100),col="cornflowerblue",
        main=paste("Fentine Hydroxyde\n(n=",
                   length(databrute$rslt_03[databrute$rslt_03!=0]),
                   "/",
                   length(databrute$rslt_03),")",sep=""),
        ylab="% germination")
abline(h=mean(as.numeric(databrute$rslt_03[databrute$rslt_03!=0])),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(databrute$rslt_03[databrute$rslt_03!=0]),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.8),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)

#export to a pdf file 3 x 7 inches (for examples)

#distribution of the % of germination at the DD
plot(as.numeric(databrute$rslt_03[order(as.numeric(databrute$rslt_03))]),
     bg="cornflowerblue",pch=21,cex=1.5,las=1,ylim=c(0,130),
     ylab="% germination",
     main="Fentine Hydroxide")

#export to a pdf file 7 x 5 inches (for examples)


##############################################################################/
#CARBENDAZIME MAP####
##############################################################################/

databrute<-databruteTOT[databruteTOT$pest_sa_id=="CARBENDAZIME" & 
                          databruteTOT$dose!=0,]

#counting the number of sampling according to their status
datacount <- cbind(
  "dep_ID" = row.names(table(
    databrute$parcel_dpt,
    databrute$rslt_03, exclude =
      ""
  )),
  "Resistant" = 0,
  "Sensible" = table(databrute$parcel_dpt,
                     as.numeric(databrute$rslt_03)>0, exclude = "")[,1],
  "Total" = rowSums(table(
    databrute$parcel_dpt,
    as.numeric(databrute$rslt_03)>0, exclude = ""
  ))
)
#export the count data by departement
write.table(
  datacount,
  file = "output/datacountCARBEN.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#extract the department coordinates
ind_list <- departe$INSEE_DEP
coorddep <- data.frame(
  "longitude" = departe@polygons[1][[1]]@labpt[1],
  "latitude" = departe@polygons[1][[1]]@labpt[2]
)
for (i in 2:length(ind_list)) {
  coorddep <- rbind(
    coorddep,
    cbind(
      "longitude" = departe@polygons[i][[1]]@labpt[1],
      "latitude" = departe@polygons[i][[1]]@labpt[2]
    )
  )
}
coorddep <- cbind("dep_ID" = ind_list, coorddep)

#combining the information of the geographic
data2map <- merge(datacount, coorddep, by = "dep_ID")
data2map <- data2map[order(as.numeric(as.character(data2map$Total)),
                           decreasing = TRUE), ]

#defining the colors of the pies
colovec <- c(brewer.pal(9, "Reds")[6])

#actual plotting
op <- par(mar = c(0, 0, 2, 0))
plot(departe,main="CARBENDAZIME")
draw.pie(
  x = data2map$longitude,
  y = data2map$latitude,
  z = cbind((as.numeric(
    as.character(data2map$Resistant)
  )),
  (as.numeric(
    as.character(data2map$Sensible)
  ))),
  col = colovec,         #colors of the pie
  lty = 1,               #line type of the pie
  border = "grey60",     #color of the border of the pie
  lwd = 0.1,             #control the width of the border
  radius = 35000, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels = NA,
  scale=FALSE # should the radius be scaled according to sample size
)

#writing the number of samples for each departement
text(
  x = data2map$longitude,
  y = data2map$latitude,
  col = "black",
  font = 2,
  labels = as.character(data2map$Total),
  cex = 1.5              #size of the text
)
par(op)

#export the map to a pdf file 7 x 7 inches (for examples)

#boxplot for the resistant populations
boxplot(as.numeric(databrute$rslt_03[databrute$rslt_03!=0]),
        boxwex=0.5,las=1,ylim=c(0,120),col="cornflowerblue",
        main=paste("Carbendazime\n(n=",
                   length(databrute$rslt_03[databrute$rslt_03!=0]),
                   "/",
                   length(databrute$rslt_03),")",sep=""),
        ylab="% germination")
abline(h=mean(as.numeric(databrute$rslt_03[databrute$rslt_03!=0])),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(databrute$rslt_03[databrute$rslt_03!=0]),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.8),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)

#export to a pdf file 3 x 7 inches (for examples)

#distribution of the % of germination at the DD
plot(as.numeric(databrute$rslt_03[order(as.numeric(databrute$rslt_03))]),
     bg="cornflowerblue",pch=21,cex=1.5,las=1,ylim=c(0,130),
     ylab="% germination",
     main="Carbendazime")

#export to a pdf file 7 x 5 inches (for examples)


##############################################################################/
#AZOXYSTROBINE RLC####
##############################################################################/

databrute<-databruteTOT[databruteTOT$pest_sa_id=="AZOXYSTROBINE" & 
                          databruteTOT$synerg_id=="SHAM" & 
                          databruteTOT$dose!=0,]

#counting the number of sampling according to their status
datacount <- cbind(
  "dep_ID" = row.names(table(
    databrute$parcel_dpt,
    databrute$rslt_03, exclude =
      ""
  )),
  "Resistant" = table(databrute$parcel_dpt,
                      as.numeric(databrute$rslt_03)>0, exclude = "")[,2],
  "Sensible" = table(databrute$parcel_dpt,
                     as.numeric(databrute$rslt_03)>0, exclude = "")[,1],
  "Total" = rowSums(table(
    databrute$parcel_dpt,
    databrute$rslt_03>0, exclude = ""
  ))
)
#export the count data by departement
write.table(
  datacount,
  file = "output/datacountAZOXSHAM.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#extract the department coordinates
ind_list <- departe$INSEE_DEP
coorddep <- data.frame(
  "longitude" = departe@polygons[1][[1]]@labpt[1],
  "latitude" = departe@polygons[1][[1]]@labpt[2]
)
for (i in 2:length(ind_list)) {
  coorddep <- rbind(
    coorddep,
    cbind(
      "longitude" = departe@polygons[i][[1]]@labpt[1],
      "latitude" = departe@polygons[i][[1]]@labpt[2]
    )
  )
}
coorddep <- cbind("dep_ID" = ind_list, coorddep)

#combining the information of the geographic
data2map <- merge(datacount, coorddep, by = "dep_ID")
data2map <- data2map[order(as.numeric(as.character(data2map$Total)),
                           decreasing = TRUE), ]

#defining the colors of the pies
colovec <- c(brewer.pal(9, "Reds")[6], brewer.pal(9, "Blues")[7])

#actual plotting
op <- par(mar = c(0, 0, 2, 0))
plot(departe,main="AZOXYSTROBINE RLC")
draw.pie(
  x = data2map$longitude,
  y = data2map$latitude,
  z = cbind((as.numeric(
    as.character(data2map$Resistant)
  )),
  (as.numeric(
    as.character(data2map$Sensible)
  ))),
  col = colovec,         #colors of the pie
  lty = 1,               #line type of the pie
  border = "grey60",     #color of the border of the pie
  lwd = 0.1,             #control the width of the border
  radius = 35000, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels = NA,
  scale=FALSE # should the radius be scaled according to sample size
)

#writing the number of samples for each departement
text(
  x = data2map$longitude,
  y = data2map$latitude,
  col = "black",
  font = 2,
  labels = as.character(data2map$Total),
  cex = 1.5              #size of the text
)
par(op)

#export the map to a pdf file 7 x 7 inches (for examples)

#boxplot for the resistant populations
boxplot(as.numeric(databrute$rslt_03[databrute$rslt_03!=0]),
        boxwex=0.5,las=1,ylim=c(0,140),col="cornflowerblue",
        main=paste("Azoxystrobine RLC\n(n=",
                   length(databrute$rslt_03[databrute$rslt_03!=0]),
                   "/",
                   length(databrute$rslt_03),")",sep=""),
        ylab="% germination")
abline(h=mean(as.numeric(databrute$rslt_03[databrute$rslt_03!=0])),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(databrute$rslt_03[databrute$rslt_03!=0]),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.8),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)

#export to a pdf file 3 x 7 inches (for examples)

#distribution of the % of germination at the DD
plot(as.numeric(databrute$rslt_03[order(as.numeric(databrute$rslt_03))]),
     bg="cornflowerblue",pch=21,cex=1.5,las=1,ylim=c(0,130),
     ylab="% germination",
     main="Azoxystrobine RLC")

#export to a pdf file 7 x 5 inches (for examples)


##############################################################################/
#AZOXYSTROBINE R TOTAL####
##############################################################################/

databrute<-databruteTOT[databruteTOT$pest_sa_id=="AZOXYSTROBINE" & 
                          databruteTOT$synerg_id=="AUCUN" &
                          databruteTOT$dose!=0,]

#counting the number of sampling according to their status
datacount <- cbind(
  "dep_ID" = row.names(table(
    databrute$parcel_dpt,
    databrute$rslt_03, exclude =
      ""
  )),
  "Resistant" = 0,
  "Sensible" = table(databrute$parcel_dpt,
                     as.numeric(databrute$rslt_03)>0, exclude = "")[,1],
  "Total" = rowSums(table(
    databrute$parcel_dpt,
    databrute$rslt_03>0, exclude = ""
  ))
)
#export the count data by departement
write.table(
  datacount,
  file = "output/datacountAZOX.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#extract the department coordinates
ind_list <- departe$INSEE_DEP
coorddep <- data.frame(
  "longitude" = departe@polygons[1][[1]]@labpt[1],
  "latitude" = departe@polygons[1][[1]]@labpt[2]
)
for (i in 2:length(ind_list)) {
  coorddep <- rbind(
    coorddep,
    cbind(
      "longitude" = departe@polygons[i][[1]]@labpt[1],
      "latitude" = departe@polygons[i][[1]]@labpt[2]
    )
  )
}
coorddep <- cbind("dep_ID" = ind_list, coorddep)

#combining the information of the geographic
data2map <- merge(datacount, coorddep, by = "dep_ID")
data2map <- data2map[order(as.numeric(as.character(data2map$Total)),
                           decreasing = TRUE), ]

#defining the colors of the pies
colovec <- c(brewer.pal(9, "Reds")[6])

#actual plotting
op <- par(mar = c(0, 0, 2, 0))
plot(departe,main="AZOXYSTROBINE")
draw.pie(
  x = data2map$longitude,
  y = data2map$latitude,
  z = cbind((as.numeric(
    as.character(data2map$Resistant)
  )),
  (as.numeric(
    as.character(data2map$Sensible)
  ))),
  col = colovec,         #colors of the pie
  lty = 1,               #line type of the pie
  border = "grey60",     #color of the border of the pie
  lwd = 0.1,             #control the width of the border
  radius = 35000, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels = NA,
  scale=FALSE # should the radius be scaled according to sample size
)

#writing the number of samples for each departement
text(
  x = data2map$longitude,
  y = data2map$latitude,
  col = "black",
  font = 2,
  labels = as.character(data2map$Total),
  cex = 1.5              #size of the text
)
par(op)

#export the map to a pdf file 7 x 7 inches (for examples)

#boxplot for the resistant populations
boxplot(as.numeric(databrute$rslt_03[databrute$rslt_03!=0]),
        boxwex=0.5,las=1,ylim=c(0,140),col="cornflowerblue",
        main=paste("Azoxystrobine\n(n=",
                   length(databrute$rslt_03[databrute$rslt_03!=0]),
                   "/",
                   length(databrute$rslt_03),")",sep=""),
        ylab="% germination")
abline(h=mean(as.numeric(databrute$rslt_03[databrute$rslt_03!=0])),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(databrute$rslt_03[databrute$rslt_03!=0]),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.8),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)

#export to a pdf file 3 x 7 inches (for examples)

#distribution of the % of germination at the DD
plot(as.numeric(databrute$rslt_03[order(as.numeric(databrute$rslt_03))]),
     bg="cornflowerblue",pch=21,cex=1.5,las=1,ylim=c(0,130),
     ylab="% germination",
     main="Azoxystrobine")

#export to a pdf file 7 x 5 inches (for examples)


##############################################################################/
#Mefentrifluconazole map 3 resistance categories####
##############################################################################/

#defining the colors of the pies
colovec <- c(
  rgb(200,200,00, max = 255, alpha = 240),
  rgb(250,150,50, max = 255, alpha = 240),
  rgb(255,50,50, max = 255, alpha = 240)
)
#changing the projection of the map
departe.wgs <- spTransform(departe,
                           CRS("+proj=longlat +datum=WGS84"))

#actual plotting
op <- par(mar = c(0, 0, 2, 0))
plot(departe.wgs,main="MEFENTRIFLUCONAZOLE")
draw.pie(
  x = as.numeric(databruteTOT$gps_long),
  y = as.numeric(databruteTOT$gps_lat),
  z = cbind(
    as.numeric(databruteTOT$FR.30),
    as.numeric(databruteTOT$FR30.100),
    as.numeric(databruteTOT$FR.100)
    ),
  col = colovec,         #colors of the pie
  lty = 1,               #line type of the pie
  border = NA,     #color of the border of the pie
  lwd = 0.01,             #control the width of the border
  radius = 0.08, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels = NA,
  scale=FALSE # should the radius be scaled according to sample size
)
par(op)


##############################################################################/
#END
##############################################################################/