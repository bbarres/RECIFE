##############################################################################/
##############################################################################/
#Code for producing a map of Cercospora beticola resistance
##############################################################################/
##############################################################################/

#loading the dataset and the necessary library
source("recif_load.R")


##############################################################################/
#New versions old active substances####
##############################################################################/

##############################################################################/
#Maps for population germination tests for old SA####
##############################################################################/

#load the resistance results for the 2019-2020 campaign
oldSA<-read.delim(
  "data/data_DC_AZ_FH_Carb_2019_2020.txt",
  header=TRUE,
  sep="\t"
)

oldSA<-oldSA[!is.na(oldSA$gps_lat),]

#turning this dataframe into a spatial dataframe (wgs84)
oldSA.wgs<-SpatialPointsDataFrame(coords=oldSA[,c("gps_long","gps_lat")],
                                  data=oldSA,
                                  proj4string=CRS("+proj=longlat +datum=WGS84")
)
oldSA<-spTransform(oldSA.wgs,CRS("+init=epsg:2154"))


##############################################################################/
#Fentine Hydroxyde MAP####
##############################################################################/

oldprod<-oldSA[oldSA$pest_sa_id=="FENTINE HYDROXYDE" & 
                 oldSA$dose!=0 & oldSA$repet==1 &
                 oldSA$milieu_cult_id!="CV8",]
oldprod$rslt_03[oldprod$rslt_03>100]<-100

#splitting the continuous percentage of germination in categories
oldprod$catgerm<-cut(oldprod$rslt_03,
                     breaks=c(0.00,0.001,5,10,20,30,40,50,
                              max(oldprod$rslt_03)),
                     include.lowest=TRUE)
nomCat<-c("0 %","]0-5] %","]5-10] %","]10-20] %","]20-30] %",
          "]30-40] %","]40-50] %",">50 %")
#defining the colors of the points
levels(oldprod$catgerm)<-brewer.pal(11,"RdYlGn")[8:1]
#defining sampling year
oldprod$year<-as.numeric(substr(oldprod$prelvt_id,1,2))+2

#actual plotting
nf<-layout(matrix(c(1,1,1,2,
                    1,1,1,2,
                    1,1,1,3,
                    4,4,4,3),4,4,byrow=TRUE))
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.1,main="",border="grey70")
plot(REG_SHP.1,lwd=2,add=TRUE)
points(
  x=as.numeric(oldprod$gps_long),
  y=as.numeric(oldprod$gps_lat),
  bg=as.character(oldprod$catgerm),  #colors of the points
  pch=oldprod$year,                  #plotting character
  cex=2                              #size of the points
)

legend(110000,7150000,title="Germination\nclasses",
       legend=nomCat,cex=1,pt.cex=1.8,
       y.intersp=0.7,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(280000,7150000,legend=c("2019","2020"),cex=1,pt.cex=1.6,
       y.intersp=0.7,x.intersp=0.8,title="Année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")
par(op)

#histogram of the distribution of the % of germination
op<-par(mar=c(6.1,4.1,2.1,2.1))
hist(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     breaks=20,bty="l",freq=FALSE,las=1,main="",xlim=c(0,100),
     col=brewer.pal(11,"RdYlGn")[c(8,6,5,5,4,4,3,3,2,2,rep(1,10))],
     xlab="Germination classes",ylab="Pourcentage")
box(bty="l")
par(op)

#boxplot for the resistant populations
op<-par(mar=c(3.1,4.1,4.1,2.1))
boxplot(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0]),
        boxwex=0.4,las=1,ylim=c(0,100),col="transparent",
        main=paste("Germination des résistants\n(n=",
                   length(oldprod$rslt_03[oldprod$rslt_03!=0]),
                   "/",
                   length(oldprod$rslt_03),")",sep=""),
        ylab="% germination",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0])),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0]),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)
par(op)

#distribution of the % of germination at the DD
op<-par(mar=c(3.1,6.1,0,2.1))
plot(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     bg=as.character(oldprod$catgerm[order(as.numeric(oldprod$rslt_03))]),
     pch=oldprod$year[order(as.numeric(oldprod$rslt_03))],
     cex=1.5,las=1,ylim=c(0,100),
     ylab="% germination",xlab="",
     main="",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0])),
       col="red",lty=2,lwd=3)
par(op)

#export to a pdf file 10 x 7 inches


##############################################################################/
#Carbendazime MAP####
##############################################################################/

oldprod<-oldSA[oldSA$pest_sa_id=="CARBENDAZIME" & 
                 oldSA$dose!=0,]
oldprod$rslt_03[oldprod$rslt_03>100]<-100

#splitting the continuous percentage of germination in categories
oldprod$catgerm<-cut(oldprod$rslt_03,
                     breaks=c(0.00,0.001,5,10,20,30,40,50,
                              max(oldprod$rslt_03)),
                     include.lowest=TRUE)
nomCat<-c("0 %","]0-5] %","]5-10] %","]10-20] %","]20-30] %",
          "]30-40] %","]40-50] %",">50 %")
#defining the colors of the points
levels(oldprod$catgerm)<-brewer.pal(11,"RdYlGn")[8:1]
#defining sampling year
oldprod$year<-as.numeric(substr(oldprod$prelvt_id,1,2))+2

#actual plotting
nf<-layout(matrix(c(1,1,1,2,
                    1,1,1,2,
                    1,1,1,3,
                    4,4,4,3),4,4,byrow=TRUE))
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.1,main="",border="grey70")
plot(REG_SHP.1,lwd=2,add=TRUE)
points(
  x=as.numeric(oldprod$gps_long),
  y=as.numeric(oldprod$gps_lat),
  bg=as.character(oldprod$catgerm),  #colors of the points
  pch=oldprod$year,                  #plotting character
  cex=2                              #size of the points
)

legend(110000,7150000,title="Germination\nclasses",
       legend=nomCat,cex=1,pt.cex=1.8,
       y.intersp=0.7,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(280000,7150000,legend=c("2019","2020"),cex=1,pt.cex=1.6,
       y.intersp=0.7,x.intersp=0.8,title="Année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")
par(op)

#histogram of the distribution of the % of germination
op<-par(mar=c(6.1,4.1,2.1,2.1))
hist(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     breaks=20,bty="l",freq=FALSE,las=1,main="",xlim=c(0,100),
     col=brewer.pal(11,"RdYlGn")[c(8,6,5,5,4,4,3,3,2,2,rep(1,10))],
     xlab="Germination classes",ylab="Pourcentage")
box(bty="l")
par(op)

#boxplot for the resistant populations
op<-par(mar=c(3.1,4.1,4.1,2.1))
boxplot(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0]),
        boxwex=0.4,las=1,ylim=c(0,100),col="transparent",
        main=paste("Germination des résistants\n(n=",
                   length(oldprod$rslt_03[oldprod$rslt_03!=0]),
                   "/",
                   length(oldprod$rslt_03),")",sep=""),
        ylab="% germination",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0])),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0]),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)
par(op)

#distribution of the % of germination at the DD
op<-par(mar=c(3.1,6.1,0,2.1))
plot(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     bg=as.character(oldprod$catgerm[order(as.numeric(oldprod$rslt_03))]),
     pch=oldprod$year[order(as.numeric(oldprod$rslt_03))],
     cex=1.5,las=1,ylim=c(0,100),
     ylab="% germination",xlab="",
     main="",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0])),
       col="red",lty=2,lwd=3)
par(op)

#export to a pdf file 10 x 7 inches


##############################################################################/
#Azoxystrobine R total MAP####
##############################################################################/

oldprod<-oldSA[oldSA$pest_sa_id=="AZOXYSTROBINE" & 
                 oldSA$synerg_id=="AUCUN" &
                 oldSA$dose!=0,]
oldprod$rslt_03[oldprod$rslt_03>100]<-100

#splitting the continuous percentage of germination in categories
oldprod$catgerm<-cut(oldprod$rslt_03,
                     breaks=c(0.00,0.001,5,10,20,30,40,50,
                              max(oldprod$rslt_03)),
                     include.lowest=TRUE)
nomCat<-c("0 %","]0-5] %","]5-10] %","]10-20] %","]20-30] %",
          "]30-40] %","]40-50] %",">50 %")
#defining the colors of the points
levels(oldprod$catgerm)<-brewer.pal(11,"RdYlGn")[8:1]
#defining sampling year
oldprod$year<-as.numeric(substr(oldprod$prelvt_id,1,2))+2

#actual plotting
nf<-layout(matrix(c(1,1,1,2,
                    1,1,1,2,
                    1,1,1,3,
                    4,4,4,3),4,4,byrow=TRUE))
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.1,main="",border="grey70")
plot(REG_SHP.1,lwd=2,add=TRUE)
points(
  x=as.numeric(oldprod$gps_long),
  y=as.numeric(oldprod$gps_lat),
  bg=as.character(oldprod$catgerm),  #colors of the points
  pch=oldprod$year,                  #plotting character
  cex=2                              #size of the points
)

legend(110000,7150000,title="Germination\nclasses",
       legend=nomCat,cex=1,pt.cex=1.8,
       y.intersp=0.7,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(280000,7150000,legend=c("2019","2020"),cex=1,pt.cex=1.6,
       y.intersp=0.7,x.intersp=0.8,title="Année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")
par(op)

#histogram of the distribution of the % of germination
op<-par(mar=c(6.1,4.1,2.1,2.1))
hist(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     breaks=20,bty="l",freq=FALSE,las=1,main="",xlim=c(0,100),
     col=brewer.pal(11,"RdYlGn")[c(8,6,5,5,4,4,3,3,2,2,rep(1,10))],
     xlab="Germination classes",ylab="Pourcentage")
box(bty="l")
par(op)

#boxplot for the resistant populations
op<-par(mar=c(3.1,4.1,4.1,2.1))
boxplot(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0]),
        boxwex=0.4,las=1,ylim=c(0,100),col="transparent",
        main=paste("Germination des résistants\n(n=",
                   length(oldprod$rslt_03[oldprod$rslt_03!=0]),
                   "/",
                   length(oldprod$rslt_03),")",sep=""),
        ylab="% germination",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0])),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0]),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)
par(op)

#distribution of the % of germination at the DD
op<-par(mar=c(3.1,6.1,0,2.1))
plot(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     bg=as.character(oldprod$catgerm[order(as.numeric(oldprod$rslt_03))]),
     pch=oldprod$year[order(as.numeric(oldprod$rslt_03))],
     cex=1.5,las=1,ylim=c(0,100),
     ylab="% germination",xlab="",
     main="",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0])),
       col="red",lty=2,lwd=3)
par(op)

#export to a pdf file 10 x 7 inches


##############################################################################/
#Azoxystrobine TSR only MAP####
##############################################################################/

oldprod<-oldSA[oldSA$pest_sa_id=="AZOXYSTROBINE" & 
                 oldSA$synerg_id=="SHAM" &
                 oldSA$dose!=0,]
oldprod$rslt_03[oldprod$rslt_03>100]<-100

#splitting the continuous percentage of germination in categories
oldprod$catgerm<-cut(oldprod$rslt_03,
                     breaks=c(0.00,0.001,5,10,20,30,40,50,
                              max(oldprod$rslt_03)),
                     include.lowest=TRUE)
nomCat<-c("0 %","]0-5] %","]5-10] %","]10-20] %","]20-30] %",
          "]30-40] %","]40-50] %",">50 %")
#defining the colors of the points
levels(oldprod$catgerm)<-brewer.pal(11,"RdYlGn")[8:1]
#defining sampling year
oldprod$year<-as.numeric(substr(oldprod$prelvt_id,1,2))+2

#actual plotting
nf<-layout(matrix(c(1,1,1,2,
                    1,1,1,2,
                    1,1,1,3,
                    4,4,4,3),4,4,byrow=TRUE))
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.1,main="",border="grey70")
plot(REG_SHP.1,lwd=2,add=TRUE)
points(
  x=as.numeric(oldprod$gps_long),
  y=as.numeric(oldprod$gps_lat),
  bg=as.character(oldprod$catgerm),  #colors of the points
  pch=oldprod$year,                  #plotting character
  cex=2                              #size of the points
)

legend(110000,7150000,title="Germination\nclasses",
       legend=nomCat,cex=1,pt.cex=1.8,
       y.intersp=0.7,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(280000,7150000,legend=c("2019","2020"),cex=1,pt.cex=1.6,
       y.intersp=0.7,x.intersp=0.8,title="Année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")
par(op)

#histogram of the distribution of the % of germination
op<-par(mar=c(6.1,4.1,2.1,2.1))
hist(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     breaks=20,bty="l",freq=FALSE,las=1,main="",xlim=c(0,100),
     col=brewer.pal(11,"RdYlGn")[c(3,3,2,2,rep(1,10))],
     xlab="Germination classes",ylab="Pourcentage")
box(bty="l")
par(op)

#boxplot for the resistant populations
op<-par(mar=c(3.1,4.1,4.1,2.1))
boxplot(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0]),
        boxwex=0.4,las=1,ylim=c(0,100),col="transparent",
        main=paste("Germination des résistants\n(n=",
                   length(oldprod$rslt_03[oldprod$rslt_03!=0]),
                   "/",
                   length(oldprod$rslt_03),")",sep=""),
        ylab="% germination",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0])),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0]),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)
par(op)

#distribution of the % of germination at the DD
op<-par(mar=c(3.1,6.1,0,2.1))
plot(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     bg=as.character(oldprod$catgerm[order(as.numeric(oldprod$rslt_03))]),
     pch=oldprod$year[order(as.numeric(oldprod$rslt_03))],
     cex=1.5,las=1,ylim=c(0,100),
     ylab="% germination",xlab="",
     main="",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03[oldprod$rslt_03!=0])),
       col="red",lty=2,lwd=3)
par(op)

#export to a pdf file 10 x 7 inches






##############################################################################/
#Maps for population germination tests for newer active substances####
##############################################################################/

#load the resistance results for the 2019-2020 campaign

#load the resistance results for the 2019-2020 campaign
newSA<-read.delim(
  "data/data_DC_classes_2019_2020.txt",
  header=TRUE,
  sep="\t"
)

#removing samples without spatial coordinates
newSA<-newSA[!is.na(newSA$gps_lat),]

#turning this dataframe into a spatial dataframe (wgs84)
newSA.wgs<-SpatialPointsDataFrame(coords=newSA[,c("gps_long","gps_lat")],
                                  data=newSA,
                                  proj4string=CRS("+proj=longlat +datum=WGS84")
)
newSA<-spTransform(newSA.wgs,CRS("+init=epsg:2154"))

#defining the colors of the pies with transparency
colovec<-brewer.pal(11,"RdYlGn")[c(5,3,1)]
#defining another color vector
colovec<-c(
  rgb(254,224,139,max=255,alpha=210),
  rgb(244,109,67,max=255,alpha=210),
  rgb(165,0,38,max=255,alpha=210)
)


##############################################################################/
#Mefentrifluconazole map 3 resistance categories####
##############################################################################/

#limiting the data set to mefentrifluconazole
newprod<-newSA[newSA$pest_sa_id=="MEFENTRIFLUCONAZOLE",]

#actual plotting
nf<-layout(matrix(c(1,1,1,2,
                    1,1,1,2,
                    1,1,1,3,
                    1,1,1,3),4,4,byrow=TRUE))
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.2,main="",border="grey70")
plot(REG_SHP.2,lwd=2,add=TRUE)
draw.pie(
  x=as.numeric(newprod$gps_long),
  y=as.numeric(newprod$gps_lat),
  z=cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100)
  ),
  col=colovec,         #colors of the pie
  lty=1,               #line type of the pie
  border="transparent",     #color of the border of the pie
  lwd=0.01,             #control the width of the border
  radius=5000, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels=NA,
  scale=FALSE # should the radius be scaled according to sample size
)
legend(467000,7100000,title="Classes de facteur\nde résistance",
       legend=c("FR<30","30<FR<100","100<FR"),
       cex=1.8,pt.cex=3.5,
       y.intersp=0.9,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=colovec,
       bg="transparent",bty="n")
par(op)

op<-par(mar=c(2.1,0,0,0))
pie(colMeans(cbind(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100))),
  labels=round(colMeans(cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100))),1),
  col=colovec,main=NA,font=1,cex=1.2,
  border="transparent",radius=0.95)
title("Fréquence globale",line=-2.5,cex.main=1.5)
par(op)

op<-par(mar=c(3.1,1.6,0,2.1))
boxplot(cbind(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),
  col=colovec,frame=FALSE,las=1,boxwex=0.4,
  names=c("FR<30","30<FR<100","100<FR"))
title("Distribution des fréquences",
      line=2,xpd=NA,cex.main=1.5)
box(bty="l")
stripchart(list(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),
  cex=1,pch=19,at=c(1:3),
  col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
  method="jitter",jitter=0.1,add=TRUE)
par(op)

#export to .pdf 10 x 7 inches


##############################################################################/
#tétraconazole map 3 resistance categories####
##############################################################################/

#limiting the data set to mefentrifluconazole
newprod<-newSA[newSA$pest_sa_id=="TETRACONAZOLE",]

#actual plotting
nf<-layout(matrix(c(1,1,1,2,
                    1,1,1,2,
                    1,1,1,3,
                    1,1,1,3),4,4,byrow=TRUE))
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.2,main="",border="grey70")
plot(REG_SHP.2,lwd=2,add=TRUE)
draw.pie(
  x=as.numeric(newprod$gps_long),
  y=as.numeric(newprod$gps_lat),
  z=cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100)
  ),
  col=colovec,         #colors of the pie
  lty=1,               #line type of the pie
  border="transparent",     #color of the border of the pie
  lwd=0.01,             #control the width of the border
  radius=5000, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels=NA,
  scale=FALSE # should the radius be scaled according to sample size
)
legend(467000,7100000,title="Classes de facteur\nde résistance",
       legend=c("FR<30","30<FR<100","100<FR"),
       cex=1.8,pt.cex=3.5,
       y.intersp=0.9,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=colovec,
       bg="transparent",bty="n")
par(op)

op<-par(mar=c(2.1,0,0,0))
pie(colMeans(cbind(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),na.rm=TRUE),
  labels=round(colMeans(cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100)),na.rm=TRUE),1),
  col=colovec,main=NA,font=1,cex=1.2,
  border="transparent",radius=0.95)
title("Fréquence globale",line=-2.5,cex.main=1.5)
par(op)

op<-par(mar=c(3.1,1.6,0,2.1))
boxplot(cbind(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),
  col=colovec,frame=FALSE,las=1,boxwex=0.4,
  names=c("FR<30","30<FR<100","100<FR"))
title("Distribution des fréquences",
      line=2,xpd=NA,cex.main=1.5)
box(bty="l")
stripchart(list(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),
  cex=1,pch=19,at=c(1:3),
  col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
  method="jitter",jitter=0.1,add=TRUE)
par(op)

#export to .pdf 10 x 7 inches


##############################################################################/
#difenoconazole map 3 resistance categories####
##############################################################################/

#limiting the data set to mefentrifluconazole
newprod<-newSA[newSA$pest_sa_id=="DIFENOCONAZOLE",]

#actual plotting
nf<-layout(matrix(c(1,1,1,2,
                    1,1,1,2,
                    1,1,1,3,
                    1,1,1,3),4,4,byrow=TRUE))
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.2,main="",border="grey70")
plot(REG_SHP.2,lwd=2,add=TRUE)
draw.pie(
  x=as.numeric(newprod$gps_long),
  y=as.numeric(newprod$gps_lat),
  z=cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100)
  ),
  col=colovec,         #colors of the pie
  lty=1,               #line type of the pie
  border="transparent",     #color of the border of the pie
  lwd=0.01,             #control the width of the border
  radius=5000, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels=NA,
  scale=FALSE # should the radius be scaled according to sample size
)
legend(467000,7100000,title="Classes de facteur\nde résistance",
       legend=c("FR<30","30<FR<100","100<FR"),
       cex=1.8,pt.cex=3.5,
       y.intersp=0.9,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=colovec,
       bg="transparent",bty="n")
par(op)

op<-par(mar=c(2.1,0,0,0))
pie(colMeans(cbind(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),na.rm=TRUE),
  labels=round(colMeans(cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100)),na.rm=TRUE),1),
  col=colovec,main=NA,font=1,cex=1.2,clockwise=FALSE,
  border="transparent",radius=0.95,init.angle=20)
title("Fréquence globale",line=-2.5,cex.main=1.5)
par(op)

op<-par(mar=c(3.1,1.6,0,2.1))
boxplot(cbind(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),
  col=colovec,frame=FALSE,las=1,boxwex=0.4,
  names=c("FR<30","30<FR<100","100<FR"))
title("Distribution des fréquences",
      line=2,xpd=NA,cex.main=1.5)
box(bty="l")
stripchart(list(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),
  cex=1,pch=19,at=c(1:3),
  col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
  method="jitter",jitter=0.1,add=TRUE)
par(op)

#export to .pdf 10 x 7 inches


##############################################################################/
#prothioconazole-desthio map 3 resistance categories####
##############################################################################/

#limiting the data set to mefentrifluconazole
newprod<-newSA[newSA$pest_sa_id=="PROTHIOCONAZOLE-DESTHIO",]

#actual plotting
nf<-layout(matrix(c(1,1,1,2,
                    1,1,1,2,
                    1,1,1,3,
                    1,1,1,3),4,4,byrow=TRUE))
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.2,main="",border="grey70")
plot(REG_SHP.2,lwd=2,add=TRUE)
draw.pie(
  x=as.numeric(newprod$gps_long),
  y=as.numeric(newprod$gps_lat),
  z=cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100)
  ),
  col=colovec,         #colors of the pie
  lty=1,               #line type of the pie
  border="transparent",     #color of the border of the pie
  lwd=0.01,             #control the width of the border
  radius=5000, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels=NA,
  scale=FALSE # should the radius be scaled according to sample size
)
legend(467000,7100000,title="Classes de facteur\nde résistance",
       legend=c("FR<30","30<FR<100","100<FR"),
       cex=1.8,pt.cex=3.5,
       y.intersp=0.9,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=colovec,
       bg="transparent",bty="n")
par(op)

op<-par(mar=c(2.1,0,0,0))
pie(colMeans(cbind(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),na.rm=TRUE),
  labels=round(colMeans(cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100)),na.rm=TRUE),1),
  col=colovec,main=NA,font=1,cex=1.2,clockwise=TRUE,
  border="transparent",radius=0.95)
title("Fréquence globale",line=-2.5,cex.main=1.5)
par(op)

op<-par(mar=c(3.1,1.6,0,2.1))
boxplot(cbind(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),
  col=colovec,frame=FALSE,las=1,boxwex=0.4,
  names=c("FR<30","30<FR<100","100<FR"))
title("Distribution des fréquences",
      line=2,xpd=NA,cex.main=1.5)
box(bty="l")
stripchart(list(
  as.numeric(newprod$FR.30),
  as.numeric(newprod$FR30.100),
  as.numeric(newprod$FR.100)),
  cex=1,pch=19,at=c(1:3),
  col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
  method="jitter",jitter=0.1,add=TRUE)
par(op)

#export to .pdf 10 x 7 inches












##############################################################################/
#END
##############################################################################/



























#load geographical data
load("data/commu.RData")
load("data/arrond.RData")
load("data/departe.RData")
load("data/regionsLight.RData")
load("data/departeLight.RData")



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

#load the resistance results for the 2019-2020 campaign
databruteTOT <- read.delim(
  "data/data_DC_AZ_FH_Carb_2019_2020.txt",
  header = TRUE,
  sep = "\t"
)

# #first we merge the resistance table with the commune info
# databrute<-merge(databrute,db_commu,by.x="parcel_cmne",by.y="NOM_COM_M")
# #in order to acces to the arrondissement ID, we create an individual ID for 
# #each arrondissement combining INSEE_DEP and INSEE_ARR
# Raox_list$DEPARR<-paste(Raox_list$INSEE_DEP,Raox_list$INSEE_ARR)
# #then we merge the resistance table with the arrondissement info
# Raox_list<-merge(Raox_list,db_arrond,by.x="DEPARR",by.y="DEPARR")





##############################################################################/
#Mefentrifluconazole map 3 resistance categories####
##############################################################################/

#load the resistance results for the 2019-2020 campaign
databruteTOT <- read.delim(
  "data/data_DC_classes_MEFENTRI_2019_2020.txt",
  header = TRUE,
  sep = "\t"
)

#changing the projection of the map
departeLight.wgs <- spTransform(departeLight,
                                CRS("+proj=longlat +datum=WGS84"))
regionsLight.wgs <- spTransform(regionsLight,
                                CRS("+proj=longlat +datum=WGS84"))
#defining the colors of the pies with transparency
colovec <- c(
  rgb(200,200,00, max = 255, alpha = 240),
  rgb(250,150,50, max = 255, alpha = 240),
  rgb(255,50,50, max = 255, alpha = 240)
)
#defining another color vector
colovec<-brewer.pal(11,"RdYlGn")[c(5,3,1)]

#actual plotting
op <- par(mar = c(0, 0, 2, 0))
plot(departeLight.wgs,main="MEFENTRIFLUCONAZOLE",border="grey70")
plot(regionsLight.wgs,lwd=2,add=TRUE)
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
  border = "black",     #color of the border of the pie
  lwd = 0.01,             #control the width of the border
  radius = 0.08, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels = NA,
  scale=FALSE # should the radius be scaled according to sample size
)
par(op)


##############################################################################/
#Prothio-destio map 3 resistance categories####
##############################################################################/

#load the resistance results for the 2019-2020 campaign
databruteTOT <- read.delim(
  "data/data_DC_classes_PROTHIO_2019_2020.txt",
  header = TRUE,
  sep = "\t"
)

#changing the projection of the map
departeLight.wgs <- spTransform(departeLight,
                                CRS("+proj=longlat +datum=WGS84"))
regionsLight.wgs <- spTransform(regionsLight,
                                CRS("+proj=longlat +datum=WGS84"))
#defining the colors of the pies with transparency
colovec <- c(
  rgb(200,200,00, max = 255, alpha = 240),
  rgb(250,150,50, max = 255, alpha = 240),
  rgb(255,50,50, max = 255, alpha = 240)
)
#defining another color vector
colovec<-brewer.pal(11,"RdYlGn")[c(5,3,1)]

#actual plotting
op <- par(mar = c(0, 0, 2, 0))
plot(departeLight.wgs,main="PROTHIO-DESTHIO",border="grey70")
plot(regionsLight.wgs,lwd=2,add=TRUE)
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
  border = "black",     #color of the border of the pie
  lwd = 0.01,             #control the width of the border
  radius = 0.08, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels = NA,
  scale=FALSE # should the radius be scaled according to sample size
)
par(op)


##############################################################################/
#Difénoconazole map 3 resistance categories####
##############################################################################/

#load the resistance results for the 2019-2020 campaign
databruteTOT <- read.delim(
  "data/data_DC_classes_DIFENO_2019_2020.txt",
  header = TRUE,
  sep = "\t"
)

#changing the projection of the map
departeLight.wgs <- spTransform(departeLight,
                                CRS("+proj=longlat +datum=WGS84"))
regionsLight.wgs <- spTransform(regionsLight,
                                CRS("+proj=longlat +datum=WGS84"))
#defining the colors of the pies with transparency
colovec <- c(
  rgb(200,200,00, max = 255, alpha = 240),
  rgb(250,150,50, max = 255, alpha = 240),
  rgb(255,50,50, max = 255, alpha = 240)
)
#defining another color vector
colovec<-brewer.pal(11,"RdYlGn")[c(5,3,1)]

#actual plotting
op <- par(mar = c(0, 0, 2, 0))
plot(departeLight.wgs,main="DIFENOCONAZOLE",border="grey70")
plot(regionsLight.wgs,lwd=2,add=TRUE)
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
  border = "black",     #color of the border of the pie
  lwd = 0.01,             #control the width of the border
  radius = 0.08, #(sqrt(as.numeric(as.character(data2map$Total))) * 16000), 
  #this number control the radius of the pies
  labels = NA,
  scale=FALSE # should the radius be scaled according to sample size
)
par(op)


##############################################################################/
#Tétraconazole map 3 resistance categories####
##############################################################################/

#load the resistance results for the 2019-2020 campaign
databruteTOT <- read.delim(
  "data/data_DC_classes_TETRACO_2019_2020.txt",
  header = TRUE,
  sep = "\t"
)

#changing the projection of the map
departeLight.wgs <- spTransform(departeLight,
                                CRS("+proj=longlat +datum=WGS84"))
regionsLight.wgs <- spTransform(regionsLight,
                                CRS("+proj=longlat +datum=WGS84"))
#defining the colors of the pies with transparency
colovec <- c(
  rgb(200,200,00, max = 255, alpha = 240),
  rgb(250,150,50, max = 255, alpha = 240),
  rgb(255,50,50, max = 255, alpha = 240)
)
#defining another color vector
colovec<-brewer.pal(11,"RdYlGn")[c(5,3,1)]

#actual plotting
op <- par(mar = c(0, 0, 2, 0))
plot(departeLight.wgs,main="TETRACONAZOLE",border="grey70")
plot(regionsLight.wgs,lwd=2,add=TRUE)
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
  border = "black",     #color of the border of the pie
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