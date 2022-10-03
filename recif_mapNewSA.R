##############################################################################/
##############################################################################/
#Maps for population germination tests for old SA
##############################################################################/
##############################################################################/

#loading the dataset and the necessary library
source("recif_load.R")

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