##############################################################################/
##############################################################################/
#Maps for population germination tests for old SA
##############################################################################/
##############################################################################/

#loading the dataset and the necessary library
source("recif_load.R")

oldSA<-oldSA[!is.na(oldSA$gps_lat),]

#turning this dataframe into a spatial dataframe (wgs84)
oldSA.wgs<-SpatialPointsDataFrame(coords=oldSA[,c("gps_long","gps_lat")],
                                  data=oldSA,
                                  proj4string=CRS("+proj=longlat +datum=WGS84")
)
oldSA.wgs<-st_as_sf(oldSA.wgs)
oldSA.wgs<-st_transform(oldSA.wgs,CRS("+init=epsg:2154"))


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

pdf(file="output/map_fentineHydrox.pdf",width=7,height=10)
#actual plotting
nf<-layout(matrix(c(1,1,1,1,
                    1,1,1,1,
                    1,1,1,1,
                    2,2,2,2,
                    3,3,3,4,
                    3,3,3,4),6,4,byrow=TRUE))
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.1,main="",border="grey70")
plot(REG_SHP.1,lwd=2,add=TRUE)
plot(oldprod,
     bg=as.character(oldprod$catgerm),col="black",
     pch=oldprod$year,
     cex=2,
     add=TRUE)
legend(160000,7150000,title="classe de fréquence\nde résistance",
       legend=nomCat,cex=1,pt.cex=1.8,
       y.intersp=1,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(360000,7150000,legend=c("2019","2020"),cex=1,pt.cex=1.6,
       y.intersp=1,x.intersp=0.8,title="année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")
text(123000,7150000,labels="A",cex=3,font=2)
par(op)

#distribution of the % of germination at the DD
op<-par(mar=c(5.1,6.1,0,2.1))
plot(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     bg=as.character(oldprod$catgerm[order(as.numeric(oldprod$rslt_03))]),
     pch=oldprod$year[order(as.numeric(oldprod$rslt_03))],
     cex=1.5,las=1,ylim=c(0,100),cex.axis=1.3,cex.lab=1.5,
     ylab="",xlab="population",
     main="",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03)),
       col="red",lty=2,lwd=3)
mtext(text="B",side=3,cex=2,at=0,font=2,las=0,adj=3.55,line=1)
mtext(text="% resistance",side=2,cex=1,font=2,line=4)
par(op)

#histogram of the distribution of the % of germination
op<-par(mar=c(5.1,6.1,3.1,5.1))
hist(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     breaks=20,bty="l",freq=FALSE,las=1,main="",xlim=c(0,100),cex.lab=1.5,
     col=brewer.pal(11,"RdYlGn")[c(8,6,5,5,4,4,3,3,2,2,rep(1,10))],
     cex.axis=1.3,
     xlab="classe de fréquence de résistance",ylab="")
box(bty="l")
mtext(text="C",side=3,cex=2,at=0,font=2,las=0,adj=3.45,line=1)
mtext(text="fréquence",side=2,cex=1,font=2,line=4)
par(op)

#boxplot for the resistant populations
op<-par(mar=c(5.1,4.1,3.1,2.1))
boxplot(as.numeric(oldprod$rslt_03),
        boxwex=0.4,las=1,ylim=c(0,100),col="transparent",
        # main=paste("Germination des résistants\n(n=",
        #            length(oldprod$rslt_03[oldprod$rslt_03!=0]),
        #            "/",
        #            length(oldprod$rslt_03),")",sep=""),
        ylab="",frame=FALSE,cex.axis=1.3,cex.lab=1.5)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03)),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(oldprod$rslt_03),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)
mtext(text="D",side=3,cex=2,at=0,font=2,las=0,adj=2,line=1)
mtext(text="% resistance",side=2,cex=1,font=2,line=4)
par(op)
#export to a pdf file 10 x 7 inches
dev.off()

min(as.numeric(oldprod$rslt_03))
max(as.numeric(oldprod$rslt_03))
mean(as.numeric(oldprod$rslt_03))
median(as.numeric(oldprod$rslt_03))
length(as.numeric(oldprod$rslt_03))


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

pdf(file="output/map_carbendaz.pdf",width=7,height=10)
#actual plotting
nf<-layout(matrix(c(1,1,1,1,
                    1,1,1,1,
                    1,1,1,1,
                    2,2,2,2,
                    3,3,3,4,
                    3,3,3,4),6,4,byrow=TRUE))
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
legend(160000,7150000,title="classe de fréquence\nde résistance",
       legend=nomCat,cex=1,pt.cex=1.8,
       y.intersp=1,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(360000,7150000,legend=c("2019","2020"),cex=1,pt.cex=1.6,
       y.intersp=1,x.intersp=0.8,title="année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")
text(123000,7150000,labels="A",cex=3,font=2)
par(op)

#distribution of the % of germination at the DD
op<-par(mar=c(5.1,6.1,0,2.1))
plot(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     bg=as.character(oldprod$catgerm[order(as.numeric(oldprod$rslt_03))]),
     pch=oldprod$year[order(as.numeric(oldprod$rslt_03))],
     cex=1.5,las=1,ylim=c(0,100),cex.axis=1.3,cex.lab=1.5,
     ylab="",xlab="population",
     main="",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03)),
       col="red",lty=2,lwd=3)
mtext(text="B",side=3,cex=2,at=0,font=2,las=0,adj=3.55,line=1)
mtext(text="% resistance",side=2,cex=1,font=2,line=4)
par(op)

#histogram of the distribution of the % of germination
op<-par(mar=c(5.1,6.1,3.1,5.1))
hist(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     breaks=20,bty="l",freq=FALSE,las=1,main="",xlim=c(0,100),cex.lab=1.5,
     col=brewer.pal(11,"RdYlGn")[c(8,6,5,5,4,4,3,3,2,2,rep(1,10))],
     cex.axis=1.3,
     xlab="classe de fréquence de résistance",ylab="")
box(bty="l")
mtext(text="C",side=3,cex=2,at=0,font=2,las=0,adj=3.45,line=1)
mtext(text="fréquence",side=2,cex=1,font=2,line=4)
par(op)

#boxplot for the resistant populations
op<-par(mar=c(5.1,4.1,3.1,2.1))
boxplot(as.numeric(oldprod$rslt_03),
        boxwex=0.4,las=1,ylim=c(0,100),col="transparent",
        # main=paste("Germination des résistants\n(n=",
        #            length(oldprod$rslt_03[oldprod$rslt_03!=0]),
        #            "/",
        #            length(oldprod$rslt_03),")",sep=""),
        ylab="",frame=FALSE,cex.axis=1.3,cex.lab=1.5)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03)),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(oldprod$rslt_03),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)
mtext(text="D",side=3,cex=2,at=0,font=2,las=0,adj=2,line=1)
mtext(text="% resistance",side=2,cex=1,font=2,line=4)
par(op)
#export to a pdf file 10 x 7 inches
dev.off()

min(as.numeric(oldprod$rslt_03))
max(as.numeric(oldprod$rslt_03))
mean(as.numeric(oldprod$rslt_03))
median(as.numeric(oldprod$rslt_03))
length(as.numeric(oldprod$rslt_03))


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

pdf(file="output/map_azoxyTot.pdf",width=7,height=10)
#actual plotting
nf<-layout(matrix(c(1,1,1,1,
                    1,1,1,1,
                    1,1,1,1,
                    2,2,2,2,
                    3,3,3,4,
                    3,3,3,4),6,4,byrow=TRUE))
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
legend(160000,7150000,title="classe de fréquence\nde résistance",
       legend=nomCat,cex=1,pt.cex=1.8,
       y.intersp=1,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(360000,7150000,legend=c("2019","2020"),cex=1,pt.cex=1.6,
       y.intersp=1,x.intersp=0.8,title="année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")
text(123000,7150000,labels="A",cex=3,font=2)
par(op)

#distribution of the % of germination at the DD
op<-par(mar=c(5.1,6.1,0,2.1))
plot(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     bg=as.character(oldprod$catgerm[order(as.numeric(oldprod$rslt_03))]),
     pch=oldprod$year[order(as.numeric(oldprod$rslt_03))],
     cex=1.5,las=1,ylim=c(0,100),cex.axis=1.3,cex.lab=1.5,
     ylab="",xlab="population",
     main="",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03)),
       col="red",lty=2,lwd=3)
mtext(text="B",side=3,cex=2,at=0,font=2,las=0,adj=3.55,line=1)
mtext(text="% resistance",side=2,cex=1,font=2,line=4)
par(op)

#histogram of the distribution of the % of germination
op<-par(mar=c(5.1,6.1,3.1,5.1))
hist(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     breaks=20,bty="l",freq=FALSE,las=1,main="",xlim=c(0,100),cex.lab=1.5,
     col=brewer.pal(11,"RdYlGn")[c(8,6,5,5,4,4,3,3,2,2,rep(1,10))],
     cex.axis=1.3,
     xlab="classe de fréquence de résistance",ylab="")
box(bty="l")
mtext(text="C",side=3,cex=2,at=0,font=2,las=0,adj=3.45,line=1)
mtext(text="fréquence",side=2,cex=1,font=2,line=4)
par(op)

#boxplot for the resistant populations
op<-par(mar=c(5.1,4.1,3.1,2.1))
boxplot(as.numeric(oldprod$rslt_03),
        boxwex=0.4,las=1,ylim=c(0,100),col="transparent",
        # main=paste("Germination des résistants\n(n=",
        #            length(oldprod$rslt_03[oldprod$rslt_03!=0]),
        #            "/",
        #            length(oldprod$rslt_03),")",sep=""),
        ylab="",frame=FALSE,cex.axis=1.3,cex.lab=1.5)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03)),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(oldprod$rslt_03),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)
mtext(text="D",side=3,cex=2,at=0,font=2,las=0,adj=2,line=1)
mtext(text="% resistance",side=2,cex=1,font=2,line=4)
par(op)
#export to a pdf file 10 x 7 inches
dev.off()

min(as.numeric(oldprod$rslt_03))
max(as.numeric(oldprod$rslt_03))
mean(as.numeric(oldprod$rslt_03))
median(as.numeric(oldprod$rslt_03))
length(as.numeric(oldprod$rslt_03))


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

pdf(file="output/map_azoxyCib.pdf",width=7,height=10)
nf<-layout(matrix(c(1,1,1,1,
                    1,1,1,1,
                    1,1,1,1,
                    2,2,2,2,
                    3,3,3,4,
                    3,3,3,4),6,4,byrow=TRUE))
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
legend(160000,7150000,title="classe de fréquence\nde résistance",
       legend=nomCat,cex=1,pt.cex=1.8,
       y.intersp=1,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(360000,7150000,legend=c("2019","2020"),cex=1,pt.cex=1.6,
       y.intersp=1,x.intersp=0.8,title="année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")
text(123000,7150000,labels="A",cex=3,font=2)
par(op)

#distribution of the % of germination at the DD
op<-par(mar=c(5.1,6.1,0,2.1))
plot(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     bg=as.character(oldprod$catgerm[order(as.numeric(oldprod$rslt_03))]),
     pch=oldprod$year[order(as.numeric(oldprod$rslt_03))],
     cex=1.5,las=1,ylim=c(0,100),cex.axis=1.3,cex.lab=1.5,
     ylab="",xlab="population",
     main="",frame=FALSE)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03)),
       col="red",lty=2,lwd=3)
mtext(text="B",side=3,cex=2,at=0,font=2,las=0,adj=3.55,line=1)
mtext(text="% resistance",side=2,cex=1,font=2,line=4)
par(op)

#histogram of the distribution of the % of germination
op<-par(mar=c(5.1,6.1,3.1,5.1))
hist(as.numeric(oldprod$rslt_03[order(as.numeric(oldprod$rslt_03))]),
     breaks=20,bty="l",freq=FALSE,las=1,main="",xlim=c(0,100),cex.lab=1.5,
     col=brewer.pal(11,"RdYlGn")[c(8,6,5,5,4,4,3,3,2,2,rep(1,10))][7:20],
     cex.axis=1.3,
     xlab="classe de fréquence de résistance",ylab="")
box(bty="l")
mtext(text="C",side=3,cex=2,at=0,font=2,las=0,adj=3.45,line=1)
mtext(text="fréquence",side=2,cex=1,font=2,line=4)
par(op)

#boxplot for the resistant populations
op<-par(mar=c(5.1,4.1,3.1,2.1))
boxplot(as.numeric(oldprod$rslt_03),
        boxwex=0.4,las=1,ylim=c(0,100),col="transparent",
        # main=paste("Germination des résistants\n(n=",
        #            length(oldprod$rslt_03[oldprod$rslt_03!=0]),
        #            "/",
        #            length(oldprod$rslt_03),")",sep=""),
        ylab="",frame=FALSE,cex.axis=1.3,cex.lab=1.5)
box(bty="l")
abline(h=mean(as.numeric(oldprod$rslt_03)),
       col="red",lty=2,lwd=3)
stripchart(as.numeric(oldprod$rslt_03),
           cex=1,pch=19,
           col=adjustcolor("grey",alpha=0.5),vertical=TRUE,
           method="jitter",jitter=0.1,add=TRUE)
mtext(text="D",side=3,cex=2,at=0,font=2,las=0,adj=2,line=1)
mtext(text="% resistance",side=2,cex=1,font=2,line=4)
par(op)
#export to a pdf file 10 x 7 inches
dev.off()

min(as.numeric(oldprod$rslt_03))
max(as.numeric(oldprod$rslt_03))
mean(as.numeric(oldprod$rslt_03))
median(as.numeric(oldprod$rslt_03))
length(as.numeric(oldprod$rslt_03))


##############################################################################/
#END
##############################################################################/