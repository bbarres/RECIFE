##############################################################################/
##############################################################################/
#Code for producing the figure of the Note commune resistance
##############################################################################/
##############################################################################/

#loading the data set and the necessary library
source("recif_load.R")


##############################################################################/
#Preparing the two data sets####
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
op<-par(mfrow=c(2,2),mar=c(0,0,0,0))
plot(DEP_SHP.1,main="",border="grey70")
title(main="Fentine Hydroxyde",font=2,line=-2)
plot(REG_SHP.1,lwd=2,add=TRUE)
points(
  x = as.numeric(oldprod$gps_long),
  y = as.numeric(oldprod$gps_lat),
  bg = as.character(oldprod$catgerm),  #colors of the points
  pch = oldprod$year,                  #plotting character
  cex = 1.7                            #size of the points
)
legend(100000,7250000,title="Germination\nclasses",
       legend=nomCat,cex=0.9,pt.cex=1.8,
       y.intersp=0.6,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(300000,7250000,legend=c("2019","2020"),cex=0.9,pt.cex=1.6,
       y.intersp=0.6,x.intersp=0.8,title="Année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")


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
plot(DEP_SHP.1,main="",border="grey70")
title(main="Azoxystrobine",font=2,line=-2)
plot(REG_SHP.1,lwd=2,add=TRUE)
points(
  x = as.numeric(oldprod$gps_long),
  y = as.numeric(oldprod$gps_lat),
  bg = as.character(oldprod$catgerm),  #colors of the points
  pch = oldprod$year,                  #plotting character
  cex = 1.7                            #size of the points
)
legend(100000,7250000,title="Germination\nclasses",
       legend=nomCat,cex=0.9,pt.cex=1.8,
       y.intersp=0.6,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=as.character(levels(oldprod$catgerm)),
       bg="transparent",bty="n")
legend(300000,7250000,legend=c("2019","2020"),cex=0.9,pt.cex=1.6,
       y.intersp=0.6,x.intersp=0.8,title="Année",title.adj=0.3,
       pch=c(21,22),col=c("black"),bg="transparent",bty="n")


##############################################################################/
#difenoconazole map 3 resistance categories####
##############################################################################/

#limiting the data set to mefentrifluconazole
newprod<-newSA[newSA$pest_sa_id=="DIFENOCONAZOLE",]

#actual plotting
plot(DEP_SHP.1,main="",border="grey70")
title(main="Difénoconazole",font=2,line=-2)
plot(REG_SHP.1,lwd=2,add=TRUE)
draw.pie(
  x = as.numeric(newprod$gps_long),
  y = as.numeric(newprod$gps_lat),
  z = cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100)
  ),
  col=colovec,           #colors of the pie
  lty=1,                 #line type of the pie
  border="transparent",  #color of the border of the pie
  lwd=0.01,              #control the width of the border
  radius=11000,          #(sqrt(as.numeric(as.character(data2map$Total))) * 
                         #16000), this number control the radius of the pies
  labels=NA,
  scale=FALSE            #should the radius be scaled according to sample size
)
legend(100000,7200000,title="Classes de facteur\nde résistance",
       legend=c("FR<30","30<FR<100","100<FR"),
       cex=0.9,pt.cex=1.8,
       y.intersp=0.6,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=colovec,
       bg="transparent",bty="n")


##############################################################################/
#tétraconazole map 3 resistance categories####
##############################################################################/

#limiting the data set to mefentrifluconazole
newprod<-newSA[newSA$pest_sa_id=="TETRACONAZOLE",]

#actual plotting
plot(DEP_SHP.1,main="",border="grey70")
title(main="Tétraconazole",font=2,line=-2)
plot(REG_SHP.1,lwd=2,add=TRUE)
draw.pie(
  x = as.numeric(newprod$gps_long),
  y = as.numeric(newprod$gps_lat),
  z = cbind(
    as.numeric(newprod$FR.30),
    as.numeric(newprod$FR30.100),
    as.numeric(newprod$FR.100)
  ),
  col=colovec,           #colors of the pie
  lty=1,                 #line type of the pie
  border="transparent",  #color of the border of the pie
  lwd=0.01,              #control the width of the border
  radius=11000,          #(sqrt(as.numeric(as.character(data2map$Total))) * 
                         #16000), this number control the radius of the pies
  labels=NA,
  scale=FALSE            #should the radius be scaled according to sample size
)
legend(100000,7200000,title="Classes de facteur\nde résistance",
       legend=c("FR<30","30<FR<100","100<FR"),
       cex=0.9,pt.cex=1.8,
       y.intersp=0.6,x.intersp=0.8,
       pch=15,title.adj=0.3,
       col=colovec,
       bg="transparent",bty="n")

par(op)

#export to a pdf file 10 x 10 inches


##############################################################################/
#loading the Myzus data sets####
##############################################################################/

#loading the data exported from PROSPER "tableau de bord"
rezMoni<-read.table("data/myzusNotCo.txt",header=TRUE,sep="\t",
                    colClasses="character")
#turning the resistance status factor into two different columns
rezMoni$rslt_RS<-factor(rezMoni$rslt_RS,levels=c("R","S"))
rezMoni$Resistant<-as.numeric(rezMoni$rslt_RS==levels(rezMoni$rslt_RS)[1])
rezMoni$Sensitive<-as.numeric(rezMoni$rslt_RS==levels(rezMoni$rslt_RS)[2])
rezMoni$Total<-rezMoni$Resistant+rezMoni$Sensitive
#grouping and counting the resistant and sensitive samples
rezMoni$SA<-as.factor(rezMoni$SA)
dataCamem<-rezMoni %>% 
  group_by(SA,dptmt,pest,host) %>% 
  summarise(Resist=sum(Resistant),Sensi=sum(Sensitive),Tot=sum(Total))

#splitting the continuous percentage of germination in categories
dataCamem$catgerm<-cut(dataCamem$Resist,
                       breaks=c(0,0.001,1,5,10,15,20),
                       include.lowest=TRUE)
nomCat<-c("[0]","[1]","]1-5]","]5-10]","]10-15]","]15-20]")
#defining the colors of the points
levels(dataCamem$catgerm)<-c(brewer.pal(4,"Paired")[4],
                             brewer.pal(11,"RdYlGn")[5:1])


##############################################################################/
#Maps of the results of the monitoring####
##############################################################################/

#mapping the results for each "programme" of the monitoring
#first

png(file=paste("output/Myzusmap",".png",sep=""),
    width=9,height=9,units="in",res=300)

op<-par(mfrow=c(2,2),mar=c(0,0,0,0))

temp<-dataCamem[dataCamem$SA=="1014",]
plot(DEP_SHP,border="grey70")
title(main="kdr",font.main=4,line=-1.4,cex.main=2)
plot(DEP_SHP[DEP_SHP$INSEE_DEP %in% temp$dptmt,],
     col=as.character(temp$catgerm),
     border="grey70",add=TRUE)
plot(REG_SHP,lwd=2,add=TRUE)
legend(70000,7150000,title="Classes de nombre de\nprélèvement(s) résistant(s)",
       legend=nomCat,cex=0.8,pt.cex=1.6,
       y.intersp=0.9,x.intersp=0.8,
       pch=15,title.adj=0.6,
       col=as.character(levels(temp$catgerm)),
       bg="transparent",bty="n")

temp<-dataCamem[dataCamem$SA=="918",]
plot(DEP_SHP,border="grey70")
title(main="skdr",font.main=4,line=-1.4,cex.main=2)
plot(DEP_SHP[DEP_SHP$INSEE_DEP %in% temp$dptmt,],
     col=as.character(temp$catgerm),
     border="grey70",add=TRUE)
plot(REG_SHP,lwd=2,add=TRUE)
legend(70000,7150000,title="Classes de nombre de\nprélèvement(s) résistant(s)",
       legend=nomCat,cex=0.8,pt.cex=1.6,
       y.intersp=0.9,x.intersp=0.8,
       pch=15,title.adj=0.6,
       col=as.character(levels(temp$catgerm)),
       bg="transparent",bty="n")

temp<-dataCamem[dataCamem$SA=="222",]
plot(DEP_SHP,border="grey70")
title(main="MACE",font.main=4,line=-1.4,cex.main=2)
plot(DEP_SHP[DEP_SHP$INSEE_DEP %in% temp$dptmt,],
     col=as.character(temp$catgerm),
     border="grey70",add=TRUE)
plot(REG_SHP,lwd=2,add=TRUE)
legend(70000,7150000,title="Classes de nombre de\nprélèvement(s) résistant(s)",
       legend=nomCat,cex=0.8,pt.cex=1.6,
       y.intersp=0.9,x.intersp=0.8,
       pch=15,title.adj=0.6,
       col=as.character(levels(temp$catgerm)),
       bg="transparent",bty="n")

par(op)
dev.off()


##############################################################################/
#END
##############################################################################/