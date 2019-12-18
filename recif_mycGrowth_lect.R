##############################################################################/
##############################################################################/
#Comparison of the measure after 14 or 20 days
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)
library(tidyr)
library(ade4)
library(factoextra)
library(RColorBrewer)

#loading the data
#datamyc<-read.table("data/cerco_mars19.txt",header=TRUE,sep=";")
datamyc2<-read.table("data/results_ec50_panels_1_2.txt",header=TRUE,sep=";")
datamyc2<-read.table("data/data_1_2_3_4.txt",header=TRUE,sep="\t")


##############################################################################/
#Regression analysis of mycelial growth experiment 14 days####
##############################################################################/

datamyc<-datamyc2[datamyc2$tps_expo=="14",]

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),read_time=factor(),
                    ED50=character(),ED95=character(),ED99=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_14D.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  
  REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    read_time=factor(),ED50=character(),
                    ED95=character(),ED99=character())
  
  for (i in 1:dim(table(data_subSA$ech_id))[1]) {
    tempdat<-data_subSA[data_subSA$ech_id==names(table(data_subSA$ech_id))[i],]
    if(tempdat[tempdat$dose==max(tempdat$dose),"rslt_03"]>50) {
      tempx<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=tempdat$ech_id[1],
                        "read_time"=data_subSA$tps_expo[1],
                        "ED50"=paste(">",max(tempdat$dose),sep=""),
                        "ED95"=paste(">",max(tempdat$dose),sep=""),
                        "ED99"=paste(">",max(tempdat$dose),sep=""))
    } else {
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   fct=LN.3())
      plot(temp.m1,ylim=c(0,110),xlim=c(0,100),
           main=paste(SAlist[j],as.character(tempdat$ech_id[1])))
      temp<-ED(temp.m1,c(50,5,1),type="absolute")
      tempx<-data.frame("Subs_Act"=SAlist[j],
                        "sample_ID"=as.character(tempdat$ech_id[1]),
                        "read_time"=data_subSA$tps_expo[1],
                        "ED50"=as.character(temp[1]),
                        "ED95"=as.character(temp[2]),
                        "ED99"=as.character(temp[3]))
    }
    
    REZSA<-rbind(REZSA,tempx)
  }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()
    

##############################################################################/
#Regression analysis of mycelial growth experiment 20 days####
##############################################################################/

datamyc<-datamyc2[datamyc2$tps_expo=="20",]

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-CompRez

#we make a subselection of the data according to the SA
pdf(file="output/plot_14D.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  
  REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    read_time=factor(),ED50=character(),
                    ED95=character(),ED99=character())
  
  for (i in 1:dim(table(data_subSA$ech_id))[1]) {
    tempdat<-data_subSA[data_subSA$ech_id==names(table(data_subSA$ech_id))[i],]
    if(tempdat[tempdat$dose==max(tempdat$dose),"rslt_03"]>50) {
      tempx<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=tempdat$ech_id[1],
                        "read_time"=data_subSA$tps_expo[1],
                        "ED50"=paste(">",max(tempdat$dose),sep=""),
                        "ED95"=paste(">",max(tempdat$dose),sep=""),
                        "ED99"=paste(">",max(tempdat$dose),sep=""))
    } else {
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   fct=LN.3())
      plot(temp.m1,ylim=c(0,110),xlim=c(0,100),
           main=paste(SAlist[j],as.character(tempdat$ech_id[1])))
      temp<-ED(temp.m1,c(50,5,1),type="absolute")
      tempx<-data.frame("Subs_Act"=SAlist[j],
                        "sample_ID"=as.character(tempdat$ech_id[1]),
                        "read_time"=data_subSA$tps_expo[1],
                        "ED50"=as.character(temp[1]),
                        "ED95"=as.character(temp[2]),
                        "ED99"=as.character(temp[3]))
    }
    
    REZSA<-rbind(REZSA,tempx)
  }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()

#exporting the result as a text file
CompRez<-CompRez[order(CompRez$Subs_Act,CompRez$sample_ID),]
write.table(CompRez, file="output/results_cerco.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#barplot to compare the ED50 of the different samples####
##############################################################################/

#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
CompRez$ED50<-as.character(CompRez$ED50)
CompRez[CompRez$ED50==">10","ED50"]<-12
CompRez[CompRez$ED50==">20","ED50"]<-22
CompRez[CompRez$ED50==">50","ED50"]<-52
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))
op<-par(mfrow=c(2,1))
par(mar=c(2,3,7,1))
barplot(as.numeric(as.character(CompRez[CompRez$read_time=="14",]$ED50)),
        ylim=c(0,52),col=CompRez[CompRez$read_time=="14",]$Subs_Act,
        main="reading at 14 Days",las=2)
par(mar=c(7,3,2,1))
barplot(as.numeric(as.character(CompRez[CompRez$read_time=="20",]$ED50)),
        ylim=c(0,52),col=CompRez[CompRez$read_time=="20",]$Subs_Act,
        names.arg=CompRez[CompRez$read_time=="20",]$sample_ID,las=2,
        main="reading at 20 Days")
par(op)
#export to pdf 32 x 8 inches

#histogramme by samples: reading at 14 days
samplelist<-as.character(names(table(CompRez$sample_ID)))
temp14<-CompRez[CompRez$read_time=="14",]

op<-par(mfrow=c(5,8))
for (i in (1:length(samplelist))) {
  barplot(temp14[temp14$sample_ID==samplelist[i],]$ED50,
          col=c(1,2,3,4,5),las=1,main=samplelist[i],
          ylim=c(0,52))
}
par(op)
#export to pdf 15 x 15

#histogramme by samples: reading at 20 days
samplelist<-as.character(names(table(CompRez$sample_ID)))
temp20<-CompRez[CompRez$read_time=="20",]

op<-par(mfrow=c(5,8))
for (i in (1:length(samplelist))) {
  barplot(temp20[temp20$sample_ID==samplelist[i],]$ED50,
          col=c(1,2,3,4,5),las=1,main=samplelist[i],
          ylim=c(0,52))
}
par(op)
#export to pdf 10 x 10


##############################################################################/
#correlation between reading after 14 and 20 days####
##############################################################################/

#preparing the dataset
outlie<-spread(CompRez[,c(1:4)],read_time,ED50)
#correlation between the 14 days and 20 days measures
cor(outlie$`14`,outlie$`20`)

plot(outlie$`14`,outlie$`20`,col=outlie$Subs_Act,pch=19,cex=1,las=1,
     xlab="IC50 after 14 days",ylab="IC50 after 20 days",
     main="Correlation between reading after 14 or 20 days")
text(outlie[outlie$`14`>20,3],outlie[outlie$`14`>20,4],
     labels=outlie[outlie$`14`>20,"sample_ID"],pos=1,xpd=NA)
abline(0,1,lty=2,lwd=3)
legend(2,45,names(table(CompRez$Subs_Act)),col=c(1,2,3,4,5),
       pch=19,bty="n")

#export to pdf 7 x 7

#same figure with log transformation of IC50
cor(log(outlie$`14`),log(outlie$`20`))

plot(log(outlie$`14`),log(outlie$`20`),col=outlie$Subs_Act,pch=19,cex=1,las=1,
     xlab="IC50 after 14 days",ylab="IC50 after 20 days",
     main="Correlation between reading after 14 or 20 days")
abline(0,1,lty=2,lwd=3)
legend(-1,-2,names(table(CompRez$Subs_Act)),col=c(1,2,3,4,5),
       pch=19,bty="n")

#export to pdf 7 x 7


##############################################################################/
#correlation between ED50 estimated for different active substances####
##############################################################################/

temp<-CompRez[,c(1:4)]
temp<-spread(temp,Subs_Act,ED50)

#a function to compute the absolute correlation between pairs of variables
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(temp[temp$read_time=="14",c(3:7)],las=1,main="14 Days",
      lower.panel=panel.smooth, upper.panel=panel.cor)

pairs(log(temp[temp$read_time=="14",c(3:7)]),las=1,main="Log 14 Days",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 8 x 8 inches

pairs(temp[temp$read_time=="20",c(3:7)],las=1,main="20 Days",
      lower.panel=panel.smooth, upper.panel=panel.cor)

pairs(log(temp[temp$read_time=="20",c(3:7)]),las=1,main="Log 20 Days",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 8 x 8 inches


##############################################################################/
#Analyzing the multisensitivity profil of the strains####
##############################################################################/

#Clusterization based on 14 days reading
temp14<-CompRez[CompRez$read_time=="14",c(1:4)]
temp14<-spread(temp14,Subs_Act,ED50)
row.names(temp14)<-temp14$sample_ID

#PCA for the reading after 14 days
truc14<-dudi.pca(temp14[,-c(1,2)],
                 scannf=FALSE,nf=3)
scatter(truc14)
#determining the optimal number of clusters
fviz_nbclust(temp14[,c(3:7)],kmeans,method="gap_stat")
clust14<-kmeans(temp14[,c(3:7)],5)
fviz_cluster(clust14,data=temp14[,c(3:7)])
plot(truc14$li[,c(1,2)],col=brewer.pal(5,"Dark2")[clust14$cluster],
     pch=19,cex=2)
hclu14<-hclust(dist(scale(temp14[,c(3:7)]),
                    method="euclidean"),
               method="ward.D2")
plot(hclu14)
fviz_dend(hclu14,k=5,cex=0.5,rect=TRUE,
          k_colors=brewer.pal(5,"Dark2"))

#Clusterization based on 20 days reading
temp20<-CompRez[CompRez$read_time=="20",c(1:4)]
temp20<-spread(temp20,Subs_Act,ED50)
row.names(temp20)<-temp20$sample_ID

#PCA for the reading after 20 days
truc20<-dudi.pca(temp20[,-c(1,2)],
                 scannf=FALSE,nf=3)
scatter(truc20)
#determining the optimal number of clusters
fviz_nbclust(temp20[,c(3:7)],kmeans,method="gap_stat")
clust20<-kmeans(temp[temp$read_time=="20",c(3:7)],5)
fviz_cluster(clust20,data=temp20[,c(3:7)])
plot(truc20$li[,c(1,2)],col=brewer.pal(5,"Dark2")[clust20$cluster],
     pch=19,cex=2)
hclu20<-hclust(dist(scale(temp20[,c(3:7)]),
                    method="euclidean"),
               method="ward.D2")
plot(hclu20)
fviz_dend(hclu20,k=5,cex=0.5,rect=TRUE,
          k_colors=brewer.pal(5,"Dark2"))


##############################################################################/
#plot of the distribution of IC50 for each active substance####
##############################################################################/

#readings after 14 days
op<-par(mfrow=c(2,3))
plot(temp14[order(c(temp14$DIFENOCONAZOLE)),"DIFENOCONAZOLE"],
     main="14D DIFENOCONAZOLE IC50",col=1,pch=19,cex=2,las=1,
     ylab="IC50")
plot(temp14[order(c(temp14$FLUTRIAFOL)),"FLUTRIAFOL"],
     main="14D FLUTRIAFOL IC50",col=2,pch=19,cex=2,las=1,
     ylab="IC50")
plot(temp14[order(c(temp14$`PROTHIOCONAZOLE-DESTHIO`)),
            "PROTHIOCONAZOLE-DESTHIO"],
     main="14D PROTHIOCONAZOLE-DESTHIO IC50",col=3,pch=19,
     cex=2,las=1,ylim=c(0,12),ylab="IC50")
plot(temp14[order(c(temp14$TEBUCONAZOLE)),"TEBUCONAZOLE"],
     main="14D TEBUCONAZOLE IC50",col=4,pch=19,cex=2,las=1,
     ylab="IC50")
plot(temp14[order(c(temp14$TETRACONAZOLE)),"TETRACONAZOLE"],
     main="14D TETRACONAZOLE IC50",col=5,pch=19,cex=2,las=1,
     ylab="IC50")
par(op)
#export to pdf 9 x 7 inches

#readings after 20 days
op<-par(mfrow=c(2,3))
plot(temp20[order(c(temp20$DIFENOCONAZOLE)),"DIFENOCONAZOLE"],
     main="20D DIFENOCONAZOLE IC50",col=1,pch=19,cex=2,las=1,
     ylab="IC50")
plot(temp20[order(c(temp20$FLUTRIAFOL)),"FLUTRIAFOL"],
     main="20D FLUTRIAFOL IC50",col=2,pch=19,cex=2,las=1,
     ylab="IC50")
plot(temp20[order(c(temp20$`PROTHIOCONAZOLE-DESTHIO`)),
            "PROTHIOCONAZOLE-DESTHIO"],
     main="20D PROTHIOCONAZOLE-DESTHIO IC50",col=3,pch=19,
     cex=2,las=1,ylim=c(0,12),ylab="IC50")
plot(temp20[order(c(temp20$TEBUCONAZOLE)),"TEBUCONAZOLE"],
     main="20D TEBUCONAZOLE IC50",col=4,pch=19,cex=2,las=1,
     ylab="IC50")
plot(temp20[order(c(temp20$TETRACONAZOLE)),"TETRACONAZOLE"],
     main="20D TETRACONAZOLE IC50",col=5,pch=19,cex=2,las=1,
     ylab="IC50")
par(op)
#export to pdf 9 x 7 inches


##############################################################################/
#END
##############################################################################/