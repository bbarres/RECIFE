##############################################################################/
##############################################################################/
#Code for production of figures with all SA, scoring after 20 days
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
#datamyc2<-read.table("data/20200309_data_temp.txt",header=TRUE,sep=";")
datamyc2<-read.table("data/20200415_data_12SA-2.txt",header=TRUE,
                     sep=";",stringAsFactors=TRUE)


##############################################################################/
#Regression analysis of mycelial growth experiment scoring 20 or 21 days####
##############################################################################/

datamyc<-datamyc2[datamyc2$tps_expo>16,]

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),read_time=factor(),
                    ED50=character(),ED95=character(),ED99=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_ASA.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  
  REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    read_time=factor(),ED50=character(),
                    ED95=character(),ED99=character())
  
  for (i in 1:dim(table(data_subSA$ech_id))[1]) {
    tempdat<-data_subSA[data_subSA$ech_id==names(table(data_subSA$ech_id))[i],]
    if(tempdat[tempdat$dose==max(tempdat$dose),"rslt_03"]>30) {
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
write.table(CompRez, file="output/ASA_results_cerco.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#barplot to compare the ED50 of the different samples####
##############################################################################/

#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
cooloor<- brewer.pal(12,"Set3")
CompRez$ED50<-as.character(CompRez$ED50)
CompRez[CompRez$ED50==">10","ED50"]<-12
CompRez[CompRez$ED50==">20","ED50"]<-22
CompRez[CompRez$ED50==">50","ED50"]<-52
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))

pdf(file="output/histo_AllInd_ASA.pdf",width=60,height=8)
op<-par(mfrow=c(1,1))
par(mar=c(8,3,3,0.5))
barplot(as.numeric(as.character(CompRez$ED50)),
        ylim=c(0,52),col=cooloor[as.numeric(CompRez$Subs_Act)],
        names.arg=CompRez$sample_ID,las=2,
        main="Comparison of the different samples by SA")
abline(h=12,lty=2)
abline(h=22,lty=2)
abline(h=52,lty=2)
legend(300,47,levels(CompRez$Subs_Act),fill=cooloor,bty="n")
par(op)
dev.off()

#histogramme by samples
samplelist<-as.character(names(table(CompRez$sample_ID)))
pdf(file="output/histo_byInd_ASA.pdf",width=9,height=20)
op<-par(mfrow=c(9,5))
for (i in (1:length(samplelist))) {
  temp<-merge(as.data.frame(levels(CompRez$Subs_Act)),
              CompRez[CompRez$sample_ID==samplelist[i],],
              by.x=1,by.y=1,
              all.x=TRUE)
  barplot(temp$ED50,col=cooloor,las=1,main=samplelist[i],
            ylim=c(0,52))
}
par(op)
dev.off()
#export to pdf 12 x 16


##############################################################################/
#correlation between ED50 estimated for different active substances####
##############################################################################/

temp<-CompRez[,c(1,2,4)]
temp<-spread(temp,Subs_Act,ED50)

#a function to compute the absolute correlation between pairs of variables
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r, 2), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(temp[,c(2:13)],las=1,main="Correlation between ActSubst",
      lower.panel=panel.smooth, upper.panel=panel.cor)

pairs(log(temp[,c(2:13)]),las=1,main="Correlation between log(ActSubst)",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 11 x 11 inches

#just for the difenoconazole and tetraconazole
pairs(log(temp[,c(3,12)]),las=1,main="Correlation between log(ActSubst)",
      lower.panel=panel.smooth, upper.panel=panel.cor)


##############################################################################/
#Analyzing the multisensitivity profil of the strains####
##############################################################################/

#Clusterization based on scoring of 10 SA
row.names(temp)<-temp$sample_ID

#PCA for the scoring on 10 SA
truc<-dudi.pca(temp[,-c(1)],
               scannf=FALSE,nf=3)
scatter(truc)
#determining the optimal number of clusters
fviz_nbclust(temp[,c(2:13)],kmeans,method="gap_stat")
clust<-kmeans(temp[,c(2:13)],5)
fviz_cluster(clust,data=temp[,c(2:13)])
plot(truc$li[,c(1,2)],col=brewer.pal(5,"Dark2")[clust$cluster],
     pch=19,cex=2)

#we remove FENTINE HYDROXYDE and TOLNAFTATE because it is of no 
#interest here as well as individuals 39 to 44 that have too many 
#missing values
hclu<-hclust(dist(scale(temp[-c(39:44),c(2:4,6:12)]),
                  method="euclidean"),
               method="ward.D2")
plot(hclu)
fviz_dend(hclu,k=5,cex=0.5,rect=TRUE,
          k_colors=brewer.pal(5,"Dark2"))
#export to pdf 9 x 6 inches


##############################################################################/
#plot of the distribution of IC50 for each active substance####
##############################################################################/

#preparing the dataset
temp<-CompRez[,c(1,2,4)]
temp<-spread(temp,Subs_Act,ED50)

#distribution of the IC50 by Active Substance
op<-par(mfrow=c(3,4))
plot(temp[order(c(temp$CYPROCONAZOLE)),"CYPROCONAZOLE"],
     main="CYPROCONAZOLE IC50",bg=cooloor[1],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$DIFENOCONAZOLE)),"DIFENOCONAZOLE"],
     main="DIFENOCONAZOLE IC50",bg=cooloor[2],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$EPOXICONAZOLE)),"EPOXICONAZOLE"],
     main="EPOXICONAZOLE IC50",bg=cooloor[3],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$`FENTINE HYDROXYDE`)),"FENTINE HYDROXYDE"],
     main="FENTINE HYDROXYDE IC50",bg=cooloor[4],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$FLUTRIAFOL)),"FLUTRIAFOL"],
     main="FLUTRIAFOL IC50",bg=cooloor[5],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$MEFENTRIFLUCONAZOLE)),"MEFENTRIFLUCONAZOLE"],
     main="MEFENTRIFLUCONAZOLE IC50",bg=cooloor[6],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$METCONAZOLE)),"METCONAZOLE"],
     main="METCONAZOLE IC50",bg=cooloor[7],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$PROCHLORAZE)),"PROCHLORAZE"],
     main="PROCHLORAZE IC50",bg=cooloor[8],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$`PROTHIOCONAZOLE-DESTHIO`)),"PROTHIOCONAZOLE-DESTHIO"],
     main="PROTHIOCONAZOLE-DESTHIO IC50",bg=cooloor[9],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$TEBUCONAZOLE)),"TEBUCONAZOLE"],
     main="TEBUCONAZOLE IC50",bg=cooloor[10],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$TETRACONAZOLE)),"TETRACONAZOLE"],
     main="TETRACONAZOLE IC50",bg=cooloor[11],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$TOLNAFTATE)),"TOLNAFTATE"],
     main="TOLNAFTATE IC50",bg=cooloor[12],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
par(op)
#export to pdf 14 x 10 inches


##############################################################################/
#END
##############################################################################/