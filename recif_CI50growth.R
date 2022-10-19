##############################################################################/
##############################################################################/
#Code for production of figures with all SA, scoring after 20 days
##############################################################################/
##############################################################################/

##loading the dataset and the necessary library
source("recif_load.R")


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
    #here we limit the dose response analyses to individuals who reach a 
    #certain threshold in inhibition
    if(mean(tempdat[tempdat$dose==max(tempdat$dose),"rslt_03"])>40) {
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
CompRez$Subs_Act<-as.factor(CompRez$Subs_Act)

pdf(file="output/histo_AllInd_ASA.pdf",width=60,height=8)
op<-par(mfrow=c(1,1))
par(mar=c(8,3,3,0.5))
barplot(as.numeric(as.character(CompRez$ED50)),
        ylim=c(0,52),col=cooloor[as.numeric(as.factor(CompRez$Subs_Act))],
        names.arg=CompRez$sample_ID,las=2,
        main="Comparison of the different samples by SA")
abline(h=12,lty=2)
abline(h=22,lty=2)
abline(h=52,lty=2)
legend(300,47,levels(CompRez$Subs_Act),fill=cooloor,bty="n")
par(op)
dev.off()

#histogram by samples
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
#Analyzing the multisensitivity profile of the strains####
##############################################################################/

temp<-CompRez[,c(1:2,4)]
temp<-spread(temp,Subs_Act,ED50)
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
#END
##############################################################################/