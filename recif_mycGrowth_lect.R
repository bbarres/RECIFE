##############################################################################/
##############################################################################/
#Comparison of the measure after 14 or 20 days
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)

#loading the data
#datamyc<-read.table("data/cerco_mars19.txt",header=TRUE,sep=";")
datamyc2<-read.table("data/panel1_2lectures.txt",header=TRUE,sep=";")


##############################################################################/
#Regression analysis of mycelial growth experiment 14 days####
##############################################################################/

datamyc<-datamyc2[datamyc2$tps_expo=="14",]

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),read_time=factor(),
                    ED50=character(),ED95=character(),ED99=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_14h.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                  & data_subSA$rslt_03>50,
                                  "ech_id"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                           read_time=factor(),ED50=character(),
                           ED95=character(),ED99=character()),
         REZSA<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "read_time"=data_subSA$tps_expo[1],
                           "ED50"=paste(">",max(data_subSA$dose),sep=""),
                           "ED95"=paste(">",max(data_subSA$dose),sep=""),
                           "ED99"=paste(">",max(data_subSA$dose),sep=""))
  )
  #we limit the dataset to the sample that reach somehow a IC of 50%
  if(dim(data_subSA[!(data_subSA$ech_id %in% SA_rez),])[1]!=0) {
    SA.dat<-data_subSA[!(data_subSA$ech_id %in% SA_rez),]
    SA.dat<-drop.levels(SA.dat)
    for (i in 1:dim(table(SA.dat$ech_id))[1]) {
      tempdat<-SA.dat[SA.dat$ech_id==names(table(SA.dat$ech_id))[i],]
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   fct=LL.3())
      plot(temp.m1,ylim=c(0,110),xlim=c(0,50),
           main=paste(SAlist[j],names(table(SA.dat$ech_id))[i]))
      temp<-ED(temp.m1,c(50,5,1),type="absolute")
      tempx<-data.frame("Subs_Act"=SAlist[j],
                        "sample_ID"=names(table(SA.dat$ech_id))[i],
                        "read_time"=data_subSA$tps_expo[1],
                        "ED50"=as.character(temp[1]),
                        "ED95"=as.character(temp[2]),
                        "ED99"=as.character(temp[3]))
      REZSA<-rbind(REZSA,tempx)}} else {
        REZSA<-REZSA
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
pdf(file="output/plot_20h.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                  & data_subSA$rslt_03>50,
                                  "ech_id"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                           read_time=factor(),ED50=character(),
                           ED95=character(),ED99=character()),
         REZSA<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "read_time"=data_subSA$tps_expo[1],
                           "ED50"=paste(">",max(data_subSA$dose),sep=""),
                           "ED95"=paste(">",max(data_subSA$dose),sep=""),
                           "ED99"=paste(">",max(data_subSA$dose),sep=""))
  )
  #we limit the dataset to the sample that reach somehow a IC of 50%
  if(dim(data_subSA[!(data_subSA$ech_id %in% SA_rez),])[1]!=0) {
    SA.dat<-data_subSA[!(data_subSA$ech_id %in% SA_rez),]
    SA.dat<-drop.levels(SA.dat)
    for (i in 1:dim(table(SA.dat$ech_id))[1]) {
      tempdat<-SA.dat[SA.dat$ech_id==names(table(SA.dat$ech_id))[i],]
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   fct=LL.3())
      plot(temp.m1,ylim=c(0,110),xlim=c(0,50),
           main=paste(SAlist[j],names(table(SA.dat$ech_id))[i]))
      temp<-ED(temp.m1,c(50,5,1),type="absolute")
      tempx<-data.frame("Subs_Act"=SAlist[j],
                        "sample_ID"=names(table(SA.dat$ech_id))[i],
                        "read_time"=data_subSA$tps_expo[1],
                        "ED50"=as.character(temp[1]),
                        "ED95"=as.character(temp[2]),
                        "ED99"=as.character(temp[3]))
      REZSA<-rbind(REZSA,tempx)}} else {
        REZSA<-REZSA
      }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()

#exporting the result as a text file
CompRez<-CompRez[order(CompRez$Subs_Act,CompRez$sample_ID),]
write.table(CompRez, file="output/results_cerco.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#Some plots####
##############################################################################/

#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with the dose 
#range used with an absurdly high value
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))
CompRez[is.na(CompRez$ED50),"ED50"]<-30
op<-par(mfrow=c(2,1))
par(mar=c(2,3,7,1))
barplot(as.numeric(as.character(CompRez[CompRez$read_time=="14",]$ED50)),
        ylim=c(0,30),col=CompRez[CompRez$read_time=="14",]$Subs_Act,
        main="reading at 14h")
par(mar=c(7,3,2,1))
barplot(as.numeric(as.character(CompRez[CompRez$read_time=="20",]$ED50)),
        ylim=c(0,30),col=CompRez[CompRez$read_time=="20",]$Subs_Act,
        names.arg=CompRez[CompRez$read_time=="20",]$sample_ID,las=2,
        main="reading at 20h")
par(op)


##############################################################################/
#END
##############################################################################/