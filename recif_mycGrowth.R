###############################################################################
###############################################################################
#R code for analyzing the ouput of mycelial growth bioassay of C. beticola
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)

#loading the data
datamyc<-read.table("data/cerco_mars19.txt",header=TRUE,sep=";")


###############################################################################
#Analysis of the 
###############################################################################

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)

#we make a subselection of the data according to the SA
data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[4],]


#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                & data_subSA$rslt_03>50,
                                "ech_id"])
REZSA<-data.frame("sample_ID"=SA_rez,"ED50"=max(data_subSA$dose))
#we limit the dataset to the sample that reach somehow a IC of 50%
SA.dat<-data_subSA[!(data_subSA$ech_id %in% SA_rez),]
SA.dat<-drop.levels(SA.dat)
for (i in 1: dim(table(SA.dat$ech_id))[1]) {
  temp.m1<-drm(rslt_03~dose,
               data=SA.dat[SA.dat$ech_id==names(table(SA.dat$ech_id))[i],],
               fct=LL.4())
  plot(temp.m1,main=names(table(SA.dat$ech_id))[i])
  temp<-ED(temp.m1,50,type="absolute")
  tempx<-data.frame("sample_ID"=names(table(SA.dat$ech_id))[i],
                    "ED50"=temp[1])
  REZSA<-rbind(REZSA,tempx)
}

#computing the mean ED50 for the reference strains
refval<-mean(REZdod[startsWith(as.character(REZdod$sample_ID),
                               pattern="13",trim=TRUE),"ED50"])

REZdod$ED50[REZdod$ED50>30]<-30
plot(REZdod$ED50[order(REZdod$ED50)]/refval,main="dodine",xlab="Souches ID",
     ylab="FR",las=1)
abline(refval/refval,0,col="green4",lwd=2)
abline((10*refval)/refval,0,col="red",lwd=2)
REZdod$FR<-REZdod$ED50/refval
#export to pdf 10 x 6 inches
write.table(REZdod,file="output/REZdodcroimyc18.txt",quote=FALSE,sep="\t",
            row.names=FALSE)

hist((REZdod$ED50[order(REZdod$ED50)])/refval,main="dodine",xlab="FR Classes",
     breaks=c(0,10,20,30,40,50,60,70,80),
     las=1,col=heat.colors(8)[8:1],ylim=c(0,30))
abline(v=10,col="red",lwd=3)
#export to pdf 4.5 x 9 inches


###############################################################################
#END
###############################################################################