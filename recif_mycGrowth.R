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
#Regression analysis of mycelial growth experiment
###############################################################################

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    ED50=numeric())
#we make a subselection of the data according to the SA
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                  & data_subSA$rslt_03>50,
                                  "ech_id"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                           ED50=numeric()),
         REZSA<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "ED50"=max(data_subSA$dose)))
  #we limit the dataset to the sample that reach somehow a IC of 50%
  SA.dat<-data_subSA[!(data_subSA$ech_id %in% SA_rez),]
  SA.dat<-drop.levels(SA.dat)
  for (i in 1:dim(table(SA.dat$ech_id))[1]) {
    temp.m1<-drm(rslt_03~dose,
                 data=SA.dat[SA.dat$ech_id==names(table(SA.dat$ech_id))[i],],
                 fct=LL.4())
    plot(temp.m1,main=names(table(SA.dat$ech_id))[i],
         ylim=c(0,110),xlim=c(0,50))
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("Subs_Act"=SAlist[j],
                      "sample_ID"=names(table(SA.dat$ech_id))[i],
                      "ED50"=temp[1])
    REZSA<-rbind(REZSA,tempx)
  }
  CompRez<-rbind(CompRez,REZSA)

}


###############################################################################
#END
###############################################################################