##############################################################################/
##############################################################################/
#DAPC analyses, results and plot
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#DAPC analysis: loading and preparing the data set####
##############################################################################/

#here is the structure of the data file, for explanation of each columns, see 
#ReadMe.txt file in the repository
head(datAgracc)
#a summary of the different variables
summary(datAgracc)
colnames(datAgracc)
#total number of individuals
dim(datAgracc)[1] #309 individuals
#we reorganize the levels of the host_corrected column, because the 
#alphabetical order doesn't fit our needs
datAgracc$host_ord<-factor(datAgracc$host,
                           levels=c("peach","oilseed_rape","tobacco",
                                    "other_crops","Aerial_trap"))
#we reorder the individuals according to the host_corrected factor
datAgracc<-datAgracc[order(datAgracc$host_ord),]
JDD<-datAgracc #name of the input file
JDD<-drop.levels(JDD)
#let's define a set of color for keeping some consistency in the plots
coloor<-c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")


##############################################################################/
#DAPC on microsatellites markers####
##############################################################################/

#converting data to a genind format, first we use only the microsatellite data
JDDmicro<-df2genind(JDD[,c("MP_27","MP_39","MP_44","MP_5","MP_7","MP_23",
                           "MP_45","MP_28","MP_9","MP_13","MP_2","MP_38",
                           "MP_4","MP_46")],
                    ncode=3,ind.names=JDD$sample_ID, 
                    pop=JDD$host,ploidy=2,NA.char="999")
#include the coordinates of the samples
JDDmicro@other$xy<-JDD[,c("longitude","latitude")]
#we can also include the resistance genotypes as supplementary information
JDDmicro@other$KDR<-JDD[,"KDR"]
JDDmicro@other$sKDR<-JDD[,"sKDR"]
JDDmicro@other$MACE<-JDD[,"MACE"]
JDDmicro@other$R81T<-JDD[,"R81T"]

#now we format the data set to analyse it with DAPC from the adegenet package
JDDade<-JDDmicro
#determination of the number of clusters
clustJDDade<-find.clusters(JDDade,max.n.clust=30)
#with 50 PCs, we lost nearly no information and after K=4, the decrease of 
#the BIC value is smaller, so we chose the maximum number of clusters to be 4 
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#We try to optimize the number of principal component (PCs) to retain to 
#perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=15)
temp<-optim.a.score(dapcJDDade)
#we chose the to keep 8 PCs in order to avoid over fitting of the 
#model. Then we do the actual DAPC analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=8)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
#scatter plot
scatter(dapcJDDade,xax=1, yax=2,col=coloor)
scatter(dapcJDDade,xax=1, yax=3,col=coloor)
scatter(dapcJDDade,xax=2, yax=3,col=coloor)

#Run the 'find.clusters' and DAPC analysis for K=3 and 4
set.seed(227)
clustJDDade3<-find.clusters(JDDade,n.pca=50,n.clust=3)
dapcJDDade3<-dapc(JDDade,clustJDDade3$grp,n.da=4,n.pca=12)
compoplot(dapcJDDade3,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
set.seed(355)
clustJDDade4<-find.clusters(JDDade,n.pca=50,n.clust=4)
dapcJDDade4<-dapc(JDDade,clustJDDade4$grp,n.da=4,n.pca=12)
compoplot(dapcJDDade4,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)

#a more beautiful scatter plot, the colors are matching colors used in 
#the structure-like plot
scatter(dapcJDDade4,xax=1,yax=2,cstar=1,cell=0,clab=0,main="Axis 1 & 2",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)
scatter(dapcJDDade4,xax=2,yax=3,cstar=1,cell=0,clab=0,main="Axis 2 & 3",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)


##############################################################################/
#Figure S3: Scatter plot for best K####
##############################################################################/

#first we plot the BIC variation to determine the "best" K
op<-par(mfrow=c(2,2))
plot(clustJDDade$Kstat,type="o",xlab="number of clusters (K)",ylab="BIC",
     col="blue",main="Detection based on BIC",las=1)
points(3,clustJDDade$Kstat[3],pch=21,lwd=3,cex=1.5,col="red")
points(4,clustJDDade$Kstat[4],pch=21,lwd=3,cex=1.5,col="red")
scatter(dapcJDDade3,xax=1,yax=2,cstar=1,cell=0,clab=0,main="K=3",
        solid=0.3,col=coloor[c(2,1,3)],pch=19,cex=3,
        scree.da=FALSE,scree.pca=TRUE,posi.pca="bottomright")
text(-100,450,labels="K=3, axes 1 and 2",cex=2,xpd=TRUE)
scatter(dapcJDDade4,xax=1,yax=2,cstar=1,cell=0,clab=0,main="K=3",
        solid=0.3,col=coloor[c(4,2,1,3)],pch=19,cex=3,
        scree.da=FALSE,scree.pca=TRUE)
text(250,450,labels="K=4, axes 1 and 2",cex=2,xpd=TRUE)
scatter(dapcJDDade4,xax=2,yax=3,cstar=1,cell=0,clab=0,main="K=3",
        solid=0.3,col=coloor[c(4,2,1,3)],pch=19,cex=3,
        scree.da=FALSE,scree.pca=TRUE,posi.pca="bottomright")
text(-100,450,labels="K=4, axes 2 and 3",cex=2,xpd=TRUE)
par(op)

#export to pdf 10 X 10 inches


##############################################################################/
#Figure S4: Structure-like plot####
##############################################################################/

#first you need to gather the number of individuals in each populations
effpop<-as.numeric(table(JDDade$pop))[c(4,2,5,3,1)]
#the names of the different populations might be useful too
poptiquet<-c("Peach","Oilseed\nrape","Tobacco","Other\nCrops","Aerial Trap")


#Now, we can easily plot several structure-like plot in the same figure
op<-par(mfrow=c(2,1),mar=c(0,4,0,0),oma=c(3,0,0,0))
structplot(t(dapcJDDade3$posterior)[c(2,1,3),],coloor,effpop,poptiquet,
           leg_y="K=3",cexy=1.2,mef=c(0,1,1,0,0),colbord=NA,spacepop=4,
           distxax=0.08)
mtext("DAPC\nK=3",side=2,line=1,cex=1.2,las=1,adj=0.5)
structplot(t(dapcJDDade4$posterior)[c(1,4,2,3),],coloor,effpop,poptiquet,
           leg_y="K=4",cexy=1.2,mef=c(0,1,1,1,0),colbord=NA,spacepop=4)
mtext("DAPC\nK=4",side=2,line=1,cex=1.2,las=1,adj=0.5)
par(op)


#in order to compare DAPC and STRUCTURE output, we prepare the STRUCTURE
#output for the same data set
strK3<-t(datAgracc[,c("K3_Q1","K3_Q2","K3_Q3")])
strK4<-t(datAgracc[,c("K4_Q1","K4_Q2","K4_Q3","K4_Q4")])
#now we plot STRUCTURE and DAPC alongside
#the plot for the different K values
layout(matrix(c(1,1,1,1,
                2,2,2,2,
                3,3,
                4,4,4,4,
                5,5,5,5),18,1,byrow=TRUE))
op<-par(mar=c(0.1,1.1,0.1,0),oma=c(4.1,6.5,1,0),font=2)
structplot(strK3,coloor,effpop,poptiquet,spacepop=4,
           leg_y="K=3",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("STRUCTURE\nK=3",side=2,line=2,cex=1.2,las=1,adj=0.5)
structplot(t(dapcJDDade3$posterior)[c(2,1,3),],coloor,effpop,poptiquet,
           leg_y="K=3",cexy=1.5,mef=c(0,1,1,1,0),colbord=NA,spacepop=4,
           distxax=0.1,cexpop=1.2)
mtext("DAPC\nK=3",side=2,line=2,cex=1.2,las=1,adj=0.5)
plot.new()
structplot(strK4,coloor[c(1,3,2,4,5)],effpop,poptiquet,spacepop=4,
           leg_y="K=3",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("STRUCTURE\nK=4",side=2,line=2,cex=1.2,las=1,adj=0.5)
structplot(t(dapcJDDade4$posterior)[c(3,2,4,1),],coloor,effpop,poptiquet,
           leg_y="K=4",cexy=1.5,mef=c(0,1,1,1,0),colbord=NA,spacepop=4,
           distxax=0.1,cexpop=1.2)
mtext("DAPC\nK=4",side=2,line=2,cex=1.2,las=1,adj=0.5)
par(op)

#export to pdf 12 X 5 inches


##############################################################################/
#END
##############################################################################/