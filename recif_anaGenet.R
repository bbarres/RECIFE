##############################################################################/
##############################################################################/
#DAPC analyses, results and plot
##############################################################################/
##############################################################################/

source("recif_load.R")


##############################################################################/
#Loading and preparing the data set####
##############################################################################/

#here is the structure of the data file, for explanation of each columns, see 
#ReadMe.txt file in the repository
head(sugmic)
#we turn the microsatellite data from factors to numeric
indcol<-colnames(sugmic)[3:14]
sugmic[indcol]<-lapply(sugmic[indcol],
                       function(x) as.numeric(as.character(x)))
#a summary of the different variables
summary(sugmic)
str(sugmic)
#we remove the individuals with genotyping issue
sugmic_clean<-sugmic[sugmic$qualmicro==1,]


#we also build a dataset consisting only of the 3 main populations
sugmic_SUG<-sugmic[sugmic$pop %in% list("A2019","A2020","A2014"),]
sugmic_SUG_clean<-sugmic_SUG[sugmic_SUG$qualmicro==1,]

#total number of individuals genotyped
dim(sugmic)[1] #1274 individuals
#total number of individuals genotyped with a good quality
dim(sugmic_clean)[1] #969 individuals

#total number of individuals included in the study
dim(sugmic_SUG)[1] #1167
#total number of individuals included in the study with a good quality
dim(sugmic_SUG_clean)[1] #890

JDD<-sugmic_SUG_clean #name of the input file
JDD<-drop.levels(JDD)
#let's define a set of color for keeping some consistency in the plots
coloor<-c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")


##############################################################################/
#Importing data and basic statistics of the microsatellites markers####
##############################################################################/

#converting data to a genind format, first we use only the microsatellite data
JDDmicro<-df2genind(JDD[,colnames(JDD)[3:14]],
                    ncode=3,ind.names=JDD$ech_id, 
                    pop=JDD$pop,ploidy=1)
#some basic information
summary(JDDmicro)
strata(JDDmicro)<-data.frame(popu=pop(JDDmicro))
#turning the data set to a genclon object
JDDmicroClon<-as.genclone(JDDmicro)
#missing data by pop and loci
info_table(JDDmicro,plot=TRUE)

#locus statistics
locus_table(JDDmicro)


##############################################################################/
#MLG analysis####
##############################################################################/

JDDstra<-strata(JDDmicroClon) %>% 
  group_by(popu) %>%
  summarize(Count=n())

# Plotting, First three arguments are necessary.
treemap(dtf = JDDstra, index = nameStrata(JDDmicroClon), vSize = "Count",
        type = "categorical", vColor = "popu", title = "Cercospora beticola")

#genotypic diversity assessment
JDDmicro_diversity<-poppr(JDDmicro)
JDDmicro_diversity

#creating clone corrected by population data set
JDDmicroCC<-clonecorrect(JDDmicro,strata =~popu,)


##############################################################################/
#DAPC analysis####
##############################################################################/

#now we format the data set to analyse it with DAPC from the adegenet package
JDDade<-JDDmicroCC
#determination of the number of clusters
clustJDDade<-find.clusters(JDDade,max.n.clust=30)
#with 30 PCs, we lost nearly no information and after K=5, the decrease of 
#the BIC value is smaller, so we chose the maximum number of clusters to be 5 
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#We try to optimize the number of principal component (PCs) to retain to 
#perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=10)
temp<-optim.a.score(dapcJDDade)
#we chose the to keep 6 PCs in order to avoid over fitting of the 
#model. Then we do the actual DAPC analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=6)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
#scatter plot
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,main="K=3",
        solid=0.3,col=coloor,
        pch=19,cex=3,
        scree.da=FALSE,scree.pca=FALSE,posi.pca="bottomright")
points(dapcJDDade$ind.coord[JDDade$pop==levels(JDDade$pop)[1],1],
       dapcJDDade$ind.coord[JDDade$pop==levels(JDDade$pop)[1],2],
       pch=20,col="black")
points(dapcJDDade$ind.coord[JDDade$pop==levels(JDDade$pop)[2],1],
       dapcJDDade$ind.coord[JDDade$pop==levels(JDDade$pop)[2],2],
       pch=4,col="black")
points(dapcJDDade$ind.coord[JDDade$pop==levels(JDDade$pop)[3],1],
       dapcJDDade$ind.coord[JDDade$pop==levels(JDDade$pop)[3],2],
       pch=17,col="black")

pairwise.WCfst(JDDade)


##############################################################################/
#Linkage disequilibrium analyses and data preparation####
##############################################################################/

coloor<-c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")

JDDmicroCC_genp<-"data/micro_genepop_CC.txt"
test_LD(JDDmicroCC_genp,"output/micro_genepop_CC.txt.DIS")
clean_workdir()

#the output file have to be edited so it can be used for representation
#warning: some space between columns are not present for the tests for
#all population. The input file must therefore be edited before importation
LDall<-readLines("output/micro_genepop_CC.txt.DIS")
#removing the first lines that are non relevant
LDcolnames<-unlist(strsplit(LDall[217],"\\s+"))
LDcolnames<-c("Loc1","Loc2",LDcolnames[3:5])
LDcolnames2<-unlist(strsplit(LDall[13],"\\s+"))
LDall<-LDall[-c(1:14)]
#for 12 markers and 3 populations there are 
12*11/2*3 #LD tests for markers pairs
#by population
LDbypops<-LDall[-c(198+1:length(LDall))]
LDbypops<-as.data.table(matrix(unlist(strsplit(LDbypops,"\\s+")),
                               nrow=length(strsplit(LDbypops,"\\s+")),
                               byrow=TRUE))
colnames(LDbypops)<-LDcolnames2
LDbypops$`P-Value`<-as.numeric(LDbypops$`P-Value`)

#for France 2019 population
LDfran19<-as.matrix(LDbypops[c(1:66),])
row.names(LDfran19)<-paste(LDfran19[,2],LDfran19[,3])
LDfran19<-rbind(LDfran19[,c(2,3,4)],LDfran19[,c(3,2,4)])
LDfran19<-spread(data.frame(LDfran19),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDfran19)<-as.character(LDfran19[,1])
#removing the first unnecessary column
LDfran19<-LDfran19[,c(-1)]
#reordering the columns and turning the object into a matrix
LDfran19<-as.matrix(LDfran19[c(1:2,11:12,3:10),
                             c(1:2,11:12,3:10)])

#for France 2020 population
LDfran20<-as.matrix(LDbypops[c(67:132),])
row.names(LDfran20)<-paste(LDfran20[,2],LDfran20[,3])
LDfran20<-rbind(LDfran20[,c(2,3,4)],LDfran20[,c(3,2,4)])
LDfran20<-spread(data.frame(LDfran20),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDfran20)<-as.character(LDfran20[,1])
#removing the first unnecessary column
LDfran20<-LDfran20[,c(-1)]
#reordering the columns and turning the object into a matrix
LDfran20<-as.matrix(LDfran20[c(1:2,11:12,3:10),
                             c(1:2,11:12,3:10)])

#for USA 2014 population
LDUSA14<-as.matrix(LDbypops[c(133:198),])
row.names(LDUSA14)<-paste(LDUSA14[,2],LDUSA14[,3])
LDUSA14<-rbind(LDUSA14[,c(2,3,4)],LDUSA14[,c(3,2,4)])
LDUSA14<-spread(data.frame(LDUSA14),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDUSA14)<-as.character(LDUSA14[,1])
#removing the first unnecessary column
LDUSA14<-LDUSA14[,c(-1)]
#reordering the columns and turning the object into a matrix
LDUSA14<-as.matrix(LDUSA14[c(1:2,11:12,3:10),
                             c(1:2,11:12,3:10)])

#for all populations
LDallpop<-LDall[c((198+7):(length(LDall)-2))]
LDallpop<-as.data.table(matrix(unlist(strsplit(LDallpop,"\\s+")),
                               nrow=length(strsplit(LDallpop,"\\s+")),
                               byrow=TRUE))
LDallpop<-as.data.table(LDallpop[,c(1,3,4:6)])
colnames(LDallpop)<-LDcolnames
LDallpop$`P-Value`<-as.numeric(sub("<","",LDallpop$`P-Value`,fixed=TRUE))
row.names(LDallpop)<-paste(LDallpop$Loc1,LDallpop$Loc2,sep="_")
LDallpop<-as.matrix(LDallpop)
LDallpop<-rbind(LDallpop[,c(1,2,5)],LDallpop[,c(2,1,5)])
LDallpop<-spread(data.frame(LDallpop),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDallpop)<-as.character(LDallpop[,1])
#removing the first unnecessary column
LDallpop<-LDallpop[,c(-1)]
#reordering the columns and turning the object into a matrix
LDallpop<-as.matrix(LDallpop[c(1:2,11:12,3:10),
                             c(1:2,11:12,3:10)])


##############################################################################/
#Figure with 4 panels
##############################################################################/

#panel figure: France 2019 population
temp<-LDfran19
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00076 & temp<0.8]<-0.6
temp[temp>0.00015 & temp<0.6]<-0.4
temp[temp>0.000015 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP1<-viewport(x=0,y=0,width=1/3,height=1,just=c("left","bottom"),
              name="vp1")
pushViewport(VP1)
LDall<-LDheatmap(temp,title=NULL,
                 add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                 color=c(rep(grey(0.8),3),
                         brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                 name="fran19",flip=TRUE,add.key=FALSE,
                 newpage=FALSE)
grid.edit(gPath("fran19","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("fran19","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.05,font=2),
          rot=0,hjust=0.87)
grid.text("France 2019",x=unit(0.5,"npc"),y=unit(0.8,"npc"),
          gp=gpar(fontsize=35),check=TRUE)
upViewport()

#panel figure: France 2020 population
temp<-LDfran20
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00076 & temp<0.8]<-0.6
temp[temp>0.00015 & temp<0.6]<-0.4
temp[temp>0.000015 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP2<-viewport(x=1/3,y=0,width=1/3,height=1,just=c("left","bottom"),
              name="vp2")
pushViewport(VP2)
LDall<-LDheatmap(temp,title=NULL,
                 add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                 color=c(rep(grey(0.8),3),
                         brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                 name="fran20",flip=TRUE,add.key=FALSE,
                 newpage=FALSE)
grid.edit(gPath("fran20","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("fran20","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.05,font=2),
          rot=0,hjust=0.87)
grid.text("France 2020",x=unit(0.5,"npc"),y=unit(0.8,"npc"),
          gp=gpar(fontsize=35),check=TRUE)
upViewport()

#panel figure: USA 2014 population
temp<-LDUSA14
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00076 & temp<0.8]<-0.6
temp[temp>0.00015 & temp<0.6]<-0.4
temp[temp>0.000015 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP3<-viewport(x=2/3,y=0,width=1/3,height=1,just=c("left","bottom"),
              name="vp3")
pushViewport(VP3)
LDall<-LDheatmap(temp,title=NULL,
                 add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                 color=c(rep(grey(0.8),3),
                         brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                 name="USA14",flip=TRUE,add.key=FALSE,
                 newpage=FALSE)
grid.edit(gPath("USA14","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("USA14","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.05,font=2),
          rot=0,hjust=0.87)
grid.text("USA 2014",x=unit(0.5,"npc"),y=unit(0.8,"npc"),
          gp=gpar(fontsize=35),check=TRUE)
upViewport()

#export to pdf landscape 4.5 x 15 inches





#extra figure for fun
#panel figure: all populations
temp<-LDallpop
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00076 & temp<0.8]<-0.6
temp[temp>0.00015 & temp<0.6]<-0.4
temp[temp>0.000015 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP1<-viewport(x=0,y=0.5,width=0.5,height=0.5,just=c("left","bottom"),
              name="vp1")
pushViewport(VP1)
LDall<-LDheatmap(temp,title=NULL,
                 add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                 color=c(rep(grey(0.8),3),
                         brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                 name="allpop",flip=TRUE,add.key=FALSE,
                 newpage=FALSE)
grid.edit(gPath("allpop","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("allpop","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("All populations",x=unit(0.5,"npc"),y=unit(0.75,"npc"),
          gp=gpar(fontsize=40),check=TRUE)
upViewport()


##############################################################################/
#END
##############################################################################/