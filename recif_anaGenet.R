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
#Linkage disequilibrium analyses####
##############################################################################/

coloor<-c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")
JDDmicroCCgenp<-genind2genpop(JDDmicroCC)

onebyhost<-"data/onebyhost.txt"
test_LD(onebyhost,"output/onebyhost.txt.DIS")
clean_workdir()

#the output file have to be edited so it can be used for representation
LDbyhost<-readLines("output/onebyhost.txt.DIS")
#removing the first lines that are non relevant
LDcolnames<-unlist(strsplit(LDbyhost[13],"\\s+"))
LDbyhost<-LDbyhost[-c(1:14)]
#for 14 markers and 4 populations there are 
14*13/2*4 #LD tests for markers pairs
LDbyhost<-LDbyhost[-c(365:length(LDbyhost))]
LDbyhost<-as.data.table(matrix(unlist(strsplit(LDbyhost,"\\s+")),
                               nrow=length(strsplit(LDbyhost,"\\s+")),
                               byrow=TRUE))
colnames(LDbyhost)<-LDcolnames
LDbyhost$`P-Value`<-as.numeric(LDbyhost$`P-Value`)

#Linkage disequilibrium in the peach population
LDpeach<-as.matrix(LDbyhost[c(1:91),])
row.names(LDpeach)<-paste(LDpeach[,2],LDpeach[,3])
LDpeach<-rbind(LDpeach[,c(2,3,4)],LDpeach[,c(3,2,4)])
LDpeach<-spread(data.frame(LDpeach),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDpeach)<-as.character(LDpeach[,1])
#removing the first unnecessary column
LDpeach<-LDpeach[,c(-1)]
#reordering the columns and turning the object into a matrix
LDpeach<-as.matrix(LDpeach[c(2,8,12,13,14,1,3:7,9:11),
                           c(2,8,12,13,14,1,3:7,9:11)])

#Linkage disequilibrium in the oil seed population
LDoils<-as.matrix(LDbyhost[c(92:182),])
row.names(LDoils)<-paste(LDoils[,2],LDoils[,3])
LDoils<-rbind(LDoils[,c(2,3,4)],LDoils[,c(3,2,4)])
LDoils<-spread(data.frame(LDoils),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDoils)<-as.character(LDoils[,1])
#removing the first unnecessary column
LDoils<-LDoils[,c(-1)]
#reordering the columns and turning the object into a matrix
LDoils<-as.matrix(LDoils[c(2,8,12,13,14,1,3:7,9:11),
                         c(2,8,12,13,14,1,3:7,9:11)])


##############################################################################/
#END
##############################################################################/