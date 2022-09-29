##############################################################################/
##############################################################################/
#Code for producing Figure 2
##############################################################################/
##############################################################################/

##loading the dataset and the necessary library
CompRez<-read.table(file="output/ASA_results_cerco.txt",header=TRUE,
                    sep="\t")
CompRez$ED50<-as.character(CompRez$ED50)
CompRez[CompRez$ED50==">10","ED50"]<-12
CompRez[CompRez$ED50==">20","ED50"]<-22
CompRez[CompRez$ED50==">50","ED50"]<-52
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))

#define a color vector
cooloor<-brewer.pal(12,"Set3")[c(1,2,3,5,6,7,8,9,10,11)]
#preparing the dataset for cross active substance analyses
multSA<-CompRez[,c(1,2,4)]
multSA<-spread(multSA,Subs_Act,ED50)
row.names(multSA)<-multSA$sample_ID
#we remove active substance that don't have enough information fenpropidine
#fluopyram and tolnaftate
multSA<-multSA[,-c(1,5,6,14)]
colnames(multSA)<-c("cyproconazole","difénoconazole","époxiconazole",
                    "flutriafol","méfentrifluconazole","metconazole",
                    "prochloraze","prothioconazole-desthio","tébuconazole",
                    "tétraconazole")


##############################################################################/
#Regression analysis of mycelial growth experiment scoring 20 or 21 days####
##############################################################################/

#ordered distribution of the EC50 of strains for the different active
#substances
plot(multSA[order(c(multSA[,1])),c(1)],
     main="EC50 distribution",bg=cooloor[1],pch=21,cex=2,las=1,
     ylab="EC50",xlab="",
     ylim=c(min(multSA,na.rm=TRUE),max(multSA,na.rm=TRUE)+10),log="y")
for (i in 2:10) {
  points(multSA[order(c(multSA[,i])),c(i)],
         bg=cooloor[i],pch=21,cex=2,las=1)
  
}

legend(70,0.4,legend=colnames(multSA),
       cex=1,pt.cex=1.3,
       y.intersp=0.7,x.intersp=0.9,
       pch=c(15),
       col=cooloor,
       bty="n")
#export to .pdf 8 x 7 inches

#boxplot of the EC50 for the different active substances
vioplot(multSA,wex=0.7,col=cooloor,las=1)
stripchart(multSA,cex=0.9,frame=FALSE,yaxt='n',xaxt='n',method="jitter",
           offset=1,add=TRUE,vertical=TRUE,pch=21,
           col=adjustcolor("black",alpha=0.5),jitter=0.25)
points(apply(multSA,MARGIN=2,FUN=mean,na.rm=TRUE),
       pch=17,cex=1.5,col="red",bg=cooloor,font=2)
#export to .pdf 8 x 7 inches

#a function to compute the absolute correlation between pairs of variables
panel.cor<-function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r, 2), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#a small modification of the 'panel.smooth' function for 'pairs' plot
panel.smooMod<-function(x,y,col=par("col"),bg=NA,pch=par("pch"), 
                        cex=1,col.smooth=2,span=2/3,iter=3,
                        lwd=par("lwd"),...) 
{
  points(x,y,pch=pch,col=col,bg=bg,cex=cex)
  ok<-is.finite(x) & is.finite(y)
  if(any(ok)) 
    lines(stats::lowess(x[ok],y[ok],f=span,iter=iter), 
          col=col.smooth,lwd=lwd*2,...)
}

diag_custom_labels<-function(x,y,labels,cex,bg,...)
{
  if (!exists('i')) i <<- 1
  text(mean(x),mean(y),
       labels=c("cyproconazole","difénoconazole","époxiconazole",
                "flutriafol","méfentrifluconazole","metconazole",
                "prochloraze","prothioconazole-desthio","tébuconazole",
                "tétraconazole")[[i]],
       cex=c(1,1,1,1,1,1,1,1,1,1)[[i]],bg=cooloor[[i]])
  i<<-i + 1
}

diag_custom_labels<-function(labels) {
  function(x,y,lbl, ...) {
    if (lbl %in% names(labels)) lbl<-labels[[lbl]]
    text(x,y,lbl,col=cooloor,...)
  }
}


pairs(log(multSA),
      text.panel=diag_custom_labels(c(CYPROCONAZOLE="cyproconazole",
                                      DIFENOCONAZOLE="difénoconazole",
                                      "époxiconazole",
                                      "flutriafol","méfentrifluconazole","metconazole",
                                      "prochloraze","prothioconazole-desthio","tébuconazole",
                                      "tétraconazole")),
      lower.panel=panel.smooMod,upper.panel=panel.cor,las=1)



https://stackoverflow.com/questions/31851537/set-backgroud-color-of-a-panel-in-pairs-function-call

##############################################################################/
#END
##############################################################################/