##############################################################################/
##############################################################################/
#Code for producing Figure 2
##############################################################################/
##############################################################################/

##loading the data set and the necessary library
source("recif_load.R")

##loading the data set and the necessary library, first run the 
#recif_CI50growth.R script
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
                    "prochloraze","prothio-desthio","tébuconazole",
                    "tétraconazole")


##############################################################################/
#correlation between IC50####
##############################################################################/

#a function to compute the absolute correlation between pairs of variables
cols=brewer.pal(9,"YlOrRd") #goes from red to white to blue
pal=colorRampPalette(cols)
cor_colors=data.frame(correlation=seq(0,1,0.01), 
                      correlation_color=pal(101)[1:101])
cor_colors$correlation_color=as.character(cor_colors$correlation_color)

panel.cor<-function(x,y,digits=2,prefix="",cex.cor, ...) {
  par(usr=c(0,1,0,1))
  u<-par('usr') 
  names(u)<-c("xleft","xright","ybottom","ytop")
  r<-abs(cor(x,y,use="pairwise.complete.obs"))
  bgcolor=cor_colors[r*100,2]    # converts correlation into a specific color
  do.call(rect,c(col=bgcolor,as.list(u))) # colors the correlation box
  txt<-format(c(r,2),digits=digits)[1]
  txt<-paste(prefix,txt,sep="")
  if(missing(cex.cor)) cex.cor<-0.3
  text(0.5,0.5,txt,cex=cex.cor*r*10,font=2)
}

#a small modification of the 'panel.smooth' function for 'pairs' plot
panel.smooMod<-function(x,y,col=par("col"),bg=NA,pch=par("pch"), 
                        cex=1,col.smooth=2,span=2/3,iter=3,
                        lwd=par("lwd"),...) {
  points(x,y,pch=pch,col=col,bg=bg,cex=cex)
  ok<-is.finite(x) & is.finite(y)
  if(any(ok)) 
    lines(stats::lowess(x[ok],y[ok],f=span,iter=iter), 
          col=col.smooth,lwd=lwd*2,...)
}

diag_custom_labels<-function(x,...) {
  par(usr=c(0,1,0,1))
  u<-par('usr')
  names(u)<-c("xleft","xright","ybottom","ytop")
  if (!exists('i')) i<<-1
  do.call(rect,c(col=cooloor[i],as.list(u)))
  text(0.5,0.5,colnames(multSA)[i],font=2,cex=1.2)
  i<<-i + 1
}

#actual plotting
pdf(file="output/correlation_plot.pdf",width=16,height=10)
i<-1
pairs(log(multSA),
      lower.panel=panel.smooMod,
      upper.panel=panel.cor,text.panel="",
      diag.panel=diag_custom_labels,las=1)
mtext(text="C",side=3,cex=4,at=0,font=2,las=1,adj=1.4,line=1.1)
#export to .pdf 16x10 inches
dev.off()


##############################################################################/
#distribution of the IC50 for the different active substances####
##############################################################################/

pdf(file="output/distribution_CI50.pdf",width=12,height=6)
op<-par(mfrow=c(1,2),mar=c(5.1,4.1,3.1,2.1))
#boxplot of the EC50 for the different active substances
vioplot(multSA,wex=0.7,col=cooloor,las=1,xaxt="n",yaxt="n",
        frame.plot=FALSE,ylog=FALSE,log="",
        font.lab=2,font=2)
stripchart(multSA,cex=0.9,frame=FALSE,yaxt='n',xaxt='n',method="jitter",
           offset=1,add=TRUE,vertical=TRUE,pch=21,
           col=adjustcolor("black",alpha=0.5),jitter=0.25)
points(apply(multSA,MARGIN=2,FUN=mean,na.rm=TRUE),
       pch=17,cex=1.5,col="red",bg=cooloor,font=2)
axis(1,lwd=2,at=c(1:10),labels=NA)
text(labels=colnames(multSA),x=1:10,y=rep(-4,10),cex=0.9,font=2,
     srt=-35,xpd=TRUE,adj=0)
axis(2,lwd=2,las=1,font=2)
legeText<-expression(paste(CI[50]," (mg.",L^-1,")",sep=""))
mtext(legeText,side=2,font=2,cex=1.3,line=2.5)
box(bty="o",lwd=2)
mtext(text="A",side=3,cex=3,at=0,font=2,las=1,adj=1.4,line=0)
#export to .pdf 8 x 7 inches

#ordered distribution of the EC50 of strains for the different active
#substances
plot(multSA[order(c(multSA[,1])),c(1)],
     main="",bg=cooloor[1],pch=21,cex=1.5,las=1,
     ylab="",xlab="",xaxt="n",yaxt="n",bty="n",
     ylim=c(min(multSA,na.rm=TRUE),max(multSA,na.rm=TRUE)+10),log="y")
for (i in 2:10) {
  points(multSA[order(c(multSA[,i])),c(i)],
         bg=cooloor[i],pch=21,cex=1.5,las=1)
}
legend(65,0.5,legend=colnames(multSA),
       cex=1,pt.cex=1.3,
       y.intersp=1,x.intersp=0.8,
       pch=c(15),
       col=cooloor,
       bty="n")
axis(1,lwd=2,at=seq(0,110,by=10),font=2)
axis(2,lwd=2,las=1,font=2,
     at=c(0.005,0.01,0.05,0.1,0.5,1,5,10,50),
     labels=c("0.005","0.01","0.05","0.1","0.5","1","5","10","50"))
box(bty="o",lwd=2)
mtext(legeText,side=2,font=2,cex=1.3,line=2.5)
mtext(text="B",side=3,cex=3,at=-4,font=2,las=1,adj=1.4,line=0)
#export to .pdf 8 x 7 inches
par(op)
#export to .pdf 12 x 6 inches
dev.off()


##############################################################################/
#END
##############################################################################/