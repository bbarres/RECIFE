##############################################################################/
##############################################################################/
#Boxplot Figure for the different cyp51 haplotypes
##############################################################################/
##############################################################################/

#loading the data set and the necessary library
source("recif_load.R")

#define a color vector
cooloor<-brewer.pal(12,"Set3")[c(1,2,3,5,6,7,8,9,10,11)]
names(cooloor)<-c("cyproconazole","difénoconazole","époxiconazole",
                  "flutriafol","méfentrifluconazole","metconazole",
                  "prochloraze","prothio-desthio","tébuconazole",
                  "tétraconazole")
#pick one color for the background
fondplot<-brewer.pal(11,"RdYlBu")[7]


##############################################################################/
#preparing the data set####
##############################################################################/

haplo51[,11:22]<-replace(haplo51[,11:22],haplo51[,11:22]==">10","12")
haplo51[,11:22]<-replace(haplo51[,11:22],haplo51[,11:22]==">20","22")
haplo51[,11:22]<-replace(haplo51[,11:22],haplo51[,11:22]==">50","52")
haplo51[,11:22]<-lapply(haplo51[,11:22],as.numeric)
haplo51$haplotype_cyp51<-factor(haplo51$haplotype_cyp51,
                                levels=c("1","3","4","7","8","9",
                                         "10","11","12","13","15",
                                         "16","17","18"))
colnames(haplo51)[11:22]<-c("cyproconazole","difénoconazole",
                            "époxiconazole","fenpropidine","fluopyram",
                            "méfentrifluconazole","metconazole","prochloraze",
                            "prothioconazole-desthio","tébuconazole",
                            "tétraconazole","tolnaftate")

#the haplotype found in the RECIFE study
haploRECIF<-haplo51[haplo51$study_pool=="RECIFE",]
haploRECIF<-haploRECIF[,c("strain_ID","haplotype_cyp51",
                          "difénoconazole",
                          "méfentrifluconazole",
                          "prothioconazole-desthio",
                          "tétraconazole")]
haploRECIF<-drop.levels(haploRECIF)
haploRECIF<-pivot_longer(haploRECIF,c("difénoconazole",
                                      "méfentrifluconazole",
                                      "prothioconazole-desthio",
                                      "tétraconazole"),
                         names_to="ActivSub",values_to="CI50")
haploRECIF$ActivSub<-as.factor(haploRECIF$ActivSub)


##############################################################################/
#boxplot plotting by haplotype####
##############################################################################/

pdf("output/haplopheno_byhaplo.pdf",width=13,height=6)
atat<-boxplot(haploRECIF$CI50~haploRECIF$ActivSub:haploRECIF$haplotype_cyp51,
              las=1,col=cooloor[c(2,5,8,10)],at=c(1:54)[-c(seq(5,54,5))],
              xaxt="n",yaxt="n",frame.plot=FALSE,ann=FALSE)
rect(0,-10,5,62,col="white")
rect(5,-10,10,62,col=fondplot)
rect(15,-10,20,62,col=fondplot)
rect(25,-10,30,62,col=fondplot)
rect(35,-10,40,62,col=fondplot)
rect(45,-10,50,62,col=fondplot)
rect(50,-10,55,62,col="white")
text(x=seq(2.5,52.5,5),y=50,font=4,
     labels=paste("n=",atat$n[(0:10)*4+1],sep=""))
axis(1,lwd=2,at=seq(2.5,52.5,5),font=2,
     labels=levels(haploRECIF$haplotype_cyp51))
mtext("haplotype",side=1,font=2,cex=1.4,line=2.5)
axis(2,lwd=2,las=1,font=2)
legeText<-expression(paste(CI[50]," (mg.",L^-1,")",sep=""))
mtext(legeText,side=2,font=2,cex=1.3,line=2.5)
legend(45.5,48,legend=levels(haploRECIF$ActivSub),
       cex=1,pt.cex=1.3,bg="white",
       y.intersp=1,x.intersp=0.8,
       pch=c(22),
       pt.bg=cooloor[c(2,5,8,10)])
box(bty="o",lwd=2)

boxplot(haploRECIF$CI50~haploRECIF$ActivSub:haploRECIF$haplotype_cyp51,
        las=1,col=cooloor[c(2,5,8,10)],at=c(1:54)[-c(seq(5,54,5))],
        xaxt="n",yaxt="n",frame.plot=FALSE,ann=FALSE,add=TRUE)

points(x=c(1:54)[-c(seq(5,54,5))],
       y=tapply(haploRECIF$CI50,
                INDEX=haploRECIF$haplotype_cyp51:haploRECIF$ActivSub,
                FUN=mean,na.rm=TRUE),
       pch=23,bg=cooloor[c(2,5,8,10)])

byhaploCI<-tapply(haploRECIF$CI50,
                  INDEX=haploRECIF$haplotype_cyp51,
                  FUN=mean,na.rm=TRUE)
segments(x0=seq(from=1,to=54,by=5)-0.5,
         y0=byhaploCI,
         x1= seq(from=4,to=54,by=5)+0.5,
         y1=byhaploCI,col="red",lty=2,lwd=3)
dev.off()
#export to .pdf 13 x 6 inches


##############################################################################/
#boxplot plotting by Active substance####
##############################################################################/

pdf("output/haplopheno_bySA.pdf",width=13,height=6)
atat<-boxplot(haploRECIF$CI50~haploRECIF$haplotype_cyp51:haploRECIF$ActivSub,
              las=1,col=rep(cooloor[c(2,5,8,10)],each=11),
              at=c(1:47)[-c(12,24,36)],
              xaxt="n",yaxt="n",frame.plot=FALSE,ann=FALSE)
rect(0,-10,12,62,col="white")
rect(12,-10,24,62,col=fondplot)
rect(36,-10,48,62,col=fondplot)
axis(1,lwd=2,at=c(6,18,30,42),font=2,
     labels=levels(haploRECIF$ActivSub))
mtext("Substance active",side=1,font=2,cex=1.4,line=2.5)
axis(2,lwd=2,las=1,font=2)
legeText<-expression(paste(CI[50]," (mg.",L^-1,")",sep=""))
mtext(legeText,side=2,font=2,cex=1.3,line=2.5)
box(bty="o",lwd=2)

boxplot(haploRECIF$CI50~haploRECIF$haplotype_cyp51:haploRECIF$ActivSub,
        las=1,col=rep(cooloor[c(2,5,8,10)],each=11),
        at=c(1:47)[-c(12,24,36)],
        xaxt="n",yaxt="n",frame.plot=FALSE,ann=FALSE,add=TRUE)

points(x=c(1:47)[-c(12,24,36)],
       y=tapply(haploRECIF$CI50,
                INDEX=haploRECIF$ActivSub:haploRECIF$haplotype_cyp51,
                FUN=mean,na.rm=TRUE),
       pch=23,bg=rep(cooloor[c(2,5,8,10)],each=11))

byactiSACI<-tapply(haploRECIF$CI50,
                  INDEX=haploRECIF$ActivSub,
                  FUN=mean,na.rm=TRUE)
segments(x0=seq(from=1,to=48,by=12)-0.5,
         y0=byactiSACI,
         x1= seq(from=11,to=48,by=12)+0.5,
         y1=byactiSACI,col="red",lty=2,lwd=3)
dev.off()
#export to .pdf 13 x 6 inches


##############################################################################/
#END
##############################################################################/