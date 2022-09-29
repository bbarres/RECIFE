##############################################################################/
##############################################################################/
#Code for producing Figure 2
##############################################################################/
##############################################################################/

##loading the dataset and the necessary library
CompRez<-read.table(file)
#define a color vector
cooloor<-brewer.pal(12,"Set3")[c(1,2,3,5,6,7,8,9,10,11)]
#preparing the dataset
temp<-CompRez[,c(1,2,4)]
temp<-spread(temp,Subs_Act,ED50)

##############################################################################/
#Regression analysis of mycelial growth experiment scoring 20 or 21 days####
##############################################################################/




#same plot but combine on one figure and with log(EC50)
plot(temp[order(c(temp$CYPROCONAZOLE)),"CYPROCONAZOLE"],
     main="EC50 distribution",bg=cooloor[1],pch=21,cex=2,las=1,
     ylab="EC50",xlab="",ylim=c(0,52))
points(temp[order(c(temp$DIFENOCONAZOLE)),"DIFENOCONAZOLE"],
       bg=cooloor[2],pch=21,cex=2,las=1)
points(temp[order(c(temp$EPOXICONAZOLE)),"EPOXICONAZOLE"],
       bg=cooloor[3],pch=21,cex=2,las=1)
points(temp[order(c(temp$FLUTRIAFOL)),"FLUTRIAFOL"],
       bg=cooloor[4],pch=21,cex=2,las=1)
points(temp[order(c(temp$MEFENTRIFLUCONAZOLE)),"MEFENTRIFLUCONAZOLE"],
       bg=cooloor[5],pch=21,cex=2,las=1)
points(temp[order(c(temp$METCONAZOLE)),"METCONAZOLE"],
       bg=cooloor[6],pch=21,cex=2,las=1)
points(temp[order(c(temp$PROCHLORAZE)),"PROCHLORAZE"],
       bg=cooloor[7],pch=21,cex=2,las=1)
points(temp[order(c(temp$`PROTHIOCONAZOLE-DESTHIO`)),
                "PROTHIOCONAZOLE-DESTHIO"],
       bg=cooloor[8],pch=21,cex=2,las=1)
points(temp[order(c(temp$TEBUCONAZOLE)),"TEBUCONAZOLE"],
       bg=cooloor[9],pch=21,cex=2,las=1,)
points(temp[order(c(temp$TETRACONAZOLE)),"TETRACONAZOLE"],
       bg=cooloor[10],pch=21,cex=2,las=1)
legend(0,55,legend=c("cyproconazole","difénoconazole","époxiconazole",
                      "flutriafol","méfentrifluconazole","metconazole",
                      "prochloraze","prothioconazole-desthio","tébuconazole",
                      "tétraconazole"),
       cex=1,pt.cex=1.3,
       y.intersp=0.7,x.intersp=0.9,
       pch=c(15),
       col=cooloor,
       bty="n")
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


pairs(log(temp[,c(2:4,7:13)]),
      text.panel=diag_custom_labels(c(CYPROCONAZOLE="cyproconazole",
                                      DIFENOCONAZOLE="difénoconazole",
                                      "époxiconazole",
                                      "flutriafol","méfentrifluconazole","metconazole",
                                      "prochloraze","prothioconazole-desthio","tébuconazole",
                                      "tétraconazole")),
      lower.panel=panel.smooMod,upper.panel=panel.cor,las=1)





##############################################################################/
#END
##############################################################################/