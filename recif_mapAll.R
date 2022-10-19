##############################################################################/
##############################################################################/
#Maps for population germination tests for old SA
##############################################################################/
##############################################################################/

#loading the dataset and the necessary library
source("recif_load.R")


##############################################################################/
#defining additional function for the mapping####
##############################################################################/

#function for a scale, found in "Auxiliary Cartographic Functions in R: 
#North Arrow, Scale Bar, and Label with a Leader Arrow", Tanimura et al 2007, 
#J of Statistical software
#The code has been slightly modified in order to convert the meter in km
scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
  if(missing(loc)) stop("loc is missing")
  if(missing(length)) stop("length is missing")
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- (x[c(1,3)]-loc[1])/1000
  labels <- append(labels,paste((x[5]-loc[1])/1000,unit))
  text(x[c(1,3,5)],y[4],labels=labels,adj=c(0.5,0),cex=division.cex)
}


##############################################################################/
#Map of all samples collected####
##############################################################################/

pdf("output/mapAll.pdf",width=12,height=8)
op<-par(mar=c(0,0,0,0))
plot(DEP_SHP.1,main="",border="grey70")
plot(REG_SHP.1,lwd=2,add=TRUE)
points(
  x=as.numeric(AllSamp$gps_long),
  y=as.numeric(AllSamp$gps_lat),
  col="darkblue",                 #colors of the points
  pch=3,                          #plotting character
  cex=1                           #size of the points
)
#adding a scale bar to the map
scalebar(c(130000,7050000),300000,"km",division.cex=1.3)
par(op)
dev.off()
#export to .pdf 12 x 8 inches


##############################################################################/
#END
##############################################################################/