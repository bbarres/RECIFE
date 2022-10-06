##############################################################################/
##############################################################################/
#Boxplot Figure for the different cyp51 haplotypes
##############################################################################/
##############################################################################/

#loading the dataset and the necessary library
source("recif_load.R")


##############################################################################/
#boxplot plotting####
##############################################################################/



CompRez$ED50<-as.character(CompRez$ED50)
CompRez[CompRez$ED50==">10","ED50"]<-12
CompRez[CompRez$ED50==">20","ED50"]<-22
CompRez[CompRez$ED50==">50","ED50"]<-52
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))









##############################################################################/
#END
##############################################################################/