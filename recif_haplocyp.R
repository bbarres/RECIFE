##############################################################################/
##############################################################################/
#Boxplot Figure for the different cyp51 haplotypes
##############################################################################/
##############################################################################/

#loading the dataset and the necessary library
source("recif_load.R")


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

haploRECIF<-haplo51[haplo51$study_pool=="RECIFE",]
haploRECIF<-haploRECIF[,c("strain_ID","haplotype_cyp51",
                          "CI50_difenoconazole",
                          "CI50_mefentrifluconazole",
                          "CI50_prothioconazole.desthio",
                          "CI50_tetraconazole")]
haploRECIF<-drop.levels(haploRECIF)


##############################################################################/
#boxplot plotting####
##############################################################################/

haploRECIF<-pivot_longer(haploRECIF,c("CI50_difenoconazole",
                                      "CI50_mefentrifluconazole",
                                      "CI50_prothioconazole.desthio",
                                      "CI50_tetraconazole"),
                         names_to="Active_subst",values_to="CI50")
haploRECIF$Active_subst<-as.factor(haploRECIF$Active_subst)

boxplot(haploRECIF$CI50~haploRECIF$Active_subst+haploRECIF$haplotype_cyp51)





##############################################################################/
#END
##############################################################################/