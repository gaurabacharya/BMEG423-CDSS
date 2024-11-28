## Clear workspace
rm(list=ls())

## Load all the data so we can quickly combine it and explore it. 
load(file.path("..","Data","training_2024-11-04.RData"))
SEPSISdat<-sepsis_data
load(file.path("..","Data","testing_2024-11-04.RData"))
nrow(SEPSISdat) # should be n=200112
SEPSISdat_test<-sepsis_data
rm("sepsis_data")

#
## Some data exploration code
id_sepsis<-unique(SEPSISdat$patient[SEPSISdat$SepsisLabel==1])
SEPSISdat$has_sepsis<-F
SEPSISdat$has_sepsis[SEPSISdat$patient %in% id_sepsis] <-T

SEPSISdat_NOsepsis <- SEPSISdat[!SEPSISdat$patient %in% unique(SEPSISdat$patient[SEPSISdat$SepsisLabel==1]),]
SEPSISdat_Sepsis <- SEPSISdat[SEPSISdat$patient %in% unique(SEPSISdat$patient[SEPSISdat$SepsisLabel==1]),]
SEPSISdat_Sepsis <- SEPSISdat_Sepsis[SEPSISdat_Sepsis$SepsisLabel==1,]

##
library('lattice')
cinc_dat_class<- rbind(SEPSISdat_NOsepsis,SEPSISdat_Sepsis)
pdf("./candidates2.pdf", width=12, height=6)
cols <- colnames(cinc_dat_class)
cols <- cols[!(cols %in% c("patient","Age","Gender","Unit1","Unit2","HospAdmTime","ICULOS","SepsisLabel"))]
cinc_dat_class$SepsisLabel <- as.factor(cinc_dat_class$SepsisLabel)
levels(cinc_dat_class$SepsisLabel)<-c(F,T)
n <- 0
for (c in cols){
  f=as.formula(paste0(c,"~SepsisLabel"))
  p<-bwplot(f,data=cinc_dat_class,horizontal=F,
            panel = function(..., box.ratio) {
              panel.violin(..., col = "lightblue",
                           varwidth = FALSE, box.ratio = box.ratio)
            })
  if (n==0){ # x0, y0, xT, yT
    print(p,position=c(0,0.5,0.25,1),more=T)
  }else if(n==1){
    print(p,position=c(0.25,0.5,0.5,1),more=T)
  }else if(n==2){
    print(p,position=c(0.5,0.5,0.75,1),more=T)
  }else if(n==3){
    print(p,position=c(0.75,0.5,1,1),more=T)
  }else if(n==4){
    print(p,position=c(0,0.0,0.25,0.5),more=T)
  }else if(n==5){
    print(p,position=c(0.25,0,0.5,0.5),more=T)
  }else if(n==6){
    print(p,position=c(0.5,0,0.75,0.5),more=T)
  }else if(n==7){
    print(p,position=c(0.75,0,1,0.5))
    n <- -1
  }
  n <- n+1
}
dev.off()
