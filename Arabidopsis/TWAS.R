setwd("TWAS/Arabidopsis/") #Please adjust the path
##From Delin Li delin.bio@gmail.com Publication: https://academic.oup.com/plphys/article/186/4/1800/6212071

#The arabidopsis FT16 in TWAS paper
#"packages for  GAPIT.Version is 2018.08.18, GAPIT 3.0
library(multtest)
library(gplots)
library(genetics)
library(compiler) #this library is already installed in R library("scatterplot3d")
library(scatterplot3d)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
library(bigmemory)
library(biganalytics)
require(compiler) #for cmpfun

#for manhattan plot
source("PlotFunction.R")

Pheno<- read.csv("study_12_values.csv") 


myGM<-read.delim("Geno.Infor.Arabidopsis.txt")
myGD<-read.delim("Arabidopsis.RNA.Quan5.Geno.Cutoff1.txt.gz")

myY<-Pheno[  !is.na(Pheno[,3]) ,c(1,3)]
myY<-myY[order(as.character(myY[,1])),]

myGD.tem<-myGD[myGD$taxa %in% myY[,1],]
myGD.tem<-myGD.tem[order(as.character(myGD.tem$taxa)),]
myGD.tem <- myGD.tem[,apply(myGD.tem,2,function(x){return(sum(is.na(x))==0)})]

myGM.tem<-myGM[myGM$Name %in% colnames(myGD.tem),]
trait<-colnames(myY)[2]

myY<-myY[ myY[,1] %in% myGD$taxa,]

#Run TWAS with CMLM
myGAPIT <- GAPIT(Y=myY,
                 GD=myGD.tem, GM=myGM.tem,
                 PCA.total=3,
                 kinship.cluster=("average"),
                 kinship.group=("Mean"),
                 Geno.View.output=F,
                 group.from=30,
                 group.to=100000,
                 group.by=10,
                 SNP.MAF=0,
                 file.output=F)
myGAPIT$GWAS$FDR <- p.adjust(myGAPIT$GWAS$P.value,method = "BH")
write.csv(myGAPIT$GWAS,"Arabidopsis.FT16.TWAS.CMLM.csv",row.names = F)
