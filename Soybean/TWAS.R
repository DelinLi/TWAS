setwd("TWAS/Soybean/")#Please adjust the path
##From Delin Li delin.bio@gmail.com Publication: https://academic.oup.com/plphys/article/186/4/1800/6212071
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

##The numeric format genotype format required by GAPIT
myGM<-read.delim("Geno.Infor.txt.gz")
myGD<-read.delim("Soybean.ALL.Quan5.Geno.txt.gz")

Pheno <-read.delim("PubescenceColor.qualitative.txt") #ignore the third column
 
myY<-Pheno[ Pheno$ID %in% myGD$taxa ,1:2] # The one having genotype and phenotype
myY<-myY[order(as.character(myY[,1])),]

myGD.tem<-myGD[myGD$taxa %in% myY[,1],]
myGD.tem <- myGD.tem[,apply(myGD.tem,2,function(x){return(sum(is.na(x))==0)})] #remove the gene with NA 
myGD.tem<-myGD.tem[order(as.character(myGD.tem$taxa)),]
myGM.tem <- myGM[myGM$Name %in%   colnames(myGD.tem),]

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
                 file.output=F
)
myGAPIT$GWAS$FDR <- p.adjust(myGAPIT$GWAS$P.value,method = "BH")
write.csv(myGAPIT$GWAS,"Soybean.PubescenceColor.TWAS.CMLM.csv",row.names = F)

png(paste0("Soybean.TWAS.FDR.png"), width = 20,height=8,units = "in",res = 200)
layout(matrix(1:1, 1, 1, byrow = TRUE))
par(mar = (c(4.2,5,2,2)+ 0.5), mgp=c(3,1.8,0))    ##
par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box
FDR.Manhattan.only(myGAPIT$GWAS[,c(2,3,10)] ,name.of.trait= "TWAS",plot.type = "Genomewise",cutOff= 0.05,
                   highliht.sig=F,cex.lab=1.8)
dev.off()
