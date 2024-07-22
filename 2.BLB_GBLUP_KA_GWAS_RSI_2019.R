#!/usr/bin/env Rscript
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=intel16|intel18
#SBATCH --array=1-1000

rm(list=ls())
setwd('/mnt/gs21/scratch/izquier7/r_bbl/year_2019')

library(BGLR)
library(SFSI)
library(GAPIT)
library(dplyr)
library(tidyverse)

load('../GB_BLB.Rdata')

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

JOBS <- expand.grid(rep=1:100, trait=1:10)

trait <- as.vector(JOBS[job,"trait"])
rep <- as.vector(JOBS[job,"rep"])

M <- as.matrix(GB_BLB$geno) # markers
Y <- GB_BLB$pheno[,-1] # phenotype
myGD <- GB_BLB$geno  # markers - GWAS
myGD$taxa = rownames(myGD)
myGD = myGD[,c(2316,1:2315)]
rownames(myGD) = NULL

myGM <- GB_BLB$positions # positions - markers
colnames(myGM)[1] = 'Name'

Y <- Y[, seq(0,24,2)]
Y <- Y[,-c(3:4)]

G <- tcrossprod(scale(M))/ncol(M) # genomic relationship
D<-as.matrix(dist(M,method="euclidean"))^2 #euclidian distance
D<-D/mean(D)

X <- as.matrix(GB_BLB$nirs19)

index <- 1:272
M1 <- M[index,]
GD <- myGD[index,]
G1 <- G[index,index]
D1 <- D[index,index]
Y <- Y[index,]
X <- scale(X[index,])


index <- apply(X,1,function(x)all(!is.na(x)))
M1 <- M1[index,]
GD <- GD[index,]
G1 <- G1[index,index]
Y <- Y[index,]
X <- X[index,]


trait_name <- colnames(Y)[trait]
y <- Y[,trait]

index <- which(!is.na(y))
M1 <- M1[index,]
GD <- GD[index,]
G1 <- G1[index,index]
D1 <- D1[index,index]
y <- y[index]
X <- X[index,]

#######Sparce selection index
n <- length(y)
cv = 10

accuracy_fold_trn = c()
accuracy_fold_tst =  c()
imax_fold_max =  c()
lambda <- c()

set.seed(rep)

### Training and Testing
tst <- sample(1:n, 0.3*n) 
trn <- seq_along(1:n)[-tst] 

trn_fold <- rep(seq(1:cv),
                ceiling((length(trn))/cv))[1:(length(trn))]
trn_fold <- sample(trn_fold)

#cross validation   
for(k in 1:cv){
  trt_ind <- which(trn_fold != k)
  tst_ind <- which(trn_fold == k)
  
  fm_fold <- getGenCov(y[trt_ind],
                       X[trt_ind,],
                       K=G1[trt_ind,trt_ind],
                       mc.cores=10)
  
  covU_fold <- fm_fold$covU 
  Xtrn_fold <- scale(X[trt_ind,])
  VAR_fold <- var(Xtrn_fold)
  
  fm2_fold <- solveEN(VAR_fold, 
                      covU_fold, 
                      X=Xtrn_fold)
  
  lambda_fold <- fm2_fold$lambda
  yHat_fold <- fm2_fold$yHat
  
  accuracy_fold <- cor(y[trt_ind],yHat_fold)[1,]
  
  imax_fold <- which.max(abs(accuracy_fold))
  
  B_fold <- coef(fm2_fold)[,imax_fold]
  
  SI_fold <- scale(as.vector(X%*%B_fold))
  accuracy_fold_trn[k] = cor(SI_fold[trt_ind],y[trt_ind])
  accuracy_fold_tst[k] = cor(SI_fold[tst_ind],y[tst_ind])
  imax_fold_max[k] = imax_fold
  lambda[k] <- lambda_fold[imax_fold]
}

lambda0 <- mean(lambda)

fm <- getGenCov(y[trn],
                X[trn,],
                K=G1[trn,trn] )

covU <- fm$covU 

Xtrn <- scale(X[trn,])
VAR <- var(Xtrn)

SVD <- eigen(VAR)
index <- SVD$values>0
invVAR <- SVD$vectors[,index]%*%diag(SVD$values[index])%*%t(SVD$vectors[,index])
B0 <- invVAR%*%covU
SI <- scale(as.vector(X%*%B0))

fm2 <- solveEN(VAR, 
               covU, 
               X=Xtrn, lambda=lambda0)

B1 <- coef(fm2)

SSI <- scale(as.vector(X%*%B1))

accuracy <- c(cor(SI[tst],y[tst]),
              cor(SSI[tst],y[tst]))

names(accuracy) <- c("SI_tst","SSI_tst")

##### Univariado GBLUP
yNA <- y
yNA[tst] <- NA

ETA <- list(list(K=G1,model="RKHS"))

nIter <- 12000
burnIn <- 2000

fm_uni <- BGLR(y=yNA, ETA=ETA,
               nIter=nIter, burnIn=burnIn)

accuracy_GBLUP = cor(y[tst], fm_uni$yHat[tst], use = 'pairwise.complete.obs')

##### Univariado -  Gaussian average kernel
h<- c(.02,1,5)
KList<-list()

for(i in 1:length(h)){
  KList[[i]]<-list(K=exp(-h[i]*D1),model='RKHS')
}

fmKA<-BGLR(y=yNA,ETA=KList,
           nIter=nIter,burnIn=burnIn,
           saveAt="KA_")

accuracy_KA = cor(y[tst],fmKA$yHat[tst],
                         use = 'pairwise.complete.obs')

### Multi-trait GBLUP
Y_SSI18 <- cbind(yNA,SSI)


fm_muv <- Multitrait(y=Y_SSI18, ETA=ETA, 
                     resCov=list(type="DIAG"), 
                     nIter=nIter, burnIn=burnIn)


accuracy_GBLUP_MT = cor(y[tst],fm_muv$ETAHat[tst])

### Multi-trait KA
fm_muv_ka <- Multitrait(y=Y_SSI18, ETA=KList, 
                     resCov=list(type="DIAG"), 
                     nIter=nIter, burnIn=burnIn)


accuracy_KA_MT = cor(y[tst],fm_muv_ka$ETAHat[tst])

### GWAS

pheno_gwas = t((rbind(rownames(G1),yNA)))
pheno_gwas = as.data.frame(pheno_gwas)
colnames(pheno_gwas) = c('taxa','y')
pheno_gwas$y = as.numeric(pheno_gwas$y)

myGAPIT <- GAPIT(
  Y = pheno_gwas, 
  GD = GD, 
  GM = myGM,
  PCA.total=3,
  model= c("FarmCPU"), file.output = F
)

bonferroni = 0.05/dim(GD)[2]
GWAS_res =filter(myGAPIT$GWAS, P.value < bonferroni)

indexQTL <- which(colnames(M) %in% GWAS_res$SNP) ## QTl for yield

accuracy_gwas = c()
accuracy_KA_gwas = c()
accuracy_KA_gwas = c()
accuracy_KA_gwas_mt = c()

if (length(indexQTL) >= 1) {
  
  # GBLUP + GWAS
  M_new <- M[,-indexQTL]
  QTL <- as.matrix(M[,indexQTL])
  G_qtl <- tcrossprod(scale(M_new))/ncol(M_new) # genomic relationship
  D_qtl <-as.matrix(dist(M_new,method="euclidean"))^2 #euclidian distance
  D_qtl<-D_qtl/mean(D_qtl)

  index <- 1:272
  M_new <- M_new[index,]
  QTL <- as.matrix(QTL[index,])
  G_qtl <- G_qtl[index,index]
  D_qtl <- D_qtl[index,index]
  
  index <- apply(X,1,function(x)all(!is.na(x)))
  M_new <- M_new[index,]
  QTL <- as.matrix(QTL[index,])
  G_qtl <- G_qtl[index,index]
  D_qtl <- D_qtl[index,index]
  
  index <- which(!is.na(y))
  M_new <- M_new[index,]
  QTL <- as.matrix(QTL[index,])
  G_qtl <- G_qtl[index,index]
  D_qtl <- D_qtl[index,index]
  
  ETA <- list(list(X=QTL,model="FIXED"),
              list(K=G_qtl,model="RKHS"))
  
  fm_gwas <- BGLR(y=yNA, ETA=ETA, nIter=nIter, 
                  burnIn=burnIn)
  
  accuracy_gwas = cor(y[tst],fm_gwas$yHat[tst],
                            use = 'pairwise.complete.obs')
  
  # GBLUP + GWAS + MT

  fm_gblup_gwas_mt <- Multitrait(y=Y_SSI18, ETA=ETA, 
                       resCov=list(type="DIAG"), 
                       nIter=nIter, burnIn=burnIn)
  
  
  accuracy_gblup_mt_gwas = cor(y[tst],fm_gblup_gwas_mt$ETAHat[tst],
                               use = 'pairwise.complete.obs')
  
  # KA + GWAS
  KList<-list()
  
  KList = list(list(X=QTL,model="FIXED"))
  
  for(i in 1:length(h)){
   KList[[i+1]]<- list(K=exp(-h[i]*D_qtl),model='RKHS')
  }
  
  fmKA_gwas<-BGLR(y=yNA,ETA=KList,
             nIter=nIter,burnIn=burnIn,
             saveAt="KA_")

  accuracy_KA_gwas = cor(y[tst],
                           fmKA_gwas$yHat[tst],
                      use = 'pairwise.complete.obs')
  
  #### KA + GWAS + MT
  fm_muv_ka_gwas_mt <- Multitrait(y=Y_SSI18, ETA=KList, 
                          resCov=list(type="DIAG"), 
                          nIter=nIter, burnIn=burnIn)
  
  
  accuracy_KA_gwas_mt = cor(y[tst],fm_muv_ka_gwas_mt$ETAHat[tst])
} else {
  accuracy_gwas = 'NA'
  accuracy_KA_gwas = 'NA'
  accuracy_gblup_mt_gwas = 'NA'
  accuracy_KA_gwas_mt = 'NA'
}

### results
gp_results <- c(accuracy, 
                accuracy_GBLUP, accuracy_KA, 
                accuracy_GBLUP_MT, accuracy_KA_MT,
                accuracy_gwas, accuracy_gblup_mt_gwas,
                accuracy_KA_gwas, accuracy_KA_gwas_mt)

names(gp_results) <- c('SI','SSI','GBLUP','KA','GBLUP_MT','KA_MT', 
                       'GBLUP_GWAS', 'GBLUP_GWAS_MT','KA_GWAS',
                       'KA_GWAS_MT')

save(gp_results,
     file=paste0("results_",trait_name,"_rep_",rep,".RData"))

save(GWAS_res,
     file=paste0("GWAS_results_",trait_name,"_rep_",rep,".RData"))
