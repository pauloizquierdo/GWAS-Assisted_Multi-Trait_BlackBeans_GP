#!/usr/bin/env Rscript
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=intel16|intel18
#SBATCH --array=1-1000

rm(list=ls())
setwd('/mnt/gs21/scratch/izquier7/r_bbl/year_2019')

library(BGLR)
library(SFSI)

load('../GB_BLB.Rdata')

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

JOBS <- expand.grid(rep=1:100, trait=1:10)

trait <- as.vector(JOBS[job,"trait"])
rep <- as.vector(JOBS[job,"rep"])

M <- as.matrix(GB_BLB$geno) # markers
Y <- GB_BLB$pheno[,-1] # phenotype

Y19 <- Y[, seq(0,24,2)]
Y19 <- Y19[,-c(3:4)]

G <- tcrossprod(scale(M))/ncol(M) # genomic relationship
D<-as.matrix(dist(M,method="euclidean"))^2 #euclidian distance
D<-D/mean(D)

index <- 1:272
M1 <- M[index,]
G1 <- G[index,index]
D1 <- D[index,index]
Y19 <- Y19[index,]

trait_name <- colnames(Y19)[trait]
y19 <- Y19[,trait]

index <- which(!is.na(y19))
M1 <- M1[index,]
G1 <- G1[index,index]
D1 <- D1[index,index]
y19 <- y19[index]

n <- length(y19)
tst <- sample(1:n, 0.3*n) 
trn <- seq_along(1:n)[-tst] 

##### GBLUP
yNA <- y19
yNA[tst] <- NA

nIter <- 12000
burnIn <- 2000

ETA <- list(list(K=G1,model="RKHS"))

fm_uni <- BGLR(y=yNA, ETA=ETA,
               nIter=nIter, burnIn=burnIn)

accuracy_gs_tst_19 = cor(y19[tst],
                         fm_uni$yHat[tst],
                         use = 'pairwise.complete.obs')

##### RKHS - Gaussian kernel - grid
h<- c(.02,1,5)
accuracy_Kernel = c()

for(i in 1:length(h)){
  # COMPUTES THE KERNEL
  K<-exp(-h[i]*D1)
  # FITS THE MODEL
  ETA<-list(list(K=K,model='RKHS'))

  fmKernel<-BGLR(y=yNA,ETA=ETA,
           nIter=nIter,burnIn=burnIn)

  accuracy_Kernel[i] = cor(y19[tst],
                        fmKernel$yHat[tst],
                        use = 'pairwise.complete.obs')
}

##### RKHS - Gaussian kernel - Average kernel
KList<-list()

for(i in 1:length(h)){
  KList[[i]]<-list(K=exp(-h[i]*D1),model='RKHS')
}

fmKA<-BGLR(y=yNA,ETA=KList,
           nIter=nIter,burnIn=burnIn,
           saveAt="KA_")

accuracy_KA_tst_19 = cor(y19[tst],
                         fmKA$yHat[tst],
                         use = 'pairwise.complete.obs')

### results
gp_results <- c(accuracy_gs_tst_19, accuracy_Kernel,
                accuracy_KA_tst_19)

names(gp_results) <- c('GBLUP', 'K1', 'K2', 'K3',
                       'KA')

save(gp_results,
     file=paste0("results_",trait_name,"_rep_",
                 rep,".RData"))
