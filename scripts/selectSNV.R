#!/usr/bin/env Rscript
argv = commandArgs(trailingOnly=TRUE)

library(devtools)
load_all(path = "/mnt/wigclust18/data/unsafe/belleau/beatHaplo/subsetMix/packages/simAnc")

sampleRef <- argv[1]
nbSim <- as.numeric(argv[2])
chr <- argv[3]
minFreq <- argv[4]

nameSel <- gsub(".rds", "", sampleRef)


PATH_1K <- "/mnt/wigclust5/data/unsafe/belleau/process1000G/samples1000gUnrelated/"
PATH_OUT <- paste0("data/simRes.", nameSel, ".", nbSim,"/", chr, "/")

fileMatFreq <- paste0(PATH_1K, "genotypeSample/", chr, "/matFreqSNV", chr,".txt.bz2")
matFreq <- read.csv2(fileMatFreq, header=FALSE)
colnames(matFreq) <- c("chr", "pos", "ref", "alt", "AF", "EAS_AF" ,"EUR_AF", "AFR_AF", "AMR_AF", "SAS_AF")

infoSNV <- parseSelMinFreq(snv=matFreq, genotype=genotype, minFreq)

saveRDS(infoSNV, paste0(PATH_OUT, "infoSNV.f", minFreq,".rds"))

q()
