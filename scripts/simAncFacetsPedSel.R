#!/usr/bin/env Rscript
argv = commandArgs(trailingOnly=TRUE)

library(devtools)
load_all(path = "/mnt/wigclust18/data/unsafe/belleau/beatHaplo/subsetMix/packages/simAnc")



sampleRef <- argv[1]
patientID <- argv[2]
fileFacets <- argv[3]
nbSim <- as.numeric(argv[4])
chr <- argv[5]
seedV <- as.numeric(argv[6])
set.seed(seedV)
nameSel <- gsub(".rds", "", sampleRef)

PATH_1K <- "/mnt/wigclust5/data/unsafe/belleau/process1000G/samples1000gUnrelated/"
PATH_OUT <- paste0("data/simRes.", nameSel, ".", nbSim,"/", chr, "/")
PATH_PED <- "data/ped1000g/"
PATH_BED <- "data/bedAll/"
PATH_FACETS <- "data/facetsFit/"

minCov <- 10
minFreq =0.01
seqError <- 0.001/3
dProp <- NA


fileBed <- paste0(PATH_BED, patientID,".bed.gz")
pedFile <- paste0(PATH_PED, sampleRef)
fileFacets <- paste0(PATH_FACETS, "facets_", fileFacets ,"_NO.rds")
filePedSel <- paste0(PATH_PED, sampleRef)
fileMatFreq <- paste0(PATH_1K, "genotypeSample/", chr, "/matFreqSNV", chr,".txt.bz2")

simulationGenotypeProfileFacets(PATH_OUT,PATH_1K,
                                patientID,
                                fileBed, fileFacets, filePedSel, fileMatFreq,
                                chr, nbSim, minCov, minFreq, seqError, dProp)

