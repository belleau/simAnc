
#' @title simulationGenotypeProfileFacets
#'
#' @description
#'
#' @param fileList a \code{list} of \code{GRanges}, the segments from multiple
#' files.
#'
#' @param snv a \code{data.frame}
#'
#' @param genotype a \code{data.frame}
#'
#' @return a \code{}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes
#' @encoding UTF-8
#' @export
simulationGenotypeProfileFacets <- function(PATH_OUT,
                                        PATH_1K,
                                        patientID,
                                        fileBed,
                                        fileFacets,
                                        filePedSel,
                                        fileMatFreq,
                                        chr,
                                        nbSim,
                                        minCov = 10,
                                        minFreq =0.01,
                                        seqError =  0.001/3,
                                        dProp = NA){
    print("Read Files")
    pedSel <- readRDS(filePedSel)
    matFreq <- read.csv2(fileMatFreq, header=FALSE)
    colnames(matFreq) <- c("chr", "pos", "ref", "alt", "AF", "EAS_AF" ,"EUR_AF", "AFR_AF", "AMR_AF", "SAS_AF")

    #Elzar
    #bedCov <- read.table(pipe(paste0("zcat ", PATH_BED, patientID,".bed.gz|grep $'", chr, "\t'")), sep="\t")[,1:3]
    #Wigclust
    bedCov <- read.table(pipe(paste0("zcat ", fileBed,"|grep -P ", chr, "'\t'")), sep="\t")[,1:3]
    colnames(bedCov) <- c("chr", "start", "end")
    facetsRes <- readRDS(fileFacets)
    mysegs <- facetsRes$cncf[which(facetsRes$cncf$chrom == as.numeric(gsub("chr", "",chr))), ]
    mysegs$lap <- rep(NA, nrow(mysegs))
    mysegs$lap[which(!(is.na(mysegs$lcn.em)) || mysegs$tcn.em != 0 )]

    print("End read files")
    infoSNV <- parseSelMinFreq(snv=matFreq, genotype=genotype, minFreq)
    print("computeBedCov")
    infoSNV$snv <- computeBedCov(bedCov, infoSNV$snv)


    snv <- NULL
    for(i in seq_len(nrow(pedSel))){
        print(paste0("Process ",patientID, ".", pedSel$sample.id[i]))
        genotype <- read.csv2(paste0(PATH_1K,
                                     "genotypeSample/",
                                     chr, "/",
                                     pedSel$sample.id[i], ".",
                                     chr, ".vcf.bz2"))
        print("End read genotype")




        resFinal <- simuleBasicGenoChr(genotype,
                                       infoSNV,
                                       bedCov,
                                       mysegs,
                                       nbSim,
                                       minCov,
                                       seqError,
                                       dProp)

        resFile <- paste0(PATH_OUT, patientID, ".", pedSel$sample.id[i], ".", chr, ".rds")
        saveRDS(resFinal, resFile)
    }

}
