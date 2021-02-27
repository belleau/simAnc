#' @title parseSelMinFreq
#'
#' @description Generate a subset of snv
#'
#' @param fileList a \code{list} of \code{GRanges}, the segments from multiple
#' files.
#'
#' @param snv a \code{data.frame} with at least the column
#' chr, pos, ref, alt, AF, EAS_AF , EUR_AF, AFR_AF, AMR_AF, SAS_AF
#'
#' @param minFreq \code{numeric} minima frequency in at least one population
#'
#' @return a \code{list} of snv = \code{data.frame} with
#' column : chr, pos, ref, alt, gtype
#' and
#' listSubset = index of the snv keep from the snv input
#'
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschenes and
#' Alexander Krasnitz
#'
#' @export


parseSelMinFreq <- function(snv, genotype, minFreq=0.01){

    #snv <- read.csv2(fileName, header = FALSE)

    #colnames(snv) <- c("chr", "pos", "ref", "alt", "AF", "EAS_AF" ,"EUR_AF", "AFR_AF", "AMR_AF", "SAS_AF")

    listPos <- which(snv$EAS_AF >= minFreq |
                         snv$EUR_AF >= minFreq |
                         snv$AFR_AF >= minFreq |
                         snv$AMR_AF >= minFreq |
                         snv$SAS_AF >= minFreq)


    snv <- snv[listPos, c("chr", "pos", "ref", "alt")]

    return(list(snv=snv, listSubset = listPos))
}

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
                                        chr,
                                        nbSim,
                                        minCov = 10,
                                        seqError =  0.001/3,
                                        dProp = NA){
    print("Read Files")
    pedSel <- readRDS(filePedSel)

    infoSNV <- readRDS(paste0(PATH_OUT, "infoSNV.rds"))

    #Elzar
    #bedCov <- read.table(pipe(paste0("zcat ", PATH_BED, patientID,".bed.gz|grep $'", chr, "\t'")), sep="\t")[,1:3]
    #Wigclust
    bedCov <- read.table(pipe(paste0("zcat ", fileBed,"|grep -P ", chr, "'\t'")), sep="\t")[,1:3]
    colnames(bedCov) <- c("chr", "start", "end")
    facetsRes <- readRDS(fileFacets)
    mysegs <- facetsRes$cncf[which(facetsRes$cncf$chrom == as.numeric(gsub("chr", "",chr))), ]
    mysegs$lap <- rep(NA, nrow(mysegs))
    selTmp <- which(!(is.na(mysegs$lcn.em)) || mysegs$tcn.em != 0 )
    mysegs$lap[selTmp] <- mysegs$lcn.em[selTmp] / mysegs$tcn.em[selTmp]

    print("End read files")

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
