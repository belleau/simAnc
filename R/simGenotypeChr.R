#' @title simuleBasicGenoChr
#'
#' @description
#'
#' @param fileList a \code{list} of \code{GRanges}, the segments from multiple
#' files.
#'
#' @param snv a \code{data.frame}
#'
#' @param genotype \code{data.frame}
#'
#' @return a \code{}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal


simuleBasicGenoChr <- function(genotype,
                       infoSNV,
                       mysegs,
                       nbSim,
                       minCov = 10,
                       seqError =  0.001/3,
                       dProp = NA){


    infoSNV$snv <- setGeno(infoSNV, genotype)

    infoSNV$snv <- parseSegLap(mysegs, infoSNV$snv, seqError, dProp)

    resSim <- simulateAllele(infoSNV$snv, minCov, nbSim)

    blockSeg <- simulateBlockSeg(infoSNV$snv[resSim$listSNV, ], nbSim)

    matGeno <- genoMatrix( infoSNV$snv, resSim, blockSeg)

    resFinal <- list(snv = infoSNV$snv[resSim$listSNV,-1*seq_len(4)],
                     listKeep = resSim$listSNV,
                     matSim = resSim$matSim,
                     blockSeg = blockSeg,
                     matGeno = matGeno)

    return(resFinal)
}


#' @title simuleCatGenoChr
#'
#' @description
#'
#' @param fileList a \code{list} of \code{GRanges}, the segments from multiple
#' files.
#'
#' @param snv a \code{data.frame}
#'
#' @param genotype \code{data.frame}
#'
#' @return a \code{}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal


simuleCatGenoChr <- function(genotype,
                             infoSNV,
                             mysegs,
                             popCur,
                             nbSim,
                             minCov = 10,
                             seqError =  0.001/3,
                             dProp = 0.5){


    infoSNV$snv <- setGeno(infoSNV, genotype)

    infoSNV$snv <- parseCatLap(mysegs, infoSNV$snv, "GName", popCur, seqError, dProp)

    resSim <- simulateAllele(infoSNV$snv, minCov, nbSim)

    blockSeg <- simulateBlockCat(infoSNV$snv[resSim$listSNV, ], nbSim)

    matGeno <- genoMatrix( infoSNV$snv, resSim, blockSeg)

    resFinal <- list(snv = infoSNV$snv[resSim$listSNV,-1*seq_len(4)],
                     listKeep = resSim$listSNV,
                     matSim = resSim$matSim,
                     blockSeg = blockSeg,
                     matGeno = matGeno)

    return(resFinal)
}

#' @title simuleSegBlockGenoChr
#'
#' @description
#'
#' @param fileList a \code{list} of \code{GRanges}, the segments from multiple
#' files.
#'
#' @param snv a \code{data.frame}
#'
#' @param genotype \code{data.frame}
#'
#' @return a \code{}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @keywords internal


simuleSegBlockGenoChr <- function(genotype,
                             infoSNV,
                             mysegs,
                             popCur,
                             nbSim,
                             minCov = 10,
                             seqError =  0.001/3,
                             dProp = 0.5){


    infoSNV$snv <- setGeno(infoSNV, genotype)

#    infoSNV$snv <- parseCatLap(mysegs, infoSNV$snv, "GName", popCur, seqError, dProp)
    infoSNV$snv <- parseSegBlockLap(mysegs, infoSNV$snv, popCur, seqError, dProp)

    resSim <- simulateAllele(infoSNV$snv, minCov, nbSim)

    blockSeg <- simulateBlockCat(infoSNV$snv[resSim$listSNV, ], nbSim)

    matGeno <- genoMatrix( infoSNV$snv, resSim, blockSeg)

    resFinal <- list(snv = infoSNV$snv[resSim$listSNV,-1*seq_len(4)],
                     listKeep = resSim$listSNV,
                     matSim = resSim$matSim,
                     blockSeg = blockSeg,
                     matGeno = matGeno)

    return(resFinal)
}
