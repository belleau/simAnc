
#' @title setGeno
#'
#' @description Generate a subset of snv with genotype
#'
#' @param fileList a \code{list} of \code{GRanges}, the segments from multiple
#' files.
#'
#' @param infoSNV a \code{list} of snv = \code{data.frame} with
#' column : chr, pos, ref, alt, gtype
#' and
#' listSubset = index of the snv keep from the snv input
#'
#' @param genotype \code{data.frame}
#' one column format 0|0, 0|1, 1|0 or 1|1
#' one row per snv same row then snv
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
#' @keywords internal

setGeno <- function(infoSNV, genotype){


    infoSNV$snv$gtype <- genotype[infoSNV$listSubset,]

    return(infoSNV$snv)
}


#' @title parseSegLap
#'
#' @description Add the column lap and seg to the snv dataframe from parseStatMinFreq
#' Here lap is lower allele proportion for the snv in the sample
#' (not minimal allele freq in the population)
#'
#' @param mysegs a \code{list} of \code{data.frame}, from one chr only
#' with at least the column start, end, lap. If the segments come from
#' Facets lap = lcn.em / tcn.em
#'
#' @param snv a \code{data.frame}
#'
#' @param seqError \code{numeric} lap for the other alleles
#' for the homozygote Default: 0.001/3
#'
#' @param dProp \code{numeric} value of lap out of the segments
#' Default: \code{NA}
#'
#' @return a \code{data.frame} the snv from the input plus column
#' lap and seg (the index of the segments)
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschenes and
#' Alexander Krasnitz
#'
#' @keywords internal


parseSegLap <- function(mysegs, snv, seqError = 0.001/3, dProp=NA){

    #for each variant position in the reference, find CN and left allele
    #fraction
    z<-cbind(c(mysegs[,"start"],mysegs[,"end"],snv[,"pos"]),
             c(seq_len(nrow(mysegs)),
               -(seq_len(nrow(mysegs))),
               rep(0,nrow(snv))))

    z<-z[order(z[,1]),,drop=F]
    mysegis<-cumsum(z[,2])[z[,2]==0]
    lap <- mysegs[,"lap"]#mysegs[,"lcn.em"]/mysegs[,"tcn.em"]
    snv[,"lap"]<-rep(dProp,nrow(snv))
    snv[,"seg"]<-rep(NA,nrow(snv))
    snv[mysegis!=0,"lap"] <- lap[mysegis[mysegis!=0]]
    snv[mysegis!=0,"seg"] <- mysegis[mysegis!=0]
    snv[which(snv$gtype %in% c("0|0", "1|1")), "lap"] <- seqError
    return(snv)

}



#' @title computeBedCov
#'
#' @description Coverage at the snv position
#'
#' @param myreads a \code{data.frame} of the start and end of reads for only
#' one chr. Column at least  start and end
#'
#' @param snv a \code{data.frame} of snv in only one chr with at
#' least the column: pos
#'
#' @return a \code{data.frame} of snv in only one chr with a new
#' column readcount
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschenes and
#' Alexander Krasnitz
#'
#' @keywords internal


computeBedCov <- function(myreads, snv){

    z<-cbind(c(myreads[,"start"],myreads[,"end"],snv[,"pos"]),
             c(rep(1,nrow(myreads)),
               rep(-1,nrow(myreads)),
               rep(0,nrow(snv))))
    z<-z[order(z[,1]),,drop=F]

    snv[,"readcount"]<-cumsum(z[,2])[z[,2]==0]
    return(snv)

}


#' @title simulateAllele
#'
#' @description Compute the readcount for lower allele proportion
#' don't need to know which allele it is at this step
#'
#' @param snvCovLap a \code{data.frame} with at least the column:
#' lap and readcount
#'
#' @param minCov a \code{numeric} the min value of readcount to keep the snv
#'
#' @param nbSim \code{integer} number of simulation
#'
#' @return a \code{list}
#' listSNV: the index of the snv keep
#  matSim: the readcount of the of the allele dim = length(listSNV) x nbSim
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschenes and
#' Alexander Krasnitz
#'
#' @keywords internal


simulateAllele <- function(snvCovLap, minCov, nbSim){

    listKeep <- which(!is.na(snvCovLap[,"lap"]) & snvCovLap$readcount >= minCov )
    snvCovLapNoNa<-snvCovLap[listKeep,,drop=F]
    snvCovLapNoNa<-snvCovLapNoNa[order(snvCovLapNoNa[,"lap"]),,drop=F]

    snvCovLapNoNa<-snvCovLapNoNa[order(snvCovLapNoNa[,"readcount"]),,drop=F]

    mcount<-as.matrix(table(snvCovLapNoNa[,c("readcount","lap")]))
    mcount<-cbind(as.numeric(dimnames(mcount)[[1]]),mcount)
    vf<-as.numeric(dimnames(mcount)[[2]])[-1]

    makelcounts<-function(v, vf, nbSim)Map(rbinom,nbSim * v[-1],rep(v[1],length(v)-1),vf)
    resSim <- t(matrix(unlist(apply(mcount,
                                    1,
                                    makelcounts,
                                    vf=vf,
                                    nbSim=nbSim)),
                       nr=nbSim,nc=length(listKeep)))

    colnames(resSim) <- paste0("s", seq_len(nbSim))

    return(list(listSNV=listKeep, matSim = resSim[order(snvCovLapNoNa$pos), , drop=FALSE]) )
}



#' @title simulateBlockSeg
#'
#' @description Select for each SNV which haplotype (left or right, father or mother) are lap
#'
#' @param snv a \code{data.frame}
#'
#' @param nbSim \code{integer} number of simulation
#'
#' @return a \code{matrix} of 0,1 coding a selection of the left or right (father or mother haplotype)
#' if seg is NA always 0
#' LAF == 1 left is lap LAF == 0 left is 1-lap because of
#' the function genoMatrix(...)
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschenes and
#' Alexander Krasnitz
#'
#' @keywords internal


simulateBlockSeg <- function(snv, nbSim){

    listPos <- which(!is.na(snv$seg))

    snvNoNa <- snv[listPos,]


    # the segment must be a sequence 1:nbSeg
    z <- cumsum(snvNoNa$seg - c(0,snvNoNa$seg[-nrow(snvNoNa)]))
    blockSeg <- matrix(sample(x = c(0,1), nbSim *(length(unique(snvNoNa$seg))), replace=TRUE),nc=nbSim)
    LAFparent <- matrix(0,nr=nrow(snv), nc=nbSim)
    LAFparent[listPos,] <- blockSeg[z,]
    return(LAFparent)
}


#' @title genoMatrix
#'
#' @description Genere twa matrix matRef and matAlt with the
#' readcount of the reference allele and readcount
#' of the allele alt respectivily
#'
#' @param snv a \code{data.frame}
#'
#' @param resSim \code{list} output of simulateAllele(...)
#'
#' @return a \code{list} of two  \code{matrix}
#' matRef: \code{matrix} of the readcount of the reference allele
#' matAlt \code{matrix} of the readcount of the alt allele
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschenes and
#' Alexander Krasnitz
#'
#' @keywords internal


genoMatrix <- function(snv, resSim, LAFparent){
    # genoL genotype left (0 for ref and 1 for alt)
    genoL <- matrix(rep(as.numeric(substr(snv$gtype[resSim$listSNV], 1,1)),
                        ncol(resSim$matSim)),
                    nc=ncol(resSim$matSim))

    # genoL genotype right (0 for ref and 1 for alt)
    genoR <- matrix(rep(as.numeric(substr(snv$gtype[resSim$listSNV], 3,3)),
                        ncol(resSim$matSim)),
                    nc=ncol(resSim$matSim))

    # matrix of the number of read
    readCount <- matrix(rep(snv$readcount[resSim$listSNV],
                            ncol(resSim$matSim)),
                        nc=ncol(resSim$matSim))

    # if heterozygote LAF == 1 left is lap LAF == 0 left is 1-lap
    # if homozygote the allele present is max
    tmp <- (genoL + ((genoL + genoR) * LAFparent ))%%2

    matAlt <- tmp*readCount + ((-1)^tmp) * resSim$matSim
    matRef <- readCount - matAlt

    return(list(matRef=matRef, matAlt=matAlt))
}


