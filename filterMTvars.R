library(MTseeker)
message("This function is now part of MTseeker; use that version instead!") 

#' sanitize PASSing mitochondrial variant calls to a moderate degree
#'
#' @param vars 	an MVRanges or MVRangesList (will be unlisted and relisted)
#' @param fp 	standard issue false positive filters? (TRUE: use Triska & RSRS)
#' @param NuMT	variants with VAF < this number will be presumed NuMTs (0.03)
#' @param covg	minimum median read coverage across chrM to be considered (20) 
#'
#' @return	a filtered set of variants 
#' 
#' @import 	GenomicRanges
#'
#' @examples
#' library(MTseekerData)
#' filterMTvars(RONKSvariants)
#' 
#' @export 
filterMTvars <- function(vars, fp=TRUE, NuMT=0.03, covg=20) {

  if (fp) { 
    data(fpFilter_Triska, package="MTseeker")
    fpFilter <- subset(gaps(fpFilter_Triska), strand == "*")
  } else { 
    fpFilter <- GRanges("chrM", IRanges(start=1, end=16569), strand="*")
  } 

  if (is(vars, "MVRanges")) {
    subset(subsetByOverlaps(vars, fpFilter), VAF >= NuMT & PASS)
  } else if (is(vars, "MVRangesList")) {
    MVRangesList(lapply(vars[genomeCoverage(vars)>=covg], 
                        filterMTvars, fp=fp, NuMT=NuMT))
  } else { 
    stop("This function is only meant for MVRanges and MVRangesList objects.")
  }

}
