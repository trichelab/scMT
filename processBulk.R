
library(MTseeker)

# pull in the reads 
Blast <- getMT("BAMs/SU353-Blast.merged.rCRS.bam") 
LSC <- getMT("BAMs/SU353-LSC.merged.rCRS.bam") 
pHSC <- getMT("BAMs/SU353-pHSC.merged.rCRS.bam")
SU353 <- MAlignmentsList(list(Blast=Blast, LSC=LSC, pHSC=pHSC))
saveRDS(SU353, file="SU353_bulkMT.rds") 

# call variants
SU353vars <- callMT(SU353)
saveRDS(SU353vars, file="SU353_bulkMTvars.rds") 

# filter variants 
source("filterMTvars.R") 
SU353_filtered <- filterMTvars(SU353vars)
saveRDS(SU353_filtered, file="SU353_bulkMTvars_filtered.rds") 

# predict consequences 
predictions <- lapply(lapply(SU353vars_filtered, predictCoding), 
                      subset, consequences != "")
saveRDS(predictions, file="SU353_bulkMTpredictions.rds")

