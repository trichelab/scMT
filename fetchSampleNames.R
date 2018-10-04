library(SRAdb)
sra_dbname <- 'SRAmetadb.sqlite'
sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
samples <- read.csv("combined_scATAC_runs.csv", row=2) # many missing!
SRRs <- paste0("('", paste(rownames(samples), collapse="','"), "')")
covs <- dbGetQuery(sra_con,
                    paste0("SELECT run.run_accession, ",
                           "       sample.sample_attribute ",
			   "  FROM sample, run, experiment ",
                           " WHERE run.run_accession IN ", SRRs,
                           "   AND experiment.experiment_accession = ",
                           "       run.experiment_accession",
                           "   AND sample.sample_accession = ",
                           "       experiment.sample_accession"))
covs$WHAT <- sub("acute myeloid leukemia, ", "", covs$sample_attribute)
covs$sample_attribute <- NULL
covs$source <- sub("source_name: ", "", 
                   sapply(strsplit(covs$WHAT, " \\|\\| "), `[`, 1))
covs$donor <- sub("don(e|o)r id: ", "", 
                  sapply(strsplit(covs$WHAT, " \\|\\| "), `[`, 3))
covs$WHAT <- NULL
covs$BAM <- paste0(covs$run_accession, ".merged.rCRS.bam")
names(covs)[1] <- "SRR"
write.table(covs, file="scATAC.mappings.txt")

