### Libraries
library(data.table)
library(stringr)

### Functions
source("src/logger.R")
source("src/globalVariables.R")
source("src/runBurdenTest.R")

# Parse arguments
args = commandArgs(trailingOnly = T)
if (length(args) != 2) stop("Invalid parameters")
varPredictionFilePath = args[1]
runFolder = args[2]

# Load measurement datasets
cachedMeasurementFilePath = sprintf(CACHED_MEASUREMENT_FILE_PATH, RAND_SEED)
cachedPermutedFilePath = sprintf(CACHED_PERMUTED_MEASUREMENT_FILE_PATH, RAND_SEED)
measurements = fread(cachedMeasurementFilePath)
permutedMeasurements = fread(cachedPermutedFilePath)

printInfo(sprintf("Started processing %s and its associated files", 
                  varPredictionFilePath))

# Load variant effect predictions for UKB variants
colsToKeep = unname(c(VARIANT_PREDICTION_FILE_COLUMNS_TO_KEEP, VARIANT_PREDICTORS))
varPredictions = fread(varPredictionFilePath, select = colsToKeep)
printInfo(sprintf("# of UKB variant predictions: %d", nrow(varPredictions)))

# Filter predictions
# 1) Remove common variants
# Fill missing AF with 0
varPredictions[is.na(gnomAD_AF), gnomAD_AF := 0]
varPredictions = varPredictions[gnomAD_AF < VARIANT_PREDICTION_FILTER_AF]
printInfo(sprintf("# of UKB variant predictions with variants passed AF threshold (<%.3f): %d",
                  VARIANT_PREDICTION_FILTER_AF, nrow(varPredictions)))
# 2) Keep only variants on canonical transcripts
# Split comma delimited columns
varPredictions = varPredictions[Consequence == VARIANT_PREDICTION_FILTER_SYN | Ensembl_transcriptid != ""]
varPredictions[, (VARIANT_PREDICTION_FILE_COLUMNS_COMMA_SPLITTED) 
               := lapply(.SD, strsplit, split=","),
               .SDcols=VARIANT_PREDICTION_FILE_COLUMNS_COMMA_SPLITTED]
varPredictions = varPredictions[, lapply(.SD, unlist), by = 1:nrow(varPredictions)]
varPredictions$nrow = NULL
# Keep variants on canonical transcript
varPredictions = varPredictions[Consequence == VARIANT_PREDICTION_FILTER_SYN | Feature == Ensembl_transcriptid]
printInfo(sprintf("# of UKB variant predictions on canonical transcripts: %d",
                  nrow(varPredictions)))

# Expand amino acid into reference and alternated
varPredictions[, c("aa_ref", "aa_alt") := tstrsplit(Amino_acids, "/", fixed = T)]

# Load VARITY scores
chr = str_extract(varPredictionFilePath, "(?<=_c).+?(?=_)")
varityPath = sprintf(CACHED_VARITY_FILE_FORMAT, chr)
varityScores = fread(
  varityPath, select = VARIANT_PREDICTOR_FILE_VARITY_COLUMNS_TO_KEEP
)

# Add VARITY scores
varPredictions$Protein_position = suppressWarnings(as.integer(varPredictions$Protein_position))
varPredictions = varPredictions[!is.na(Protein_position)]
varPredictions = merge(varPredictions, varityScores,
                       by.x = c("Uniprot_acc", "Protein_position", "aa_ref", "aa_alt"),
                       by.y = c("p_vid", "aa_pos", "aa_ref", "aa_alt"), all.x = T)
rm(varityScores)
invisible(gc(full = T))
printInfo("Added VARITY scores")

# Load filtered mutations in UKB
filteredVariantFilePath = str_replace(
  varPredictionFilePath, VARIANT_PREDICTION_FILE_SUFFIX,
  FILTERED_VARIANT_FILE_SUFFIX
)
filteredVariants = fread(filteredVariantFilePath)
filteredVariants[, eid := as.numeric(eid)]
filteredVariants = filteredVariants[!is.na(eid)]
printInfo(sprintf("Loaded %d filtered UKB vairants", nrow(filteredVariants)))

# Here we make sure the allele freq is below zero in both gnomad and UKB AF
# This rule is added on 2022-01-04 while investigating strange burden test outputs
# We found out that we cannot just rely on gnomad AF because there're cases 
# where gnomad v2 (which is the pipeline is based on) has no variant but UKB AF is > 0.5.
# For example, variant X-48058816-A-G has a AF = 0.47 in gnomad v3
# but the variant is not in gnomad v2.
filteredVariants = filteredVariants[AF < VARIANT_PREDICTION_FILTER_AF]
printInfo(sprintf("# of filtered UKB variants past allele frequency filter: %d",
                  nrow(filteredVariants)))

# Assign variant type
filteredVariants[, type := "SNP"]
filteredVariants[nchar(ref) > 1 | nchar(alt) > 1, type := "INDEL"]

# Filter variants
# 1) Remove variants with shallow read depth
filteredVariants = filteredVariants[(DP > 7 & type == "SNP") |
                                      (DP > 10 & type == "INDEL")]
printInfo(sprintf("# of filtered UKB variants past read depth filter: %d",
                  nrow(filteredVariants)))

# Set up output file path
burdenTestOutputFilePath = str_replace(
  varPredictionFilePath, VARIANT_PREDICTION_FILE_SUFFIX,
  BURDEN_TEST_PVAL_FILE_SUFFIX
)
burdenTestOutputFilePath = str_replace(
  burdenTestOutputFilePath, fixed("ukb_variants"), runFolder
)
# If file exists, remove file
if (file.exists(burdenTestOutputFilePath)) file.remove(burdenTestOutputFilePath)

# Run burden test
phenotypes = setdiff(colnames(measurements), "eid")
invisible(runBurdenTest(varPredictions, filteredVariants, phenotypes, measurements,
              permutedMeasurements, burdenTestOutputFilePath))

printInfo(sprintf("Finished processing %s and its associated files", 
                  varPredictionFilePath))
