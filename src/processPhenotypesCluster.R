### Libraries
library(data.table)
library(stringr)

### Functions
source("src/logger.R")
source("src/globalVariables.R")
source("src/validator.R")

###
# === Note ===
# This is a simplified version of the phenotype processing script.
# The complete script (processPhenotypes.R) should be run locally first
# to prepare phenotypes.
# Because the complete script takes ~12 hours to run and only needs to be run
# once whenever new data is released, we generate a cached file
# which will be used as the input for this simplified script.
# ============
###

# Load processed and filtered UKB phenotype measurements
measurements = fread(PHENOTYPE_MEASUREMENT_PROCESSED_FILE_PATH)
printInfo(sprintf("# of UKB participants: %d", nrow(measurements)))

# Load phenotype selections
selection = fread(PHENOTYPE_SELECTION_FILE_PATH)
validatePhenotypeSelection(selection, PHENOTYPE_SELECTION_FILE_PATH, throwError)

# Remove phenotypes with less than participant cutoff
phenotypes = selection[type == "continuous",
                       CJ(field_id, side = c("top", "bottom"))]
phenotypes = phenotypes[, paste(field_id, side, sep = "_")]
phenotypes = c(
  phenotypes,
  selection[type == "binary", 
            ifelse(subset == "", field_id, paste(field_id, subset, sep = "_"))]
)
colsToRemove = measurements[, lapply(.SD, sum, na.rm = T), .SDcols = phenotypes]
colsToRemove = names(which(unlist(colsToRemove) < PHENOTYPE_MIN_PARTICIPANT_CUTOFF))
measurements[, (colsToRemove) := NULL]
printInfo(sprintf("Removed phenotypes %s because they did not meet min. participant number cutoff (%d)",
                  paste0(colsToRemove, collapse = ", "), PHENOTYPE_MIN_PARTICIPANT_CUTOFF))
phenotypes = setdiff(phenotypes, colsToRemove)

# Randomize measurements to generate permuted data
set.seed(RAND_SEED)
permutedMeasurements = copy(measurements)
# Separately permute male (field 31 = 1) and female (field 31 = 0)
permutedMeasurements[`31` == 1, eid := sample(eid)] # Male
permutedMeasurements[`31` == 0, eid := sample(eid)]
printInfo(sprintf("Permuted measurements with random seed: %d", RAND_SEED))

# Remove basic phenotypes (e.g. sex)
measurements[, (setdiff(colnames(measurements), c("eid", phenotypes))) := NULL]
permutedMeasurements[, (setdiff(colnames(permutedMeasurements), c("eid", phenotypes))) := NULL]

# Save real and permuted measurements
cachedMeasurementFilePath = sprintf(CACHED_MEASUREMENT_FILE_PATH, RAND_SEED)
cachedPermutedFilePath = sprintf(CACHED_PERMUTED_MEASUREMENT_FILE_PATH, RAND_SEED)
fwrite(measurements, cachedMeasurementFilePath)
fwrite(permutedMeasurements, cachedPermutedFilePath)
printInfo(sprintf("Cached real and permuted measurements to %s and %s, respectively", 
                  cachedMeasurementFilePath,
                  cachedPermutedFilePath))
printInfo("Completed processing phenotypes!")
