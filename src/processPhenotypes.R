### Libraries
library(data.table)
library(stringr)

### Functions
source("src/logger.R")
source("src/globalVariables.R")
source("src/validator.R")

# Load UKB phenotype measurements
measurements = fread(PHENOTYPE_MEASUREMENT_FILE_PATH)
printInfo(sprintf("# of UKB participants: %d", nrow(measurements)))

# For phenotypes with multiple visits, only keep the first visit with data
phenotypes = colnames(measurements)[-1] # Remove the eid column
phenotypes = str_split(phenotypes, fixed("_"), simplify = T)
phenotypes = as.data.table(phenotypes)
phenotypes = phenotypes[, .(phenotype = str_remove(V1, "p"),
                            instance = str_remove(V2, "i"),
                            array = str_remove(V3, "a"))]
res = lapply(unique(phenotypes$phenotype), function(pId) {
 cat("Process phenotype ID", pId, "\n")
  subMeasurements = measurements[, which(phenotypes$phenotype == pId) + 1, with = F]
  
  # If the phenotype contains multiple measurements within one visit
  # combine data from all arrays 
  if (nrow(phenotypes[phenotype == pId & array != ""]) > 0) {
    inds = phenotypes[phenotype == pId, unique(sprintf("p%s_i%s", pId, instance))]
    res = lapply(inds, function(ind) {
      res = apply(
        subMeasurements[, colnames(subMeasurements)[str_starts(colnames(subMeasurements), ind)], with = F], 1,
        function(row) {
          if (all(is.na(row))) return(NA)
          return(paste0(row, collapse = ","))
        })
      res = data.table(res)
      colnames(res) = ind
      return(res)
    })
    subMeasurements = do.call(cbind, res)
  }
  # Select first visit with data
  res = data.table(apply(subMeasurements, 1, function(row) na.omit(row)[1]))
  colnames(res) = pId
  return(res)
})
saveRDS(res, "cache/measurements_merged.rds")
eids = measurements[, .(eid)]
measurements = do.call(cbind, append(list(eids), res))
rm(res)
gc()
fwrite(measurements, "cache/measurements_merged.csv")

# Filter measurements
# 1) To keep only participants with exome sequencing data
eids = readLines(EIDS_FILE_PATH)
measurements = measurements[eid %in% eids]
printInfo(sprintf("# of participants with exome sequencing data: %d",
                  nrow(measurements)))
# 2) Remove withdrawn participants
withdrawnParticipants = readLines(WITHDRAWN_PARTICIPANT_FILE_PATH, warn = F)
measurements = measurements[!eid %in% withdrawnParticipants]
printInfo(sprintf("# of participants after removing withdrawn participants: %d",
                  nrow(measurements)))
# 3) To keep only participants with specified ethnicity
toKeep = startsWith(measurements[, as.character(get(ETHNICITY_FIELD[["id"]]))],
                    ETHNICITY_FIELD[["prefix"]])
measurements = measurements[toKeep]
printInfo(sprintf("# of participants with %s ethnicity: %d",
                  ETHNICITY_FIELD[["description"]], nrow(measurements)))
fwrite(measurements, "cache/measurements_filtered.csv")

# Select only relevant phenotypes
# Load phenotype selections
selection = fread(PHENOTYPE_SELECTION_FILE_PATH)
validatePhenotypeSelection(selection, PHENOTYPE_SELECTION_FILE_PATH, throwError)

# Keep only columns representing relevant phenotypes
colsToKeep = c("eid", unique(selection$field_id))
measurements = measurements[, colsToKeep, with = F]
fwrite(measurements, "cache/measurements_filtered.csv")

# Process phenotype measurements
# 1) Convert quantitative phenotypes to qualitative phenotypes
colsToRemove = apply(selection[type == "continuous"], 1, function(row) {
  # Extract relevant info
  fieldId = str_trim(as.character(row[["field_id"]]))
  subset = row[["subset"]]

  # Subset measurements
  subMeasurements = measurements[, get(fieldId)]
  needsSubsetField = any(str_detect(subMeasurements, fixed(",")))
  
  # Aggregate measurements based on value set in the subset field
  if (!is.na(needsSubsetField) & needsSubsetField) {
    if (subset == "")
      stop(sprintf("missing subset field for field %s that needs to be aggregated", fieldId))
    if (subset == "average") {
      subMeasurements = sapply(str_split(subMeasurements, fixed(",")), function(m) {
        m = suppressWarnings(na.omit(as.numeric(m)))
        if (length(m) < 1) return(NA)
        mean(m)
      })
    }
  }
  
  # Split into two binary phenotypes
  cutoffs = quantile(subMeasurements, probs = c(0.1, 0.9), na.rm = T, names = F)
  newFieldIds = paste(fieldId, c("top", "bottom"), sep = "_")
  # 1) Phenotype 1 (_top): convert top 10% to cases, middle 80% as control
  #    and discard the remaining bottom 10%
  measurements[, (newFieldIds[1]) := 
                 ifelse(subMeasurements < cutoffs[1], NA, 
                        as.integer(subMeasurements > cutoffs[2]))]
  # 2) Phenotype 2 (_bottom): convert bottom 10% to cases, middle 80% as control
  #    and discard the top 10%
  measurements[, (newFieldIds[2]) := 
                 ifelse(subMeasurements > cutoffs[2], NA, 
                        as.integer(subMeasurements < cutoffs[1]))]
  
  return(fieldId)
})
colsToRemove = unique(unlist(colsToRemove))
measurements[, (colsToRemove) := NULL]
# 2) Expand binary phenotypes
colsToRemove = apply(selection[type == "binary"], 1, function(row) {
  # Extract relevant info
  fieldId = str_trim(as.character(row[["field_id"]]))
  subset = str_trim(as.character(row[["subset"]]))
  
  # If has subset, create a new column
  if (subset != "") {
    newFieldId = paste0(fieldId, "_", subset)
    measurements[, (newFieldId) := as.integer(str_detect(get(fieldId), subset))]
  } else {
    measurements[, (fieldId) := as.integer(!is.na(get(fieldId)))]
  }
  
  return(fieldId)
})
colsToRemove = unique(unlist(colsToRemove))
measurements[, (colsToRemove) := NULL]
fwrite(measurements, PHENOTYPE_MEASUREMENT_PROCESSED_FILE_PATH)

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
