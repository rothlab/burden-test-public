### Libraries
library(data.table)

### Functions
source("src/logger.R")
source("src/globalVariables.R")

# Parse arguments
args = commandArgs(trailingOnly = T)
if (length(args) != 2) stop("Invalid parameters")
runFolder = args[1]

# Load p value files
pValPaths = list.files(runFolder, BURDEN_TEST_PVAL_FILE_SUFFIX, full.names = T)

# Get row names
colNames = c("mode_inheritance", "type_variants", "qualifying_variants",
             "is_real_data", "p_val")

# Set up conditions
predictorCols = unname(c(VARIANT_PREDICTORS, VARIANT_PREDICTOR_VARITY))
conditions = expand.grid(
  mode_inheritance = MODE_INHERITANCE,
  variant_type = VARIANT_TYPE,
  variant_effect = c("raw_count", predictorCols)
)

# Submit job for each condition
invisible(apply(conditions, 1, function(row) {
  varCol = as.character(row[["variant_effect"]])
  varType = as.character(row[["variant_type"]])
  modeInheritance = as.character(row[["mode_inheritance"]])
  printInfo(paste("Process condition:", varCol, varType, modeInheritance))
  
  # Submit jobs for each trait prevalence and burden combination
  combs = expand.grid(trait_bin = 1:10, burden_bin = 1:10)
  invisible(mapply(function(trait, burden) {
    pFilePath = sprintf(
      "%s/%s_%s_%s_%d-%d_p_cache.txt", runFolder, varCol, varType, 
      modeInheritance, trait, burden)
    permutedPFilePath = sprintf(
      "%s/%s_%s_%s_%d-%d_permuted_p_cache.txt", runFolder, varCol, varType, 
      modeInheritance, trait, burden)
    jobName = paste(varCol, varType, modeInheritance, sep = "_")
    jobScript = sprintf("Rscript src/mergeParseOutputResultsWithSplittedNull.R %s %s", pFilePath, permutedPFilePath)
    system2("submitjob.sh",
            args = c("-n", jobName, "-t", "24:00:00", "-m", "20G", "-c", "16",
                     "-e", sprintf("%s/%s.log", runFolder, jobName),
                     jobScript), wait = T)
  }, combs$trait_bin, combs$burden_bin))
}))