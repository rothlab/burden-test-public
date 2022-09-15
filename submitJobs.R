### Libraries
library(stringr)
library(stringi)

### Functions
source("src/globalVariables.R")
source("src/logger.R")

# Parse arguments
args = commandArgs(trailingOnly = T)
if (length(args) > 0) {
  runId = args[1]
  runFolder = paste("ukb_output/run", runId, sep = "_")
  if (!dir.exists(runFolder)) throwError("The associated run does not exist")
} else {
  # Generate run ID and make sure it's not used yet
  runId = stri_rand_strings(1, 10)
  runFolder = paste("ukb_output/run", runId, sep = "_")
  while(dir.exists(runFolder)) {
    runId = stri_rand_strings(1, 10)
    runFolder = paste("ukb_output/run", runId, sep = "_")
  }
  dir.create(runFolder)
}

# List all variant effect prediction files
# and use them to guide analysis
varPredictionFilePaths = list.files(VARIANT_INPUT_FOLDER, VARIANT_PREDICTION_FILE_SUFFIX, full.names = T)

# If no measurements cached, create them
cachedMeasurementFilePath = sprintf(CACHED_MEASUREMENT_FILE_PATH, RAND_SEED)
cachedPermutedFilePath = sprintf(CACHED_PERMUTED_MEASUREMENT_FILE_PATH, RAND_SEED)
if (!file.exists(cachedMeasurementFilePath) |
    !file.exists(cachedPermutedFilePath))
  source("src/processPhenotypes.R")

# Select jobs to submit
jobsToSubmit = basename(varPredictionFilePaths)
jobsToSubmit = str_remove(jobsToSubmit, VARIANT_PREDICTION_FILE_SUFFIX)

# List completed jobs and remove them from jobs to submit
completedJobs = list.files("ukb_output", BURDEN_TEST_PVAL_FILE_SUFFIX)
if (length(completedJobs) > 0) {
  completedJobs = basename(completedJobs)
  completedJobs = str_remove(completedJobs, BURDEN_TEST_PVAL_FILE_SUFFIX)
  printInfo(sprintf("%d complete jobs detected", length(completedJobs)))
  jobsToSubmit = setdiff(jobsToSubmit, completedJobs)
}

# Submit jobs
lapply(jobsToSubmit, function(jobName) {
  varPredictionFilePath = sprintf(
    "%s/%s%s", VARIANT_INPUT_FOLDER, jobName, VARIANT_PREDICTION_FILE_SUFFIX)
  jobScript = sprintf("Rscript src/main.R %s %s", varPredictionFilePath, runFolder)
  system2("submitjob.sh",
          args = c("-n", jobName, "-t", "6:00:00", "-m", "8G",
                   "-e", sprintf("%s/%s.log", runFolder, jobName),
                   jobScript), wait = T)
})
