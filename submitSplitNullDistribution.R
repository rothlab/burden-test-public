### Libraries
library(data.table)

### Functions
source("src/logger.R")
source("src/globalVariables.R")

# Parse arguments
args = commandArgs(trailingOnly = T)
if (length(args) != 1) stop("Invalid parameters")
runFolder = args[1]

# Set up conditions
predictorCols = unname(c(VARIANT_PREDICTORS, VARIANT_PREDICTOR_VARITY))
conditions = expand.grid(
  mode_inheritance = MODE_INHERITANCE,
  variant_type = VARIANT_TYPE,
  variant_effect = c("raw_count", predictorCols)
)

# Submit job for each condition
invisible(sapply(1:nrow(conditions), 1, function(index) {
  row = conditions[index,]
  varCol = as.character(row[["variant_effect"]])
  varType = as.character(row[["variant_type"]])
  modeInheritance = as.character(row[["mode_inheritance"]])
  printInfo(paste("Process condition:", varCol, varType, modeInheritance))
  
  jobName = paste(varCol, varType, modeInheritance, sep = "_")
  jobScript = sprintf("Rscript src/splitNullDistribution.R %s %s", runFolder, index)
  system2("submitjob.sh",
          args = c("-n", jobName, "-t", "24:00:00", "-m", "50G", "-c", "1",
                   "-e", sprintf("%s/%s.log", runFolder, jobName),
                   jobScript), wait = T)
}))