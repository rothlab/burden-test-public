### Libraries
library(data.table)

### Functions
source("src/logger.R")
source("src/globalVariables.R")

# Parse arguments
args = commandArgs(trailingOnly = T)
if (length(args) != 2) stop("Invalid parameters")
runFolder = args[1]
index = args[2]

# Load p value files
pValPaths = list.files(runFolder, BURDEN_TEST_PVAL_FILE_SUFFIX, full.names = T)

# Set up conditions
predictorCols = unname(c(VARIANT_PREDICTORS, VARIANT_PREDICTOR_VARITY))
conditions = expand.grid(
  mode_inheritance = MODE_INHERITANCE,
  variant_type = VARIANT_TYPE,
  variant_effect = c("raw_count", predictorCols)
)

# Submit job for each condition
row = conditions[index,]
varCol = as.character(row[["variant_effect"]])
varType = as.character(row[["variant_type"]])
modeInheritance = as.character(row[["mode_inheritance"]])
printInfo(paste("Process condition:", varCol, varType, modeInheritance))

# Load subset p values
subPVals = lapply(pValPaths, function(filePath) {
  printInfo(paste("===> Loading p values file:", filePath))
  pVals = fread(filePath)
  # Subset p values
  pVals[
    mode_inheritance == modeInheritance
    & type_variants == varType 
    & qualifying_variants == varCol
    ]
})
subPVals = rbindlist(subPVals)
gc()
printInfo("===> Finished loading P values")
if (nrow(subPVals) < 1) {
  printInfo("===> No P values. Exited.")
  invokeRestart("abort") # Exit the program
}

# Compute trait prevalence and burden (cases and controls with variants)
subPVals[, num_cases := num_cases_with_variants + num_cases_without_variants]
subPVals[, gene_burden := num_cases_with_variants + num_controls_with_variants]

# Assign bins to trait prevalence and burden
subPVals[, num_cases_bin := 
           as.numeric(cut(num_cases, quantile(num_cases, probs = 0:10/10)))]
subPVals[, gene_burden_bin :=
           as.numeric(cut(gene_burden, quantile(gene_burden, probs = 0:10/10)))]

# Calculate FDRs
combs = expand.grid(trait_bin = 1:10, burden_bin = 1:10)
calculateFDR = function(inds, ps, permutedPs) {
  fdr = sapply(inds, function(index) {
    p = ps[index]
    nPs = sum(ps <= p)
    nPermutedPs = sum(permutedPs <= p)
    ratio = min(1, max(nPermutedPs / nPs, 1/length(permutedPs))) # Set nominal FDR to the number of permuted values
    
    return(c(p = p, num_p = nPs, num_permuted_p = nPermutedPs, fdr = ratio))
  })
  fdr = as.data.table(t(fdr))
  
  return(fdr)
}
subPVals = mapply(function(trait, burden) {
  pVals = subPVals[num_cases_bin == trait & gene_burden_bin == burden]
  ps = pVals[is_real_data == T]$p_val
  permutedPs = pVals[is_real_data == F]$p_val
  
  # Select significant p values (< 0.05)
  subSignificantPIndices = which(ps < 0.05)
  fdr = calculateFDR(subSignificantPIndices, ps, permutedPs)

  # Merge FDRs into the p-value table
  pVals = cbind(pVals[is_real_data == T][subSignificantPIndices, ], fdr)
  pVals$is_real_data = NULL
  
  # Check if output p values matching input p values
  if (nrow(pVals[p != p_val]) > 0) stop("====> Mismatching p values and FDRs. Stopped.")
  pVals$p = NULL
  
  return(pVals)
}, combs$trait_bin, combs$burden_bin, SIMPLIFY = F)
subPVals = rbindlist(subPVals)
gc()

# Determine P value cutoff
# Here, warning (no p_vals with FDR less than cutoff) is suppressed and -Inf is returned
# because we will just use -Inf as the cutoff.
# Nothing will pass the cutoff, as expected.
subPVals[, past_fdr_cutoff := fdr < FDR_CUTOFF]
subPVals[, fdr := unlist(fdr)]

# Save all gene-trait combinations and only the significant ones
outputPath = paste(varCol, varType, modeInheritance, FDR_ALL_FILE_PATH, sep = "_")
outputPath = paste(runFolder, outputPath, sep = "/")
fwrite(subPVals, outputPath)
outputPath = paste(varCol, varType, modeInheritance, FDR_SIGNIFICANT_FILE_PATH, sep = "_")
outputPath = paste(runFolder, outputPath, sep = "/")
outputPath = sprintf(outputPath, as.character(FDR_CUTOFF))
fwrite(subPVals[past_fdr_cutoff == T], outputPath)
printInfo("===> Saved all gene-trait combinations and only the significant ones to files")

printInfo("All done!")