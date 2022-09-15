### Libraries
library(data.table)
library(ggplot2)
library(stringr)

### Functions
source("src/logger.R")
source("src/globalVariables.R")

# Parse arguments
args = commandArgs(trailingOnly = T)
if (length(args) != 2) stop("Invalid parameters")
runFolder = args[1]
index = as.numeric(args[2])

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
    & qualifying_variants == varCol]
})
subPVals = rbindlist(subPVals)
printInfo("===> Finished loading P values")

if (nrow(subPVals) < 1) {
  printInfo("===> No P values. Exited.")
  return(NULL)
}

# Compute trait prevalence and burden (cases and controls with variants)
subPVals[, num_cases := num_cases_with_variants + num_cases_without_variants]
subPVals[, gene_burden := num_cases_with_variants + num_controls_with_variants]

# Assign bins to trait prevalence and burden
subPVals[, num_cases_bin := 
           as.numeric(cut(num_cases, quantile(num_cases, probs = 0:10/10)))]
subPVals[, gene_burden_bin :=
           as.numeric(cut(gene_burden, quantile(gene_burden, probs = 0:10/10)))]
getBreaks = function(values, levels) {
  breaks = levels(cut(values, levels))
  breaks = str_split(breaks, ",", simplify = T)[, 2]
  breaks = str_sub(breaks, end = -2)
  return(as.integer(breaks))
}
numCaseBreaks = getBreaks(subPVals$num_cases, 
                          quantile(subPVals$num_cases, probs = 0:10/10))
geneBurdenBreaks = getBreaks(subPVals$gene_burden,
                             quantile(subPVals$gene_burden, probs = 0:10/10))

# Cache files
combs = expand.grid(trait_bin = 1:10, burden_bin = 1:10)
summary = mapply(function(trait, burden) {
  pVals = subPVals[num_cases_bin == trait & gene_burden_bin == burden]
  pFilePath = sprintf(
    "%s/%s_%s_%s_%d-%d_p_cache.txt", runFolder, varCol, varType,
    modeInheritance, trait, burden)
  writeLines(as.character(pVals[is_real_data == TRUE, p_val]), pFilePath)
  permutedPFilePath = sprintf(
    "%s/%s_%s_%s_%d-%d_permuted_p_cache.txt", runFolder, varCol, varType,
    modeInheritance, trait, burden)
  writeLines(as.character(pVals[is_real_data == FALSE, p_val]), permutedPFilePath)
  
  return(data.table(trait_bin = trait, burden_bin = burden,
                    num_p = nrow(pVals[is_real_data == TRUE]),
                    num_permuted_p = nrow(pVals[is_real_data == FALSE])))
}, combs$trait_bin, combs$burden_bin, SIMPLIFY = F)
summary = rbindlist(summary)
printInfo("===> Finished caching files")

# Plot distribution
summary$log_num_permuted_p = log10(summary$num_permuted_p)
plot = ggplot(summary, aes(x = trait_bin, y = burden_bin)) +
  geom_tile(aes(fill = log_num_permuted_p)) +
  geom_text(aes(label = round(log_num_permuted_p, 1))) +
  scale_x_discrete(name = "# of participants w/ trait", limits = factor(1:10),
                   labels = numCaseBreaks) +
  scale_y_discrete(name = "# of participants w/ qualifying variants",
                   limits = factor(1:10),
                   labels = geneBurdenBreaks) +
  scale_fill_gradientn(colours = topo.colors(2), limits = c(0, 7)) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", legend.title.align = 1) +
  labs(title = "P values in null distribution amongst trait prevalence and gene burden",
       subtitle = sprintf("Predictor: %s; Variant Type: %s; Mode: %s", 
                          varCol, varType, modeInheritance), 
       fill = "# of gene-trait combinations\n(log10 scale)")

# Save plot
plotFile = sprintf("%s/%s_%s_%s_p_distribution.png", 
                   runFolder, varCol, varType, modeInheritance)
ggsave(plotFile, plot, width = 10, height = 10, units = "in", dpi = 300)
