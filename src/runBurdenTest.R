### Libraries
library(data.table)
library(stringr)

### Functions
source("src/globalVariables.R")

runBurdenTest = function(varPredictions, filteredVariants, phenotypes,
                         measurements, permutedMeasurements, burdenTestOutputFilePath) {
  # Get unique genes
  uniqueGenes = unique(varPredictions$Gene)
  predictorCols = unname(c(VARIANT_PREDICTORS, VARIANT_PREDICTOR_VARITY))
  
  # Create burden test conditions
  conditions = expand.grid(
    phenotype = phenotypes,
    mode_inheritance = MODE_INHERITANCE,
    variant_type = VARIANT_TYPE,
    variant_effect = c("raw_count", predictorCols)
  )
  apply(conditions[1:5,], 1, function(row) {
    phenoCol = as.character(row[["phenotype"]])
    varCol = as.character(row[["variant_effect"]])
    varType = as.character(row[["variant_type"]])
    modeInheritance = as.character(row[["mode_inheritance"]])
    
    # Associate variant predictions to filtered UKB variants
    variants = merge(
      varPredictions, filteredVariants, by.x = "#Uploaded_variation", 
      by.y = "variant_name"
    )
    
    # Filter variants by variant types
    if (varType == "synonymous")
      variants = variants[Consequence == VARIANT_PREDICTION_FILTER_SYN]
    if (varType == "missense")
      variants = variants[Consequence %in% VARIANT_PREDICTION_FILTER_MISSENSE]
    if (varType == "protein_truncating")
      variants = variants[Consequence %in% VARIANT_PREDICTION_FILTER_PTV]
    if (varType == "missense_and_protein_truncating")
      variants = variants[Consequence %in% c(VARIANT_PREDICTION_FILTER_MISSENSE,
                                             VARIANT_PREDICTION_FILTER_PTV)]
    
    # Merge variants with phenotype measurements
    variantMeasurements = merge(
      variants,
      measurements[, c("eid", phenoCol), with = F],
      by = "eid"
    )
    variantPermutedMeasurements = merge(
      variants,
      permutedMeasurements[, c("eid", phenoCol), with = F],
      by = "eid"
    )
    
    if (varCol != "raw_count") {
      if (varType %in% c("missense", "missense_and_protein_truncating")) {
        # Convert variant predictions into PASS/FAIL qualification
        # Picking the cutoff as the 95th percentile of the control distribution
        predScores = variantMeasurements[
          Consequence %in% VARIANT_PREDICTION_FILTER_MISSENSE,
          c(phenoCol, varCol), with = F
        ]
        predScores = predScores[get(phenoCol) == 0, na.omit(get(varCol))]
        
        # If number of control variants less than the threshold, return NULL
        if (length(predScores) < CONTROL_NUM_VARIANT_CUTOFF) return(NULL)
        cutoff = quantile(predScores, c(0.95))
        cutoff = cutoff[[1]]
        variantMeasurements[, (varCol) := as.integer(get(varCol) > cutoff)]
        variantPermutedMeasurements[, (varCol) := as.integer(get(varCol) > cutoff)]
      }
      
      if (varType %in% c("protein_truncating", "missense_and_protein_truncating")) {
        # Convert protein truncating variants to full weight
        variantMeasurements[
          Consequence %in% VARIANT_PREDICTION_FILTER_PTV,
          (varCol) := 1
        ]
        variantPermutedMeasurements[
          Consequence %in% VARIANT_PREDICTION_FILTER_PTV,
          (varCol) := 1
        ]
      }
    }
    
    # Remove NA from measurements
    # NA removal only needs to be done once for each phenotype
    # It's moved to here to save time
    mReal = measurements[!is.na(get(phenoCol))]
    mPermuted = permutedMeasurements[!is.na(get(phenoCol))]
    
    # Run burden test on all unique genes
    pVals = lapply(uniqueGenes, function(gene) {
      # Get relevant variants and collapse into gene
      # Rule: keep participants with at least one qualifying variant
      if (varCol == "raw_count") {
        subVariantMeasurements = variantMeasurements[
          Gene == gene,
          c("eid", phenoCol),
          with = F
        ]
        subVariantPermutedMeasurements = variantPermutedMeasurements[
          Gene == gene,
          c("eid", phenoCol),
          with = F
        ]
        # If recessive model, only consider participants with at least two participants
        if (modeInheritance == "recessive") {
          subVariantMeasurements = subVariantMeasurements[duplicated(subVariantMeasurements)]
          subVariantPermutedMeasurements = subVariantPermutedMeasurements[duplicated(subVariantPermutedMeasurements)]
        }
        subVariantMeasurements = unique(subVariantMeasurements)
        subVariantPermutedMeasurements = unique(subVariantPermutedMeasurements)
      } else {
        subVariantMeasurements = variantMeasurements[
          Gene == gene & get(varCol) == 1,
          c("eid", phenoCol, varCol),
          with = F
        ]
        subVariantPermutedMeasurements = variantPermutedMeasurements[
          Gene == gene & get(varCol) == 1,
          c("eid", phenoCol, varCol),
          with = F
        ]
        # If recessive model, only consider participants with at least two participants
        if (modeInheritance == "recessive") {
          subVariantMeasurements = subVariantMeasurements[duplicated(subVariantMeasurements)]
          subVariantPermutedMeasurements = subVariantPermutedMeasurements[duplicated(subVariantPermutedMeasurements)]
        }
        subVariantMeasurements = unique(subVariantMeasurements)
        subVariantPermutedMeasurements = unique(subVariantPermutedMeasurements)
      }
            
      # Compute p value for each gene, phenotype and variant effect predictor
      computePVal = function(useRealData) {
        # Determine which dataset to use
        if (useRealData) {
          m = mReal
          sub = subVariantMeasurements[!is.na(get(phenoCol))]
        } else {
          m = mPermuted
          sub = subVariantPermutedMeasurements[!is.na(get(phenoCol))]
        }
        
        # Compute p value for raw variant count
        sumSub = sub[, sum(get(phenoCol))]
        sumM = m[, sum(get(phenoCol))]
        numCasesWithVariants = sumSub
        numCasesWithoutariants = sumM - sumSub
        numControlsWithVariants = nrow(sub) - sumSub
        numControlsWithoutVariants = nrow(m) - sumM - numControlsWithVariants
        
        # If the number of cases and controls does not meet the cutoff, return NULL
        if ((numCasesWithVariants + numCasesWithoutariants) < PHENOTYPE_NUM_PARTICIPANT_CUTOFF |
            (numControlsWithVariants + numControlsWithoutVariants) < PHENOTYPE_NUM_PARTICIPANT_CUTOFF) {
          return(NULL)
        }
        
        # If the number of participants with qualifying variant is less than the cutoff
        # return NULL
        if ((numCasesWithVariants + numControlsWithVariants) < PHENOTYPE_NUM_PARTICIPANT_CUTOFF) {
          return(NULL)
        }
        
        # Compute fisher's exact P value
        pVal = data.table(
          gene = gene, phenotype = phenoCol, 
          type_variants = varType, qualifying_variants = varCol,
          mode_inheritance = modeInheritance, is_real_data = useRealData,
          num_cases_with_variants = numCasesWithVariants,
          num_cases_without_variants = numCasesWithoutariants,
          num_controls_with_variants = numControlsWithVariants,
          num_controls_without_variants = numControlsWithoutVariants
        )
        p = fisher.test(matrix(c(numCasesWithVariants, numControlsWithVariants,
                                 numCasesWithoutariants, numControlsWithoutVariants), 
                               nrow = 2, byrow = T),
                        alternative = "two.sided")
        pVal[, c("p_val", "odds_ratio") := list(p$p.value, p$estimate)]
        
        return(pVal) 
      }
      
      # Compute p values for real and permuted data
      pVals = rbind(computePVal(T), computePVal(F))
      
      return(pVals)
    })
    pVals = rbindlist(pVals)
    
    # Append p values to output file
    fwrite(pVals, burdenTestOutputFilePath,
           append = file.exists(burdenTestOutputFilePath))
    
    invisible()
  }, simplify = F)
}
