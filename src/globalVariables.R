### Variables in main.R
FILTERED_VARIANT_FILE_SUFFIX = "_filtered_mut.csv"
VARIANT_PREDICTION_FILE_SUFFIX = "_all_weights.csv"
VARIANT_PREDICTION_FILE_COLUMNS_TO_KEEP = c(
  "#Uploaded_variation", "Gene", "Feature", "Feature_type", "Protein_position",
  "Amino_acids", "Codons", "Ensembl_geneid", "Ensembl_proteinid",
  "Ensembl_transcriptid", "genename", "Uniprot_acc", "gnomAD_AF", "Consequence",
  "HGVSp_VEP"
)
VARIANT_PREDICTION_FILE_COLUMNS_COMMA_SPLITTED = c(
  "Ensembl_geneid", "Ensembl_proteinid", "Ensembl_transcriptid", "genename",
  "Uniprot_acc", "HGVSp_VEP"
)
VARIANT_PREDICTION_FILTER_SYN = "synonymous_variant"
VARIANT_PREDICTION_FILTER_MISSENSE = c("missense_variant,splice_region_variant",
                                       "missense_variant")
VARIANT_PREDICTION_FILTER_PTV = c("stop_gained", "stop_lost", "start_lost",
                                  "stop_gained,splice_region_variant")
VARIANT_PREDICTION_FILTER_AF = 0.001
VARIANT_PREDICTORS = c("REVEL" = "REVEL_score")
VARIANT_PREDICTOR_FILE_VARITY = "data/varity_all_predictions.txt"
VARIANT_PREDICTOR_FILE_VARITY_COLUMNS_TO_KEEP = c(
  "p_vid", "aa_pos", "aa_ref", "aa_alt", "VARITY_R"
)
VARIANT_PREDICTOR_VARITY = "VARITY_R"
BURDEN_TEST_PVAL_FILE_SUFFIX = "_pvals.csv"

### Variables in processPhenotypes.R
EIDS_FILE_PATH = "data/450k_eids.txt"
EXOME_SEQ_FIELD_ID = "23141"
PHENOTYPE_MEASUREMENT_FILE_PATH = "data/phenotypes_all_450k.csv"
PHENOTYPE_MEASUREMENT_PROCESSED_FILE_PATH = "cache/measurements_filtered_processed.csv"
PHENOTYPE_SELECTION_FILE_PATH = "data/phenotype_description.csv"
PHENOTYPE_MIN_PARTICIPANT_CUTOFF = 1000
WITHDRAWN_PARTICIPANT_FILE_PATH = "data/w51135_all.csv"
ETHNICITY_FIELD = c("id" = "21000", "prefix" = "1", "description" = "Caucasian")
RAND_SEED = 117119
CACHED_MEASUREMENT_FILE_PATH = "cache/measurements_rand_seed_%d.csv"
CACHED_PERMUTED_MEASUREMENT_FILE_PATH = "cache/permuted_measurements_rand_seed_%d.csv"

### Variables in runBurdenTest.R
PHENOTYPE_NUM_PARTICIPANT_CUTOFF = 10
CONTROL_NUM_VARIANT_CUTOFF = 50
MODE_INHERITANCE = c("dominant", "recessive")
VARIANT_TYPE = c("missense", "protein_truncating", "synonymous",
                 "missense_and_protein_truncating")

### Variable in processVarityScores.R
CACHED_VARITY_FILE_FORMAT = "cache/varity_predictions_chr_%s.txt"

### Variable in parseOutput.R
FDR_ALL_FILE_PATH = "fdr_all_combined.csv"
FDR_SIGNIFICANT_FILE_PATH = "fdr_significant_cutoff-%s_combined.csv"
FDR_SIGNIFICANT_COLLAPSED_FILE_PATH = "fdr_significant_collapsed_cutoff-%s_combined.csv"
FDR_CUTOFF = 0.05
SCATTERPLOT_FILE_PATH = "pval_fdr_scatterplot.png"
HISTOGRAM_FILE_PATH = "pval_fdr_histogram.png"
VENN_DIAGRAM_FILE_PATH = "gene-trait_venn_cutoff-%s.png"
VENN_DIAGRAM_BY_CONDITION_FILE_PATH = "gene-trait_by_condition_venn_cutoff-%s.png"
