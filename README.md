# UK Biobank Burden Scan Pipeline

This Github repository contains the burden scan pipeline used in the manuscript "Empowering rare variant burden-based gene-trait association studies via optimized computational predictor choice" (in preparation).

**Please note**, this repository does NOT contain raw data from the UK Biobank as they are available upon application. Please visit the [UK Biobank website](https://www.ukbiobank.ac.uk/) to apply for data access.

## Install

To install the pipeline, first make sure your R version is at least R 4.0. You can check by typing the following into your R console:

```{r}
R.version()$major
```

Next, clone the repository by typing the following into your command line interface:

```{bash}
git clone https://github.com/rothlab/ukb_burden_test
```

## Required input files

### 1. VARITY predictions

You need to download VARITY predictions (http://varity.varianteffect.org/downloads/varity_all_predictions.tar.gz), unzip the file (varity_all_predictions.txt) and move it to the `data` folder.

### 2. IDs for participants with exome sequences

You need to write a list of EIDs for participants with exome sequences into a txt file (`450k_eids.txt`) and store it in the `data` folder. The text file should contain one EID per line. A template file is included in the `data` folder for your reference.

### 3. Phenotype measurements

You need to write all phenotype measurements into a CSV file (`phenotypes_all_450k.csv`) and store it in the `data` folder. 

The first column (`eid`) stores the EID of the participant. The remaining columns store phenotype measurements. The first row contains column names. Each row thereafter represents one UKB participant. The column arrangement can be determined using the template file included in the `data` folder as a refernece.

### 4. UKB variants

UK Biobank variants must be pre-processed and filtered (see Method of the manuscript for details) on the UKB Research Analysis Platform. Filtered variatns are securely transferred to an UKB-approved local cluster, together with phenotype measurements mentioned above. Use the template file (`example_filtereed_mut.csv`) included in the `ukb_variants` folder as a reference.

Variant effect predictor scores must also be curated. Use the template file (`example_all_weights.csv`) included in the `ukb_variants` folder as a reference.

## Confiugre the running parameters

The `src/globalVariables.R` file lists all the running parameters that you can configure.

Here we document all parameters supported by the pipeline and their functions.

| Parameter | Description | Default Value |
| --- | --- | --- |
| FILTERED_VARIANT_FILE_SUFFIX | The suffix of filtered variant files | _filtered_mut.csv |
| VARIANT_PREDICTION_FILE_SUFFIX | The suffix of variant prediction files | _all_weights.csv |
| VARIANT_PREDICTION_FILE_COLUMNS_TO_KEEP | Columns in variant prediction files to keep | *!incomplete list!* "#Uploaded_variation", "Gene", "Feature", ... |
| VARIANT_PREDICTION_FILE_COLUMNS_COMMA_SPLITTED | Columns in variant prediction files that contain comma splitted data | *!incomplete list!* "Ensembl_geneid", "Ensembl_proteinid", "Ensembl_transcriptid", ... |
| VARIANT_PREDICTION_FILTER_SYN | Name of synonymous variants in variant prediction file | synonymous_variant |
| VARIANT_PREDICTION_FILTER_MISSENSE | Name of missense variants in variant prediction file | "missense_variant,splice_region_variant", "missense_variant" |
| VARIANT_PREDICTION_FILTER_PTV | Name of protein-truncating variants in variant prediction file | *!incomplete list!* "stop_gained", "stop_lost", ... |
| VARIANT_PREDICTION_FILTER_AF | The maximum allele frequency that we will allow when filtering variants | 0.001 |
| VARIANT_PREDICTORS | The list of variant effect predictors (except VARITY) to use for the burden scan | "REVEL" = "REVEL_score" |
| VARIANT_PREDICTOR_FILE_VARITY | The path to VARITY scores | data/varity_all_predictions.txt |
| VARIANT_PREDICTOR_FILE_VARITY_COLUMNS_TO_KEEP | Columns to keep in the variant prediction files | *!incomplete list!* "p_vid", "aa_pos", ... |
| VARIANT_PREDICTOR_VARITY | The VARITY flavour to use | VARITY_R |
| BURDEN_TEST_PVAL_FILE_SUFFIX | The suffix of p value output files | _pvals.csv |
| EIDS_FILE_PATH | The file containing all participant IDs (EIDs) with exome sequences available | data/450k_eids.txt |
| EXOME_SEQ_FIELD_ID | The field ID of exomes | 23141 |
| PHENOTYPE_MEASUREMENT_FILE_PATH | The file of phenotype measurements | data/phenotypes_all_450k.csv |
| PHENOTYPE_MEASUREMENT_PROCESSED_FILE_PATH | The processed file of phenotype measurements | cache/measurements_filtered_processed.csv |
| PHENOTYPE_SELECTION_FILE_PATH | The file containing relevant phenotypes | data/phenotype_description.csv |
| PHENOTYPE_MIN_PARTICIPANT_CUTOFF | The number of participants that each phenotype needs to have in order to be considered | 1000 |
| WITHDRAWN_PARTICIPANT_FILE_PATH | The file of participants that have withdrawn their consent to participate in UK Biobank | data/withdraws.csv |
| ETHNICITY_FIELD | Only participants of this ethnitiy are considered | "id" = "21000", "prefix" = "1", "description" = "Caucasian" |
| RAND_SEED | The random seed for permutation | 117119 |
| CACHED_MEASUREMENT_FILE_PATH | The pattern of cached measurement file | cache/measurements_rand_seed_%d.csv |
| CACHED_PERMUTED_MEASUREMENT_FILE_PATH | The pattern of cached permuted measurement file | cache/permuted_measurements_rand_seed_%d.csv |
| PHENOTYPE_NUM_PARTICIPANT_CUTOFF | The minimum number of participant that each phenotype must have that also carries a qualifying variant | 10 |
| CONTROL_NUM_VARIANT_CUTOFF | The minimum number of control variants for a gene-trait combination to be considered | 50 |
| MODE_INHERITANCE | The mode of inheritance models used | "dominant", "recessive" |
| VARIANT_TYPE | The types of variants to consider | "missense", "protein_truncating", "synonymous", "missense_and_protein_truncating" |
| CACHED_VARITY_FILE_FORMAT | The pattern of cached VARITY scores separated by chromosome | cache/varity_predictions_chr_%s.txt |
| FDR_ALL_FILE_PATH | The file to save all FDRs | fdr_all_combined.csv |
| FDR_SIGNIFICANT_FILE_PATH | The file to save *significant* FDRs | fdr_significant_cutoff-%s_combined.csv |
| FDR_SIGNIFICANT_COLLAPSED_FILE_PATH | The file to save FDRs collapsed by gene-trait combination | fdr_significant_collapsed_cutoff-%s_combined.csv |
| FDR_CUTOFF | The FDR threshold | 0.05 |
| SCATTERPLOT_FILE_PATH | The P value vs. FDR scatterplot | pval_fdr_scatterplot.png |
| HISTOGRAM_FILE_PATH | The histograms of P values and FDRs | pval_fdr_histogram.png |
| VENN_DIAGRAM_FILE_PATH | The overall Venn diagram of the number of unique gene-trait combinations detected | gene-trait_venn_cutoff-%s.png |
| VENN_DIAGRAM_BY_CONDITION_FILE_PATH | The Venn diagram of the number of unique gene-trait combinations detected in each condition | gene-trait_by_condition_venn_cutoff-%s.png |

## Run the pipeline

The analysis pipeline is designed to be run on a cluster with the SLURM scheduler. Job scheduling relies on the [clusterUtil](https://github.com/rothlab/clusterUtil) toolkit. Please make sure you have the compatible scheduler and toolkit installed before proceeding. **This pipelien will NOT run on a local machine.**

### Computing P values

To run ths part of the pipeline, type the following line into your command line interface:

```
Rscript submitJobs.R
```

If, for some reason, the submission script (`submitJobs.R`) was terminated before all jobs were complete, you can resume the pipeline using this command:

```
Rscript <run-id>
```

You can find the run ID by looking at the most recently created folder in the `output` folder.

Next, we need to split the null distribution into smaller quaduants by running the following line:

```
Rscript submitSplitNullDistribution.R <output-folder>
```

The program requires one arguments:

- Output folder: the output folder from the previous step.

### Computing FDRs

After P values were calculated, we next need to convert P values to FDRs.

To run ths part of the pipeline, type the following line into your command line interface:

```
Rscript submitParseOutputJobsWithSplittedNull.R <output-folder>
```

The program requires one arguments:

- Output folder: the output folder from the previous step.

### Finalizing results

Lastly, we need to finalize results and generate plots by running the following line:

```
Rscript finalizeResults <output-folder>
```

The program requires one arguments:

- Output folder: the output folder from the previous step.

## Contact us

If you have any feedback, suggestions or questions, please reach out via [email](mailto:kvn.kuang@mail.utoronto.ca) or open an issue on github.