library(data.table)
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(stringr)

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

# Collect FDRs from all conditions
fdrs = apply(conditions, 1, function(row) {
  varCol = as.character(row[["variant_effect"]])
  varType = as.character(row[["variant_type"]])
  modeInheritance = as.character(row[["mode_inheritance"]])
  printInfo(paste("Process condition:", varCol, varType, modeInheritance))
  
  # Load significant gene-trait combinations
  fdrPath = paste(varCol, varType, modeInheritance, FDR_SIGNIFICANT_FILE_PATH, sep = "_")
  fdrPath = paste(runFolder, fdrPath, sep = "/")
  fdrPath = sprintf(fdrPath, as.character(FDR_CUTOFF))
  if (!file.exists(fdrPath)) return(NULL)
  fdr = fread(fdrPath)
  
  # If existing gene-trait combinations from the same condition with multiple records,
  # keep the one with smallest P value
  fdr = fdr[, .SD[which.min(p_val)], 
            by = c("gene", "phenotype", "type_variants", "qualifying_variants", "mode_inheritance")]
  
  return(fdr)
}, simplify = F)
fdrs = rbindlist(fdrs)

# Adjust the cutoff
fdrs = fdrs[fdr < FDR_CUTOFF]

# Save combined FDRs to file
outputPath = paste0(runFolder, "/all_conditions_", FDR_SIGNIFICANT_FILE_PATH)
outputPath = sprintf(outputPath, as.character(FDR_CUTOFF))
fwrite(fdrs, outputPath)

# Determine overlapping gene-trait combinations
plotTable = fdrs[
  mode_inheritance == "dominant" &
    type_variants == "missense_and_protein_truncating", .(
  comb = paste(gene, phenotype, sep = "-"), mode_inheritance, 
  type_variants, qualifying_variants, p_val
  )]
rc_unique = plotTable[
   qualifying_variants == "raw_count", unique(str_remove(comb, "_top|_bottom"))
  ]
r_unique = plotTable[
  qualifying_variants == "REVEL_score", unique(str_remove(comb, "_top|_bottom"))
  ]
v_unique = plotTable[
  qualifying_variants == "VARITY_R", unique(str_remove(comb, "_top|_bottom"))
]

# Plot Venn diagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn.diagram(
  x = list(rc_unique, r_unique, v_unique),
  category.names = c("Raw Count", "REVEL", "VARITY"),
  main = "Overlap of gene-trait combinations from different conditions",
  sub = sprintf("Variant type: Missense + PTV; Mode: Dom; FDR < 0.05"),
  filename = "missense_ptv_dom_venn_fdr_0.05.png", output = T,
  lwd = 2, lty = "blank", fill = brewer.pal(3, "Pastel2"),
  cex = 1.5, fontface = "bold", fontfamily = "sans",
  cat.cex = 1, cat.fontface = "bold", 
  cat.default.pos = "outer", cat.fontfamily = "sans",
  imagetype = "png", height = 5.5, width = 5.5, units = "in", resolution = 300
  )

# Plot P value distribution
plot = ggplot(plotTable) +
  geom_histogram(aes(x = p_val, fill = qualifying_variants), binwidth = 4) +
  geom_vline(xintercept = max(plotTable$p_val), linetype = 8, colour = "white") +
  annotate(geom = "label", x = max(plotTable$p_val), 
           y = Inf, vjust = 1, hjust = 0.75,
           label = sprintf("Max P: %s", format(max(plotTable$p_val), 
                                                scientific = T, digits = 2)),
           label.padding = unit(2, "mm")) +
  scale_x_log10() +
  ggtitle(label = "P value distribution of significant gene-trait combinations",
          subtitle = "Mode of Inheritance: dominant; Variant Types: Missense + PTV; FDR < 0.05") +
  labs(x = "P value", y = "# of Gene-trait combinations",
       fill = "Method to filter variants") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom")
ggsave("p_distribution_sig_combs_missense+ptv_dom_fdr_0.05.png", plot,
       width = 10, height = 5, units = "in", dpi = 300)
