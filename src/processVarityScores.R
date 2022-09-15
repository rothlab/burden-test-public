### Libraries
library(data.table)

### Functions
source("src/globalVariables.R")

# Load VARITY scores
varityScores = fread(VARIANT_PREDICTOR_FILE_VARITY)

# Split VARITY scores by chromosome
lapply(unique(varityScores$chr), function(c) {
  sub = varityScores[chr == c,
                     VARIANT_PREDICTOR_FILE_VARITY_COLUMNS_TO_KEEP, with = F]
  fwrite(sub, sprintf(CACHED_VARITY_FILE_FORMAT, c))
  
  invisible()
})
