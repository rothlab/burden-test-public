validatePhenotypeSelection = function(selection, filePath, error) {
  # Global variables
  SELECTION_COLUMNS = c(
    "field_id", "subset", "phenotype_description", "type"
    )
  SELECTION_COLUMN_TYPE = c("continuous", "binary", "basic")

  # Check if selection has the right columns
  if (!identical(sort(colnames(selection)), sort(SELECTION_COLUMNS)))
    throwError(sprintf("Column mismatch: %s expected; %s in file %s",
               paste0(SELECTION_COLUMNS, collapse = ", "),
               paste0(colnames(selection), collapse = ", "), filePath))
  
  # Check if the type column has only valid values
  colValues = unique(selection$type)
  if (!identical(sort(colValues), sort(SELECTION_COLUMN_TYPE)))
    throwError(sprintf("Type column invalid value: %s expected; %s in file %s",
                       paste0(SELECTION_COLUMN_TYPE, collapse = ", "),
                       paste0(colValues, collapse = ", "), filePath))
}
