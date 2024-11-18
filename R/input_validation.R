#' Validate plasmid data structure and contents
#'
#' @description
#' Internal function that performs comprehensive validation of plasmid data structure
#' and contents. Checks for required columns, appropriate data types, and logical
#' consistency of values.
#'
#' @param data A data frame containing plasmid annotation data
#'
#' @return The validated data frame, unchanged if all checks pass
#'
#' @details
#' Performs the following validations:
#' * Checks for presence of required columns: "start", "end", "length"
#' * Validates that 'start' contains non-negative numeric values
#' * Validates that 'end' values are numeric and greater than corresponding start values
#' * Ensures all entries have the same length value
#' * Verifies that length is greater than or equal to maximum end position
#' @noRd
#' @keywords internal
validate_plasmid_data <- function(data) {
  required_cols <- c("start", "end", "length")
  # Check for required columns
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    rlang::abort(
      message = sprintf(
        "Missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call = NULL
    )
  }

  if (!all(is.numeric(data$start)) || !all(data$start >= 0)) {
    rlang::abort("'start' column must contain non-negative numeric values.", call = NULL)
  }

  if (!all(is.numeric(data$end)) || !all(data$end > data$start)) {
    rlang::abort("'end' column must contain numeric values greater than corresponding start values.", call = NULL)
  }

  # Check length consistency
  unique_lengths <- unique(data$length)
  if (length(unique_lengths) != 1) {
    rlang::abort("All entries must have the same length value.", call = NULL)
  }

  if (unique_lengths < max(data$end)) {
    rlang::abort("Length must be greater than or equal to the maximum end position.", call = NULL)
  }

  data
}

#' Read and validate plasmid data from various input sources
#'
#' @description
#' Internal function that handles reading plasmid data from different input sources
#' and validates the data structure and contents.
#'
#' @param input Either a file path (character) to a CSV/Excel file or a data frame
#'
#' @return A validated data frame containing plasmid data
#'
#' @details
#' Supports the following input formats:
#' * CSV files (.csv extension)
#' * Excel files (.xlsx or .xls extension, requires readxl package)
#' * Data frames (passed directly)
#'
#' After reading the data, calls `validate_plasmid_data()` to ensure data integrity.
#'
#' @section Dependencies:
#' * readxl package is required for reading Excel files
#' * utils package for reading CSV files
#'
#' Throws error if:
#' * Input format is unsupported
#' * Required package (readxl) is missing for Excel files
#' * File reading fails
#' * Data validation fails
#' @noRd
#' @keywords internal
read_plasmid_data <- function(input) {
  if (is.character(input)) {
    # Check file extension
    ext <- tolower(tools::file_ext(input))

    if (ext == "csv") {
      data <- utils::read.csv(input, stringsAsFactors = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        rlang::abort("Package 'readxl' needed to read Excel files. Please install it.", call = NULL)
      }
      data <- readxl::read_excel(input)
      data <- as.data.frame(data)
    } else {
      rlang::abort("Unsupported file format. Please provide CSV or Excel file.", call = NULL)
    }
  } else if (is.data.frame(input)) {
    data <- input
  } else {
    rlang::abort("Input must be a data frame or a file path to CSV/Excel file.", call = NULL)
  }

  validated_data <- validate_plasmid_data(data)

  return(validated_data)
}
