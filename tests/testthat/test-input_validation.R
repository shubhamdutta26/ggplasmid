# Create sample data for tests
valid_data <- data.frame(
  start = c(1, 100, 500),
  end = c(50, 200, 800),
  length = c(1000, 1000, 1000),
  stringsAsFactors = FALSE
)

test_that("validate_plasmid_data accepts valid data", {
  expect_no_error(validate_plasmid_data(valid_data))
  expect_identical(validate_plasmid_data(valid_data), valid_data)
})

test_that("validate_plasmid_data checks for required columns", {
  # Missing 'start' column
  invalid_data <- valid_data[, c("end", "length")]
  expect_error(
    validate_plasmid_data(invalid_data),
    "Missing required columns: start"
  )

  # Missing multiple columns
  invalid_data <- valid_data[, "length", drop = FALSE]
  expect_error(
    validate_plasmid_data(invalid_data),
    "Missing required columns: start, end"
  )
})

test_that("validate_plasmid_data validates start values", {
  # Negative start values
  invalid_data <- valid_data
  invalid_data$start[1] <- -5
  expect_error(
    validate_plasmid_data(invalid_data),
    "'start' column must contain non-negative numeric values"
  )

  # Non-numeric start values
  invalid_data <- valid_data
  invalid_data$start <- as.character(invalid_data$start)
  expect_error(
    validate_plasmid_data(invalid_data),
    "'start' column must contain non-negative numeric values"
  )
})

test_that("validate_plasmid_data validates end values", {
  # End values less than start values
  invalid_data <- valid_data
  invalid_data$end[1] <- invalid_data$start[1] - 1
  expect_error(
    validate_plasmid_data(invalid_data),
    "'end' column must contain numeric values greater than corresponding start values"
  )

  # Non-numeric end values
  invalid_data <- valid_data
  invalid_data$end <- as.character(invalid_data$end)
  expect_error(
    validate_plasmid_data(invalid_data),
    "'end' column must contain numeric values greater than corresponding start values"
  )
})

test_that("validate_plasmid_data validates length consistency", {
  # Different length values
  invalid_data <- valid_data
  invalid_data$length[1] <- 2000
  expect_error(
    validate_plasmid_data(invalid_data),
    "All entries must have the same length value"
  )

  # Length less than maximum end position
  invalid_data <- valid_data
  invalid_data$length <- 700 # Less than max end position of 800
  expect_error(
    validate_plasmid_data(invalid_data),
    "Length must be greater than or equal to the maximum end position"
  )
})

# Tests for read_plasmid_data function
test_that("read_plasmid_data handles data frame input", {
  expect_no_error(read_plasmid_data(valid_data))
  expect_identical(read_plasmid_data(valid_data), valid_data)
})

test_that("read_plasmid_data validates input type", {
  expect_error(
    read_plasmid_data(list(a = 1)),
    "Input must be a data frame or a file path to CSV/Excel file"
  )

  expect_error(
    read_plasmid_data(42),
    "Input must be a data frame or a file path to CSV/Excel file"
  )
})

test_that("read_plasmid_data handles file extensions", {
  expect_error(
    read_plasmid_data("data.txt"),
    "Unsupported file format. Please provide CSV or Excel file"
  )
})

# File-based tests using temporary files
test_that("read_plasmid_data can read CSV files", {
  # Create temporary CSV file
  temp_csv <- tempfile(fileext = ".csv")
  write.csv(valid_data, temp_csv, row.names = FALSE)

  # Test reading
  result <- read_plasmid_data(temp_csv)
  expect_equal(result, valid_data)

  # Clean up
  unlink(temp_csv)
})

# Mock test for Excel file handling
test_that("read_plasmid_data handles Excel file requirements", {
  skip_if_not_installed("readxl")

  # Test with non-existent Excel file to check error message
  expect_error(
    read_plasmid_data("nonexistent.xlsx"),
    "cannot open|cannot read|does not exist"
  )
})

# Test data validation after reading
test_that("read_plasmid_data validates data after reading", {
  # Create temporary CSV with invalid data
  invalid_data <- valid_data
  invalid_data$start[1] <- -5
  temp_csv <- tempfile(fileext = ".csv")
  write.csv(invalid_data, temp_csv, row.names = FALSE)

  # Test that validation occurs after reading
  expect_error(
    read_plasmid_data(temp_csv),
    "'start' column must contain non-negative numeric values"
  )

  # Clean up
  unlink(temp_csv)
})
