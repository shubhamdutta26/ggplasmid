# Test helper functions
create_test_features <- function() {
  data.frame(
    start = c(1, 100, 200, 300),
    end = c(50, 150, 250, 350),
    length = rep(1000, 4),
    direction = c(1, -1, 0, 1),
    name = c("gene1", "gene2", "gene3", "gene4"),
    color = c("red", "blue", "green", "yellow"),
    stringsAsFactors = FALSE
  )
}

# Tests for text_pos function
test_that("text_pos returns correct positions", {
  # Test outer positions
  expect_equal(text_pos(0), "right")                # 0° - right
  expect_equal(text_pos(pi/2), "t_center")          # 90° - top center
  expect_equal(text_pos(pi), "left")                # 180° - left
  expect_equal(text_pos(3*pi/2), "b_center")        # 270° - bottom center

  # Test quarter positions
  expect_equal(text_pos(pi/3), "t_center")          # 60° - start of top center
  expect_equal(text_pos(2*pi/3), "t_center")        # 120° - end of top center
  expect_equal(text_pos(4*pi/3), "b_center")        # 240° - start of bottom center
  expect_equal(text_pos(5*pi/3), "b_center")        # 300° - end of bottom center

  # Test inner positions (shifted by π radians)
  expect_equal(text_pos(0, "inner"), "left")        # Inner: shifts to 180° - left
  expect_equal(text_pos(pi/2, "inner"), "b_center") # Inner: shifts to 270° - bottom center
  expect_equal(text_pos(pi, "inner"), "right")      # Inner: shifts to 0° - right
  expect_equal(text_pos(3*pi/2, "inner"), "t_center") # Inner: shifts to 90° - top center

  # Test edge cases
  expect_equal(text_pos(2*pi), "right")             # Full circle (same as 0°)
  expect_equal(text_pos(-pi/2), "b_center")         # Negative angle (normalizes to 270°)
})

# Tests for calc_glyphs function
test_that("calc_glyphs generates correct coordinates", {
  test_row <- data.frame(
    level = 0,
    direction = 1,
    rstart = 0,
    rend = pi/2,
    line_color = "black",
    fill_color = "white",
    stringsAsFactors = FALSE
  )

  result <- calc_glyphs(test_row)

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("x", "y", "Lx1", "Ly1", "anno_line_color",
                         "lineX", "lineY", "theta", "anno_pos"))

  # Check coordinates are numeric vectors
  expect_type(result$x, "double")
  expect_type(result$y, "double")
  expect_type(result$Lx1, "double")
  expect_type(result$Ly1, "double")

  # Check for valid coordinates
  expect_true(all(!is.na(result$x)))
  expect_true(all(!is.na(result$y)))
  expect_true(all(abs(result$x) <= 1))
  expect_true(all(abs(result$y) <= 1))
})

# Tests for calc_num_markers function
test_that("calc_num_markers generates appropriate markers", {
  # Test with different plasmid lengths
  small_result <- calc_num_markers(1000)
  medium_result <- calc_num_markers(5000)
  large_result <- calc_num_markers(10000)

  # Check structure
  expect_s3_class(small_result, "data.frame")
  expect_named(small_result, c("lineX", "lineY", "theta", "text_align",
                               "bp", "Lx1", "Ly1", "size"))

  # Check marker spacing
  expect_true(nrow(medium_result) > nrow(small_result))
  expect_false(all(diff(small_result$bp) >= 500))
  expect_true(all(small_result$bp <= 1000))
  expect_true(all(small_result$bp > 0))

  # Check text alignment values
  expect_true(all(small_result$text_align %in% c("right", "left", "t_center", "b_center")))
})

# Tests for calc_level function
test_that("calc_level assigns non-overlapping levels", {
  # Test with non-overlapping features
  non_overlapping <- data.frame(
    start = c(1, 101, 201),
    end = c(50, 150, 250),
    length = rep(300, 3),
    stringsAsFactors = FALSE
  )

  result <- calc_level(non_overlapping)
  expect_equal(max(result$level), 0)

  # Test with overlapping features
  overlapping <- data.frame(
    start = c(1, 25, 75),
    end = c(50, 75, 125),
    length = rep(300, 3),
    stringsAsFactors = FALSE
  )

  result <- calc_level(overlapping)
  expect_true(max(result$level) > 0)

  # Test with wrapping features
  wrapping <- data.frame(
    start = c(250, 50),
    end = c(50, 150),
    length = rep(300, 2),
    stringsAsFactors = FALSE
  )

  result <- calc_level(wrapping)
  expect_true(all(!is.na(result$level)))

  # Test empty input
  empty <- data.frame(
    start = numeric(0),
    end = numeric(0),
    length = numeric(0)
  )

  expect_equal(nrow(calc_level(empty)), 0)
})
