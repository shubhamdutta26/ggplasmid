# Tests for plot_plasmid function

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

test_that("plot_plasmid generates valid plots", {
  test_data <- create_test_features()

  # Test basic plot
  p <- plot_plasmid(test_data)
  expect_s3_class(p, "ggplot")

  # Test with custom fill
  p <- plot_plasmid(test_data, fill = "color")
  expect_s3_class(p, "ggplot")

  # Test with disabled features
  p <- plot_plasmid(test_data,
                    show.center.label = FALSE,
                    show.labels = FALSE,
                    show.ticks = FALSE)
  expect_s3_class(p, "ggplot")

  # Test with custom center label
  p <- plot_plasmid(test_data,
                    center.label.text = "Test Plasmid",
                    center.label.color = "red",
                    center.label.size = 6)
  expect_s3_class(p, "ggplot")

  # Test with empty data
  empty_data <- test_data[0, ]
  p <- plot_plasmid(empty_data)
  expect_s3_class(p, "ggplot")

  # Test error handling
  expect_error(plot_plasmid(test_data, fill = "nonexistent_column"))
})

# Integration tests
test_that("full pipeline works correctly", {
  test_data <- create_test_features()

  # Test complete workflow
  test_data$level <- calc_level(test_data)$level
  test_data$rstart <- (test_data$start / test_data$length[1]) * 2 * pi
  test_data$rend <- (test_data$end / test_data$length[1]) * 2 * pi

  # Generate glyphs for each feature
  glyphs <- lapply(seq_len(nrow(test_data)), function(i) calc_glyphs(test_data[i, ]))

  # Check glyph generation
  expect_equal(length(glyphs), nrow(test_data))
  expect_true(all(sapply(glyphs, function(g) all(c("x", "y") %in% names(g)))))

  # Create final plot
  p <- plot_plasmid(test_data, fill = "name")
  expect_s3_class(p, "ggplot")

  # Check plot layers
  expect_true(length(p$layers) > 0)
  expect_true(any(sapply(p$layers, function(l) "GeomPolygon" %in% class(l$geom))))
})

# Visual tests (requires manual inspection)
test_that("visual elements are correctly positioned", {
  skip_on_ci()
  test_data <- create_test_features()

  # Generate plot with all visual elements
  p <- plot_plasmid(test_data,
                    fill = "name",
                    feature.outline.color = "black",
                    show.center.label = TRUE,
                    show.labels = TRUE,
                    show.ticks = TRUE)

  # Save plot for manual inspection
  # ggsave("test_plasmid.pdf", p, width = 8, height = 8)

  # Basic checks that can be automated
  expect_true(length(p$layers) >= 4)  # Should have multiple layers
  expect_equal(p$coordinates$limits$x, c(-0.35, 0.35))
  expect_equal(p$coordinates$limits$y, c(-0.35, 0.35))
})

# Performance tests
test_that("functions handle large datasets efficiently", {
  skip_on_ci()

  # Create large dataset
  n_features <- 1000
  large_data <- data.frame(
    start = seq(1, n_features * 10, by = 10),
    end = seq(9, n_features * 10, by = 10),
    length = rep(n_features * 10, n_features),
    direction = sample(c(-1, 0, 1), n_features, replace = TRUE),
    name = paste0("gene", 1:n_features),
    stringsAsFactors = FALSE
  )

  # Test calc_level performance
  time_start <- Sys.time()
  result <- calc_level(large_data)
  time_end <- Sys.time()
  expect_true(as.numeric(difftime(time_end, time_start, units = "secs")) < 5)

  # Test plot_plasmid performance
  time_start <- Sys.time()
  p <- plot_plasmid(large_data[1:100, ])  # Test with first 100 features
  time_end <- Sys.time()
  expect_true(as.numeric(difftime(time_end, time_start, units = "secs")) < 5)
})

# Error handling tests
test_that("functions handle invalid inputs gracefully", {
  # Test invalid angles
  expect_error(text_pos("invalid"))
  expect_error(text_pos(NA))

  # Test invalid data frames
  invalid_data <- data.frame(
    wrong_column = 1:3
  )
  expect_error(calc_level(invalid_data))
  expect_error(plot_plasmid(invalid_data))

  # Test invalid directions
  invalid_direction <- create_test_features()
  invalid_direction$direction <- 2
  expect_error(plot_plasmid(invalid_direction))

  # Test negative positions
  negative_positions <- create_test_features()
  negative_positions$start <- -1
  expect_error(plot_plasmid(negative_positions))
})
