BASE_RADIUS <- 0.18
THICKNESS <- 0.017 # Moved constant out of function

# Optimize text position calculation with pre-calculated constants
TEXT_POSITIONS <- list(
  TR_QUARTER = pi / 3,
  TL_QUARTER = 2 * pi / 3,
  BL_QUARTER = 4 * pi / 3,
  BR_QUARTER = 5 * pi / 3
)

text_pos <- function(theta, pos = "outer") {
  # Vectorized position calculation
  theta <- if (pos == "inner") (theta - pi) %% (2 * pi) else theta %% (2 * pi)

  # Vectorized position determination
  ifelse(
    theta >= TEXT_POSITIONS$BL_QUARTER & theta <= TEXT_POSITIONS$BR_QUARTER,
    "b_center",
    ifelse(
      theta >= TEXT_POSITIONS$TR_QUARTER & theta <= TEXT_POSITIONS$TL_QUARTER,
      "t_center",
      ifelse(
        theta <= TEXT_POSITIONS$TR_QUARTER | theta >= TEXT_POSITIONS$BR_QUARTER,
        "right",
        "left"
      )
    )
  )
}

calc_glyphs <- function(row) {
  # Pre-calculate common values
  feat_radius <- BASE_RADIUS + THICKNESS * 2.3 * row$level
  shift <- pi / 2

  # Determine segment direction
  seg_len <- abs(row$rend - row$rstart)
  if (row$sframe == 1) {
    r_start <- row$rend
    r_end <- row$rstart
  } else {
    r_start <- row$rstart
    r_end <- row$rend
  }

  # Optimize segment calculation
  n_segments <- max(as.integer(25 * seg_len) + 3, 10) # Added minimum segments
  theta <- seq(shift - r_start, shift - r_end, length.out = n_segments)

  # Vectorized coordinate calculations
  radius_outer <- feat_radius + THICKNESS
  radius_inner <- feat_radius - THICKNESS

  x1 <- radius_outer * cos(theta)
  y1 <- radius_outer * sin(theta)
  x2 <- radius_inner * cos(rev(theta))
  y2 <- radius_inner * sin(rev(theta))

  # Handle feature orientation
  has_orientation <- !is.null(row$sframe) &&
    abs(row$sframe) == 1 &&
    !is.null(row$type) &&
    row$type %in% c("CDS", "promoter")

  if (has_orientation) {
    # Efficient vector operations for arrow
    x1 <- x1[1:(length(x1) - 2)]
    y1 <- y1[1:(length(y1) - 2)]
    x1 <- c(x1, feat_radius * cos(shift - r_end))
    y1 <- c(y1, feat_radius * sin(shift - r_end))
    x2 <- x2[3:length(x2)]
    y2 <- y2[3:length(y2)]
  }

  # Calculate text and line positions efficiently
  line_theta_avg <- mean(c(r_start, r_end))
  theta <- (pi / 2) - line_theta_avg

  # Pre-calculate common trig values
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)

  # Calculate positions using pre-computed values
  text_radius <- feat_radius * 1.6
  text_x <- cos_theta * text_radius
  text_y <- sin_theta * text_radius

  # Efficient line coordinate calculation
  long_radius <- feat_radius * 1.5
  lineX <- c(
    cos_theta * (feat_radius + THICKNESS),
    cos_theta * long_radius
  )
  lineY <- c(
    sin_theta * (feat_radius + THICKNESS),
    sin_theta * long_radius
  )

  # Determine annotation color efficiently
  anno_line_color <- if (is.null(row$fill_color) || row$fill_color == "#ffffff") {
    row$line_color
  } else {
    row$fill_color
  }

  # Return results as a list
  list(
    x = c(x1, x2),
    y = c(y1, y2),
    Lx1 = text_x,
    Ly1 = text_y,
    anno_line_color = anno_line_color,
    lineX = lineX,
    lineY = lineY,
    theta = theta,
    anno_pos = text_pos(theta)
  )
}

calc_num_markers <- function(plas_len) {
  # Optimize chunk size calculation
  chunk_size <- max(round((plas_len %/% 5) / 500) * 500, 500)

  # Vectorized sequence generation
  chunks <- seq(0, plas_len - as.integer(chunk_size / 2),
    by = as.integer(chunk_size)
  )
  chunks <- chunks[chunks < plas_len]
  chunks[chunks == 0] <- 1

  # Vectorized calculations
  chunksR <- (chunks / plas_len) * 2 * pi
  theta <- (pi / 2) - chunksR

  # Pre-calculate constants
  offset <- 0.155
  longRadius <- offset / 1.08
  text_radius <- longRadius * 0.85

  # Vectorized coordinate calculations
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)

  # Create data frame efficiently
  data.frame(
    lineX = I(Map(c, cos_theta * offset, cos_theta * longRadius)),
    lineY = I(Map(c, sin_theta * offset, sin_theta * longRadius)),
    theta = theta,
    text_align = vapply(theta, text_pos, character(1), pos = "inner"),
    bp = chunks,
    Lx1 = cos_theta * text_radius,
    Ly1 = sin_theta * text_radius,
    size = "12px"
  )
}

calc_level <- function(annotations) {
  if (nrow(annotations) == 0) {
    return(annotations)
  }

  # Calculate feature sizes efficiently
  feature_sizes <- annotations$qend - annotations$qstart
  wrap_idx <- feature_sizes < 0
  feature_sizes[wrap_idx] <- feature_sizes[wrap_idx] + annotations$qlen[wrap_idx]

  # Sort features by size once
  annotations <- annotations[order(-feature_sizes), ]

  # Initialize levels vector
  n_features <- nrow(annotations)
  levels <- integer(n_features)

  # First feature always at level 0
  if (n_features > 1) {
    # Process remaining features
    for (i in 2:n_features) {
      s <- annotations$qstart[i]
      e <- annotations$qend[i]

      if (s >= e) {
        s <- s - annotations$qlen[i]
      }

      # Vectorized overlap check
      s_j <- annotations$qstart[1:(i - 1)]
      e_j <- annotations$qend[1:(i - 1)]
      idx <- s_j >= e_j
      s_j[idx] <- s_j[idx] - annotations$qlen[idx]

      overlaps <- !(e <= s_j | s >= e_j)

      if (!any(overlaps)) {
        levels[i] <- 0
      } else {
        occupied <- levels[1:(i - 1)][overlaps]
        levels[i] <- min(setdiff(0:max(occupied), occupied))
      }
    }
  }

  annotations$level <- levels
  annotations
}
