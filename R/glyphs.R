#' @title Constants for Plasmid Visualization
#' @description Define basic geometric constants used in plasmid visualization
#' @noRd
#' @keywords internal
BASE_RADIUS <- 0.18
THICKNESS <- 0.017
TEXT_POSITIONS <- list(
  TR_QUARTER = pi / 3,
  TL_QUARTER = 2 * pi / 3,
  BL_QUARTER = 4 * pi / 3,
  BR_QUARTER = 5 * pi / 3
)

#' Calculate Text Position Based on Angle
#'
#' @description
#' Visual explanation of angles and expected returns:
#                    t_center (60° to 120°)
#                         |
#                    left |  right
#                         |
#      left  -------------+------------- right
#                         |
#                    left |  right
#                         |
#                    b_center (240° to 300°)
#'
#' @param theta Numeric angle in radians
#' @param pos Character string indicating position type ("outer" or "inner")
#'
#' @return Character string indicating text alignment position ("b_center", "t_center", "right", or "left")
#' @noRd
#' @details
#' Determines the optimal text positioning based on the angle around the plasmid circle.
#' Used for both feature labels and numerical markers.
#' @keywords internal
text_pos <- function(theta, pos = "outer") {
  # Simplified text position calculation
  theta <- (if (pos == "inner") (theta - pi) else theta) %% (2 * pi)

  if (theta >= TEXT_POSITIONS$BL_QUARTER && theta <= TEXT_POSITIONS$BR_QUARTER) {
    return("b_center")
  }
  if (theta >= TEXT_POSITIONS$TR_QUARTER && theta <= TEXT_POSITIONS$TL_QUARTER) {
    return("t_center")
  }
  if (theta <= TEXT_POSITIONS$TR_QUARTER || theta >= TEXT_POSITIONS$BR_QUARTER) {
    return("right")
  }
  return("left")
}

#' Calculate Glyph Coordinates for Plasmid Features
#'
#' @param row A data frame row containing feature information with the following required columns:
#'   \describe{
#'     \item{level}{Integer indicating the feature's level (ring position)}
#'     \item{direction}{Numeric: 1 for clockwise, -1 for anti-clockwise, 0 or NULL for no direction}
#'     \item{rstart}{Numeric start position in radians}
#'     \item{rend}{Numeric end position in radians}
#'     \item{line_color}{Character string specifying line color (optional)}
#'     \item{fill_color}{Character string specifying fill color (optional)}
#'   }
#'
#' @return A list containing coordinates and styling information:
#'   \describe{
#'     \item{x,y}{Vectors of x,y coordinates for the feature polygon}
#'     \item{Lx1,Ly1}{Text label coordinates}
#'     \item{anno_line_color}{Color for annotation lines}
#'     \item{lineX,lineY}{Coordinates for annotation lines}
#'     \item{theta}{Angle for text placement}
#'     \item{anno_pos}{Text position indicator}
#'   }
#'
#' @details
#' Calculates all coordinates needed to draw a single plasmid feature, including
#' the main glyph shape, annotation lines, and text placement positions.
#' Handles directional arrows if direction is specified.
#' @noRd
#' @keywords internal
calc_glyphs <- function(row) {
  feat_radius <- BASE_RADIUS + THICKNESS * 2.3 * row$level
  shift <- pi / 2
  has_direction <- !is.null(row$direction) && abs(row$direction) == 1

  # Calculate start and end positions
  r_start <- row$rstart
  r_end <- row$rend

  # Calculate segment points
  seg_len <- abs(r_end - r_start)
  n_segments <- max(as.integer(25 * seg_len) + 3, 10)

  # For clockwise (direction = 1), we go from start to end
  # For anti-clockwise (direction = -1), we go from end to start
  if (!has_direction || row$direction == 1) {
    theta <- seq(shift - r_start, shift - r_end, length.out = n_segments)
  } else {
    # For anti-clockwise, reverse the sequence
    theta <- rev(seq(shift - r_start, shift - r_end, length.out = n_segments))
  }

  # Calculate radii points
  radius_outer <- feat_radius + THICKNESS
  radius_inner <- feat_radius - THICKNESS

  # Generate coordinates
  x1 <- radius_outer * cos(theta)
  y1 <- radius_outer * sin(theta)
  x2 <- radius_inner * cos(rev(theta))
  y2 <- radius_inner * sin(rev(theta))

  # Handle arrows if needed
  if (has_direction) {
    if (row$direction == 1) { # Clockwise
      # Arrow points toward r_end
      arrow_theta <- shift - r_end
      x1 <- c(x1[1:(length(x1) - 2)], feat_radius * cos(arrow_theta))
      y1 <- c(y1[1:(length(y1) - 2)], feat_radius * sin(arrow_theta))
    } else { # Anti-clockwise
      # Arrow points toward r_start
      arrow_theta <- shift - r_start
      x1 <- c(x1[1:(length(x1) - 2)], feat_radius * cos(arrow_theta))
      y1 <- c(y1[1:(length(y1) - 2)], feat_radius * sin(arrow_theta))
    }
    x2 <- x2[3:length(x2)]
    y2 <- y2[3:length(y2)]
  }

  # Calculate text and line positions
  theta_mean <- pi / 2 - mean(c(r_start, r_end))
  cos_theta <- cos(theta_mean)
  sin_theta <- sin(theta_mean)
  text_radius <- feat_radius * 1.6
  line_radius <- feat_radius * 1.5

  # Determine colors
  line_color <- if (!is.null(row$line_color)) row$line_color else "black"
  fill_color <- if (!is.null(row$fill_color)) row$fill_color else "#ffffff"
  anno_line_color <- if (is.null(fill_color) || fill_color == "#ffffff") line_color else fill_color

  list(
    x = c(x1, x2),
    y = c(y1, y2),
    Lx1 = cos_theta * text_radius,
    Ly1 = sin_theta * text_radius,
    anno_line_color = anno_line_color,
    lineX = c(cos_theta * (feat_radius + THICKNESS), cos_theta * line_radius),
    lineY = c(sin_theta * (feat_radius + THICKNESS), sin_theta * line_radius),
    theta = theta_mean,
    anno_pos = text_pos(theta_mean)
  )
}

#' Calculate Numerical Marker Positions
#'
#' @param plas_len Numeric value indicating total plasmid length in base pairs
#'
#' @return A data frame containing marker information:
#'   \describe{
#'     \item{lineX,lineY}{List columns with coordinates for marker lines}
#'     \item{theta}{Numeric vector of angles for text placement}
#'     \item{text_align}{Character vector of text alignment positions}
#'     \item{bp}{Numeric vector of base pair positions}
#'     \item{Lx1,Ly1}{Numeric vectors of text label coordinates}
#'     \item{size}{Character vector of text sizes}
#'   }
#'
#' @details
#' Calculates positions for numerical markers around the plasmid circle.
#' Markers are placed at regular intervals, with spacing adjusted based on
#' plasmid size.
#' @noRd
#' @keywords internal
calc_num_markers <- function(plas_len) {
  chunk_size <- max(round((plas_len %/% 5) / 500) * 500, 500)
  chunks <- seq(0, plas_len - chunk_size / 2, by = chunk_size)
  chunks <- chunks[chunks < plas_len]
  chunks[chunks == 0] <- 1

  # Calculate positions
  chunksR <- (chunks / plas_len) * 2 * pi
  theta <- pi / 2 - chunksR

  # Calculate radii
  offset <- 0.155
  long_radius <- offset / 1.08
  text_radius <- long_radius * 0.85

  # Generate coordinates
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)

  data.frame(
    lineX = I(Map(c, cos_theta * offset, cos_theta * long_radius)),
    lineY = I(Map(c, sin_theta * offset, sin_theta * long_radius)),
    theta = theta,
    text_align = vapply(theta, text_pos, character(1), pos = "inner"),
    bp = chunks,
    Lx1 = cos_theta * text_radius,
    Ly1 = sin_theta * text_radius,
    size = "12px"
  )
}

#' Calculate Non-overlapping Feature Levels
#'
#' @param annotations Data frame containing feature annotations with columns:
#'   \describe{
#'     \item{start}{Numeric start position}
#'     \item{end}{Numeric end position}
#'     \item{length}{Numeric total length of the plasmid}
#'   }
#'
#' @return Data frame with original annotations plus a new 'level' column
#'   indicating the assigned level (ring) for each feature
#'
#' @details
#' Assigns features to different levels (concentric rings) to prevent overlap.
#' Features are sorted by size and assigned to the innermost available level
#' where they don't overlap with existing features.
#' @noRd
#' @keywords internal
calc_level <- function(annotations) {
  if (nrow(annotations) == 0) {
    return(annotations)
  }

  # Calculate feature sizes and handle wrapping
  feature_sizes <- annotations$end - annotations$start
  wrap_idx <- feature_sizes < 0
  feature_sizes[wrap_idx] <- feature_sizes[wrap_idx] + annotations$length[wrap_idx]

  # Sort by feature size
  annotations <- annotations[order(-feature_sizes), ]

  # Initialize levels
  n_features <- nrow(annotations)
  levels <- integer(n_features)

  if (n_features > 1) {
    for (i in 2:n_features) {
      s <- annotations$start[i]
      e <- annotations$end[i]

      if (s >= e) s <- s - annotations$length[i]

      s_j <- annotations$start[1:(i - 1)]
      e_j <- annotations$end[1:(i - 1)]
      idx <- s_j >= e_j
      s_j[idx] <- s_j[idx] - annotations$length[idx]

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
