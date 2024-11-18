#' Plot Circular Plasmid Diagram
#'
#' @param data Data frame containing feature information with required columns:
#'   \describe{
#'     \item{start}{Numeric start position}
#'     \item{end}{Numeric end position}
#'     \item{length}{Numeric total length of the plasmid}
#'     \item{direction}{Numeric: 1 for forward, -1 for reverse, 0 for no direction (optional)}
#'   }
#' @param fill Character name of column to use for feature colors (optional)
#' @param feature.outline.color Character color for feature outlines (default: "black")
#' @param feature.outline.linewidth Numeric width for feature outlines (default: 1)
#' @param show.center.label Logical whether to show center label (default: TRUE)
#' @param center.label.text Character override text for center label (optional)
#' @param center.label.color Character color for center label (default: "blue")
#' @param center.label.size Numeric size for center label (default: 4)
#' @param show.labels Logical whether to show feature labels (default: TRUE)
#' @param show.ticks Logical whether to show position markers (default: TRUE)
#'
#' @return A ggplot2 object containing the plasmid visualization
#'
#' @details
#' Creates a circular visualization of plasmid features using ggplot2.
#' Features are drawn as curved polygons arranged in concentric rings to prevent
#' overlap. Optional elements include numerical markers, feature labels, and
#' a central base pair count.
#'
#' @examples
#' file <- system.file("extdata", "plasmid_1.csv", package = "ggplasmid")
#' plasmid_data <- read.csv(file)
#' plot_plasmid(plasmid_data, fill = "feature")
#'
#' @import ggplot2
#' @importFrom rlang inform abort
#' @export
plot_plasmid <- function(data,
                         fill = NULL,
                         feature.outline.color = "black",
                         feature.outline.linewidth = 1,
                         show.center.label = TRUE,
                         center.label.text = NULL,
                         center.label.color = "blue",
                         center.label.size = 4,
                         show.labels = TRUE,
                         show.ticks = TRUE) {

  data <- read_plasmid_data(input = data)

  if (nrow(data) == 0) {
    return(ggplot2::ggplot())
  }

  # Validate direction column if present
  if ("direction" %in% names(data)) {
    invalid_directions <- !data$direction %in% c(-1, 0, 1)
    if (any(invalid_directions)) {
      invalid_values <- unique(data$direction[invalid_directions])
      rlang::abort(
        c(
          "Direction must be one of: -1 (reverse), 0 (none), 1 (forward).",
          "x" = paste("Invalid values found in direction column:", paste(invalid_values, collapse = ", "))
        )
      )
    }
  } else {
    data$direction <- 0
    rlang::inform("No direction column found in data. Features will be drawn without arrows.")
  }

  # Handle fill colors
  if (is.null(fill)) {
    data$default_fill <- "feature"
    fill <- "default_fill"
  } else if (!fill %in% names(data)) {
    rlang::abort(paste("Column", fill, "not found in data frame"))
  }

  # Calculate positions and levels
  plas_len <- data$length[1]
  data <- calc_level(data)
  data$rstart <- (data$start / plas_len) * 2 * pi
  data$rend <- (data$end / plas_len) * 2 * pi

  # Handle wrapping
  wrap_idx <- data$rstart < 0
  data$rstart[wrap_idx] <- data$rstart[wrap_idx] + 2 * pi
  wrap_idx <- data$rend < 0
  data$rend[wrap_idx] <- data$rend[wrap_idx] + 2 * pi
  wrap_idx <- data$rend < data$rstart
  data$rend[wrap_idx] <- data$rend[wrap_idx] + 2 * pi

  # Set fill colors
  data$fill_color <- if (is.null(fill) || fill == "default_fill") "gray50" else data[[fill]]

  # Calculate glyphs
  glyphs <- lapply(seq_len(nrow(data)), function(i) calc_glyphs(data[i, ]))

  # Create base plot
  p <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )

  # Add backbone circle
  theta <- seq(0, 2 * pi, length.out = 100)
  p <- p + ggplot2::geom_path(
    data = data.frame(
      x = BASE_RADIUS * cos(theta),
      y = BASE_RADIUS * sin(theta)
    ),
    ggplot2::aes(x, y),
    color = "black",
    linewidth = 1
  )

  # Add features and annotations
  for (i in seq_along(glyphs)) {
    glyph <- glyphs[[i]]

    # Add feature polygon
    p <- p + ggplot2::geom_polygon(
      data = data.frame(
        x = glyph$x,
        y = glyph$y,
        fill = if (fill == "default_fill") "feature" else data[[fill]][i]
      ),
      ggplot2::aes(x, y, fill = fill),
      color = feature.outline.color,
      linewidth = feature.outline.linewidth
    )

    # Add labels if requested
    if (show.labels) {
      p <- p + ggplot2::geom_path(
        data = data.frame(x = glyph$lineX, y = glyph$lineY),
        ggplot2::aes(x, y),
        color = feature.outline.color,
        alpha = 0.5,
        linewidth = 1
      )

      hjust <- switch(glyph$anno_pos,
        "right" = 0,
        "left" = 1,
        0.5
      )
      vjust <- if (glyph$anno_pos %in% c("t_center", "b_center")) 0.5 else 0

      p <- p + ggplot2::geom_text(
        data = data.frame(
          x = glyph$Lx1,
          y = glyph$Ly1,
          label = data[[fill]][i]
        ),
        ggplot2::aes(x, y, label = label),
        hjust = hjust,
        vjust = vjust,
        size = 3
      )
    }
  }

  # Add tick marks if requested
  if (show.ticks) {
    ticks <- calc_num_markers(plas_len)

    for (i in seq_len(nrow(ticks))) {
      p <- p + ggplot2::geom_path(
        data = data.frame(
          x = unlist(ticks$lineX[i]),
          y = unlist(ticks$lineY[i])
        ),
        ggplot2::aes(x, y),
        color = "black",
        alpha = 0.5,
        linewidth = 0.5
      )

      hjust <- switch(ticks$text_align[i],
        "right" = 0,
        "left" = 1,
        0.5
      )

      p <- p + ggplot2::geom_text(
        data = data.frame(x = ticks$Lx1[i], y = ticks$Ly1[i]),
        ggplot2::aes(x, y),
        label = format(ticks$bp[i], big.mark = ","),
        size = 2.5,
        alpha = 0.5,
        hjust = hjust
      )
    }
  }

  # Add center label if requested
  if (show.center.label) {
    p <- p + ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = if (is.null(center.label.text)) {
        paste0(format(plas_len, big.mark = ","), " bp")
      } else {
        center.label.text
      },
      color = center.label.color,
      size = center.label.size
    )
  }

  # Set plot coordinates and fill scale
  p <- p + ggplot2::coord_fixed(
    xlim = c(-0.35, 0.35),
    ylim = c(-0.35, 0.35)
  )

  if (fill == "default_fill") {
    p <- p + ggplot2::scale_fill_manual(
      values = c("feature" = "gray50"),
      guide = "none"
    )
  }

  p
}
