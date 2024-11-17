plot_plasmid <- function(
    data,
    fill = NULL,
    show.center.label = TRUE,
    center.label.text = NULL,
    center.label.color = "blue",
    center.label.size = 4,
    show.labels = TRUE,
    show.ticks = TRUE) { # Added show.ticks parameter

  if (nrow(data) == 0) {
    return(ggplot())
  }

  # Handle the case where fill is NULL
  if (is.null(fill)) {
    data$default_fill <- "feature" # Add a dummy column for filling
    fill <- "default_fill" # Point fill to the dummy column
  } else if (!fill %in% names(data)) {
    rlang::abort(paste("Column", fill, "not found in data frame"))
  }

  # Calculate levels and positions efficiently
  plas_len <- data$qlen[1]
  data <- calc_level(data)

  # Vectorized position calculations
  data$rstart <- (data$qstart / plas_len) * 2 * pi
  data$rend <- (data$qend / plas_len) * 2 * pi

  # Handle wrapping efficiently
  wrap_idx <- data$rstart < 0
  data$rstart[wrap_idx] <- data$rstart[wrap_idx] + 2 * pi
  wrap_idx <- data$rend < 0
  data$rend[wrap_idx] <- data$rend[wrap_idx] + 2 * pi
  wrap_idx <- data$rend < data$rstart
  data$rend[wrap_idx] <- data$rend[wrap_idx] + 2 * pi

  # Calculate colors efficiently
  if (is.null(fill) || fill == "default_fill") {
    # If fill is NULL or using default_fill, use gray for all features
    data$fill_color <- "gray50"
  } else {
    unique_fills <- unique(data[[fill]])
    colors <- scale_fill_discrete(aesthetics = "fill")$palette(length(unique_fills))
    names(colors) <- unique_fills
    data$fill_color <- colors[data[[fill]]]
  }

  # Calculate glyphs
  glyphs <- lapply(seq_len(nrow(data)), function(i) calc_glyphs(data[i, ]))

  # Create base plot
  p <- ggplot() +
    theme_void() +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank()
    )

  # Add backbone circle efficiently
  theta <- seq(0, 2 * pi, length.out = 100)
  circle_points <- data.frame(
    x = BASE_RADIUS * cos(theta),
    y = BASE_RADIUS * sin(theta)
  )
  p <- p + geom_path(data = circle_points, aes(x, y), color = "black", size = 1)

  # Add features and annotations efficiently
  for (i in seq_along(glyphs)) {
    glyph <- glyphs[[i]]
    # Add feature polygon
    p <- p + geom_polygon(
      data = data.frame(
        x = glyph$x,
        y = glyph$y,
        fill_value = if (fill == "default_fill") "feature" else data[[fill]][i]
      ),
      aes(x, y, fill = fill_value),
      color = data$line_color[i],
      size = 1
    )

    # Add annotation line and label only if show.labels is TRUE
    if (show.labels) {
      # Add annotation line
      p <- p + geom_path(
        data = data.frame(x = glyph$lineX, y = glyph$lineY),
        aes(x, y),
        color = glyph$anno_line_color,
        alpha = 0.5,
        size = 1
      )

      # Add label with optimized positioning
      hjust <- switch(glyph$anno_pos,
                      "right" = 0,
                      "left" = 1,
                      0.5
      )
      vjust <- if (glyph$anno_pos %in% c("t_center", "b_center")) 0.5 else 0
      p <- p + geom_text(
        data = data.frame(
          x = glyph$Lx1,
          y = glyph$Ly1,
          label = data[[fill]][i]
        ),
        aes(x, y, label = label),
        hjust = hjust,
        vjust = vjust,
        size = 3
      )
    }
  }

  # Add BP markers efficiently, only if show.ticks is TRUE
  if (show.ticks) {
    ticks <- calc_num_markers(plas_len)

    # Add marker lines and labels
    for (i in seq_len(nrow(ticks))) {
      p <- p + geom_path(
        data = data.frame(
          x = unlist(ticks$lineX[i]),
          y = unlist(ticks$lineY[i])
        ),
        aes(x, y),
        color = "black",
        alpha = 0.5,
        size = 0.5
      )
      hjust <- switch(ticks$text_align[i],
                      "right" = 0,
                      "left" = 1,
                      0.5
      )
      p <- p + geom_text(
        data = data.frame(x = ticks$Lx1[i], y = ticks$Ly1[i]),
        aes(x, y),
        label = format(ticks$bp[i], big.mark = ","),
        size = 2.5,
        alpha = 0.5,
        hjust = hjust
      )
    }
  }

  # Add center label and finalize plot
  if (show.center.label == TRUE) {
    p <- p +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = ifelse(is.null(center.label.text),
                       paste0(format(plas_len, big.mark = ","), " bp"),
                       center.label.text
        ),
        color = center.label.color,
        size = center.label.size
      )
  }

  # Set fill scale based on whether we're using default gray or custom colors
  if (fill == "default_fill") {
    p <- p + scale_fill_manual(values = c("feature" = "gray50"), guide = "none")
  } else {
    p <- p + scale_fill_discrete(name = fill)
  }

  # Set fixed coordinate limits
  p <- p + coord_fixed(
    xlim = c(-0.35, 0.35),
    ylim = c(-0.35, 0.35)
  )

  p
}
