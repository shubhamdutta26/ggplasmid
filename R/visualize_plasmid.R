BASE_RADIUS <- 0.18

text_pos <- function(theta, pos = "outer") {
  if (pos == "inner") {
    theta <- theta - pi
  }
  if (theta < 0) {
    theta <- theta + 2 * pi
  }

  trQ <- pi / 3
  tlQ <- 2 * pi / 3
  blQ <- 4 * pi / 3
  brQ <- 5 * pi / 3

  if (theta >= blQ && theta <= brQ) {
    return("b_center")
  } else if (theta >= trQ && theta <= tlQ) {
    return("t_center")
  } else if (theta <= trQ || theta >= brQ) {
    return("right")
  } else {
    return("left")
  }
}

calc_glyphs <- function(row) {
  r1 <- row$rend
  r2 <- row$rstart
  frame <- row$sframe
  level <- row$level

  THICKNESS <- 0.017
  feat_radius <- BASE_RADIUS
  feat_radius <- feat_radius + THICKNESS * 2.3 * level

  seg_len <- r1 - r2
  if (frame == 1) {  # reverses for direction
    tmp <- r1
    r1 <- r2
    r2 <- tmp
  }

  shift <- pi / 2  # corrects for starting at the correct polar space
  n_segments <- as.integer(25 * seg_len) + 3  # number of lines/sampling size
  theta <- seq(shift - r1, shift - r2, length.out = n_segments)

  x1 <- (feat_radius + THICKNESS) * cos(theta)
  y1 <- (feat_radius + THICKNESS) * sin(theta)

  line_theta_avg <- mean(c(r1, r2))

  has_orientation <- !is.null(row$sframe) && abs(row$sframe) == 1

  if (has_orientation) {
    x1 <- x1[1:(length(x1)-2)]
    y1 <- y1[1:(length(y1)-2)]

    line_theta_avg <- atan2(mean(x1), mean(y1))

    x1 <- c(x1, feat_radius * cos(shift - r2))
    y1 <- c(y1, feat_radius * sin(shift - r2))
  }

  x2 <- (feat_radius - THICKNESS) * cos(rev(theta))
  y2 <- (feat_radius - THICKNESS) * sin(rev(theta))

  if (has_orientation) {
    x2 <- x2[3:length(x2)]
    y2 <- y2[3:length(y2)]
  }

  x <- c(x1, x2)
  y <- c(y1, y2)

  theta <- (pi / 2) - line_theta_avg

  # Calculate line start position
  Lx0 <- cos(theta) * (feat_radius + THICKNESS)
  Ly0 <- sin(theta) * (feat_radius + THICKNESS)

  # Calculate line end position
  long_radius <- feat_radius * 1.5

  # Calculate text position even further out
  text_radius <- feat_radius * 1.6

  # Line endpoint
  Lx1 <- cos(theta) * long_radius
  Ly1 <- sin(theta) * long_radius

  # Text position
  text_x <- cos(theta) * text_radius
  text_y <- sin(theta) * text_radius

  lineX <- c(Lx0, Lx1)
  lineY <- c(Ly0, Ly1)

  anno_line_color <- row$fill_color
  if (is.null(anno_line_color) || anno_line_color == "#ffffff") {
    anno_line_color <- row$line_color
  }

  anno_pos <- text_pos(theta)

  list(
    x = x,
    y = y,
    Lx1 = text_x,
    Ly1 = text_y,
    anno_line_color = anno_line_color,
    lineX = lineX,
    lineY = lineY,
    theta = theta,
    anno_pos = anno_pos
  )
}

calc_num_markers <- function(plas_len) {
  chunk_size <- round((plas_len %/% 5) / 500) * 500
  if (chunk_size == 0) {
    chunk_size <- 500
  }

  chunks <- seq(0, plas_len - as.integer(chunk_size / 2), by = as.integer(chunk_size))
  chunks <- chunks[chunks < plas_len]
  chunks[chunks == 0] <- 1

  chunksR <- (chunks / plas_len) * 2 * pi
  theta <- (pi / 2) - chunksR

  offset <- 0.155
  Lx0 <- cos(theta) * offset
  Ly0 <- sin(theta) * offset
  longRadius <- offset / 1.08
  text_radius <- longRadius * 0.85

  Lx1 <- cos(theta) * longRadius
  Ly1 <- sin(theta) * longRadius

  text_x <- cos(theta) * text_radius
  text_y <- sin(theta) * text_radius

  lineX <- Map(c, Lx0, Lx1)
  lineY <- Map(c, Ly0, Ly1)

  ticks <- data.frame(
    lineX = I(lineX),
    lineY = I(lineY),
    theta = theta
  )

  ticks$text_align <- sapply(theta, text_pos, pos = "inner")
  ticks$bp <- chunks
  ticks$Lx1 <- text_x
  ticks$Ly1 <- text_y
  ticks$size <- "12px"

  return(ticks)
}

calc_level <- function(annotations) {
  if (nrow(annotations) == 0) return(annotations)

  # Sort features by size (larger features get priority)
  feature_sizes <- annotations$qend - annotations$qstart
  feature_sizes[feature_sizes < 0] <- feature_sizes[feature_sizes < 0] + annotations$qlen[feature_sizes < 0]
  annotations <- annotations[order(-feature_sizes),]

  annotations$level <- NA
  annotations$level[1] <- 0  # First feature always at level 0

  if (nrow(annotations) > 1) {
    for (i in 2:nrow(annotations)) {
      s <- annotations$qstart[i]
      e <- annotations$qend[i]

      if (s >= e) {
        s <- s - annotations$qlen[i]
      }

      # Find overlapping features
      overlaps <- sapply(1:(i-1), function(j) {
        s_j <- annotations$qstart[j]
        e_j <- annotations$qend[j]

        if (s_j >= e_j) {
          s_j <- s_j - annotations$qlen[j]
        }

        !(e <= s_j || s >= e_j)
      })

      if (!any(overlaps)) {
        annotations$level[i] <- 0
      } else {
        occupied_levels <- annotations$level[overlaps]
        new_level <- 0
        while (new_level %in% occupied_levels) {
          new_level <- new_level + 1
        }
        annotations$level[i] <- new_level
      }
    }
  }

  annotations$level[is.na(annotations$level)] <- 3
  return(annotations)
}

visualize_plasmid <- function(df, linear = FALSE) {
  if (nrow(df) == 0) {
    return(ggplot())
  }

  # Get plasmid length from first row
  plas_len <- df$qlen[1]

  # Calculate levels for features
  df <- calc_level(df)

  # Convert positions to radians
  df$rstart <- (df$qstart / plas_len) * 2 * pi
  df$rend <- (df$qend / plas_len) * 2 * pi
  df$rstart[df$rstart < 0] <- df$rstart[df$rstart < 0] + 2 * pi
  df$rend[df$rend < 0] <- df$rend[df$rend < 0] + 2 * pi
  df$rend[df$rend < df$rstart] <- df$rend[df$rend < df$rstart] + 2 * pi

  # Calculate glyphs for each feature
  glyphs <- lapply(1:nrow(df), function(i) calc_glyphs(df[i,]))

  # Create base plot
  p <- ggplot() +
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.position = "bottom"
    )

  # Add backbone circle
  circle_points <- data.frame(
    x = BASE_RADIUS * cos(seq(0, 2*pi, length.out = 100)),
    y = BASE_RADIUS * sin(seq(0, 2*pi, length.out = 100))
  )
  p <- p + geom_path(data = circle_points, aes(x, y), color = "black", size = 1)

  # Add features
  for (i in 1:length(glyphs)) {
    glyph <- glyphs[[i]]
    feature_data <- data.frame(
      x = glyph$x,
      y = glyph$y
    )

    p <- p + geom_polygon(
      data = feature_data,
      aes(x, y),
      fill = df$fill_color[i],
      color = df$line_color[i],
      size = 1
    )

    # Add annotation lines
    line_data <- data.frame(
      x = glyph$lineX,
      y = glyph$lineY
    )
    p <- p + geom_path(
      data = line_data,
      aes(x, y),
      color = glyph$anno_line_color,
      alpha = 0.5,
      size = 1
    )

    # Add labels
    text_data <- data.frame(
      x = glyph$Lx1,
      y = glyph$Ly1,
      label = df$Feature[i]
    )

    hjust <- switch(glyph$anno_pos,
                    "right" = 0,
                    "left" = 1,
                    0.5
    )
    vjust <- if(glyph$anno_pos %in% c("t_center", "b_center")) 0.5 else 0

    p <- p + geom_text(
      data = text_data,
      aes(x, y, label = label),
      hjust = hjust,
      vjust = vjust,
      size = 3
    )
  }

  # Calculate and add bp markers
  ticks <- calc_num_markers(plas_len)

  for (i in 1:nrow(ticks)) {
    # Add marker lines
    marker_line <- data.frame(
      x = unlist(ticks$lineX[i]),
      y = unlist(ticks$lineY[i])
    )
    p <- p + geom_path(
      data = marker_line,
      aes(x, y),
      color = "black",
      alpha = 0.5,
      size = 0.5
    )

    # Add bp labels with correct position value
    p <- p + geom_text(
      data = data.frame(x = ticks$Lx1[i], y = ticks$Ly1[i]),
      aes(x, y),
      label = format(ticks$bp[i], big.mark = ","),
      size = 2.5,
      alpha = 0.5,
      hjust = switch(ticks$text_align[i],
                     "right" = 0,
                     "left" = 1,
                     0.5
      )
    )
  }

  # Add total bp label in center
  p <- p + annotate(
    "text",
    x = 0,
    y = -0.02,
    label = paste0(format(plas_len, big.mark = ","), " bp"),
    color = "#7b7b7b",
    size = 4
  )

  # Set plot limits
  plot_size <- 0.35
  p <- p + coord_fixed(
    xlim = c(-plot_size, plot_size),
    ylim = c(-plot_size, plot_size)
  )

  return(p)
}
