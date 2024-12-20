
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggplasmid

<!-- badges: start -->

[![R-CMD-check](https://github.com/shubhamdutta26/ggplasmid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/shubhamdutta26/ggplasmid/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/shubhamdutta26/ggplasmid/graph/badge.svg)](https://app.codecov.io/gh/shubhamdutta26/ggplasmid)
<!-- badges: end -->

ggplasmid is an R package that provides tools for visualizing plasmid
maps from various data sources. It creates publication-ready, circular
plasmid diagrams with customizable features, labels, and annotations.

## Features

- Supports multiple input formats:
  - Data frames
  - CSV files
  - Excel files (xlsx/xls)
- Generates circular plasmid visualizations with:
  - Customizable feature arrows and colors
  - Automatic feature level calculation to prevent overlapping
  - Position markers with base pair annotations
  - Customizable center labels
  - Flexible label positioning
- Robust data validation and error handling
- Built on ggplot2 for high-quality graphics

## Installation

You can install the development version of ggplasmid from GitHub with:

``` r
# install.packages("pak")
pak::pak("shubhamdutta26/ggplasmid")
```

## Dependencies

Required packages: - ggplot2 - readxl (for Excel file support)

## Usage

### Basic Usage

``` r
library(ggplasmid)
library(ggplot2)

file <- system.file("extdata", "plasmid_eg.csv", package = "ggplasmid")
plasmid_data <- read.csv(file)
plot_plasmid(plasmid_data, fill = "feature")
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
plot_plasmid(plasmid_data, fill = "type", show.labels = FALSE) +
  scale_fill_brewer(palette = "Set3")
```

<img src="man/figures/README-example-2.png" width="100%" />

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
For major changes, please open an issue first to discuss what you would
like to change.

## License

MIT License

## Citation

If you use this package in your research, please cite:

    To be added

## Acknowledgments

This package was inspired by various plasmid visualization tools and
built to provide a flexible, R-based solution for creating
publication-quality plasmid maps.
