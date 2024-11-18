library(ggplot2)
df <- read.csv("inst/extdata/plasmid_eg.csv") |>
  dplyr::select(start:length)
plot_plasmid(data = df,
             show.labels = T,
             show.ticks = T,
             show.center.label = T,
             center.label.text = "XXX",
             center.label.size = 10)

# Create visualization
p <- plot_plasmid("inst/extdata/plasmid_eg.csv",
             fill = "type",
             show.labels = F,
             center.label.text = "",
             center.label.size = 10,
             feature.outline.linewidth = 0.5,
             show.ticks = F) +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")
p
# Save the plot if desired
ggsave("plasmid_visualization.png", p, dpi = 300,
       height = 5, width = 5, device = ragg::agg_png)
