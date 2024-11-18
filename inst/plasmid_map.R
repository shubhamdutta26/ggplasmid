library(ggplot2)
df <- read.csv("inst/extdata/plasmid_1.csv") |>
  dplyr::select(start:length)
plot_plasmid(data = df,
             show.labels = T,
             show.ticks = T,
             show.center.label = T,
             center.label.text = "XXX",
             center.label.size = 10)

# Create visualization
plot_plasmid("inst/extdata/pSC101.csv",
             fill = "feature",
             show.labels = T)
p
# Save the plot if desired
ggsave("plasmid_visualization.pdf", p)
