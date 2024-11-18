library(ggplot2)
df <- read.csv("inst/extdata/plasmid_1.csv") |>
  dplyr::select(-direction)
plot_plasmid(data = df,
             show.labels = T,
             show.ticks = T,
             show.center.label = T,
             center.label.text = "XXX",
             center.label.size = 10)

# Create visualization
plot_plasmid(df, fill = "type", show.labels = T) +
  #ggsci::scale_fill_npg()+
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Y", fill = "XXX")
p
# Save the plot if desired
ggsave("plasmid_visualization.pdf", p)
