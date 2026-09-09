# Install and load required packages
library(usmap)
library(ggplot2)

# 1. Define the IRI Region Mapping
iri_data <- data.frame(
  state = c(
    "CT", "MA", "ME", "NH", "NJ", "NY", "PA", "RI", "VT",        # Northeast
    "IL", "IN", "MI", "OH", "WI",                               # Great Lakes
    "DE", "DC", "KY", "MD", "NC", "SC", "TN", "VA", "WV",        # Mid-South
    "AL", "FL", "GA", "MS",                                     # Southeast
    "IA", "KS", "MN", "MO", "NE", "ND", "SD",                   # Plains
    "AR", "LA", "OK", "TX",                                     # South Central
    "AZ", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY", # West
    "CA"                                                        # California
  ),
  Region = c(
    rep("Northeast", 9),
    rep("Great Lakes", 5),
    rep("Mid-South", 9),
    rep("Southeast", 4),
    rep("Plains", 7),
    rep("South Central", 4),
    rep("West", 10),
    "California"
  )
)
iri_data$Region = factor(iri_data$Region, levels = c("Northeast", "Mid-South", "West", "South Central", "Plains", "California", "Southeast", "Great Lakes"))
color_region = brewer.pal(8, "Set3")
color_region[6] = "#FFED6F"
# 2. Plot the Map
svg("Figures/IRI_Region_Map_May072026.svg", width = 4, height = 3)
plot_usmap(data = iri_data, values = "Region", color = "white", labels = F, exclude = c("AK","HI")) +
  scale_fill_manual(values=color_region) + 
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    panel.background = element_rect(fill = "white")
  ) 
dev.off()