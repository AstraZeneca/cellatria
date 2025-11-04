# ===================================
# =================================== 
getBaseTheme <- function() {
  # Define universal plot settings
  base_theme <- theme_bw() + 
    theme(text=element_text(size=10),
          plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), units="line")) +
    theme(plot.title=element_text(size=12, hjust=0.5)) +
    theme(panel.border = element_rect(size = 1.25)) +
    theme(plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
    ) +
    #theme(strip.background=element_rect(fill='white')) + 
    theme(strip.background=element_blank(),
          strip.text=element_text(size=11, face = "bold")) +
    theme(axis.title=element_text(size=11, vjust=0.25, face = "bold"),
          axis.text=element_text(size=10, face = "bold")) +
    theme(legend.margin=margin(t=-0.1, r=0, b=0, l=0, unit="line"),
          legend.position='bottom',
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          legend.key.height=grid::unit(10, "points"), 
          legend.key.width=grid::unit(15, "points"),
          legend.spacing=grid::unit(2, "points"))
  
  return(base_theme)
}
# ===================================
# =================================== 
color_helper <- function(set.names, panel="Set1") {
  x <- brewer.pal.info[brewer.pal.info$category == "qual",]
  panel.size <- x$maxcolors[rownames(x) == panel]
  getPalette <- colorRampPalette(brewer.pal(panel.size, panel))  
  colrs <- setNames(getPalette(length(set.names)), 
                    set.names)
  return(colrs)
}
# ===================================
# ===================================
# Helper function to create the projection plot
plot_projection <- function(df, x_col, y_col, feature, clr_map, geomLabelRepel = FALSE, plot_alpha = 0.7) {

# Convert feature to a sorted factor (ensures ordering like 0,1,2, S1, S2, ...)
df <- df %>%
  dplyr::mutate(
    feature = as.character(!!sym(feature)),  # Ensure it's a character before conversion
    feature = factor(feature, levels = mixedsort(unique(feature)), ordered = TRUE)
  ) %>%
  dplyr::arrange(feature)  # Arrange data based on ordered feature

  # Randomize row order to reduce overplotting bias
  df <- df %>% dplyr::slice_sample(n = nrow(df))

  # Compute centroids AFTER ordering feature
  centroids <- df %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(
      x_center = median(!!sym(x_col), na.rm = TRUE),
      y_center = median(!!sym(y_col), na.rm = TRUE)
    )

  p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col), color = feature)) +
    getBaseTheme() +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, face = "bold")
    ) +
    geom_point(size = 0.35, alpha = plot_alpha) +
    scale_color_manual(values = clr_map) +
    guides(
      fill = guide_legend(override.aes = list(size = 4), ncol = 1),
      color = guide_legend(override.aes = list(size = 4), ncol = 1)
    )

  if (geomLabelRepel) {
    p <- p + geom_label_repel(data = centroids, aes(x = x_center, y = y_center, label = feature),
                     size = 3, color = "black")  # Add labels at centroid positions
  } else {
    p <- p + geom_text_repel(data = centroids, aes(x = x_center, y = y_center, label = feature),
                     size = 3, color = "black")  # Add labels at centroid positions
  }

return(p)
}
# ===================================
# ===================================
# Define density plot function
plot_density_ggplot <- function(df, x_col, y_col, feature) {
  ggplot(df, aes_string(x = x_col, y = y_col, color = feature)) +
    getBaseTheme() +
    theme(
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, face = "bold"),
      legend.title = element_blank(),
      legend.position = "right",
      legend.key.height = grid::unit(15, "points"), 
      legend.key.width = grid::unit(5, "points"),
      legend.spacing = grid::unit(2, "points"),
      legend.text = element_text(size = 7, face = "bold")) +
    geom_point(size = 0.35, alpha = 0.8) +  # Adjust size and transparency as needed
    scale_color_viridis_c()  # Continuous color scale
}
# ===================================
# ===================================