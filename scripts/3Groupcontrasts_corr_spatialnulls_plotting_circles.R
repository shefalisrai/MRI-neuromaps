# Neuromaps Spatial Correlations -plotting circles

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(scales)

# ── 0. Paths ──────────────────────────────────────────────────────────────────
out_dir <- "/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/"
df      <- read.csv(file.path(out_dir, "neuromaps_results.csv"), stringsAsFactors = FALSE)
cat(sprintf("Loaded %d rows\n", nrow(df)))

# ── 1. Metadata ───────────────────────────────────────────────────────────────
contrast_order <- c("EDR_HC", "EDBP_HC", "EDR_EDBP")

contrast_labels <- c(
  EDR_HC   = "Restricting vs. Healthy Controls",
  EDBP_HC  = "Binge/Purging vs. Healthy Controls",
  EDR_EDBP = "Restricting vs. Binge/Purging"
)

measure_order <- c("T1w_contrast", "sulc", "thickness", "area",
                   "RND-gm", "RNI-gm", "RND-wm", "RNI-wm")

measure_labels <- c(
  "T1w_contrast" = "T1w Contrast",
  "sulc"         = "Sulcal Depth",
  "thickness"    = "Cortical Thickness",
  "area"         = "Surface Area",
  "RND-gm"       = "RND (GM)",
  "RNI-gm"       = "RNI (GM)",
  "RND-wm"       = "RND (WM)",
  "RNI-wm"       = "RNI (WM)"
)

# ── 2. Detect annotations & build long data ───────────────────────────────────
r_cols      <- grep("^r_", names(df), value = TRUE)
annot_names <- sub("^r_", "", r_cols)

# Remove near-zero annotation columns
annot_names <- annot_names[!annot_names %in% c("sert", "5ht2a", "5ht1b", "5ht4", "gabaa", "cbf","cbv", "fcgradient01","fcgradient02","fcgradient03","FChomology")]

plot_rows <- lapply(annot_names, function(name) {
  data.frame(
    annotation = name,
    measure    = df$measure,
    contrast   = df$contrast,
    r          = df[[paste0("r_", name)]],
    p          = df[[paste0("p_", name)]],
    stringsAsFactors = FALSE
  )
})

plot_df <- bind_rows(plot_rows) %>%
  mutate(
    significant = p < 0.05,
    measure     = factor(measure,    levels = rev(measure_order)),  # top → bottom
    contrast    = factor(contrast,   levels = contrast_order),
    annotation  = factor(annotation, levels = annot_names)
  )

# ── 3. Shared theme ───────────────────────────────────────────────────────────
pub_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid          = element_blank(),
    panel.border        = element_rect(color = "grey30", fill = NA, linewidth = 0.5),
    panel.background    = element_rect(fill = "white", color = NA),
    plot.background     = element_rect(fill = "white", color = NA),
    # X-axis: bigger, bold, angled
    axis.text.x         = element_text(angle = 45, hjust = 1, vjust = 1,
                                       size = 8.5, color = "grey10",
                                       face = "bold"),
    # Y-axis: bigger
    axis.text.y         = element_text(size = 8.5, color = "grey10"),
    axis.ticks          = element_line(color = "grey60", linewidth = 0.3),
    axis.title          = element_blank(),
    strip.text          = element_text(face = "bold", size = 11),
    # Legend
    legend.title        = element_text(size = 9, face = "bold"),
    legend.text         = element_text(size = 8.5),
    legend.key.width    = unit(0.5, "cm"),
    legend.key.height   = unit(2.5, "cm"),
    legend.position     = "right",
    # Panel title: centered, larger
    plot.title          = element_text(face = "bold", size = 11, hjust = 0.5,
                                       margin = margin(b = 6)),
    plot.title.position = "plot",
    plot.margin         = margin(8, 6, 8, 6)
  )

# COLOR
r_limit   <- 0.5
pur_scale <- scale_fill_gradientn(
  colors = c("#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#67A9CF", "#2166AC"),
  limits = c(-r_limit, r_limit),
  oob    = scales::squish,
  name   = "Pearson r",
  breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
  labels = c("−0.5", "−0.25", "0", "0.25", "0.5")
)

# ── 5. Panel builder ──────────────────────────────────────────────────────────
make_panel <- function(contrast_id, show_y = TRUE, show_legend = FALSE) {
  
  sub <- plot_df %>% filter(contrast == contrast_id)
  
  ggplot(sub, aes(x = annotation, y = measure, fill = r, color = r)) +
    geom_tile(fill = "white", color = "grey90", linewidth = 0.4) +
    # Non-significant: just bold r text, colored by r value, get rid of close to r=0
    geom_text(data = filter(sub, !significant & abs(r) >= 0.15),
              aes(label = sprintf("%.2f", r), color = r),
              size = 3.0, fontface = "bold") +
    # Significant: filled circle sized by |r|, colored by r value  
    geom_point(data = filter(sub, significant),
               aes(fill = r, size = abs(r)),
               shape = 21, color = "grey20", stroke = 0.4) +
    # r value text on top of circles
    geom_text(data = filter(sub, significant),
              aes(label = sprintf("%.2f", r)),
              size = 2.5, color = "white", fontface = "bold") +
    scale_color_gradientn(
      colors = c("#B2182B","#EF8A62","#FDDBC7","#F7F7F7","#D1E5F0","#67A9CF","#2166AC"),
      limits = c(-r_limit, r_limit), guide = "none") +
    scale_size_continuous(range = c(3, 12), guide = "none") +
    pur_scale +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(
      expand = c(0, 0),
      labels = if (show_y) measure_labels else NULL
    ) +
    labs(title = contrast_labels[[contrast_id]]) +
    pub_theme +
    theme(legend.position = if (show_legend) "right" else "none")
}

# ── 6. Build panels ───────────────────────────────────────────────────────────
p1 <- make_panel("EDR_HC",   show_y = TRUE, show_legend = FALSE)
p2 <- make_panel("EDBP_HC",  show_y = TRUE, show_legend = FALSE)
p3 <- make_panel("EDR_EDBP", show_y = TRUE, show_legend = TRUE)

# ── 7. Extract legend, then strip from all panels ────────────────────────────
legend_only <- get_legend(
  p3 + theme(
    legend.position   = "right",
    legend.key.height = unit(5, "cm"),
    legend.key.width  = unit(0.55, "cm"),
    legend.title      = element_text(size = 10, face = "bold"),
    legend.text       = element_text(size = 9)
  )
)

p1_nl <- p1 + theme(legend.position = "none")
p2_nl <- p2 + theme(legend.position = "none")
p3_nl <- p3 + theme(legend.position = "none")

# ── 8. Stack panels vertically ───────────────────────────────────────────────
stacked <- plot_grid(
  p1_nl, p2_nl, p3_nl,
  ncol        = 1,
  rel_heights = c(1, 1, 1),
  align       = "v",
  axis        = "lr"
)

# ── 9. Attach legend to far right of full stack ───────────────────────────────
combined <- plot_grid(
  stacked, legend_only,
  ncol       = 2,
  rel_widths = c(1, 0.07)
)

# ── 10. Add overall title + subtitle (subtitle now visible & larger) ──────────
title_row <- ggdraw() +
  draw_label("Spatial Correlations",
             fontface = "bold", size = 15, hjust = 0.5, x = 0.5)

subtitle_row <- ggdraw() +
  draw_label("* p < 0.05, spin test",
             size = 10.5, color = "grey35", hjust = 0.5, x = 0.5,
             fontface = "italic")

final <- plot_grid(
  title_row, subtitle_row, combined,
  ncol        = 1,
  rel_heights = c(0.045, 0.03, 1)
)

# ── 11. Save ──────────────────────────────────────────────────────────────────
n_annot   <- length(annot_names)
n_measure <- length(measure_order)

fig_w <- n_annot * 0.52 + 3.5
fig_h <- n_measure * 0.52 * 3 + 3.0

ggsave(
  filename = file.path(out_dir, "unified_heatmap_neuromaps_ggplot_circles.png"),
  plot     = final,
  width    = fig_w,
  height   = fig_h,
  units    = "in",
  dpi      = 600
)
