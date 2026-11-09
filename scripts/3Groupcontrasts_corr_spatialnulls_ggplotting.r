# Neuromaps Spatial Correlations-plotting heatmap unified

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(scales)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")  # ensure UTF-8 so the dagger (†) renders, not "..."

out_dir <- "/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/"
df<- read.csv(file.path(out_dir, "neuromaps_results.csv"), stringsAsFactors = FALSE)
cat(sprintf("Loaded %d rows\n", nrow(df)))

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

r_cols      <- grep("^r_", names(df), value = TRUE)
annot_names <- c("genepc1", "cbf", "cbv", "cmrglc_raichle", "cmro2",
                 "SAaxis", "evoexp_xu", "devexp")
annot_labels <- c(
  "genepc1" = "Gene Expression PC1",
  "cbf" = "Cerebral Blood Flow",
  "cbv" = "Cerebral Blood Volume",
  "cmrglc_raichle" = "Glucose Metabolism",
  "cmro2"= "Oxygen Metabolism",
  "SAaxis" = "Sensorimotor-Association Axis",
  "evoexp_xu" = "Evolutionary Expansion",
  "devexp" = "Developmental Expansion(*RH)"
)

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

# # Original no FDR
# plot_df <- bind_rows(plot_rows) %>%
#   mutate(
#     significant = p < 0.05,
#     measure     = factor(measure,    levels = rev(measure_order)),
#     contrast    = factor(contrast,   levels = contrast_order),
#     annotation  = factor(annotation, levels = annot_names)
#   )


# FDR correction - add nominal significance tier
plot_df <- bind_rows(plot_rows) %>%
  group_by(measure) %>%
  mutate(p_fdr = p.adjust(p, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    sig_fdr     = p_fdr < 0.05,
    sig_nominal = p < 0.05 & !sig_fdr,          # nominal only
    measure     = factor(measure,    levels = rev(measure_order)),
    contrast    = factor(contrast,   levels = contrast_order),
    annotation  = factor(annotation, levels = annot_names)
  )

pub_theme <- theme_minimal(base_size = 11) +
  theme(
        panel.grid       = element_blank(),
        panel.border     = element_rect(color = "grey30", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1,
                                        size = 11, color = "black", face = "bold"),
        axis.text.y      = element_text(size = 11, color = "black", face = "bold"),
        axis.ticks       = element_line(color = "grey60", linewidth = 0.3),
        axis.title       = element_blank(),
        legend.title     = element_text(size = 9, face = "bold"),
        legend.text      = element_text(size = 8.5),
        legend.key.width  = unit(0.5, "cm"),
        legend.key.height = unit(2.5, "cm"),
        legend.position  = "right",
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9),
        plot.margin      = margin(8, 6, 8, 6)
  )

r_limit   <- 0.5
pur_scale <- scale_fill_gradientn(
  colors = brewer.pal(11, "PuOr")[2:10],
  limits = c(-r_limit, r_limit),
  oob    = scales::squish,
  name   = "Pearson r",
  breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
  labels = c("−0.5", "−0.25", "0", "0.25", "0.5")
)


make_panel <- function(contrast_id, show_y = TRUE, show_legend = FALSE, show_x = TRUE) {
  
  sub <- plot_df %>% filter(contrast == contrast_id)
  
  ggplot(sub, aes(x = annotation, y = measure)) +
    geom_tile(aes(fill = r), color = "white", linewidth = 0.5) +
    
    # ── FDR-significant: white bg tile + filled colored circle ──
    geom_tile(data = filter(sub, sig_fdr),
              fill = "#F4F4F6", color = "white", linewidth = 0.5) +
    geom_point(data = filter(sub, sig_fdr),
               aes(fill = r),
               shape = 21, size = 12, color = "black", stroke = 1.2) +
    
    # ── Nominal only: white bg tile + dashed-outline circle ──
    geom_tile(data = filter(sub, sig_nominal),
              fill = "#F4F4F6", color = "white", linewidth = 0.5) +
    geom_point(data = filter(sub, sig_nominal),
               aes(fill = r),
               shape = 21, size = 12,
               color  = "white",   # softer outline
               stroke = 0.7) +  
    
    # Non-significant (neither): grey plain text
    geom_text(data = filter(sub, !sig_fdr & !sig_nominal),
              aes(label = sprintf("%.2f", r)),
              size = 3.0, color = "black", fontface = "plain", vjust = 0.5) +
    # Nominal only: italic text with a dagger †
    geom_text(data = filter(sub, sig_nominal),
              aes(label = sprintf("%.2f†", r)),
              size = 3.0, color = "black", fontface = "italic", vjust = 0.5) +
    # FDR significant: bold text with asterisk
    geom_text(data = filter(sub, sig_fdr),
              aes(label = sprintf("%.2f*", r)),
              size = 3.0, color = "black", fontface = "bold", vjust = 0.5) +
    
    geom_hline(yintercept = seq(0.5, length(measure_order) + 0.5, 1), color = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = seq(0.5, length(annot_names) + 0.5, 1), color = "grey85", linewidth = 0.3) +
    pur_scale +
    scale_x_discrete(expand = c(0, 0), labels = if (show_x) annot_labels else NULL) +
    scale_y_discrete(expand = c(0, 0), labels = if (show_y) measure_labels else NULL) +
    labs(title = contrast_labels[[contrast_id]]) +
    pub_theme +
    theme(legend.position = if (show_legend) "right" else "none")
}

p1 <- make_panel("EDR_HC",   show_y = TRUE, show_legend = FALSE, show_x = FALSE)
p2 <- make_panel("EDBP_HC",  show_y = TRUE, show_legend = FALSE, show_x = FALSE)
p3 <- make_panel("EDR_EDBP", show_y = TRUE, show_legend = TRUE,  show_x = TRUE)


legend_only <- get_legend(
  p3 + theme(
    legend.position   = "right",
    legend.key.height = unit(5, "cm"),
    legend.key.width  = unit(0.55, "cm"),
    legend.title      = element_text(size = 11, face = "bold"),
    legend.text       = element_text(size = 11)
  )
)

p1_nl <- p1 + theme(legend.position = "none")
p2_nl <- p2 + theme(legend.position = "none")
p3_nl <- p3 + theme(legend.position = "none")

stacked <- plot_grid(
  p1_nl, p2_nl, p3_nl,
  ncol        = 1,
  rel_heights = c(1, 1, 1.5), 
  align       = "v",
  axis        = "lr"
)


combined <- plot_grid(
  stacked, legend_only,
  ncol       = 2,
  rel_widths = c(1, 0.12)
)


final <- plot_grid(
  combined,
  ncol        = 1
)

n_annot   <- length(annot_names)
n_measure <- length(measure_order)

fig_w <- n_annot * 0.52 + 3.5
fig_h <- n_measure * 0.52 * 3 + 3.0

ggsave(
  filename = file.path(out_dir, "sigcircles_heatmap_neuromaps_acrossmeasuresFDRggplot.png"),
  plot     = final,
  width    = fig_w,
  height   = fig_h,
  units    = "in",
  dpi      = 600
)
