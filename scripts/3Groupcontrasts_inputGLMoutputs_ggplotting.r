library(png)
library(grid)
library(gridExtra)

img_dir <- "/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/GLM_output/vertexwise_plots"

files <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
fnames <- basename(files)

# Keep only the baseline files (no allsubs, no postpsubs, no x_Age)
fnames_filtered <- fnames[
  !grepl("allsubs", fnames) &
    !grepl("postpsubs", fnames) &
    !grepl("x_Age", fnames)
]

cat("Files after filtering:\n")
print(fnames_filtered)

# Parse from file
# Ex: "FEMA_wrapper_results_vertex_RNI-gm_sm256.matblack bg vertical z=|6| limits RNI-gm_sm256 ~ EDR_HC [zmat] fmax=..."
parse_fname <- function(fname) {
  # Extract measure from "limits MEASURE ~ CONTRAST"
  measure  <- sub(".*limits ([^ ]+) ~ .*", "\\1", fname)
  contrast <- sub(".*limits [^ ]+ ~ ([^ ]+) .*", "\\1", fname)
  
  data.frame(
    filepath = file.path(img_dir, fname),
    measure  = measure,
    contrast = contrast,
    stringsAsFactors = FALSE
  )
}

df <- do.call(rbind, lapply(fnames_filtered, parse_fname))

cat("\nParsed measures:", unique(df$measure), "\n")
cat("Parsed contrasts:", unique(df$contrast), "\n")

# If duplicates exist keep theone with higher fmax
df$fmax <- as.numeric(sub(".*fmax=([0-9e.+]+)\\.png", "\\1", df$filepath))
df <- df[order(df$measure, df$contrast, -df$fmax), ]
df <- df[!duplicated(paste(df$measure, df$contrast)), ]

# ORDER
measure_order <- c("T1w_contrast_sm256","sulc_sm256", "thickness_sm256","area_sm256", 
                   "RND-gm_sm256","RNI-gm_sm256",
                   "RND-wm_sm256","RNI-wm_sm256" )
contrast_order <- c("EDR_HC", "EDBP_HC", "EDR_EDBP")

measure_order <- measure_order[measure_order %in% df$measure]
df$measure  <- factor(df$measure,  levels = measure_order)
df$contrast <- factor(df$contrast, levels = contrast_order)

n_rows <- length(measure_order)
n_cols <- length(contrast_order)


grob_list <- list()
for (m in measure_order) {
  for (c in contrast_order) {
    row <- df[df$measure == m & df$contrast == c, ]
    if (nrow(row) >= 1) {
      img <- readPNG(as.character(row$filepath[1]))
      grob_list[[length(grob_list)+1]] <- rasterGrob(img, interpolate=TRUE)
    } else {
      grob_list[[length(grob_list)+1]] <- textGrob("missing",
                                                   gp=gpar(fontsize=10, col="red"))
    }
  }
}


full_list <- list()
full_list[[1]] <- textGrob("")


for (c in contrast_order) {
  full_list[[length(full_list)+1]] <- textGrob(c, gp=gpar(fontsize=14, fontface="bold"))
}


img_idx <- 1
for (i in seq_along(measure_order)) {
  # remove label cleaning
  label <- gsub("_sm256", "", measure_order[i])
  full_list[[length(full_list)+1]] <- textGrob(label, gp=gpar(fontsize=11, fontface="bold"), rot=90)
  for (j in seq_along(contrast_order)) {
    full_list[[length(full_list)+1]] <- grob_list[[img_idx]]
    img_idx <- img_idx + 1
  }
}

png("/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/GLM_output/vertexwise_plots/GLMoutputs_vertexwise_combined_zmatfigure.png",
    width  = (n_cols + 1) * 600,
    height = (n_rows + 1) * 500,
    res = 150)

grid.arrange(
  grobs   = full_list,
  nrow    = n_rows + 1,
  ncol    = n_cols + 1,
  widths  = c(0.5, rep(4, n_cols)),
  heights = c(0.4, rep(4, n_rows))
)

dev.off()
