# ----- Loadings Heatmap Comparison (2x2 layout) -----
# Panels: Original | Mean
#         1x OOB   | 2x OOB
# Rows: traits (aboveground lowercase, belowground uppercase)
# Cols: PC1..PC4
# Cell: loading value; color by value; bold = largest |loading| in that trait/method
#
# Output: Figures/Loadings_heatmap.png
#
# Prerequisites: scripts 06, 061, 062, 063 (aggregated loadings summaries exist)

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
PCA_original <- readRDS("data/PCATotal_ImputedObs.rds")
PCA_mean     <- readRDS("data/PCATotal_mean_imputation.rds")
loadings_1x  <- readRDS("data/imputed_bootstrap/PCA_loadings_boot_summary.rds")
loadings_2x  <- readRDS("data/imputed_bootstrap_2x/PCA_loadings_boot_summary.rds")

# Trait order (above lowercase, below uppercase)
trait_order <- c("ph", "ssd", "sm", "la", "ln", "sla", "SRL", "D", "RTD", "N")

# Helper to extract loadings matrix
extract_load <- function(mat) {
  m <- as.matrix(mat[trait_order, 1:4, drop = FALSE])
  colnames(m) <- paste0("PC", 1:4)
  m
}

load_original <- extract_load(PCA_original$PCA$loadings)
load_mean     <- extract_load(PCA_mean$PCA$loadings)
load_1x       <- extract_load(loadings_1x$mean)
load_2x       <- extract_load(loadings_2x$mean)

# Build long data frame
make_df <- function(mat, method) {
  as.data.frame(mat, optional = TRUE, stringsAsFactors = FALSE) %>%
    setNames(paste0("PC", 1:4)) %>%
    mutate(Trait = rownames(mat)) %>%
    pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "value") %>%
    mutate(Method = method)
}

long_df <- bind_rows(
  make_df(load_original, "Original"),
  make_df(load_mean,     "Mean"),
  make_df(load_1x,       "1x OOB"),
  make_df(load_2x,       "2x OOB")
)

# Flag the largest absolute value per Trait & Method
long_df <- long_df %>%
  group_by(Method, Trait) %>%
  mutate(is_max = abs(value) == max(abs(value))) %>%
  ungroup()

# Order factors for layout
long_df$Trait  <- factor(long_df$Trait, levels = trait_order)
long_df$Method <- factor(long_df$Method, levels = c("Original", "Mean", "1x OOB", "2x OOB"))
long_df$PC     <- factor(long_df$PC, levels = paste0("PC", 1:4))

# Labels
long_df$label <- sprintf("%.2f", long_df$value)

# Colors: diverging
col_scale <- scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0)

# Plot
p <- ggplot(long_df, aes(PC, Trait, fill = value)) +
  geom_tile(color = "white", size = 0.3) +
  geom_text(aes(label = label, fontface = ifelse(is_max, "bold", "plain")), size = 3) +
  col_scale +
  scale_y_discrete(limits = rev(trait_order)) +
  facet_wrap(~ Method, ncol = 2) +
  labs(title = "PCA Loadings by Imputation Method",
       subtitle = "Bold = largest |loading| per trait within method",
       x = "Principal Component",
       y = "Trait",
       fill = "Loading") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0)
  )

# Save
dir.create("Figures", showWarnings = FALSE)
ggsave("Figures/Loadings_heatmap.png", p, width = 18, height = 16, units = "cm", dpi = 300)

cat("Saved: Figures/Loadings_heatmap.png\n")
