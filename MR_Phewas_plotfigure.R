suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(stringr)
  library(scales)     # pseudo_log_trans() is here
})

## 1) Prepare data_processed (reuse your existing columns; no recomputation)
data_processed <- data_processed_ukb %>%
  mutate(
    broad_group = factor(broad_group, levels = desired_order)
  ) %>%
  arrange(broad_group, Chromosome)

## 2) Right boundary & group centers per broad_group
grp_right <- data_processed %>%
  group_by(broad_group) %>%
  summarise(xmax = max(Chromosome), .groups = "drop")

centers <- data_processed %>%
  group_by(broad_group) %>%
  summarise(center = mean(Chromosome), .groups = "drop")

## 3) Colors
broad_group_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(
  nlevels(data_processed$broad_group)
)
names(broad_group_colors) <- levels(data_processed$broad_group)

## 4) Use only “robustly significant” as label set (align with main data by outcome + broad_group)
label_key <- "outcome"  # If not using outcome, change this to your key column name

label_set <- MR_AGT_signif %>%
  transmute(
    outcome = .data[[label_key]],
    broad_group = factor(broad_group, levels = levels(data_processed$broad_group))
  ) %>%
  distinct()

label_data <- data_processed %>%
  semi_join(label_set, by = c("outcome", "broad_group")) %>%   # Note: no !! needed here
  left_join(grp_right, by = "broad_group") %>%
  mutate(
    label_x = pmax(Chromosome, xmax - 5L),  # Fixed right-side column within group (increase 5L to push farther right)
    seg_x   = Chromosome,                   # Segment start (original x)
    seg_xe  = label_x,                      # Segment end (right-side label column)
    seg_y   = neg_log10_qpvalue,
    seg_ye  = neg_log10_qpvalue
  )

# ===== 0) Choose exposure to plot (if only one, keep as-is) =====
this_exp  <- unique(data_processed_ukb$exposure)[1]
label_key <- "outcome"     # Column name used for labels (change here if needed)
sigma_use <- 1.5           # Pseudo-log compression strength: smaller -> more log-like, stronger compression (0.5–2 common)

# ===== 1) Data prep (keep your “fixed right-side” logic) =====
data_processed <- data_processed_ukb %>%
  filter(exposure == this_exp) %>%
  mutate(
    # If you already have qval/neg_log10_qpvalue, comment out these two lines
    qval = p.adjust(pval, method = "fdr"),
    neg_log10_qpvalue = -log10(qval),
    
    broad_group = factor(broad_group, levels = desired_order),
    shape_dir   = ifelse(b > 0, "up", "down")
  ) %>%
  arrange(broad_group, Chromosome)

# Right boundary & group centers per broad_group (kept exactly as before)
grp_right <- data_processed %>%
  group_by(broad_group) %>%
  summarise(xmax = max(Chromosome), .groups = "drop")

centers <- data_processed %>%
  group_by(broad_group) %>%
  summarise(center = mean(Chromosome), .groups = "drop")

# Colors (assign only to groups present; also unchanged)
broad_group_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(
  nlevels(droplevels(data_processed$broad_group))
)
names(broad_group_colors) <- levels(droplevels(data_processed$broad_group))

# Only label “robustly significant” (from MR_AGT_signif) and cap per-group labels at N to avoid crowding
max_labels_per_group <- 12

label_data <- data_processed %>%
  inner_join(
    MR_AGT_signif %>%
      filter(exposure == this_exp) %>%
      distinct(exposure, outcome, id.outcome),
    by = c("exposure","outcome","id.outcome")
  ) %>%
  group_by(broad_group) %>%
  slice_max(order_by = neg_log10_qpvalue, n = max_labels_per_group, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(grp_right, by = "broad_group") %>%
  mutate(
    # —— Key for “fixed right-side without overflow”: anchor label column near the group's right edge ——
    label_x = pmax(Chromosome, xmax - 5L),   # Fixed right-side column within group
    seg_x   = Chromosome, seg_xe = label_x,  # Segment start/end
    seg_y   = neg_log10_qpvalue, seg_ye = neg_log10_qpvalue,
    label   = str_wrap(.data[[label_key]], width = 32)
  )

# Generate neat ticks for pseudo-log
y_max <- max(data_processed$neg_log10_qpvalue, na.rm = TRUE)
pad_top <- max(0.8, 0.04 * y_max)
y_lim_top <- y_max + pad_top
candidate_breaks <- c(0,1,2,3,5,8,12,20,30,40,60,80,100)
y_breaks <- candidate_breaks[candidate_breaks <= y_lim_top]

manhattan_plot <- ggplot(
  data_processed,
  aes(x = Chromosome, y = neg_log10_qpvalue, fill = broad_group)
) +
  geom_point(aes(color = broad_group, shape = shape_dir), size = 2.5, alpha = 0.95) +
  geom_segment(
    data = label_data,
    aes(x = seg_x, xend = seg_xe, y = seg_y, yend = seg_ye),
    inherit.aes = FALSE, linewidth = 0.25, alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(x = label_x, y = neg_log10_qpvalue, label = .data[[label_key]]),
    inherit.aes = FALSE,
    direction = "y", hjust = 0,
    nudge_x = 0, nudge_y = 0,
    size = 2.8, box.padding = 0.5, point.padding = 0.5,
    min.segment.length = 0, max.overlaps = Inf,
    force = 1.5, max.iter = 3000, max.time = 2,
    segment.alpha = 0.7, seed = 42,
    ylim = c(0, y_lim_top * 0.985)       # Keep labels within plotting area
  ) +
  geom_hline(yintercept = -log10(0.05), color = "gray50",
             linetype = "dashed", linewidth = 0.5) +
  
  # —— Use pseudo-log + top padding ——
  scale_y_continuous(
    trans  = scales::pseudo_log_trans(sigma = sigma_use),  # e.g., 0.8
    breaks = y_breaks,
    limits = c(0, y_lim_top),
    expand = c(0, 0)
  ) +
  
  # Keep x-axis within bounds
  scale_x_continuous(
    breaks = centers$center,
    labels = levels(droplevels(data_processed$broad_group)),
    expand = c(0.02, 0.18)
  ) +
  scale_color_manual(values = broad_group_colors) +
  scale_fill_manual(values = broad_group_colors) +
  scale_shape_manual(values = c("up" = 24, "down" = 25)) +
  
  labs(x = NULL, y = "-log10(FDR)") +
  annotate(
    "text",
    x = min(data_processed$Chromosome, na.rm = TRUE) + 1,
    y = y_lim_top - 0.02 * y_lim_top,
    label = "▲ Beta estimate > 0\n▼ Beta estimate < 0",
    hjust = 0, vjust = 1, size = 3.2
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 8),
    axis.title  = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

manhattan_plot
library(export)

graph2ppt(x = manhattan_plot, 
          file = 'manhattan_plotu_agt.pptx',
          width = 12, height = 9)
