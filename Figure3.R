############################################################
## Combined script: Ins/Del bars + homopolymer bars +
##   HK/Fly error curves + accuracy/number (fly & HK2025)
############################################################

rm(list = ls())

############################################################
## Libraries (all at once)
############################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(forcats)

############################################################
## Working directory — all paths below are relative to this
############################################################
setwd("~/Desktop/benchmark/")

# ==============================================================================
# PART 1 — Mean insertion / deletion per read (Observed_information.txt)
# ==============================================================================

hk2025_sim <- read.table(
  "giraffe/hk2025_simplex/2_Observed_quality/Observed_information.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
hk2025_dup <- read.table(
  "giraffe/hk2025_duplex/2_Observed_quality/Observed_information.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
fly_sim_obs <- read.table(
  "giraffe/D.melanogaster/simplex/reference_sim/2_Observed_quality/Observed_information.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
fly_dup_obs <- read.table(
  "giraffe/D.melanogaster/duplex/reference_dux/2_Observed_quality/Observed_information.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

hk2025_sim$Dataset <- "Dirofilaria Hongkongensis_Simplex"
hk2025_dup$Dataset <- "Dirofilaria Hongkongensis_Duplex"
fly_sim_obs$Dataset <- "Drosophila melanogaster_Simplex"
fly_dup_obs$Dataset <- "Drosophila melanogaster_Duplex"

df_obs <- bind_rows(hk2025_sim, hk2025_dup, fly_sim_obs, fly_dup_obs) %>%
  mutate(
    Ins = as.numeric(Ins),
    Del = as.numeric(Del)
  )

df_obs$Dataset <- factor(
  df_obs$Dataset,
  levels = c(
    "Dirofilaria Hongkongensis_Simplex",
    "Dirofilaria Hongkongensis_Duplex",
    "Drosophila melanogaster_Simplex",
    "Drosophila melanogaster_Duplex"
  )
)

df_ins <- df_obs %>% group_by(Dataset) %>% summarise(mean_ins = mean(Ins), .groups = "drop")
df_del <- df_obs %>% group_by(Dataset) %>% summarise(mean_del = mean(Del), .groups = "drop")

dataset_colors <- c(
  "Drosophila melanogaster_Simplex"   = "#7ADAA5",
  "Drosophila melanogaster_Duplex"    = "#239BA7",
  "Dirofilaria Hongkongensis_Simplex" = "#ECECBB",
  "Dirofilaria Hongkongensis_Duplex"  = "#E1AA36"
)

y_max_ins_del <- max(df_ins$mean_ins, df_del$mean_del)

p_ins <- ggplot(df_ins, aes(x = Dataset, y = mean_ins, fill = Dataset)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = dataset_colors) +
  scale_y_continuous(limits = c(0, y_max_ins_del)) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "Mean insertion count per read", title = "Insertion errors")

p_del <- ggplot(df_del, aes(x = Dataset, y = mean_del, fill = Dataset)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = dataset_colors) +
  scale_y_continuous(limits = c(0, y_max_ins_del)) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "Mean deletion count per read", title = "Deletion errors")

print(p_ins)
print(p_del)

# ==============================================================================
# PART 2 — Homopolymer insertion/deletion counts (stacked bars, hardcoded)
# ==============================================================================

homopolymer_simplex <- data.frame(
  type = c(
    "1 bp Homopolymer Deletion",
    "1 bp Homopolymer Insertion",
    "2 bp Homopolymer Deletion",
    "3 bp Homopolymer Deletion",
    "4 bp Homopolymer Deletion",
    "2 bp Homopolymer Insertion",
    "5-7 bp Homopolymer Insertion"
  ),
  count = c(14, 3, 5, 1, 1, 0, 0),
  data = "D. hongkongensis_simplex"
)

homopolymer_duplex <- data.frame(
  type = c(
    "1 bp Homopolymer Deletion",
    "1 bp Homopolymer Insertion",
    "2 bp Homopolymer Deletion",
    "3 bp Homopolymer Deletion",
    "4 bp Homopolymer Deletion",
    "2 bp Homopolymer Insertion",
    "5-7 bp Homopolymer Insertion"
  ),
  count = c(10, 3, 5, 1, 0, 0, 0),
  data = "D. hongkongensis_duplex"
)

hp_fly_sim_df <- data.frame(
  type = c(
    "1 bp Homopolymer Deletion",
    "1 bp Homopolymer Insertion",
    "2 bp Homopolymer Deletion",
    "3 bp Homopolymer Deletion",
    "4 bp Homopolymer Deletion",
    "2 bp Homopolymer Insertion",
    "5-7 bp Homopolymer Insertion"
  ),
  count = c(28, 13, 1, 0, 3, 2, 1),
  data = "D. melanogaster_simplex"
)

hp_fly_dup_df <- data.frame(
  type = c(
    "1 bp Homopolymer Deletion",
    "1 bp Homopolymer Insertion",
    "2 bp Homopolymer Deletion",
    "3 bp Homopolymer Deletion",
    "4 bp Homopolymer Deletion",
    "2 bp Homopolymer Insertion",
    "5-7 bp Homopolymer Insertion"
  ),
  count = c(21, 17, 2, 1, 0, 1, 0),
  data = "D. melanogaster_duplex"
)

hp_merged <- bind_rows(hp_fly_sim_df, hp_fly_dup_df, homopolymer_simplex, homopolymer_duplex)

hp_merged$data <- factor(
  hp_merged$data,
  levels = c(
    "D. hongkongensis_simplex",
    "D. hongkongensis_duplex",
    "D. melanogaster_simplex",
    "D. melanogaster_duplex"
  )
)

hp_merged$category <- ifelse(grepl("Insertion", hp_merged$type), "Insertion", "Deletion")

hp_merged <- hp_merged %>%
  group_by(category, type) %>%
  mutate(total_type_count = sum(count)) %>%
  ungroup()

hp_merged$type <- fct_reorder(hp_merged$type, hp_merged$total_type_count, .desc = FALSE)

insertion_data <- filter(hp_merged, category == "Insertion")
deletion_data  <- filter(hp_merged, category == "Deletion")

insertion_summary <- insertion_data %>% group_by(data) %>% summarise(total = sum(count), .groups = "drop")
deletion_summary  <- deletion_data  %>% group_by(data) %>% summarise(total = sum(count), .groups = "drop")

global_ymax_hp <- max(
  max(insertion_summary$total),
  max(deletion_summary$total)
) * 1.15

insertion_colors <- c(
  "1 bp Homopolymer Insertion" = "#8DD3C7",
  "2 bp Homopolymer Insertion" = "#BEBADA",
  "5-7 bp Homopolymer Insertion" = "#B3DE69"
)

deletion_colors <- c(
  "1 bp Homopolymer Deletion" = "#8DD3C7",
  "2 bp Homopolymer Deletion" = "#BEBADA",
  "3 bp Homopolymer Deletion" = "#FB8072",
  "4 bp Homopolymer Deletion" = "#80B1D3"
)

bar_insertion <- ggplot(insertion_data, aes(x = data, y = count, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = insertion_summary,
    aes(x = data, y = total, label = paste0("n = ", total)),
    inherit.aes = FALSE,
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = insertion_colors) +
  scale_y_continuous(limits = c(0, global_ymax_hp)) +
  labs(y = "Insertion Count", title = "Homopolymer Insertion Errors") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

bar_deletion <- ggplot(deletion_data, aes(x = data, y = count, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = deletion_summary,
    aes(x = data, y = total, label = paste0("n = ", total)),
    inherit.aes = FALSE,
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = deletion_colors) +
  scale_y_continuous(limits = c(0, global_ymax_hp)) +
  labs(y = "Deletion Count", title = "Homopolymer Deletion Errors") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(bar_insertion)
print(bar_deletion)

# ==============================================================================
# PART 3 — HK2025 + Fly homopolymer error rate (%) line plots
# ==============================================================================

hk2025 <- data.frame(
  homopolymer = c("4T", "5T", "6T", "7T", "8T", "9T", "10T", "11T", "12T", "13T", "14T", "15T", "18T"),
  count = c(265, 159, 95, 63, 35, 24, 12, 7, 4, 6, 3, 3, 1),
  error_simplex = c(265, 159, 95, 62, 33, 24, 11, 3, 0, 0, 1, 0, 0),
  error_duplex  = c(265, 159, 95, 62, 33, 24, 12, 4, 1, 0, 1, 1, 0)
) %>%
  mutate(
    hongkong_simplex = error_simplex / count * 100,
    hongkong_duplex  = error_duplex  / count * 100
  )

hk_levels <- c(
  "4T", "5T", "6T", "7T", "8T", "9T",
  "10T", "11T", "12T", "13T", "14T", "15T", "18T"
)

long_df_hk <- hk2025 %>%
  select(homopolymer, hongkong_simplex, hongkong_duplex) %>%
  pivot_longer(cols = c(hongkong_simplex, hongkong_duplex), names_to = "label", values_to = "error") %>%
  mutate(homopolymer = factor(homopolymer, levels = hk_levels))

hk_colors <- c(
  "hongkong_simplex" = "#ECECBB",
  "hongkong_duplex"  = "#E1AA36"
)

hk2025_figure <- ggplot(long_df_hk, aes(homopolymer, error, color = label, group = label)) +
  geom_point() +
  geom_line(alpha = 0.5, linewidth = 1) +
  scale_color_manual(values = hk_colors) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(labels = number_format(accuracy = 1)) +
  labs(x = "Homopolymer", y = NULL, title = "Dirofilaria hongkongensis") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

fly_a <- data.frame(
  homopolymer = c("4A", "5A", "6A", "7A", "8A", "9A", "10A", "11A", "21A"),
  count = c(185, 87, 54, 17, 18, 9, 3, 1, 1),
  error_simplex = c(185, 87, 54, 17, 14, 7, 2, 0, 0),
  error_duplex  = c(185, 87, 54, 17, 14, 6, 3, 0, 0)
) %>%
  mutate(
    fly_a_simplex = error_simplex / count * 100,
    fly_a_duplex  = error_duplex  / count * 100
  )

fly_t <- data.frame(
  homopolymer = c("4T", "5T", "6T", "7T", "8T", "9T", "10T", "11T", "20T"),
  count = c(197, 83, 48, 10, 10, 16, 3, 2, 1),
  error_simplex = c(197, 83, 48, 9, 8, 12, 3, 0, 0),
  error_duplex  = c(197, 83, 48, 9, 8, 12, 3, 0, 0)
) %>%
  mutate(
    fly_t_simplex = error_simplex / count * 100,
    fly_t_duplex  = error_duplex  / count * 100
  )

fly_a_levels <- c("4A", "5A", "6A", "7A", "8A", "9A", "10A", "11A", "21A")
fly_t_levels <- c("4T", "5T", "6T", "7T", "8T", "9T", "10T", "11T", "20T")

long_df_a <- fly_a %>%
  select(homopolymer, fly_a_simplex, fly_a_duplex) %>%
  pivot_longer(cols = c(fly_a_simplex, fly_a_duplex), names_to = "label", values_to = "error") %>%
  mutate(homopolymer = factor(homopolymer, levels = fly_a_levels))

long_df_t <- fly_t %>%
  select(homopolymer, fly_t_simplex, fly_t_duplex) %>%
  pivot_longer(cols = c(fly_t_simplex, fly_t_duplex), names_to = "label", values_to = "error") %>%
  mutate(homopolymer = factor(homopolymer, levels = fly_t_levels))

fly_a_colors <- c("fly_a_simplex" = "#7ADAA5", "fly_a_duplex" = "#239BA7")
fly_t_colors <- c("fly_t_simplex" = "#7ADAA5", "fly_t_duplex" = "#239BA7")

fly_a_figure <- ggplot(long_df_a, aes(homopolymer, error, color = label, group = label)) +
  geom_point() +
  geom_line(alpha = 0.5, linewidth = 1) +
  scale_color_manual(values = fly_a_colors) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(labels = number_format(accuracy = 1)) +
  labs(x = "Homopolymer", y = NULL, title = "D. melanogaster — A") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

fly_t_figure <- ggplot(long_df_t, aes(homopolymer, error, color = label, group = label)) +
  geom_point() +
  geom_line(alpha = 0.5, linewidth = 1) +
  scale_color_manual(values = fly_t_colors) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(labels = number_format(accuracy = 1)) +
  labs(x = "Homopolymer", y = NULL, title = "D. melanogaster — T") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(hk2025_figure)
print(fly_a_figure)
print(fly_t_figure)

# Optional layout (uncomment to use):
# figure_part3 <- hk2025_figure | (fly_a_figure / fly_t_figure)
# print(figure_part3)

# ==============================================================================
# PART 4 — Accuracy & Number (log2) vs homopolymer length (from Giraffe txt)
# ==============================================================================

acc_fly_simplex <- read.table(
  "giraffe/D.melanogaster/simplex/reference_sim/2_Observed_quality/fly_sim.txt",
  header = TRUE
)
acc_fly_simplex$label <- "Drosophila melanogaster_Simplex"

acc_fly_duplex <- read.table(
  "giraffe/D.melanogaster/duplex/reference_dux/2_Observed_quality/fly_dux.txt",
  header = TRUE
)
acc_fly_duplex$label <- "Drosophila melanogaster_Duplex"

acc_hk_simplex <- read.table(
  "giraffe/hk2025_simplex/2_Observed_quality/hk2025_sim.txt",
  header = TRUE
)
acc_hk_simplex$label <- "Dirofilaria Hongkongensis_Simplex"

acc_hk_duplex <- read.table(
  "giraffe/hk2025_duplex/2_Observed_quality/hk2025_dux.txt",
  header = TRUE
)
acc_hk_duplex$label <- "Dirofilaria Hongkongensis_Duplex"

acc_fly_simplex <- acc_fly_simplex[acc_fly_simplex$base %in% c("A", "T"), ]
acc_fly_duplex  <- acc_fly_duplex[acc_fly_duplex$base %in% c("A", "T"), ]
acc_hk_simplex  <- acc_hk_simplex[acc_hk_simplex$base %in% c("T"), ]
acc_hk_duplex   <- acc_hk_duplex[acc_hk_duplex$base %in% c("T"), ]

merge_acc_fly <- rbind(acc_fly_simplex, acc_fly_duplex)
merge_acc_hk  <- rbind(acc_hk_simplex, acc_hk_duplex)

merge_acc_fly$Number <- log2(merge_acc_fly$Number)
merge_acc_hk$Number  <- log2(merge_acc_hk$Number)

merge_acc_fly$length <- factor(merge_acc_fly$length, levels = sort(unique(merge_acc_fly$length)))
merge_acc_hk$length  <- factor(merge_acc_hk$length, levels = sort(unique(merge_acc_hk$length)))

group_colors <- c(
  "Drosophila melanogaster_Simplex"   = "#7ADAA5",
  "Drosophila melanogaster_Duplex"    = "#239BA7",
  "Dirofilaria Hongkongensis_Simplex" = "#ECECBB",
  "Dirofilaria Hongkongensis_Duplex"  = "#E1AA36"
)

plot_accuracy <- function(df) {
  ggplot(df, aes(length, rate, group = label, color = label)) +
    geom_line(linewidth = 1, alpha = 0.7) +
    geom_point(size = 2, alpha = 0.7) +
    facet_wrap(~base) +
    scale_color_manual(values = group_colors) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(x = "Homopolymer length", y = "Accuracy (%)") +
    theme(
      text = element_text(size = 12),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

plot_number <- function(df) {
  ggplot(df, aes(length, Number, group = label, color = label)) +
    geom_line(linewidth = 1, alpha = 0.7) +
    geom_point(size = 2, alpha = 0.7) +
    facet_wrap(~base) +
    scale_color_manual(values = group_colors) +
    labs(x = "Homopolymer length", y = "Number (log2)") +
    theme(
      legend.position = "none",
      text = element_text(size = 12),
      axis.title = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

plot1 <- plot_accuracy(merge_acc_fly)
plot2 <- plot_number(merge_acc_fly)
plot3 <- plot_accuracy(merge_acc_hk)
plot4 <- plot_number(merge_acc_hk)

print(plot1)
print(plot2)
print(plot3)
print(plot4)

############################################################
## End
############################################################
