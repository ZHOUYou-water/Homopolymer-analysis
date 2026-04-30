rm(list = ls())
setwd("~/Desktop/benchmark/")

library(tidyverse)
library(ggdist)

################################################################################
## 一、Read accuracy：Observed vs Estimated（箱线图）
################################################################################

read_est_obs <- function(obs_path, est_path, group_name) {
  tmp <- read.table(obs_path, header = TRUE)
  tmp <- tmp[, c("ID", "Acc")]
  colnames(tmp) <- c("ReadID", "Accuracy")
  tmp$Group <- group_name
  tmp$group <- "Observed"

  est <- read.table(est_path, header = TRUE)
  est <- est[, c("ReadID", "Accuracy")]
  est$ReadID <- gsub("@", "", est$ReadID)
  est$Group <- group_name
  est$group <- "Estimated"

  rbind(tmp, est)
}

human_qc <- read_est_obs(
  "giraffe/human/Giraffe_Results/2_Observed_quality/Observed_information.txt",
  "giraffe/human/Giraffe_Results/1_Estimated_quality/Estimated_information.txt",
  "H.sapiens"
)

danio_qc <- read_est_obs(
  "giraffe/Danio/Giraffe_Results/2_Observed_quality/Observed_information.txt",
  "giraffe/Danio/Giraffe_Results/1_Estimated_quality/Estimated_information.txt",
  "D.rerio"
)

K_qc <- read_est_obs(
  "giraffe/K.azureus/Giraffe_Results/2_Observed_quality/Observed_information.txt",
  "giraffe/K.azureus/Giraffe_Results/1_Estimated_quality/Estimated_information.txt",
  "K.azureus"
)

cgc1_qc <- read_est_obs(
  "giraffe/CGC1/2_Observed_quality/Observed_information.txt",
  "giraffe/CGC1/1_Estimated_quality/Estimated_information.txt",
  "C.elegans"
)

fly_sim_qc <- read_est_obs(
  "giraffe/D.melanogaster/simplex/reference_sim/2_Observed_quality/Observed_information.txt",
  "giraffe/D.melanogaster/simplex/reference_sim/1_Estimated_quality/Estimated_information.txt",
  "D.melanogaster_simplex"
)

fly_duplex_qc <- read_est_obs(
  "giraffe/D.melanogaster/duplex/reference_dux/2_Observed_quality/Observed_information.txt",
  "giraffe/D.melanogaster/duplex/reference_dux/1_Estimated_quality/Estimated_information.txt",
  "D.melanogaster_duplex"
)

hk2025_sim_qc <- read_est_obs(
  "giraffe/hk2025_simplex/2_Observed_quality/Observed_information.txt",
  "giraffe/hk2025_simplex/1_Estimated_quality/Estimated_information.txt",
  "D.hongkongensis_simplex"
)

hk2025_duplex_qc <- read_est_obs(
  "giraffe/hk2025_duplex/2_Observed_quality/Observed_information.txt",
  "giraffe/hk2025_duplex/1_Estimated_quality/Estimated_information.txt",
  "D.hongkongensis_duplex"
)

df_qc <- bind_rows(
  danio_qc,
  human_qc,
  K_qc,
  cgc1_qc,
  fly_sim_qc,
  fly_duplex_qc,
  hk2025_sim_qc,
  hk2025_duplex_qc
)

df_qc$Accuracy <- df_qc$Accuracy * 100
df_qc$group_type <- factor(df_qc$group, levels = c("Estimated", "Observed"))

df_qc$Group <- factor(
  df_qc$Group,
  levels = c(
    "D.rerio",
    "H.sapiens",
    "K.azureus",
    "C.elegans",
    "D.melanogaster_simplex",
    "D.melanogaster_duplex",
    "D.hongkongensis_simplex",
    "D.hongkongensis_duplex"
  )
)

group_colors <- c(
  "H.sapiens" = "#7ADAA5",
  "C.elegans" = "#239BA7",
  "K.azureus" = "#ECECBB",
  "D.rerio" = "#E1AA36",
  "D.melanogaster_simplex" = "#B07AA1",
  "D.melanogaster_duplex" = "#8F5AA8",
  "D.hongkongensis_simplex" = "#A6CEE3",
  "D.hongkongensis_duplex" = "#1F78B4"
)

box1 <- ggplot(
  df_qc,
  aes(x = Group, y = Accuracy, color = Group, fill = group_type)
) +
  geom_boxplot(
    position = position_dodge2(width = 0.7, preserve = "single"),
    width = 0.6,
    outlier.shape = NA
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(
    values = c("Estimated" = "grey90", "Observed" = "white")
  ) +
  coord_cartesian(ylim = c(70, 100)) +
  labs(x = "Species", y = "Read accuracy (%)") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(box1)
# ggsave("figure1_box1_all.pdf", box1, dpi = 300, width = 10, height = 6)

################################################################################
## 二、Homopolymer accuracy（散点图）
################################################################################

human_hom <- read.table(
  "giraffe/human/Giraffe_Results/2_Observed_quality/Homoploymer_summary.txt",
  header = TRUE
)

danio_hom <- read.table(
  "giraffe/Danio/Giraffe_Results/2_Observed_quality/Homoploymer_summary.txt",
  header = TRUE
)

K_hom <- read.table(
  "giraffe/K.azureus/Giraffe_Results/2_Observed_quality/Homoploymer_summary.txt",
  header = TRUE
)

cgc1_hom <- read.table(
  "giraffe/CGC1/2_Observed_quality/Homoploymer_summary.txt",
  header = TRUE
)

fly_sim_hom <- read.table(
  "giraffe/D.melanogaster/simplex/reference_sim/2_Observed_quality/Homoploymer_summary.txt",
  header = TRUE
)

fly_duplex_hom <- read.table(
  "giraffe/D.melanogaster/duplex/reference_dux/2_Observed_quality/Homoploymer_summary.txt",
  header = TRUE
)

hk_sim_hom <- read.table(
  "giraffe/hk2025_simplex/2_Observed_quality/Homoploymer_summary.txt",
  header = TRUE
)

hk_duplex_hom <- read.table(
  "giraffe/hk2025_duplex/2_Observed_quality/Homoploymer_summary.txt",
  header = TRUE
)

df_hom <- bind_rows(
  human_hom,
  danio_hom,
  K_hom,
  cgc1_hom,
  fly_sim_hom,
  fly_duplex_hom,
  hk_sim_hom,
  hk_duplex_hom
) %>%
  mutate(
    Group = recode(
      Group,
      "human_all" = "H.sapiens",
      "Danio_all" = "D.rerio",
      "k.azureus_all" = "K.azureus",
      "CGC1" = "C.elegans",
      "fly_simplex" = "D.melanogaster_simplex",
      "fly_duplex" = "D.melanogaster_duplex",
      "simplex" = "D.hongkongensis_simplex",
      "duplex" = "D.hongkongensis_duplex"
    ),
    Accuracy = Accuracy * 100,
    Group = factor(
      Group,
      levels = c(
        "D.rerio",
        "H.sapiens",
        "K.azureus",
        "C.elegans",
        "D.melanogaster_simplex",
        "D.melanogaster_duplex",
        "D.hongkongensis_simplex",
        "D.hongkongensis_duplex"
      )
    )
  )

p_hom <- ggplot(df_hom, aes(x = Base, y = Accuracy, color = Group)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  labs(
    x = "Base",
    y = "Overall accuracy of homopolymer (%)"
  ) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(p_hom)
# ggsave("homopolymer_accuracy_all_species.pdf", p_hom, width = 8, height = 5, dpi = 300)
