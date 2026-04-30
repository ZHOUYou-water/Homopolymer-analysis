library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

rm(list=ls())

setwd("~/Desktop/")
taxon <- read.csv("ncbi_lineages_2025-09-18.csv")
homo <- read_xlsx("homopolymer_stats.xlsx")
taxon <- taxon[taxon$domain %in% "Eukaryota",]
taxon <- taxon[,c(2:8)]

taxon <- taxon %>%
  filter(if_all(everything(), ~ . != ""))

# ---- Step 1: 先尝试 Species vs species ----
match_species <- homo %>%
  inner_join(taxon, by = c("Species" = "species")) %>%
  mutate(match_type = "species")

# ---- Step 2: 对没有匹配到 species 的再尝试 Species vs genus ----
species_matched_ids <- match_species$Species

match_genus <- homo %>%
  filter(!(Species %in% species_matched_ids)) %>%
  inner_join(taxon, by = c("Species" = "genus")) %>%
  mutate(match_type = "genus")

# ---- Step 3: 合并结果 ----
final_df <- bind_rows(match_species, match_genus) %>%
  distinct()   # 去重（防止 taxon 里同一物种多次出现）

final_df <- final_df %>%
  mutate(phylum2 = ifelse(phylum %in% names(which(table(phylum) < 100)), 
                          "Others", 
                          phylum))

# 辅助函数：从单个字符串（如 "4A:127, 5A:62, 11A:1"）提取出左侧的长度数字（4,5,11）
extract_lengths_from_string <- function(s) {
  s <- as.character(s)
  if (is.na(s) || str_trim(s) == "") return(integer(0))
  parts <- unlist(strsplit(s, ","))
  parts <- str_trim(parts)
  # 只提取每个部分开头的数字（长度），比如 "4A:127" -> "4"
  nums <- str_extract(parts, "^[0-9]+")
  nums <- suppressWarnings(as.numeric(nums))
  nums[!is.na(nums)]
}

# 检查列名是否存在（根据你实际的列名调整）
cols <- c("A-homopolymers", "T-homopolymers", "C-homopolymers", "G-homopolymers")
missing_cols <- setdiff(cols, names(final_df))
if (length(missing_cols) > 0) stop("缺少列: ", paste(missing_cols, collapse = ", "))

# 方法1：用 apply（一行一行处理，返回 0/1）
final_df$homo_over10 <- apply(final_df[cols], 1, function(row) {
  nums <- unlist(lapply(as.list(row), extract_lengths_from_string))
  if (length(nums) == 0) return(0L)
  as.integer(any(nums > 10))
})

final_df$long_bases <- apply(final_df[cols], 1, function(row) {
  bases <- c("A","T","C","G")
  has_long <- mapply(function(cell, base) {
    lengths <- extract_lengths_from_string(cell)
    any(lengths > 10)
  }, as.list(row), bases)
  res <- paste(bases[has_long], collapse = ",")
  if (res == "") NA_character_ else res
})

library(dplyr)
library(ggplot2)

# 先汇总得到 df_bar - 计算比例
df_bar <- final_df %>%
  group_by(phylum2, homo_over10) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(phylum2) %>%
  mutate(
    total = sum(count),
    proportion = count / total,
    percentage_label = paste0(count, "\n(", round(proportion * 100, 1), "%)"),
    total_label = paste0("Total: ", total)
  ) %>%
  ungroup()

# 排序 phylum2
phylum_order <- df_bar %>%
  group_by(phylum2) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(phylum2)

df_bar$phylum2 <- factor(df_bar$phylum2, levels = phylum_order)

# 绘图 - 比例堆叠条形图
bar <- ggplot(df_bar, aes(x = phylum2, y = proportion, fill = factor(homo_over10))) +
  geom_bar(stat = "identity", position = "stack") +
  # 添加数量标签
  geom_text(aes(label = percentage_label, y = proportion),
            position = position_stack(vjust = 0.5),
            size = 3, ) +
  # 在顶部添加总数标签
  geom_text(data = df_bar %>% distinct(phylum2, total_label),
            aes(x = phylum2, y = 1.05, label = total_label),
            inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = c("0" = "#7ADAA5", "1" = "#239BA7"),
                    labels = c("≤10", ">10"),
                    name = "Homopolymer length") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +  # 为顶部标签留空间
  labs(x = "Phylum", y = "Percentage") +
  theme_bw() +
  theme(
    text = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(),
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right","top")
  )

library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

rm(list=ls())

setwd("~/Desktop/")
taxon <- read.csv("ncbi_lineages_2025-09-18.csv")
homo <- read_xlsx("homopolymer_stats.xlsx")
taxon <- taxon[taxon$domain %in% "Eukaryota",]
taxon <- taxon[,c(2:8)]

taxon <- taxon %>%
  filter(if_all(everything(), ~ . != ""))

# ---- Step 1: 先尝试 Species vs species ----
match_species <- homo %>%
  inner_join(taxon, by = c("Species" = "species")) %>%
  mutate(match_type = "species")

# ---- Step 2: 对没有匹配到 species 的再尝试 Species vs genus ----
species_matched_ids <- match_species$Species

match_genus <- homo %>%
  filter(!(Species %in% species_matched_ids)) %>%
  inner_join(taxon, by = c("Species" = "genus")) %>%
  mutate(match_type = "genus")

# ---- Step 3: 合并结果 ----
final_df <- bind_rows(match_species, match_genus) %>%
  distinct()   # 去重（防止 taxon 里同一物种多次出现）

final_df <- final_df %>%
  mutate(phylum2 = ifelse(phylum %in% names(which(table(phylum) < 100)), 
                          "Others", 
                          phylum))

# 辅助函数：从单个字符串（如 "4A:127, 5A:62, 11A:1"）提取出左侧的长度数字（4,5,11）
extract_lengths_from_string <- function(s) {
  s <- as.character(s)
  if (is.na(s) || str_trim(s) == "") return(integer(0))
  parts <- unlist(strsplit(s, ","))
  parts <- str_trim(parts)
  # 只提取每个部分开头的数字（长度），比如 "4A:127" -> "4"
  nums <- str_extract(parts, "^[0-9]+")
  nums <- suppressWarnings(as.numeric(nums))
  nums[!is.na(nums)]
}

# 检查列名是否存在（根据你实际的列名调整）
cols <- c("A-homopolymers", "T-homopolymers", "C-homopolymers", "G-homopolymers")
missing_cols <- setdiff(cols, names(final_df))
if (length(missing_cols) > 0) stop("缺少列: ", paste(missing_cols, collapse = ", "))

# 方法1：用 apply（一行一行处理，返回 0/1）
final_df$homo_over10 <- apply(final_df[cols], 1, function(row) {
  nums <- unlist(lapply(as.list(row), extract_lengths_from_string))
  if (length(nums) == 0) return(0L)
  as.integer(any(nums > 10))
})

final_df$long_bases <- apply(final_df[cols], 1, function(row) {
  bases <- c("A","T","C","G")
  has_long <- mapply(function(cell, base) {
    lengths <- extract_lengths_from_string(cell)
    any(lengths > 10)
  }, as.list(row), bases)
  res <- paste(bases[has_long], collapse = ",")
  if (res == "") NA_character_ else res
})

library(dplyr)
library(ggplot2)

# 先汇总得到 df_bar - 计算比例
df_bar <- final_df %>%
  group_by(phylum2, homo_over10) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(phylum2) %>%
  mutate(
    total = sum(count),
    proportion = count / total,
    percentage_label = paste0(count, "\n(", round(proportion * 100, 1), "%)"),
    total_label = paste0("Total: ", total)
  ) %>%
  ungroup()

# 排序 phylum2
phylum_order <- df_bar %>%
  group_by(phylum2) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(phylum2)

df_bar$phylum2 <- factor(df_bar$phylum2, levels = phylum_order)

# 绘图 - 比例堆叠条形图
bar <- ggplot(df_bar, aes(x = phylum2, y = proportion, fill = factor(homo_over10))) +
  geom_bar(stat = "identity", position = "stack") +
  # 添加数量标签
  geom_text(aes(label = percentage_label, y = proportion),
            position = position_stack(vjust = 0.5),
            size = 3, family = "Times New Roman") +
  # 在顶部添加总数标签
  geom_text(data = df_bar %>% distinct(phylum2, total_label),
            aes(x = phylum2, y = 1.05, label = total_label),
            inherit.aes = FALSE, size = 3, family = "Times New Roman") +
  scale_fill_manual(values = c("0" = "#7ADAA5", "1" = "#239BA7"),
                    labels = c("≤10", ">10"),
                    name = "Homopolymer length") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +  # 为顶部标签留空间
  labs(x = "Phylum", y = "Percentage") +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times New Roman"),
    axis.text.y = element_text(family = "Times New Roman"),
    axis.title.x = element_text(family = "Times New Roman"),
    axis.title.y = element_text(family = "Times New Roman"),
    legend.text = element_text(family = "Times New Roman",size = 12),
    legend.title = element_text(family = "Times New Roman", size = 13),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right","top")
  )

bar

# 分析 Nematoda 门的数据
nematode <- final_df[final_df$phylum2 %in% "Nematoda",]
nematode <- nematode[,-c(14,15,16)]

# 汇总 Nematoda 的数据 - 计算比例
df_bar2 <- nematode %>%
  group_by(order, homo_over10) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(order) %>%
  mutate(
    total = sum(count),
    proportion = count / total,
    percentage_label = paste0(count, "\n(", round(proportion * 100, 1), "%)"),
    total_label = paste0("Total: ", total)
  ) %>%
  ungroup()

# 排序 order
phylum_order2 <- df_bar2 %>%
  group_by(order) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(order)

df_bar2$order <- factor(df_bar2$order, levels = phylum_order2)

# 绘图 - 比例堆叠条形图
bar2 <- ggplot(df_bar2, aes(x = order, y = proportion, fill = factor(homo_over10))) +
  geom_bar(stat = "identity", position = "stack") +
  # 添加数量标签
  geom_text(aes(label = percentage_label, y = proportion),
            position = position_stack(vjust = 0.5),
            size = 3, family = "Times New Roman") +
  # 在顶部添加总数标签
  geom_text(data = df_bar2 %>% distinct(order, total_label),
            aes(x = order, y = 1.05, label = total_label),
            inherit.aes = FALSE, size = 3, family = "Times New Roman") +
  scale_fill_manual(values = c("0" = "#7ADAA5", "1" = "#239BA7"),
                    labels = c("≤10", ">10"),
                    name = "Homopolymer length") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +  # 为顶部标签留空间
  labs(x = "Order", y = "Percentage") +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times New Roman"),
    axis.text.y = element_text(family = "Times New Roman"),
    axis.title.x = element_text(family = "Times New Roman"),
    axis.title.y = element_text(family = "Times New Roman"),
    legend.text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_text(family = "Times New Roman", size = 13),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right","top")
  )

# ggsave("figure5_bar.pdf", dpi = 300, plot = bar, width = 8, height = 6)

# 分析 Nematoda 门的数据
nematode <- final_df[final_df$phylum2 %in% "Nematoda",]
nematode <- nematode[,-c(14,15,16)]

# 汇总 Nematoda 的数据 - 计算比例
df_bar2 <- nematode %>%
  group_by(order, homo_over10) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(order) %>%
  mutate(
    total = sum(count),
    proportion = count / total,
    percentage_label = paste0(count, "\n(", round(proportion * 100, 1), "%)"),
    total_label = paste0("Total: ", total)
  ) %>%
  ungroup()

# 排序 order
phylum_order2 <- df_bar2 %>%
  group_by(order) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(order)

df_bar2$order <- factor(df_bar2$order, levels = phylum_order2)

# 绘图 - 比例堆叠条形图
bar2 <- ggplot(df_bar2, aes(x = order, y = proportion, fill = factor(homo_over10))) +
  geom_bar(stat = "identity", position = "stack") +
  # 添加数量标签
  geom_text(aes(label = percentage_label, y = proportion),
            position = position_stack(vjust = 0.5),
            size = 3) +
  # 在顶部添加总数标签
  geom_text(data = df_bar2 %>% distinct(order, total_label),
            aes(x = order, y = 1.05, label = total_label),
            inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = c("0" = "#7ADAA5", "1" = "#239BA7"),
                    labels = c("≤10", ">10"),
                    name = "Homopolymer length") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +  # 为顶部标签留空间
  labs(x = "Order", y = "Percentage") +
  theme_bw() +
  theme(
    text = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(),
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right","top")
  )

# ggsave("figure5_bar2.pdf", dpi = 300, plot = bar2, width = 8, height = 6)
