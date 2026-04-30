############################################################
## Merged script: (1) feature homopolymer bar plot
##                  (2) duplex T-homopolymer accuracy vs position
##                  (3) HK2025 homopolish bed/sites + stacked bar
############################################################

rm(list = ls())

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(readr)
library(dplyr)

setwd("/Users/yzhou799/Desktop/")

# ==============================================================================
# PART 1 — Homopolymer vs GFF features (PX108868 / nc031365)
# ==============================================================================

# -------------------------
# 1.1 读取 homopolymer 数据
# -------------------------
hp <- read_table(
  "nc031365.txt",
  col_names = c("base", "start", "end")
) %>%
  filter(
    grepl("^[0-9]+$", start),
    grepl("^[0-9]+$", end)
  ) %>%
  mutate(
    start = as.integer(start),
    end   = as.integer(end),
    length = end - start + 1,
    size_class = ifelse(length > 10, ">10", "<=10"),
    hp_id = row_number()
  )

# -------------------------
# 1.2 拆分跨基因 Arg 18T
# -------------------------
arg18_idx <- which(hp$length == 18 & grepl("T", hp$base))
if (length(arg18_idx) == 1) {
  new_homos <- tibble(
    base = c("T", "T", "T"),
    start = c(
      hp$start[arg18_idx],
      hp$start[arg18_idx] + 12,
      hp$start[arg18_idx] + 16
    ),
    end = c(
      hp$start[arg18_idx] + 11,
      hp$start[arg18_idx] + 15,
      hp$start[arg18_idx] + 17
    ),
    length = end - start + 1,
    size_class = ifelse(length > 10, ">10", "<=10"),
    hp_id = max(hp$hp_id) + 1:3
  )
  hp <- hp[-arg18_idx, ] %>%
    bind_rows(new_homos) %>%
    arrange(start) %>%
    mutate(hp_id = row_number())
}

# -------------------------
# 1.3 转 GRanges
# -------------------------
hp_gr <- GRanges(
  seqnames = "PX108868",
  ranges = IRanges(hp$start, hp$end),
  hp_id = hp$hp_id,
  length = hp$length,
  size_class = hp$size_class
)

# -------------------------
# 1.4 读取 GFF 并筛选 feature
# -------------------------
gff <- import("Phylogenetic tree/new-whole genome tree/comparison/PX108868.gff3")
feature_df <- gff %>%
  as.data.frame() %>%
  filter(type %in% c("CDS", "tRNA", "rRNA")) %>%
  mutate(
    feature = case_when(
      type == "CDS" ~ Name,
      TRUE ~ ifelse(is.na(Name), product, Name)
    ),
    feature_type = type
  ) %>%
  mutate(
    feature = case_when(
      type == "tRNA" ~ str_remove(feature, "^tRNA-"),
      type == "rRNA" & grepl("16S", feature) ~ "16S",
      type == "rRNA" & grepl("12S", feature) ~ "12S",
      TRUE ~ feature
    )
  )

feature_gr <- GRanges(
  seqnames = feature_df$seqnames,
  ranges = IRanges(feature_df$start, feature_df$end),
  feature = feature_df$feature,
  feature_type = feature_df$feature_type
)

# -------------------------
# 1.5 找 overlaps；未重叠的 hp 归为 non-coding
# -------------------------
hits <- findOverlaps(hp_gr, feature_gr)
hp_in_feature <- tibble(
  hp_id = queryHits(hits),
  feature = mcols(feature_gr)$feature[subjectHits(hits)],
  feature_type = mcols(feature_gr)$feature_type[subjectHits(hits)]
)
hp_non_coding <- hp %>%
  filter(!hp_id %in% hp_in_feature$hp_id) %>%
  mutate(feature = "non-coding", feature_type = "non-coding")
hp_feature <- hp_in_feature %>%
  left_join(hp, by = "hp_id") %>%
  bind_rows(hp_non_coding)

# -------------------------
# 1.6 统计 feature 的 homopolymer
# -------------------------
feature_summary <- hp_feature %>%
  count(feature, feature_type, size_class)

if (!"non-coding" %in% feature_summary$feature) {
  feature_summary <- feature_summary %>%
    bind_rows(tibble(
      feature = "non-coding",
      feature_type = "non-coding",
      size_class = c("<=10", ">10"),
      n = c(0, 2)
    ))
}

# 指定基因显示顺序（Ala 在上方）：从下到上 = non-coding ... Ala
feature_order <- c(
  "non-coding", "COX3", "Leu1", "CYTB", "Gln", "Arg", "ND6", "Trp", "COX1", "ND4",
  "Thr", "ND2", "Ser1", "Glu", "ND5", "Val", "Asp", "Pro", "Ser2", "Cys",
  "ND3", "16S", "His", "COX2", "Gly", "Ile", "ATP6", "Phe", "ND1", "Tyr",
  "12S", "ND4L", "Lys", "Met", "Asn", "Leu2", "Ala"
)

feature_len_df <- feature_df %>%
  mutate(feature = case_when(
    type == "tRNA" ~ str_remove(feature, "^tRNA-"),
    type == "rRNA" & grepl("16S", feature) ~ "16S",
    type == "rRNA" & grepl("12S", feature) ~ "12S",
    TRUE ~ feature
  )) %>%
  distinct(feature, .keep_all = TRUE) %>%
  mutate(feature_length = end - start + 1) %>%
  select(feature, feature_length) %>%
  filter(feature != "non-coding") %>%
  bind_rows(tibble(feature = "non-coding", feature_length = 406))

hp_gt10_summary <- hp_feature %>%
  filter(size_class == ">10") %>%
  group_by(feature) %>%
  summarise(
    total_hp_length = sum(length),
    hp_count = n(),
    .groups = "drop"
  ) %>%
  filter(feature != "non-coding") %>%
  bind_rows(tibble(feature = "non-coding", total_hp_length = 24, hp_count = 2))

feature_prop <- feature_summary %>%
  group_by(feature) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

hp_gt10_pos <- feature_prop %>%
  filter(size_class == ">10") %>%
  group_by(feature) %>%
  summarise(y_pos = sum(prop), .groups = "drop") %>%
  left_join(hp_gt10_summary, by = "feature")

all_features <- unique(c(feature_summary$feature, feature_len_df$feature, hp_gt10_pos$feature))
feature_order_actual <- c(intersect(feature_order, all_features), setdiff(all_features, feature_order))
feature_summary <- feature_summary %>% mutate(feature = factor(feature, levels = feature_order_actual))
feature_len_df <- feature_len_df %>% mutate(feature = factor(feature, levels = feature_order_actual))
hp_gt10_pos <- hp_gt10_pos %>% mutate(feature = factor(feature, levels = feature_order_actual))

# -------------------------
# 1.7 绘图 Part 1
# -------------------------
p1_feature_homo <- ggplot(feature_summary, aes(x = feature, y = n, fill = size_class)) +
  geom_col(position = "fill") +
  coord_flip(clip = "off") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = c(">10" = "#239BA7", "<=10" = "#7ADAA5")) +
  geom_text(
    data = feature_len_df,
    aes(x = feature, y = 1.05, label = paste0(feature_length, " bp")),
    inherit.aes = FALSE, hjust = 0, size = 3
  ) +
  geom_text(
    data = hp_gt10_pos,
    aes(x = feature, y = y_pos + 0.02, label = paste0(total_hp_length, " bp (n=", hp_count, ")")),
    inherit.aes = FALSE, hjust = 0, size = 3, color = "#1B6F75"
  ) +
  theme_classic() +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5)) +
  labs(x = "Feature", y = "Proportion of homopolymers", fill = "Homopolymer length")

print(p1_feature_homo)

cat("长度>10的homopolymer总数:", sum(hp_feature$size_class == ">10"), "\n")

feature_homo_lengths <- hp_feature %>%
  group_by(feature, feature_type) %>%
  summarise(
    homo_lengths = paste(length, collapse = ", "),
    n_homo = n(),
    .groups = "drop"
  )
print(feature_homo_lengths)

# ==============================================================================
# PART 2 — Duplex: T-homopolymer accuracy vs genomic position
# ==============================================================================

px108868 <- read.table(
  "duplex.homopolymer_in_reference.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

df_T <- px108868 %>%
  filter(str_detect(type, "T$")) %>%
  mutate(
    start = as.numeric(str_split(pos, "_", simplify = TRUE)[, 2]),
    end   = as.numeric(str_split(pos, "_", simplify = TRUE)[, 3]),
    mid_pos = (start + end) / 2,
    hp_len = as.numeric(str_remove(type, "T"))
  ) %>%
  mutate(
    accuracy = num_of_mat / depth,
    hp_group = ifelse(hp_len > 10, ">10T", "≤10T")
  )

p2_duplex_T_acc <- ggplot(df_T, aes(x = mid_pos, y = accuracy)) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) +
  geom_point(aes(color = hp_group), size = 1, alpha = 0.6) +
  scale_color_manual(
    values = c(
      "≤10T" = "#7ADAA5",
      ">10T" = "#239BA7"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Genomic position",
    y = "T-homopolymer accuracy",
    color = "Homopolymer length"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

print(p2_duplex_T_acc)

# ==============================================================================
# PART 3 — HK2025: bed + homopolish sites + genome map + category stack
# ==============================================================================

# 若 bed 为制表符分隔，可将 sep 改为 "\t"
bed <- read.table(
  "Phylogenetic tree/mt/whole mt-genome/HK2025.bed",
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE
)

bed$name <- trimws(as.character(bed$name))

gr_bed <- GRanges(
  seqnames = bed$contig,
  ranges   = IRanges(start = bed$start + 1, end = bed$end),
  gene     = bed$name,
  type     = bed$type
)

sites_raw <- read_csv(
  "Phylogenetic tree/mt/whole mt-genome/homopolish_sites.csv",
  show_col_types = FALSE
)

sites <- sites_raw %>%
  filter(!is.na(`After Homopolish`) & !is.na(Gene)) %>%
  select(
    Gene        = Gene,
    Before_Pos  = `Before Homopolish`,
    Before_Base = `...3`,
    After_Pos   = `After Homopolish`,
    After_Base  = `...5`
  ) %>%
  mutate(
    Before_Pos = as.numeric(Before_Pos),
    After_Pos  = as.numeric(After_Pos),
    pos1 = After_Pos
  )

sites$Gene <- trimws(as.character(sites$Gene))
sites$After_Base <- as.character(sites$After_Base)
sites$After_Base[is.na(sites$After_Base)] <- "."

gr_sites <- GRanges(
  seqnames = "HK2025_polish",
  ranges   = IRanges(start = sites$pos1, end = sites$pos1),
  gene_csv = sites$Gene,
  base     = sites$After_Base
)

ov <- findOverlaps(gr_sites, gr_bed, ignore.strand = TRUE)
mapped_gene <- rep(NA_character_, length(gr_sites))
mapped_gene[queryHits(ov)] <- gr_bed$gene[subjectHits(ov)]

mapped_gene_na_idx <- which(is.na(mapped_gene))
if (length(mapped_gene_na_idx) > 0) {
  mapped_gene[mapped_gene_na_idx] <- sites$Gene[mapped_gene_na_idx] %>% trimws()
  mapped_gene[is.na(mapped_gene)] <- "Intergenic"
}

bed_df <- as_tibble(gr_bed) %>% mutate(gene = trimws(as.character(gene)))
sites_df <- as_tibble(gr_sites) %>%
  mutate(
    gene_csv = trimws(as.character(gene_csv)),
    base = as.character(base),
    mapped_gene = mapped_gene
  )

bed_df <- bed_df %>% arrange(start)
gene_levels <- unique(bed_df$gene) %>% rev()
bed_df$gene <- factor(bed_df$gene, levels = gene_levels)

sites_df$mapped_gene <- factor(
  sites_df$mapped_gene,
  levels = c(gene_levels, setdiff(unique(sites_df$mapped_gene), gene_levels))
)
gene_y_lookup <- setNames(seq_along(gene_levels), gene_levels)
sites_df$y_idx <- gene_y_lookup[as.character(sites_df$mapped_gene)]
sites_df$y_idx[is.na(sites_df$y_idx)] <- 1

stem_height <- 3

p3_homopolish_map <- ggplot() +
  geom_rect(
    data = bed_df,
    aes(
      xmin = start, xmax = end,
      ymin = as.numeric(gene) - 0.3,
      ymax = as.numeric(gene) + 0.3,
      fill = type
    ),
    color = "black", linewidth = 0.3, alpha = 0.7
  ) +
  geom_segment(
    data = sites_df,
    aes(
      x = start, xend = start,
      y = y_idx + 0.3,
      yend = y_idx + 0.3 + stem_height
    ),
    color = "grey70"
  ) +
  geom_point(
    data = sites_df,
    aes(x = start, y = y_idx + 0.3 + stem_height, color = base),
    size = 2, alpha = 0.7
  ) +
  scale_y_continuous(
    breaks = seq_along(gene_levels),
    labels = gene_levels,
    limits = c(0, length(gene_levels) + 0.3)
  ) +
  scale_fill_manual(values = c("CDS" = "#6AECE1", "rRNA" = "#FFF57E", "tRNA" = "#FFB76C")) +
  scale_color_manual(values = c("." = "tomato", "T" = "deepskyblue"), na.value = "black") +
  theme_bw() +
  labs(
    title = "Mitochondrial Genome Homopolish Sites",
    x = "Position",
    y = "Genes",
    fill = "Type",
    color = "Base"
  ) +
  theme(axis.text.y = element_text(size = 8))

print(p3_homopolish_map)

sites2 <- sites %>%
  mutate(
    Category = case_when(
      grepl("^tRNA", Gene) ~ "tRNA",
      Gene %in% c("12S", "16S") ~ "rRNA",
      Gene == "Non-coding region" ~ "Non-coding",
      TRUE ~ "CDS"
    )
  )

sites_unique <- sites2 %>%
  distinct(Gene, Before_Pos, .keep_all = TRUE)

count_df <- sites_unique %>%
  count(Category)

count_df$Category <- factor(
  count_df$Category,
  levels = c("CDS", "tRNA", "rRNA", "Non-coding")
)

stack_df2 <- sites_unique %>%
  mutate(
    Gene_clean = trimws(Gene),
    Category = case_when(
      grepl("^tRNA", Gene_clean) ~ "tRNA",
      Gene_clean %in% c("12S", "16S") ~ "rRNA",
      Gene_clean == "Non-coding region" ~ "Non-coding",
      TRUE ~ "CDS"
    )
  ) %>%
  count(Category, Gene_clean, name = "n")

stack_df2$Category <- factor(
  stack_df2$Category,
  levels = c("CDS", "tRNA", "rRNA", "Non-coding")
)

stack_df2 <- stack_df2 %>%
  group_by(Category) %>%
  arrange(Category, Gene_clean) %>%
  mutate(
    y_mid = cumsum(n) - n / 2,
    Label = paste0(Gene_clean, " (", n, ")")
  ) %>%
  ungroup()

p4_category_stack <- ggplot(stack_df2, aes(x = Category, y = n, fill = Category)) +
  geom_col(width = 0.8, color = "grey60", linewidth = 0.3) +
  geom_text(aes(y = y_mid, label = Label), size = 3) +
  scale_fill_manual(
    values = c(
      "CDS"        = "#6AECE1",
      "rRNA"       = "#FFF57E",
      "tRNA"       = "#FFB76C",
      "Non-coding" = "grey70"
    )
  ) +
  labs(
    x = "Gene category",
    y = "Number of unique Homopolish-corrected sites"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

print(p4_category_stack)

library(tidyverse)
library(viridis)

# 数据准备
df <- read.table(text = "
Gene Homopolymer
tRNA-Ser1 8T
ND2 13T
Non-coding 13T
ND4 13T
ND4 12T
ND4 11T
Non-coding 7T
ND6 14T
Non-coding 4A
Non-coding 18T
Non-coding 11T
ND4L 12T
12S 8T
ND1 15T
16S 15T
16S 13T
16S 11T
ND3 12T
ND5 13T
ND5 14T
ND5 13T
", header = TRUE, stringsAsFactors = FALSE)

df$Length <- as.numeric(gsub("[A-Z]", "", df$Homopolymer))

# 点图
ggplot(df, aes(x = Gene, y = Length, color = Gene)) +
  geom_jitter(
    width = 0.15,
    height = 0,
    size = 3,
    alpha = 0.85
  ) +
  geom_hline(
    yintercept = 10,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey40"
  ) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = "Gene / region",
    y = "Homopolymer length"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "none"
  )

# 可选：保存图形
# ggsave("part1_feature_homo.pdf", p1_feature_homo, width = 10, height = 8)
# ggsave("part2_duplex_T_acc.pdf", p2_duplex_T_acc, width = 10, height = 5)
# ggsave("part3_homopolish_map.pdf", p3_homopolish_map, width = 12, height = 8)
# ggsave("part4_category_stack.pdf", p4_category_stack, width = 8, height = 6)
