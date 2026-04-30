# =============================================================================
# 线粒体基因组 Circos：多物种合并脚本 (danio, human, K.azureus, cgc1, fly, hk2025)
# =============================================================================
rm(list = ls())

library(circlize)
library(dplyr)

setwd("~/Desktop/benchmark/")

# -----------------------------------------------------------------------------
# 公共函数：同源聚合物三层轨道
# -----------------------------------------------------------------------------
add_homopolymer_tracks <- function(homopolymer) {
  hp47 <- homopolymer[homopolymer$length >= 4 & homopolymer$length <= 7, ]
  hp810 <- homopolymer[homopolymer$length >= 8 & homopolymer$length <= 10, ]
  hp10 <- homopolymer[homopolymer$length > 10, ]

  circos.genomicTrack(
    hp47, ylim = c(0, 1), track.height = 0.13,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = "#63C8FF", border = NA)
    }
  )
  circos.genomicTrack(
    hp810, ylim = c(0, 1), track.height = 0.13,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = "#E14434", border = NA)
    }
  )
  circos.genomicTrack(
    hp10, ylim = c(0, 1), track.height = 0.13,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = "#FF2DD1", border = NA)
    }
  )
}

# 基因条颜色：cgc1 需先判 intergenic-region
resolve_gene_color <- function(g, gene_colors, cgc1_mode = FALSE) {
  if (grepl("trn", g, ignore.case = TRUE)) {
    return(unname(gene_colors["tRNA"]))
  }
  if (grepl("rRNA", g, ignore.case = TRUE)) {
    return(unname(gene_colors["rRNA"]))
  }
  if (grepl("CDS", g, ignore.case = TRUE)) {
    return(unname(gene_colors["CDS"]))
  }
  if (cgc1_mode && grepl("intergenic-region", g, ignore.case = TRUE)) {
    return(unname(gene_colors["intergenic-region"]))
  }
  rest <- setdiff(
    names(gene_colors),
    c("tRNA", "rRNA", "CDS", "intergenic-region")
  )
  if (length(rest) >= 1) {
    return(unname(gene_colors[rest[1]]))
  }
  unname(gene_colors["CDS"])
}

add_gene_track <- function(genome_df, gene_colors, cgc1_mode = FALSE) {
  circos.genomicTrack(
    genome_df, ylim = c(0, 1), track.height = 0.15,
    panel.fun = function(region, value, ...) {
      for (i in seq_len(nrow(region))) {
        g <- value$gene[i]
        col <- resolve_gene_color(g, gene_colors, cgc1_mode = cgc1_mode)
        circos.genomicRect(
          region[i, , drop = FALSE], value[i, , drop = FALSE],
          ybottom = 0, ytop = 1, col = col, border = NA
        )
      }
    }
  )
}

# 读 homopolymer 并统一 length
read_homopolymer <- function(path, chr_name) {
  hp <- read.table(file = path, header = TRUE)
  hp$chr <- chr_name
  hp <- hp[, c("chr", "start", "end", "base")]
  hp$length <- as.numeric(gsub("[^0-9]", "", hp$base))
  hp
}

# 一次完整出图
run_one_circos <- function(
  genome_df,
  homopolymer_df,
  gene_colors,
  pdf_file,
  legend_cex = 0.6,
  legend_homopolymer = FALSE,
  cgc1_mode = FALSE
) {
  pdf(pdf_file, width = 8, height = 6)
  on.exit(
    {
      circos.clear()
      par(family = "")
      dev.off()
    },
    add = TRUE
  )

  par(family = "serif")
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4)
  circos.genomicInitialize(genome_df)

  add_gene_track(genome_df, gene_colors, cgc1_mode = cgc1_mode)
  add_homopolymer_tracks(homopolymer_df)

  if (legend_homopolymer) {
    legend(
      "bottomright",
      legend = c(
        names(gene_colors),
        "Homopolymer 4-7",
        "Homopolymer 8-10",
        "Homopolymer >10"
      ),
      fill = c(
        unname(gene_colors),
        "#63C8FF",
        "#E14434",
        "#FF2DD1"
      ),
      title = "Annotation",
      cex = legend_cex,
      box.lty = 0,
      bg = "white"
    )
  } else {
    legend(
      "bottomright",
      legend = names(gene_colors),
      fill = unname(gene_colors),
      title = "Gene Types",
      cex = legend_cex,
      box.lty = 0,
      bg = "white"
    )
  }
}

# =============================================================================
# 1. Danio (figure2_danio.pdf)
# =============================================================================
danio <- read.csv("circos/danio.bed", sep = "", header = FALSE)
danio <- danio[, -c(2, 3)]
colnames(danio)[colnames(danio) == "V1"] <- "chr"
colnames(danio)[colnames(danio) == "V4"] <- "start"
colnames(danio)[colnames(danio) == "V5"] <- "end"
colnames(danio)[colnames(danio) == "V7"] <- "gene"
colnames(danio)[colnames(danio) == "V6"] <- "strand"
danio <- danio[, c("chr", "start", "end", "gene", "strand")]

target_values <- c("tRNA", "rRNA", "D-loop")
danio$gene <- ifelse(!(danio$gene %in% target_values), paste0("CDS", danio$gene), danio$gene)

hp_danio <- read_homopolymer("circos/danio_homopolymers_output.txt", "NC_002333.2")

gene_colors_danio <- c(
  "tRNA" = "#E4004B",
  "rRNA" = "#ED775A",
  "CDS" = "#FAD691",
  "D-loop" = "#FCC61D"
)

run_one_circos(danio, hp_danio, gene_colors_danio, "figure2_danio.pdf", legend_cex = 0.8)

# =============================================================================
# 2. Human (figure2_human.pdf)
# =============================================================================
human <- read.csv("circos/human.bed", sep = "", header = FALSE)
human <- human[, -c(2, 3)]
colnames(human)[colnames(human) == "V1"] <- "chr"
colnames(human)[colnames(human) == "V4"] <- "start"
colnames(human)[colnames(human) == "V5"] <- "end"
colnames(human)[colnames(human) == "V6"] <- "gene"
colnames(human)[colnames(human) == "V7"] <- "strand"

target_values_h <- c("trn", "rRNA", "D-loop")
human$gene <- ifelse(!(human$gene %in% target_values_h), paste0("CDS", human$gene), human$gene)

hp_human <- read_homopolymer("circos/human_homopolymers_output.txt", "NC_012920.1")

gene_colors_human <- c(
  "tRNA" = "#E4004B",
  "rRNA" = "#ED775A",
  "CDS" = "#FAD691",
  "D-loop" = "#FCC61D"
)

run_one_circos(human, hp_human, gene_colors_human, "figure2_human.pdf", legend_cex = 0.6)

# =============================================================================
# 3. K.azureus (figure2_k.pdf) — 图例含 homopolymer
# =============================================================================
k <- read.csv("circos/K.bed", sep = "", header = FALSE)
k <- k[, -5]
colnames(k)[colnames(k) == "V1"] <- "chr"
colnames(k)[colnames(k) == "V2"] <- "start"
colnames(k)[colnames(k) == "V3"] <- "end"
colnames(k)[colnames(k) == "V4"] <- "gene"
colnames(k)[colnames(k) == "V6"] <- "strand"

target_values_k <- c("trn", "rRNA", "D-loop")
k$gene <- ifelse(!(k$gene %in% target_values_k), paste0("CDS", k$gene), k$gene)

hp_k <- read_homopolymer("circos/PP437478_homopolymers_output.txt", "PP437478_1")

gene_colors_k <- c(
  "tRNA" = "#E4004B",
  "rRNA" = "#ED775A",
  "CDS" = "#FAD691",
  "D-loop" = "#FCC61D"
)

run_one_circos(
  k, hp_k, gene_colors_k, "figure2_k.pdf",
  legend_cex = 0.6,
  legend_homopolymer = TRUE
)

# =============================================================================
# 4. CGC1 (figure2_cgc1.pdf)
# =============================================================================
cgc1 <- read.csv("circos/cgc1.bed", sep = "", header = FALSE)
cgc1 <- cgc1[, -c(5)]
colnames(cgc1)[colnames(cgc1) == "V1"] <- "chr"
colnames(cgc1)[colnames(cgc1) == "V2"] <- "start"
colnames(cgc1)[colnames(cgc1) == "V3"] <- "end"
colnames(cgc1)[colnames(cgc1) == "V4"] <- "gene"
colnames(cgc1)[colnames(cgc1) == "V6"] <- "strand"

target_values_c <- c("tRNA", "rRNA", "intergenic-region", "AT-region")
cgc1$gene <- ifelse(!(cgc1$gene %in% target_values_c), paste0("CDS", cgc1$gene), cgc1$gene)

hp_cgc <- read_homopolymer("circos/cgc_homopolymers_output.txt", "AP038780_1")

gene_colors_cgc <- c(
  "tRNA" = "#E4004B",
  "rRNA" = "#ED775A",
  "CDS" = "#FAD691",
  "AT-region" = "#CADCAE",
  "intergenic-region" = "#AE75DA"
)

run_one_circos(
  cgc1, hp_cgc, gene_colors_cgc, "figure2_cgc1.pdf",
  legend_cex = 0.8,
  cgc1_mode = TRUE
)

# =============================================================================
# 5. Fly (figure4_fly.pdf)
# =============================================================================
fly <- read.csv("circos/fly.bed", sep = "", header = FALSE)
fly <- fly[, -c(2, 3)]
colnames(fly)[colnames(fly) == "V1"] <- "chr"
colnames(fly)[colnames(fly) == "V4"] <- "start"
colnames(fly)[colnames(fly) == "V5"] <- "end"
colnames(fly)[colnames(fly) == "V7"] <- "gene"
colnames(fly)[colnames(fly) == "V6"] <- "strand"

target_values_f <- c("tRNA", "rRNA", "A+T-region")
fly$gene <- ifelse(!(fly$gene %in% target_values_f), paste0("CDS", fly$gene), fly$gene)

hp_fly <- read_homopolymer("circos/fly_homopolymers_output.txt", "NC_024511.2")

gene_colors_fly <- c(
  "tRNA" = "#E4004B",
  "rRNA" = "#ED775A",
  "CDS" = "#FAD691",
  "A+T-region" = "#FCC61D"
)

run_one_circos(fly, hp_fly, gene_colors_fly, "figure4_fly.pdf", legend_cex = 0.6)

# =============================================================================
# 6. HK2025 (figure4_hk2025.pdf)
# =============================================================================
hk2025 <- read.csv("circos/HK2025.bed", sep = "", header = FALSE)
colnames(hk2025)[colnames(hk2025) == "V1"] <- "chr"
colnames(hk2025)[colnames(hk2025) == "V2"] <- "start"
colnames(hk2025)[colnames(hk2025) == "V3"] <- "end"
colnames(hk2025)[colnames(hk2025) == "V4"] <- "gene"
hk2025$strand <- "+"

target_values_hk <- c("tRNA", "rRNA")
hk2025$gene <- ifelse(!(hk2025$gene %in% target_values_hk), paste0("CDS", hk2025$gene), hk2025$gene)

hp_hk <- read_homopolymer("circos/hk2025_homopolymers_output.txt", "HK2025")

gene_colors_hk <- c(
  "tRNA" = "#E4004B",
  "rRNA" = "#ED775A",
  "CDS" = "#FAD691"
)

run_one_circos(hk2025, hp_hk, gene_colors_hk, "figure4_hk2025.pdf", legend_cex = 0.6)

message("完成：已生成 figure2_danio.pdf, figure2_human.pdf, figure2_k.pdf, figure2_cgc1.pdf, figure4_fly.pdf, figure4_hk2025.pdf")
