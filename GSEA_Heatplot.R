# GSEA 通路富集的 Core enrichment 基因与通路的热图
# 用于展示哪些基因影响了更多的通路
# 哪些通路有许多的共享基因
# 可以选择按照 ES 正负值分组画图
# 基因展示前 50 个
# 通路最多展示 25 条
# 默认差异基因是 DESeq2 的结果，如果不是要修改表头
# 默认 GSEA 结果是 clusterProfiler 产生的

# 需要以下 R 包支持：
# argparse, tidyverse, ComplexHeatmap, BuenColors, circlize

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(BuenColors))
suppressPackageStartupMessages(library(circlize))

scriptDescription="做 GSEA 通路富集和差异基因的热图。使用 clusterProfiler 包进行 GSEA 的结果的 core_enrichment 基因，而不是通路的所有基因。"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--GSEA", dest="GSEA", help="csv 格式 GSEA 结果", required=TRUE)
parser$add_argument("--pvalue-cutoff", dest="PVAL", help="P 值阈值，只分析符合条件的 GSEA 通路，默认 0.25", default=0.25)
parser$add_argument("--DEGs", dest="DEGs", help="csv 格式 DESeq2 差异基因分析结果", required=TRUE)
parser$add_argument("--output-dir", dest="OUT", help="输出目录，默认当前目录", default=".")
parser$add_argument("--basename", dest="BASENAME", help="输出文件名，默认\"Unknown\"", default="Unknown")
parser$add_argument("--plot-title", dest="TITLE", help="热图标题，默认 \"GSEA Pathway Core Enrichment Genes\"", default="GSEA Pathway Core Enrichment Genes")
parser$add_argument("--split", action="store_true", dest="SPLIT", help="是否区分通路方向")

# 取得频次最高的基因列表
getGeneList <- function(pathway_data) {
  geneList <- stringr::str_c(pathway_data$core_enrichment, collapse = "/") %>% 
    strsplit(split = "/", fixed = TRUE) %>% unlist() %>% table() %>% unlist()
  rankGeneList <- geneList[order(geneList, decreasing = TRUE)] %>% names()
  if (length(rankGeneList) > 50) {
    rankGeneList <- rankGeneList[1:50]
  }
  return(rankGeneList)
}
 
# 整理通路基因的差异数据，通路不包含该基因则为 NA
# 按照基因频率排序
getPathwayDEGs <- function(pathway_data, degs_FC) {
  pathwayList <- list()
  pathwayGene <- pathway_data$core_enrichment
  pathwayDescription <- pathway_data$Description
  names(pathwayGene) <- pathwayDescription
  rankGeneList <- getGeneList(pathway_data)
  for (i in 1:length(pathwayGene)) {
    # 按照 rankGeneList 的基因顺序找到相应基因差异倍数，如果相应基因不在 core enrichment 里，那么就赋值NA
    # 然后每一条通路的结果存到 list
    pathwayName <- pathwayDescription[i]
    coreGene <- strsplit(pathwayGene[i], split = "/", fixed = TRUE) %>% unlist()
    geneMatch <- match(rankGeneList, coreGene)
    matchCore <- coreGene[geneMatch]
    geneFC <- degs_FC[match(matchCore, names(degs_FC))]
    names(geneFC) <- rankGeneList
    pathwayList[[pathwayName]] <- geneFC
  }
  return(pathwayList)
}

# 将数据转换到画图矩阵
getPlotData <- function(pathway_data, degs_FC, degs_data) {
  pathwayList <- getPathwayDEGs(pathway_data, degs_FC)
  rankGeneList <- getGeneList(pathway_data)
  geneMap <- tibble(entrezgene_id=as.double(rankGeneList)) %>% left_join(degs_data, by="entrezgene_id") %>% 
    dplyr::select(entrezgene_id, hgnc_symbol)
  plotData <- as.data.frame(pathwayList)
  rownames(plotData) <- geneMap$hgnc_symbol
  colnames(plotData) <- pathway_data$Description
  plotData1 <- plotData[, colSums(!is.na(plotData)) > 0]
  plotData2 <- t(plotData1)
  head(plotData2, n=c(5, 5)) %>% print()
  return(plotData2)
}

# 平均行名长度
meanCharLength <- function(row_names) {
  rowNum <- length(row_names)
  rowLength <- sapply(row_names, nchar)
  totalLength <- sum(rowLength)
  meanLength <- totalLength / rowNum %>% round()
  return(meanLength)
}

# 配置画图参数
# pdf 高度由行数和行名平均长度决定
getPlotConfig <- function(plot_data) {
  plotConfig <- list()
  w <- 27
  rowNames <- rownames(plot_data)
  rowNum <- length(rowNames)
  meanRowLength <- meanCharLength(rowNames)
  heightF1 <- meanRowLength / 40 %>% ceiling()
  heightF2 <- rowNum * 0.03
  h <- rowNum * heightF1 * heightF2 * 0.4 %>% ceiling()
  if (h < 7) {
    h = 7
  }
  minL <- min(plot_data, na.rm = TRUE)
  maxL <- max(plot_data, na.rm = TRUE)
  color <- jdb_palette("ocean_brick", type="continuous")
  if (all(abs(minL) <= 1, abs(maxL) <= 1)) {
    fromL <- -1
    toL <- 1
    breaksL <- c(-1, 0, 1)
    labelsL <- c("-1", "0", "1")
    heightL <- 8
  } else if (all(abs(minL) <= 2, abs(maxL) <= 2)) {
    fromL <- -2
    toL <- 2
    breaksL <- c(-2, 0, 2)
    labelsL <- c("-2", "0", "2")
    heightL <- 8
  } else {
    fromL <- -4
    toL <- 4
    breaksL <- c(-4, -2, 0, 2, 4)
    labelsL <- c("<= -4", "-2", "0", "2", ">= 4")
    heightL <- 12
  }
  plotConfig[["w"]] <- w
  plotConfig[["h"]] <- h
  plotConfig[["color"]] <- color
  plotConfig[["fromL"]] <- fromL
  plotConfig[["toL"]] <- toL
  plotConfig[["breaksL"]] <- breaksL
  plotConfig[["labelsL"]] <- labelsL
  plotConfig[["heightL"]] <- heightL
  return(plotConfig)
}

# 修改行名长度，40 字符换行显示
wrapName <- function(name_str) {
  new_str <- strwrap(name_str, 40) %>% paste(sep = "\n", collapse = "\n")
  return(new_str)
}

# 画图
hmPlot <- function(plot_data, plot_config) {
  colorFun <- colorRamp2(seq(from = plot_config[["fromL"]], to = plot_config[["toL"]], 
                             length.out = length(plot_config[["color"]])), plot_config[["color"]])
  rowLabels <- sapply(rownames(plot_data), wrapName)
  hm <- Heatmap(plot_data, name="FoldChange(Log2)", col=colorFun, column_title = titleText, 
                column_title_side = "top", cluster_rows = FALSE, cluster_columns = FALSE, 
                show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = TRUE, 
                show_column_names = TRUE, row_names_side = "left", column_names_side = "bottom", 
                column_names_rot = 45, row_labels=rowLabels, row_names_max_width=unit(10, "cm"), na_col = "white", column_order = order(colSums(is.na(plot_data))), 
                row_order = order(rowSums(!is.na(plot_data))), border = TRUE, rect_gp = gpar(col="white"), 
                heatmap_legend_param = list(grid_height = unit(plot_config[["heightL"]], "mm"), grid_width = unit(6, "mm"), at = plot_config[["breaksL"]], labels = plot_config[["labelsL"]]))
  return(hm)
}

# 保存图片
drawPlot <- function(pathway_data, plot_path, degs_FC, degs_data) {
  plotData <- getPlotData(pathway_data, degs_FC, degs_data)
  plotConfig <- getPlotConfig(plotData)
  p <- hmPlot(plotData, plotConfig)
  pdf(plot_path, width = plotConfig[["w"]], height = plotConfig[["h"]])
  draw(p)
  dev.off()
}

# 选取前 25 通路
subsetPathways <- function(pathway_data) {
  if (nrow(pathway_data) > 25) {
    pathwayData <- dplyr::slice_head(pathway_data, n=25)
  } else if (nrow(pathway_data) < 5) { 
    writeLines("\no(≧口≦)o \n通路数目小于 5 请检查数据和参数！")
    q(save="no")
  } else
    {
    pathwayData <- pathway_data
  }
  return(pathwayData)
}

argvs <- parser$parse_args()
gseaPath <- file.path(argvs$GSEA)
pCutoff <- as.double(argvs$PVAL)
degsPath <- file.path(argvs$DEGs)
outputDir <- file.path(argvs$OUT)
baseName <- argvs$BASENAME
titleText <- argvs$TITLE
splitPathway <- argvs$SPLIT

# 默认用 entrezgene_id 进行 GSEA
# 用 SYMBOL 进行倒序排序，这样重复 entrezgene_id 时能优先选择有 symbol 的
# 用字符串倒序排序时，空字符串在非空字符串下，NA 在最下
degsData <- read_csv(degsPath) %>% dplyr::arrange(desc(hgnc_symbol)) %>% 
  dplyr::filter(!is.na(entrezgene_id)) %>% dplyr::distinct(entrezgene_id, .keep_all=TRUE)
degsFC <- degsData$log2FoldChange
names(degsFC) <- degsData$entrezgene_id

pathwayData <- read_csv(gseaPath) %>% dplyr::filter(`p.adjust` < pCutoff)
stopifnot(nrow(pathwayData) > 1)
if (splitPathway) {
  # 按照 ES 值分组画图
  pathwayData1 <- dplyr::filter(pathwayData, enrichmentScore > 0) %>% subsetPathways()
  pathwayData2 <- dplyr::filter(pathwayData, enrichmentScore < 0) %>% subsetPathways()
  fileName1 <- paste(baseName, "GSEA_Heatplot", "1", "pdf", sep=".")
  fileName2 <- paste(baseName, "GSEA_Heatplot", "2", "pdf", sep=".")
  plotPath1 <- file.path(outputDir, fileName1)
  plotPath2 <- file.path(outputDir, fileName2)
  drawPlot(pathwayData1, plotPath1, degsFC, degsData)
  drawPlot(pathwayData2, plotPath2, degsFC, degsData)
} else {
  fileName <- paste(baseName, "GSEA_Heatplot", "pdf", sep=".")
  plotPath <- file.path(outputDir, fileName)
  drawPlot(subsetPathways(pathwayData), plotPath, degsFC, degsData)
}

writeLines("ヽ(✿ﾟ▽ﾟ)ノ")