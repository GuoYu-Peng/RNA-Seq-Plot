# 通路富集 hit 到的基因与通路热图
# 用于展示哪些基因影响了更多的通路
# 哪些通路有许多的共享基因
# 基因展示前 50 个
# 通路最多展示 25 条
# 默认差异基因是 DESeq2 的结果，如果不是要修改表头
# 默认 clusterProfiler 作为通路富集工具，并且使用 entrezgene_id 富集
# 生成的热图会使用 hgnc_symbol 方便阅读

# 需要以下 R 包支持：
# argparse, tidyverse, ComplexHeatmap, BuenColors, circlize

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(BuenColors))
suppressPackageStartupMessages(library(circlize))

scriptDescription="通路富集和差异基因的热图。使用 clusterProfiler 包进行通路富集的结果 hit 到基因，而不是通路的所有基因。"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--Enrichment", dest="ENRICHMENT", help="csv 格式通路富集结果，默认使用 entrezgene_id 进行富集分析", required=TRUE)
parser$add_argument("--pvalue-cutoff", dest="PVAL", help="P 值阈值，只分析符合条件的通路，默认 0.25", default=0.25)
parser$add_argument("--DEGs", dest="DEGs", help="csv 格式 DESeq2 差异基因分析结果", required=TRUE)
parser$add_argument("--output-dir", dest="OUTDIR", help="输出目录，默认当前目录", default=".")
parser$add_argument("--basename", dest="BASENAME", help="输出文件名，默认\"Unknown\"", default="Unknown")
parser$add_argument("--plot-title", dest="TITLE", help="热图标题，默认 \"Genes of Enrichment Pathways\"", default="Genes of Enrichment Pathways")

argvs <- parser$parse_args()
inputPath <- file.path(argvs$ENRICHMENT)
outputDir <- file.path(argvs$OUTDIR)
degsPath <- file.path(argvs$DEGs)
pvalueCutoff <- as.double(argvs$PVAL)
baseName <- as.character(argvs$BASENAM)
plotTitle <- as.character(argvs$TITLE)

# 将富集结果 hit 基因拆分成向量
enrichGenes <- function(hit_gene) {
  genes <- strsplit(hit_gene, split = "/", fixed = TRUE) %>% unlist()
  return(genes)
}

# 修改行名长度，40 字符换行显示
wrapName <- function(name_str) {
  new_str <- strwrap(name_str, 40) %>% paste(sep = "\n", collapse = "\n")
  return(new_str)
}

# 平均行名长度
meanCharLength <- function(row_names) {
  rowNum <- length(row_names)
  rowLength <- sapply(row_names, nchar)
  totalLength <- sum(rowLength)
  meanLength <- totalLength / rowNum %>% round()
  return(meanLength)
}

# 只保留 25 条通路记录，slice 函数越界不会出错
pathwayData <- read_csv(inputPath) %>% dplyr::filter(`p.adjust` < pvalueCutoff) %>% 
  arrange(`p.adjust`) %>% slice_head(n = 25)
if (nrow(pathwayData) < 5) {
  writeLines("\n通路数目过少，请检查数据和参数，脚本将退出")
  q(save = "no")
}

degsData <- read_csv(degsPath) %>% dplyr::filter(!is.na(entrezgene_id)) %>% dplyr::select(entrezgene_id, hgnc_symbol, log2FoldChange) %>% 
  arrange(desc(hgnc_symbol)) %>% distinct(entrezgene_id, .keep_all=TRUE)

# 取得每个通路 Hit 基因并保存到列表
pathwayNames <- pathwayData$Description
pathwayHits <- pathwayData$geneID
hitGenes <- list()
for (i in 1:length(pathwayNames)) {
  pathwayName <- pathwayNames[i]
  pathwayHit <- pathwayHits[i]
  pathwayGenes <- enrichGenes(pathwayHit)
  hitGenes[[pathwayName]] <- pathwayGenes
}
str(hitGenes)

allHits <- unlist(hitGenes)
hitFreq <- table(allHits)
hitFreqRank <- hitFreq[order(hitFreq, decreasing = TRUE)]
# 只保留频率最高的 50 基因
if (length(hitFreqRank) > 50) {
  hitFreqRank <- hitFreqRank[1:50]
}
# 按照频率进行排序
degsData2 <- dplyr::filter(degsData, entrezgene_id %in% names(hitFreqRank)) %>% 
  mutate(entrezgene_id = factor(entrezgene_id, levels = names(hitFreqRank))) %>% arrange(entrezgene_id)
geneNum <- nrow(degsData2)
pathwayNum <- length(pathwayNames)

# 先构建 NA 矩阵再填充基因数据
hitMatrix <- matrix(data = NA, nrow = pathwayNum, ncol = geneNum)
degsFc <- degsData2$log2FoldChange
degsSymbol <- degsData2$hgnc_symbol
degsEntrez <- degsData2$entrezgene_id
for (i in 1:pathwayNum) {
  ithHitGenes <- hitGenes[[i]]
  for (j in 1:geneNum) {
    jthLogFC <- degsFc[j]
    geneID <- degsEntrez[j]
    if (geneID %in% ithHitGenes) {
      hitMatrix[i, j] <- jthLogFC
    }
  }
}
rownames(hitMatrix) <- pathwayNames
colnames(hitMatrix) <- degsSymbol

# 调整图例
maxAbsFC <- max(abs(hitMatrix), na.rm = TRUE)
minFC <- NA
maxFC <- NA
if (maxAbsFC <= 1) {
  minFC <- -1
  maxFC <- 1
} else if (maxAbsFC <=2) {
  minFC <- -4
  maxFC <- 4
} else {
  minFC <- -2
  maxFC <- 2
}

rowLabels <- sapply(pathwayNames, wrapName)
colSeq <- jdb_palette("ocean_brick", type="continuous")
colFun <- colorRamp2(breaks = seq(from = minFC, to = maxFC, length.out=length(colSeq)), colors = colSeq)
heatMap <- Heatmap(matrix = hitMatrix, col = colFun, name = "FoldChange(Log2)", column_title = plotTitle, 
                   column_title_side = "top", cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, 
                   show_column_names = TRUE, row_names_side = "left", column_names_side = "bottom", 
                   column_names_rot = 45, row_labels = rowLabels, row_order = order(rowSums(!is.na(hitMatrix))),
                   column_order = order(colSums(!is.na(hitMatrix)), decreasing = TRUE),row_names_max_width=unit(10, "cm"), 
                   na_col = "white", border = TRUE, rect_gp = gpar(col="white"), 
                   heatmap_legend_param =list(at=c(minFC, 0, maxFC), labels=c(minFC, 0, maxFC)))

# 调整合适的 PDF 高度
meanRowLength <- meanCharLength(pathwayNames)
f1 <- meanRowLength / 40 %>% ceiling()
pdfHeight <- pathwayNum * f1 * 0.012 %>% ceiling()
if (pdfHeight < 7) {
  pdfHeight = 7
}

pdfName <- paste(baseName, "ORA_Heatplot", "pdf", sep = ".")
pdfPath <- file.path(outputDir, pdfName)
pdf(pdfPath, width = 27, height = pdfHeight)
draw(heatMap)
dev.off()

writeLines("\nヽ(✿ﾟ▽ﾟ)ノ")