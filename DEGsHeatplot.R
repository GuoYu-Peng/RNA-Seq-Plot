# DESeq2 差异基因分析结果的表达热图
# 假定输入的矩阵都包含 ensembl_gene_id, entrezgene_id, hgnc_symbol 三列
# 默认使用 ensembl_gene_id

# 需要以下的 R 包依赖：
# argparse, tidyverse, ComplexHeatmap, BuenColors

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(BuenColors))


scriptDescription <- "DESeq2 差异基因分析结果的表达热图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--expression", dest="EXPR", help="csv 格式表达文件", required=TRUE)
parser$add_argument("--DEGs", dest="DEGs", help="csv 格式差异基因文件，要求是筛选后的显著差异基因", required=TRUE)
parser$add_argument("--output-dir", dest="OUT", help="输出目录，默认为当前目录", default=".")
parser$add_argument("--basename", dest="BASENAME", help="输出文件名，默认 \"Unknown\"", default="Unknown")
parser$add_argument("--z-score", dest="ZSCORE", help="是否 Z Score 处理，默认 true; 可选 [true, false]", default="true", choices=c("true", "false"))
parser$add_argument("--plot-title", dest="TITLE", help="图像标题，默认 \"DEGs Expression\"", default="DEGs Expression")
parser$add_argument("--color-palette", dest="PALETTE", help="使用 BuenColors 调色板名字，默认 Zissou", default="Zissou")

argvs <- parser$parse_args()
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASENAME
exprPath <- file.path(argvs$EXPR)
degsPath <- file.path(argvs$DEGs)
plotTitle <- argvs$TITLE
zScore <- argvs$ZSCORE
colorPalette <- argvs$PALETTE

# 使用 ensembl_gene_id
degsData <- read_csv(degsPath)
exprData <- read_csv(exprPath) %>% dplyr::select(-entrezgene_id, -hgnc_symbol) %>% dplyr::filter(!is.na(ensembl_gene_id)) %>% 
  dplyr::filter(ensembl_gene_id %in% degsData$ensembl_gene_id) %>% distinct(ensembl_gene_id, .keep_all=TRUE) %>% as.data.frame()
rownames(exprData) <- exprData$ensembl_gene_id
exprData$ensembl_gene_id <- NULL

if (zScore=="true") {
  plotData <- as.matrix(exprData) %>% t() %>% scale() %>% t()
} else {
  plotData <- as.matrix(exprData)
}

colorFun <- jdb_palette(name=colorPalette, type="continuous")
hmPlot <- Heatmap(matrix=plotData, col=colorFun, name="Normalized Expression", cluster_columns=FALSE, column_title=plotTitle, 
                  column_title_side="top", show_column_names=TRUE, column_names_side="top", column_names_rot=45, 
                  cluster_rows=TRUE, show_row_names=FALSE, row_dend_side="left")
pdfHeight <- (7 + (nrow(plotData) / 1000 )) %>% round()
plotFile <- paste(baseName, "DEGsHeatplot", "pdf", sep=".")
plotPath <- file.path(outDir, plotFile)
pdf(plotPath, width=7, height=pdfHeight)
draw(hmPlot)
dev.off()