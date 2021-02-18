# 以 DESeq2 的差异基因和表达数据，分析样品之间关系
# 差异基因结果指过滤后的数据，只包含差异基因
# 假定输入文件都包含 ensembl_gene_id, entrezgene_id, hgnc_symbol
# 将选定 ensembl_gene_id 作为基因 ID 进行分析
# 假定SampleGroup 有 Sample 和 Group 2列

# 需要的依赖包：
# argparse, tidyverse, BuenColors, ComplexHeatmap, ggrepel

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(BuenColors))
suppressPackageStartupMessages(library(ggrepel))

scriptDescription="出图展示 RNA-seq 样品分布"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--expression", dest="EXPR", help="csv 格式的表达数据", required=TRUE)
parser$add_argument("--sample-group", dest="GROUP", help="csv 格式样品分组文件，要求包含 Sample, Group 2 列", required=TRUE)
parser$add_argument("--DEGs", dest="DEGs", help="csv 格式 DESeq2 差异基因分析结果，要求是筛选后认为显著差异的，而不是所有的基因", required=TRUE)
parser$add_argument("--output-dir", dest="OUT", help="输出目录，默认当前目录", default=".")
parser$add_argument("--basename", dest="BASENAME", help="输出文件名，默认 \"Unknown\"", default="Unknown")

argvs <- parser$parse_args()
expressionPath <- file.path(argvs$EXPR)
groupPath <- file.path(argvs$GROUP)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASENAME
degsPath <- file.path(argvs$DEGs)

sampleGroup <- read_csv(groupPath)
degsData <- read_csv(degsPath)
exprData <- read_csv(expressionPath) %>% dplyr::select(-entrezgene_id, -hgnc_symbol) %>% 
  dplyr::filter(ensembl_gene_id %in% degsData$ensembl_gene_id) %>% 
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>% as.data.frame()
rownames(exprData) <- exprData$ensembl_gene_id
exprData$ensembl_gene_id <- NULL

dataPca <- prcomp(t(exprData))
pcaSummary <- summary(dataPca)
print(pcaSummary)
# 取得 PC1,PC2 解释占比
summaryTable <- pcaSummary$importance
pc1Por <- summaryTable[[2, 1]] %>% round(3)
pc2Por <- summaryTable[[2, 2]] %>% round(3)
pcaData <- dataPca$x %>% as_tibble(rownames="Sample") %>% dplyr::select(Sample, PC1, PC2) %>% 
  left_join(sampleGroup, by="Sample")

xLab <- str_glue("PC1({pc1Por})")
yLab <- str_glue("PC2({pc2Por})")
pcaPlot <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
           geom_point(aes(shape=Group)) +
           ggrepel::geom_text_repel(aes(label=Sample)) +
           labs(title = "Sample PCA", x = xLab, y = yLab) +
           theme_bw()

dataDist <- dist(t(exprData), method = "euclidean")
dataHc <- hclust(dataDist, method = "average")

colFun <- jdb_palette("brewer_fire", type = "continuous") %>% rev()
hmData <- as.matrix(dataDist)
hm <- Heatmap(hmData, name = "Euclidean", col = colFun, cluster_rows = TRUE, 
              cluster_columns = TRUE, show_row_names = TRUE, row_names_side = "right", 
              show_row_dend = FALSE, show_column_dend = FALSE, column_title = "Sample Distance", 
              column_title_side = "top")


pcaName <- paste(baseName, "Sample_PCA", "pdf", sep=".")
treeName <- paste(baseName, "Sample_Clust", "pdf", sep=".")
hmName <- paste(baseName, "Sample_HeatPlot", "pdf", sep=".")
treePath <- file.path(outDir, treeName)
hmPath <- file.path(outDir, hmName)

ggsave(filename=pcaName, plot=pcaPlot, device="pdf", path=outDir)
pdf(treePath)
plot(dataHc, xlab = "Sample")
dev.off()

pdf(hmPath)
draw(hm)
dev.off()