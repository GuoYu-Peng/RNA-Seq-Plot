# 取 clusterProfiler GSEA 结果画条形图展示，默认展示 P < 0.25 的通路

# 需要的依赖包：
# argparse, tidyverse

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

scriptDescription <- "对 clusterProfiler GSEA 分析结果画条形图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--GSEA", dest="GSEA", help="csv 格式 GSEA 结果", required=TRUE)
parser$add_argument("--output-dir", dest="OUT", help="输出目录，默认当前目录", default=".")
parser$add_argument("--basename", dest="BASENAME", help="输出文件名，默认 \"Unknown\"", default="Unknown")
parser$add_argument("--plot-title", dest="TITLE", help="图像标题，默认 \"GSEA\"", default="GSEA")
parser$add_argument("--pvalue-cutoff", dest="PVAL", help="P 值阈值，默认 0.25", default=0.25)
parser$add_argument("--color1", dest="COL1", help="最大 P 值颜色，默认 \"#F9886B\"", default="#F9886B")
parser$add_argument("--color2", dest="COL2", help="最小 P 值颜色，默认 \"#9E0013\"", default="#9E0013")
parser$add_argument("--plot-format", dest="FORMAT", help="保存图片格式，默认 pdf; 可选 [pdf, png, svg, eps]", 
                    default="pdf", choices=c("pdf", "png", "svg", "eps"))

argvs <- parser$parse_args()
gseaPath <- file.path(argvs$GSEA)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASENAME
plotTitle <- argvs$TITLE
pVal <- as.double(argvs$PVAL)
color1 <- argvs$COL1
color2 <- argvs$COL2
plotFormat <- argvs$FORMAT

plotData <- read_csv(gseaPath) %>% dplyr::filter(`p.adjust` < pVal) %>% dplyr::arrange(desc(enrichmentScore)) %>% 
  mutate(Pathway=factor(Description, levels=Description))
if (nrow(plotData) < 5) {
  writeLines("X﹏X")
  writeLines("通路数目太少，请检查数据和参数！")
  q(save="no")
}

labelWrap <- function(label_text) {
  return(str_wrap(label_text, width = 35))
}

barPlot <- ggplot(plotData, aes(Pathway, enrichmentScore)) +
  geom_bar(aes(fill = `p.adjust`), stat = "identity") +
  scale_fill_gradient(high=color1, low=color2) +
  labs(y="Enrichment Score", title=plotTitle, x="Pathway Name", fill="P value") +
  scale_x_discrete(labels = labelWrap) +
  theme(panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  coord_flip()

plotFileName <- paste(baseName, "GSEA_Barplot", plotFormat, sep=".")
plotHeight <- nrow(plotData) * 6
if (plotHeight < 100) {
  plotHeight = 100
}
ggsave(filename = plotFileName, plot = barPlot, device = plotFormat, path = outDir, width = 220, 
       height = plotHeight, units = "mm", limitsize = FALSE)

writeLines("\n╰(*°▽°*)╯")