# 从 DESeq2 差异基因结果生成火山图

# 需要的依赖包：
# argparse, tidyverse

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

scriptDescription <- "展示差异基因的火山图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--DEGs", dest="DEGs", help="DESeq2 差异基因分析结果", required=TRUE)
parser$add_argument("--output-dir", dest="OUT", help="输出目录，默认 当前目录", default=".")
parser$add_argument("--basename", dest="BASENAME", help="输出文件名，默认 \"Unknown\"", default="Unknown")
parser$add_argument("--plot-title", dest="TITLE", help="图像标题，默认 \"Differential Expression Genes\"", default="Differential Expression Genes")
parser$add_argument("--pvalue-cutoff", dest="PVAL", help="P 值阈值，默认 0.05", default=0.05)
parser$add_argument("--ratio-cutoff", dest="RATIO", help="差异倍数（log2）绝对值阈值，默认 1", default=1)
parser$add_argument("--y-max", dest="YMAX", help="Y 轴最大值，超过这个值的点将被往下调到此值，默认 20", default=20)
parser$add_argument("--x-max", dest="XMAX", help="X 轴最大值，注意脚本将输出对称 X 轴，所以这个设置将同时限制最大和最小值，默认 6", default=6)
parser$add_argument("--color1", dest="COL1", help="显著差异基因的颜色，默认 \"#CD6155\"", default="#CD6155")
parser$add_argument("--color2", dest="COL2", help="其他基因的颜色，默认 \"#566573\"", default="#566573")
parser$add_argument("--plot-format", dest="FORMAT", help="保存图片格式，默认 pdf; 可选 [pdf, png, svg, eps]", 
                    default="pdf", choices=c("pdf", "png", "svg", "eps"))

argvs <- parser$parse_args()
inPath <- file.path(argvs$DEGs)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASENAME
titleText <- argvs$TITLE
yMax <- as.double(argvs$YMAX)
xMax <- as.double(argvs$XMAX)
pVal <- as.double(argvs$PVAL)
ratioCut <- as.double(argvs$RATIO)
color1 <- argvs$COL1
color2 <- argvs$COL2
plotFormat <- argvs$FORMAT

# 被压缩的点用三角形，以示区别
addShape <- function(p_value, log2_foldchange, max_x, max_y) {
  if (-log10(p_value) <= max_y & (log2_foldchange >= -max_x & log2_foldchange <= max_x)) {
    shape <- "circle"
  } else {
    shape <- "triangle"
  }
  return(shape)
}

# 调整 Y 值
modifyYvalue <- function(p_value, max_y) {
  y <- -log10(p_value)
  if (y > max_y) {
    y <- max_y
  }
  return(y)
}

# 调整 X 值
modifyXvalue <- function(log2_foldchange, max_x) {
    if (log2_foldchange > max_x) {
    new_log2 <- max_x
  } else if (log2_foldchange < (-max_x)) {
    new_log2 <- -max_x
  } else {
    new_log2 <- log2_foldchange
  }
  return(new_log2)
}

degsData <- read_csv(inPath) %>% dplyr::filter(!is.na(padj)) %>% 
  mutate(x=map2_dbl(log2FoldChange, xMax, modifyXvalue), y=map2_dbl(padj, yMax, modifyYvalue), 
         dot_shape=pmap_chr(list(padj, log2FoldChange, xMax, yMax), addShape))
glimpse(degsData)

# 不显示 Legend
# 设置 expand 让图像框不覆盖点
volcanoPlot <- ggplot(degsData, aes(x, y)) + 
      geom_point(aes(colour = (abs(x) >= ratioCut & padj < pVal), shape=dot_shape), alpha = 0.5, show.legend = FALSE) + 
	  scale_colour_manual(values = c("TRUE" = color1, "FALSE" = color2)) + 
    scale_shape_manual(values = c("circle" = "circle", "triangle" = "triangle")) +
	  geom_vline(xintercept = c(-ratioCut, ratioCut), alpha = 0.8, linetype = "dashed") + 
	  geom_hline(yintercept = -log10(pVal), alpha = 0.8, linetype = "dashed") + 
	  labs(x = "Fold Change(log2)", y = "P value(-log10)", title = titleText) + 
    scale_x_continuous(limits = c(-xMax, xMax), expand=expansion(mult = c(0.005))) +
    scale_y_continuous(limits = c(0, yMax), expand=expansion(c(0, 0.005))) +
	  theme_bw() + 
	  theme(panel.grid = element_blank())

plotName <- paste(baseName, "Volcano", plotFormat, sep=".")
ggsave(filename = plotName, plot = volcanoPlot, dpi = 600, device = plotFormat, path=outDir)

writeLines("\no(^▽^)o")