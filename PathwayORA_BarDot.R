# 从 clusterProfiler 通路富集结果输出泡泡图和柱状图
# 可以选择展示的通路数目，默认按P值排序后筛选相应数目

# 需要以下包依赖：
# argparse, tidyverse, BuenColors


suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(BuenColors))
suppressPackageStartupMessages(library(tidyverse))

scriptDescription <- "生成 clusterProfiler 通路富集结果的泡泡图和柱状图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--enrichment", dest="ENRICHMENT", help="csv 格式通路富集结果", required=TRUE)
parser$add_argument("--output-dir", dest="OUT", help="输出目录，默认为当前目录", default=".")
parser$add_argument("--basename", dest="BASENAME", help="输出文件名，默认 \"Unknown\"", default="Unknown")
parser$add_argument("--show-number", dest="NUMBER", help="展示的通路数目，默认 20", default=20)
parser$add_argument("--plot-title", dest="TITLE", help="图像标题，默认 \"Pathway Enrichment\"", default="Pathway Enrichment")
parser$add_argument("--plot-format", dest="FORMAT", help="保存图片格式，默认 pdf; 可选 [pdf, png, svg, eps]", 
                    default="pdf", choices=c("pdf", "png", "svg", "eps"))

argvs <- parser$parse_args()
inPath <- file.path(argvs$ENRICHMENT)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASENAME
showNum <- as.integer(argvs$NUMBER) 
plotTitle <- argvs$TITLE
plotFormat <- argvs$FORMAT

pathway <- read_csv(inPath) %>% arrange(`p.adjust`) %>% separate(GeneRatio, into=c("k", "n"), sep="/") %>% 
  separate(BgRatio, into=c("M", "N"), sep="/") %>% mutate(RichFactor=as.numeric(k)/as.numeric(M))
if (nrow(pathway) > showNum) {
  pathway <- dplyr::slice(pathway, 1:showNum)
}
glimpse(pathway)

dotPathway <- arrange(pathway, RichFactor) %>% mutate(Description=factor(Description, levels = Description))
barPathway <- arrange(pathway, desc(Count)) %>% mutate(Description=factor(Description, levels = Description))

brewerRed <- rev(jdb_palette("brewer_red"))
lowerCut <- as.integer(length(brewerRed) * 0.2)
upperCut <- as.integer(length(brewerRed) * 0.8)
brewerRed <- brewerRed[lowerCut:upperCut]

dotPlot <- ggplot(dotPathway, aes(x=Description, y=RichFactor)) +
  geom_point(aes(size=Count, colour=`p.adjust`)) +
  scale_colour_gradientn(colors=brewerRed) +
  scale_x_discrete(labels=scales::wrap_format(35)) +
  labs(y="RichFactor", title=plotTitle, size="Gene count", colour="P value") +
  theme(axis.title.y=element_blank(), panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), 
        panel.grid.major=element_line(color="#DCDCDC", linetype="solid")) +
  coord_flip()


barPlot <- ggplot(barPathway, aes(x=Description, y=Count)) +
  geom_bar(aes(fill=`p.adjust`), stat="identity") +
  scale_fill_gradientn(colors=brewerRed) +
  scale_x_discrete(labels=scales::wrap_format(35)) +
  labs(y="Gene count", title=plotTitle, x="Pathway ID", fill="P value") +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())

dotFileName <- paste(baseName, "ORA_Dot", plotFormat, sep=".")
barFileName <- paste(baseName, "ORA_Bar", plotFormat, sep=".")

dotHeight <- showNum * 10 + 30
if (dotHeight < 150) {
  dotHeight<- 150
}
barWidth <- showNum * 10 + 20
if (barWidth < 120) {
  barWidth <- 120
}

ggsave(filename=dotFileName, dpi=600, plot=dotPlot, device=plotFormat, width=150, height=dotHeight, unit = "mm", path=outDir)
ggsave(filename=barFileName, dpi=600, plot=barPlot, device=plotFormat, width=barWidth, height=150, unit = "mm", path=outDir)