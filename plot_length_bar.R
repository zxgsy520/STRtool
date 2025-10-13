#!/usr/bin/env Rscript

library("ggplot2") #绘图
library("argparse") #传参

options(bitmapType="cairo") #关闭服务器与界面的互动响应

add_help_args <- function(){

    version <- "Version: v1.0.0\nAuthor: Xingguo Zhang\tEmail: invicoun@foxmail.com\n"
    desc <-  "Function:Plot a length bar chart distribution\n"
    usage <- "plot_length_bar.R stat_length.tsv -p prefix\n
Input file format:
    Length\tNumber
    17\t20
    18\t45\n"

    parser <- ArgumentParser(description=desc, usage=usage)
    parser$add_argument("-v", "--version", action="version", version=version,
        help="Print version information.")
    parser$add_argument("input", type="character",
        help="Input length statistics table.")
    parser$add_argument("-p", "--prefix", type="character", default="out",
        help="Set the prefix of the figure, default=out.")
    parser$add_argument("--width", type="integer", default=150,
        help="Set the figure width, default=150")
    parser$add_argument("--height", type="integer", default=80,
        help="Set the figure height,default=80")
    parser$add_argument("--maxx", type="integer", default=100,
        help="Set the maximum X range for display, default=100")
    parser$add_argument("--xlab", type="character", default="TRF Length(mer)",
        help="Set the maximum X range for display, default=TRF Length(mer)")

    args <- parser$parse_args()
    return(args)
}


plot_length_bar <- function(data, prefix="out", width=100, height=250, maxx=100, xlab="mer"){

    data <- read.delim(data, sep="\t", head=TRUE, check.names=FALSE, quote="", stringsAsFactors=F)
    colnames(data) <- c("Length", "Count")

    p <- ggplot(data, aes(Length, Count)) + geom_bar(stat="identity", fill="#E64B35",
         color="#E64B35", width=0.5, position=position_dodge(0.1)) + scale_y_log10(expand=c(0, 0))

    if(maxx >=1){
        p <- p + scale_x_continuous(limits=c(min(data$Length), maxx), expand=c(0, 0))
    }else{
        p <- p + scale_x_continuous(expand=c(0, 0))
    }
    p <- p + xlab(xlab) + theme_classic()+
        theme(axis.text=element_text(color="black", size=12),)+
        theme(plot.title=element_text(vjust=1), legend.key=element_blank())

    ggsave(p, file=paste(prefix, ".png", sep=""), units="mm", width=width, height=height, dpi=300)
    ggsave(p, file=paste(prefix, ".pdf", sep=""), units="mm", width=width, height=height)
1
}


args <- add_help_args()
plot_length_bar(data=args$input, prefix=args$prefix, width=args$width, height=args$height, maxx=args$maxx, xlab=args$xlab)
