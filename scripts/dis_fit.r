#!/usr/bin/env Rscript

library(ggplot2)
library(readr)
library(plyr)
library("optparse")
suppressMessages(library(dplyr))

option_list = list(
    make_option(c("-f","--file"), type = "character", default = NULL,
    help = "tsv file name", metavar = "character"),
    make_option(c("-o","--out"), type = "character", default = "output.pdf",
    help = "output pdf script in category name [default = %default]", metavar = "character"),
    make_option(c("-s","--split"), action = "store_true", default = FALSE,
    help = "split by name [default]", metavar = "character"),
    make_option(c("-b","--binwidth"), type = "numeric", default = NULL,
    help = "output histogram bin width", metavar = "numeric"),
    make_option(c("-t","--type"), type = "character", default = NULL,
    help = "output according to different type", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
    print_help(opt_parser)
    stop ("At least one argument must be supplied (input file) .n", call.=FALSE)
}
all <- read_tsv(opt$file, col_names = T, show_col_types = FALSE)

if (opt$type == "type"){
    wild <- subset(all, exist == "wild")
    wildfit <- ddply(wild, "type", summarise, grp.median = median(fit))
    if (is.null(opt$binwidth)){
        pfit <- ggplot(wild, aes(x = fit, fill = type)) +
                geom_histogram(alpha = 0.5, position = "identity") +
                geom_vline(data = wildfit, aes(xintercept = grp.median, color = type), linetype = "dashed")
    } else {
        pfit <- ggplot(wild, aes(x = fit, fill = type)) +
                geom_histogram(binwidth = opt$binwidth, alpha = 0.5, position = "identity") +
                geom_vline(data = wildfit, aes(xintercept = grp.median, color = type), linetype = "dashed")
    }
} else if (opt$type == "exist") {
    existfit <- ddply(all, "exist", summarise, grp.median = median(fit))
    if (is.null(opt$binwidth)){
        pfit <- ggplot(all, aes(x = fit, fill = exist)) +
                geom_histogram(alpha = 0.5, position = "identity") +
                geom_vline(data = existfit, aes(xintercept = grp.median, color = exist), linetype = "dashed")
    } else {
        pfit <- ggplot(all, aes(x = fit, fill = exist)) +
                geom_histogram(binwidth = opt$binwidth, alpha = 0.5, position = "identity") +
                geom_vline(data = existfit, aes(xintercept = grp.median, color = exist), linetype = "dashed")
    }
}

if (opt$split) {
    pfit <- pfit + facet_wrap(~type)
}

pdf(opt$out, width = 6, height = 6)
plot(pfit)
dev.off()
