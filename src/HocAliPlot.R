#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Cairo))

# Arguments
option_list <- list(
	make_option(c("-r", "--ref"), dest = "ref", default = ""),
	make_option(c("-q", "--query"), dest = "query", default = ""),
#    make_option(c("-c", "--ref"), dest = "refs", default = ""),
#	make_option(c("-n", "--querys"), dest = "querys", default = ""),
	make_option(c("-o", "--out"), dest = "out", default = "out")
)

parser <- OptionParser(usage = "Rscript HocAliPlot.R [options] -r [ref] -q [query]", option_list = option_list)
arguments <- parse_args(parser, positional_arguments = c(0,Inf))

ref <- arguments$options$ref
query <- arguments$options$query
#refs <- arguments$options$refs
#querys <- arguments$options$querys
out <- arguments$options$out

if (ref == "" || query == "") {
  cat("Please provide both ref and query parameters.\n")
  cat("Rscript HocAliPlot.R [options] -r [ref] -q [query]\n")
  quit(status = 1)
}

name_y <- strsplit(basename(ref),".",fixed = T)[[1]][1]
name_x <- strsplit(basename(query),".",fixed = T)[[1]][1]

ord1 <- read.table(ref, col.names = c("hoc", "chr1", "pos1"))
ord2 <- read.table(query, col.names = c("hoc", "chr2", "pos2"))

ord <- left_join(ord1, ord2)
# [,c("hoc","chr1","pos1","chr2","pos2")]
data <- na.omit(ord)

# data$chr1 <- factor(data$chr1, levels = unlist(strsplit(refs,",")))
# data$chr2 <- factor(data$chr2, levels = unlist(strsplit(querys,",")))

changetoM <- function(position) {
	position <- position / 1000000
	paste(position, "M", sep = "")
}

p <- ggplot(data = data, aes(x = pos2, y = pos1)) +
	geom_point(size = 0.5) +
	facet_grid(chr1 ~ chr2, scales = "free", space = "free") +
	theme_grey(base_size = 30) +
	labs(x = name_x, y = name_y) +
	scale_x_continuous(labels = changetoM) +
	scale_y_continuous(labels = changetoM) +
	theme(
		axis.line = element_blank(),
		panel.background = element_blank(),
		panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
		axis.text.y = element_text(colour = "black"),
		legend.position = "none",
		axis.text.x = element_text(angle = 300, hjust = 0, vjust = 1, colour = "black")
	)
ggsave(paste0(out,".pdf"),width = 12,height = 12)
