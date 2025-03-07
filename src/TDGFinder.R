#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--one2many"),
    dest = "one2many", default = "",
    help = "one2many file path"
  ),
  make_option(c("-b", "--bed"),
    dest = "bed", default = "",
    help = "bed file path"
  ),
  make_option(c("-o", "--output"),
    dest = "output", default = "results",
    help = "output file path"
  ),
  make_option(c("-t", "--threads"),
    dest = "threads", default = 14,
    help = "threads to use in para"
  ),
  make_option(c("-d", "--distance"),
    dest = "distance", default = 5, 
    help = "the max distance of tow tandem genes"
  ),
  make_option(c("-s", "--subgenome"),
    dest = "subgenome", default = "",
    help = "the file of subgenome"
  )
)

parser <- OptionParser(
  usage = "goldminer TDGFinder [options] file",
  option_list = option_list, description = "Description:
        Step2: Identify clusters in each genome (N = 1) \
        Example: \
        goldminer TDGFinder -f one2many_path -b bed_path -o output_path -s subgenome [ -t 14 ] [ -d 5 ]"
)

arguments <- parse_args(parser, positional_arguments = c(0, Inf))
path_homo <- arguments$options$one2many
path_bed <- arguments$options$bed
path_output <- arguments$options$output
t <- arguments$options$threads
d <- arguments$options$distance
subG.tb <- arguments$options$subgenome

# 检查subgenome参数是否为空
if (arguments$options$subgenome == "") {
  stop("Error: subgenome file path cannot be empty. Please provide a valid file path.")
}else if (!file.exists(subG.tb)) {
  stop("Error: subgenome file does not exist. Please provide a valid file path.")
}else if (arguments$options$one2many == "" | arguments$options$bed == "" ) {
  stop("Error: one2many and bed file path cannot be empty. Please provide a valid file path.")
}

# 检查环境依赖
if (!require(igraph)) {
  stop("Error: package 'igraph' is not installed. Please install it first.")
}
if (!require(reshape2)) {
  stop("Error: package'reshape2' is not installed. Please install it first.")
}
if (!require(dplyr)) {
  stop("Error: package 'dplyr' is not installed. Please install it first.")
}


TDGFinder <- function(i){
  tryCatch({
    ##加载包
    suppressPackageStartupMessages(library(igraph))
    suppressPackageStartupMessages(library(reshape2))
    suppressPackageStartupMessages(library(dplyr))
    
    genome_from <- genomes[i]
    genome_to <- genome_from
    bedfile <- paste0(path_bed,"/",genome_from,".bed")
    singleton <-  paste0(path_homo,"/",genome_from,"_",genome_to,"/",genome_from,"_",genome_to,"itself.singleton")
    homologsfile <- paste0(path_homo,"/",genome_from,"_",genome_to,"/",genome_from,"_",genome_to,"itself.one2many")

    if (file.exists(homologsfile) & file.exists(bedfile) & file.exists(singleton)) {
      ##读入数据
      gene_one2many <- read.table(homologsfile,header = F,col.names = c("from","to","score"))
      gene_singleton <- read.table(singleton,header = F,col.names = c("from"))
      gene_singleton$to <- gene_singleton$from
      gene_singleton$score <- rep(100,nrow(gene_singleton)) 
      gene_pair_data <- rbind(gene_one2many, gene_singleton) 
      gene_bed <- read.table(bedfile,header = F,col.names = c("chr","start","end","id","none","chain"))

      ##筛选同染色体基因对
      names(gene_pair_data) <- c("id","to","score")
      gene_pair_data_bef <- dplyr::left_join(gene_pair_data,gene_bed,by = "id")[,c(1,2,4)]
      
      names(gene_pair_data_bef) <- c("from","id","chr")
      gene_pair_data_aft <- dplyr::left_join(gene_pair_data_bef,gene_bed,by = "id")[,c(1:4)]
      
      gene_pair_data <- gene_pair_data_aft
      names(gene_pair_data) <- c("from","to","from_chr","to_chr")
      gene_pair_chr_match <- droplevels(subset(gene_pair_data,gene_pair_data$from_chr == gene_pair_data$to_chr))
      
      ##清理数据
      rm(gene_pair_data_bef,gene_pair_data_aft,gene_one2many,gene_singleton)

      ##建立基因的index
      gene_num_by_chr <- table(gene_bed$chr)
      index <- list()
      for (i in (1:length(gene_num_by_chr))) {
      index[[i]] <- c(1:gene_num_by_chr[i])
      }
      gene_bed$index <- unlist(index)

      ##计算基因间距离
      gene_bed$cen <- round((gene_bed$start + gene_bed$end)/2)
      gene_bed <- gene_bed[c(4,1:3,7:8,6)]

      ##利用dplyr快速合并基因对与bed信息
      from_data <- data.frame(from = gene_bed$id,f.start = gene_bed$start,f.end = gene_bed$end,
                            f.index = gene_bed$index,f.cen = gene_bed$cen)
      to_data <- data.frame(to = gene_bed$id,t.start = gene_bed$start,t.end = gene_bed$end,
                          t.index = gene_bed$index,t.cen = gene_bed$cen)

      gene_pair_data_pre <- left_join(gene_pair_chr_match,from_data,by = "from")
      gene_pair_data_new <- left_join(gene_pair_data_pre,to_data,by = "to")

      #利用距离筛选数据，d < 5ind
      gene_pair_data_new$dis_cen <- abs(gene_pair_data_new$f.cen - gene_pair_data_new$t.cen)
      gene_pair_data_new$dis_ind <- abs(gene_pair_data_new$f.index - gene_pair_data_new$t.index)
      data_for_igraph <- droplevels(subset(gene_pair_data_new,gene_pair_data_new$dis_ind <= d))

      ##清理数据
      rm(gene_pair_data_new,gene_pair_chr_match,gene_pair_data,from_data,to_data)

      ##创建igraph对象g
      g <- graph_from_data_frame(data_for_igraph, vertices = gene_bed)  

      #利用节点属性计算边属性
      ends_g <- ends(g,E(g))

      E(g)$dis_cen <- abs(V(g)[ends_g[,1]]$cen - V(g)[ends_g[,2]]$cen)
      E(g)$dis_ind <- abs(V(g)[ends_g[,1]]$index - V(g)[ends_g[,2]]$index)
      E(g)$dis <- (abs(V(g)[ends_g[,2]]$start - V(g)[ends_g[,1]]$end) + abs(V(g)[ends_g[,1]]$start - V(g)[ends_g[,2]]$end))/2
      E(g)$weight <- 1/E(g)$dis

      #The edge weights. Larger edge weights increase the probability that an edge is selected by the random walker. 
      #In other words, larger edge weights correspond to stronger connections.

      #利用随机游走模型通过距离权重发现基因簇
      wcg <- walktrap.community(g, weights = E(g)$weight, steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE)
      clu_from_g <- wcg[sizes(wcg) >= 1]

      #存成易读的结构
      cluster <- data.frame(gene = wcg$names,membership = wcg$membership)
      df_clu <- cluster[order(cluster$membership),] 

      names(gene_bed) <- c("gene","chr","start","end","index","cen","chain")
      new_df_clu <- dplyr::left_join(df_clu,gene_bed,by = "gene")[,c(3,4,5,1,2,8,6)]
      last_df_clu <- new_df_clu[order(new_df_clu$chr, new_df_clu$index),]
      clusize <- table(new_df_clu$membership)
      tdif <- ifelse(clusize >= 2,"td","ntd")
      last_df_clu$type <- tdif[last_df_clu$membership]
      
      write.table(last_df_clu,file = paste0(path_output,"/",genome_from,".cludb"),quote = F,row.names = F,col.names = F)
      cleanlist <- function(k,wcg) {
        gene <- wcg[[k]][1]
        len <- length(wcg[[k]])
        gene_end <- wcg[[k]][len]
        cluster <- paste(wcg[[k]], collapse = ",")
        df <- data.frame(gene = gene, gene.e = gene_end, length = len, cluster = cluster)
        return(df)
      }

      new_li <- lapply(1:length(clu_from_g), cleanlist, wcg = clu_from_g)
      l <- melt(new_li, id.vars = c("gene","gene.e","cluster"))[,c(1,2,3,5)]
      clu_li <- dplyr::left_join(l,gene_bed,by = "gene")[,c(1:5)]
      last_df_cluli <- clu_li[order(clu_li$chr, -clu_li$value),]
      write.table(last_df_cluli,file = paste0(path_output,"/",genome_from,".clu"),quote = F,row.names = F,col.names = F)
      
      subg <- dplyr::filter(subG,Gm == as.character(genome_from))
      
      for (i in 1:length(subg$subG)) {
        chrlist <- unlist(strsplit(as.character(subg$CHRlist[i]),split = ","))
        
        Clu <- dplyr::filter(last_df_clu,(chr %in% chrlist))
        Cluli <- dplyr::filter(last_df_cluli,(chr %in% chrlist))
        write.table(Clu,file = paste0(path_output,"/",genome_from,"_",subg$subG[i],".cludb"),quote = F,row.names = F,col.names = F)
        write.table(Cluli,file = paste0(path_output,"/",genome_from,"_",subg$subG[i],".clu"),quote = F,row.names = F,col.names = F)
      }
    }else{
      message("There is no Adequate Data for Tandem Duplicated Genes Determine! Exiting ......")
    }
  }, error = function(e) {
    message(paste("Error in genome:", genome_from, " - ", e$message))
  })
}

# 创建结果存储目录
system(paste0("mkdir -p ",path_output))

## 读入基因组列表
subG <- read.table(subG.tb,col.names = c("Gm","subG","CHRlist"))
genomes <- unique(subG$Gm)

## R 并行
suppressPackageStartupMessages(library(parallel))
l <- length(genomes)
x <- 1:l
if (l >= t) {
  cl <- makeCluster(t)
}else{
  cl <- makeCluster(l)
}

clusterExport(cl, c("path_homo", "path_bed", "path_output", "d", "subG", "genomes"), environment())
parSapply(cl, x, TDGFinder)
stopCluster(cl)
