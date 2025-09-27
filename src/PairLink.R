#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
option_list <- list(
	make_option(c("-l", "--first"),
							dest = "first", default = "",
							help = "Prefix name of first file"
	),
	make_option(c("-f", "--second"),
							dest = "second", default = "",
							help = "Prefix name of second file"
	),
	make_option(c("-t", "--threads"),
							dest = "threads", default = 4,
							help = "threads to use in para"
	),
	make_option(c("-o", "--out"),
							dest = "out", default = "",
							help = "link file output path"
	),
	make_option(c("-c", "--clu"),
							dest = "clu", default = "/data2/user2/xiexm/projs/GeneODL/pipeline/00.PrepData/ClusterDB_20210821/",
							help = "clu file path"
	),
	make_option(c("-i", "--iblocks"),
							dest = "iblocks", default = "/data/user/chenym/data_fordatabase/microcollinearity_file1/",
							help = "iblocks file path"
	),
	make_option(c("-b", "--blocks"),
							dest = "blocks", default = "/data/user/chenym/shiny_TGT_database/singleBest_colblocks/",
							help = "blocks file path"
	),
	make_option(c("-m", "--tosub"),
							dest = "tosub", default = "",
							help = "Path to tosub.txt file containing genome and chromosome info"
	),
	make_option(c("-w", "--selfweight"),
							dest = "selfweight", default = 1,
							help = "Weight value for self-links"
	)
)

parser <- OptionParser(
	usage = "goldminer PairLink [options] file",
	option_list = option_list, description = "Description:
        Step3: Connect clusters between pairwise genomes (N = 2) \
        Example: \
        goldminer PairLink -l maize_N -f rice_N -o out_path -c clu_path -i iblocks_path -b blocks_path [ -t 14 ] \
        OR use batch mode: \
        goldminer PairLink -m tosub.txt -o out_path -c clu_path -i iblocks_path -b blocks_path [ -t 14 -w 1 ]"
)

arguments <- parse_args(parser, positional_arguments = c(0, Inf))

gm1 <- arguments$options$first
gm2 <- arguments$options$second
cl_num <- arguments$options$threads
outpath <- arguments$options$out
cludb_path <- arguments$options$clu
blocks_path <- arguments$options$iblocks
blocks_path2 <- arguments$options$blocks

# outpath <- "/data2/user2/xiexm/projs/GeneODL/pipeline/00.PrepData/CollinearPairsDB_20220122/"
# cludb_path <- "/data2/user2/xiexm/projs/GeneODL/pipeline/00.PrepData/ClusterDB_20210821/"
# blocks_path <- "/data/user/chenym/data_fordatabase/microcollinearity_file1/"
# blocks_path2 <- "/data/user/chenym/shiny_TGT_database/singleBest_colblocks/"

matchClu <- function(i) {
	blockCheck <- function(TDGclu, from.cludb, to.cludb, blocks, block) {
		TDclu <- unlist(strsplit(as.character(TDGclu), ","))
		SG <- c()

		for (l in 1:length(TDclu)) {
			from.gene <- as.character(TDclu[l])
			if (! from.gene %in% blocks$gm1) {break}
			from.ind <- blocks$ind[blocks$gm1 == from.gene]
			to.gene <- as.character(blocks$gm2[from.ind])

			if (to.gene != ".") {
				from.c <- as.character(from.cludb$Chr[from.cludb$Gene == from.gene])
				from.s <- from.cludb$start[from.cludb$Gene == from.gene]
				from.e <- from.cludb$end[from.cludb$Gene == from.gene]

				to.c <- as.character(to.cludb$Chr[to.cludb$Gene == to.gene])
				to.s <- to.cludb$start[to.cludb$Gene == to.gene]
				to.e <- to.cludb$end[to.cludb$Gene == to.gene]

				from.block <- dplyr::filter(block, from.chr == from.c & from.start <= from.s & from.end >= from.e)
				if (length(from.block$score) != 0) {
					toChr <- unique(from.block$to.chr)
					toStart <- min(from.block$to.start)
					toEnd <- max(from.block$to.end)
					if (toChr == to.c & toStart <= to.s & toEnd >= to.e) {
						for (i in 1:5) {
							ind <- from.ind - i
							up.gene <- ifelse(ind > 0,
																as.character(blocks$gm2[ind]),
																as.character(blocks$gm2[from.ind])
							)
							if (up.gene != ".") { break }
						}
						for (j in 1:5) {
							ind <- from.ind + j
							down.gene <- ifelse(ind < length(blocks$gm2),
																	as.character(blocks$gm2[ind]),
																	as.character(blocks$gm2[from.ind])
							)
							if (down.gene != ".") { break }
						}

						if (up.gene != "." & down.gene != ".") {
							to.ind <- to.cludb$ind2[to.cludb$Gene == as.character(to.gene)]
							up.ind <- to.cludb$ind2[to.cludb$Gene == as.character(up.gene)]
							down.ind <- to.cludb$ind2[to.cludb$Gene == as.character(down.gene)]

							up.dv <- abs(up.ind - to.ind)
							down.dv <- abs(down.ind - to.ind)
							if (up.dv <= 20 & down.dv <= 20) {
								# return(as.character(to.gene))
								sg <- as.character(to.gene)
							} else {
								# return(".")
								sg <- "."
							}
						} else {
							# return(".")
							sg <- "."
						}
					} else {
						# return(".")
						sg <- "."
					}
				} else {
					# return(".")
					sg <- "."
				}
			} else {
				# return(".")
				sg <- "."
			}
			m <- length(sg)
			if (m != 0 & sg != ".") {
				SG <- append(SG, sg)
			}
		}
		return(SG)
	}

	# match Cluid by .clu file [id = row order]
	getInd <- function(SG, Cludb, tdClu) {
		if (is.null(SG)) {
			IDdf <- data.frame(ID = c("."), Freq = c(0), sg = c("."), ind = c(0))
		} else {
			ID <- c()
			for (m in 1:length(SG)) {
				id <- Cludb$cluid[Cludb$Gene == as.character(SG[m])]
				ID <- append(ID, id)
			}
			IDdf <- data.frame(table(ID))

			IDdf$sg <- integer(nrow(IDdf))
			IDdf$ind <- integer(nrow(IDdf))
			for (n in 1:nrow(IDdf)) {
				IDdf$sg[n] <- as.character(Cludb$Gene[Cludb$ind2 == min(Cludb$ind2[Cludb$cluid == IDdf$ID[n]])])
				IDdf$ind[n] <- tdClu$ind[tdClu$firstG == as.character(IDdf$sg[n])]
			}
		}
		return(IDdf)
	}


	TDGclu <- tdClu1$TDGClu[i]
	SG <- blockCheck(
		TDGclu = TDGclu,
		from.cludb = tdCludb1,
		to.cludb = tdCludb2,
		blocks = blockTwo,
		block = blockTwo2
	)
	ind <- getInd(
		SG = SG,
		Cludb = tdCludb2,
		tdClu = tdClu2
	)

	# 当一个簇有多个对应簇时呈现多行####
	for (k in 1:nrow(ind)) {
		out <- paste(c(gm1, tdClu1$ind[i], gm2, ind$ind[k], ind$Freq[k]), collapse = "  ")
		file <- paste0(outpath, "/", gm1, ".", gm2, ".link")
		write.table(out, file, row.names = FALSE, col.names = FALSE, append = T, quote = F)
	}
}

matchClu_old <- function(i) {
	# match Collinearity gene list only MicroCollinearity
	blockCheck <- function(TDGclu, from.cludb, to.cludb, blocks) {
		TDclu <- unlist(strsplit(as.character(TDGclu), ","))
		SG <- c()

		for (l in 1:length(TDclu)) {
			# l = 1
			from.gene <- as.character(TDclu[l])
			if (! from.gene %in% blocks$gm1) {break}
			from.ind <- blocks$ind[blocks$gm1 == from.gene]
			to.gene <- as.character(blocks$gm2[from.ind])

			if (to.gene != ".") {
				for (i in 1:5) {
					ind <- from.ind - i
					up.gene <- ifelse(ind > 0,
														as.character(blocks$gm2[ind]),
														as.character(blocks$gm2[from.ind])
					)
					if (up.gene != ".") { break }
				}
				for (j in 1:5) {
					ind <- from.ind + j
					down.gene <- ifelse(ind < length(blocks$gm2),
															as.character(blocks$gm2[ind]),
															as.character(blocks$gm2[from.ind])
					)
					if (down.gene != ".") { break }
				}

				if (up.gene != "." & down.gene != ".") {
					to.ind <- to.cludb$ind2[to.cludb$Gene == as.character(to.gene)]
					up.ind <- to.cludb$ind2[to.cludb$Gene == as.character(up.gene)]
					down.ind <- to.cludb$ind2[to.cludb$Gene == as.character(down.gene)]

					up.dv <- abs(up.ind - to.ind)
					down.dv <- abs(down.ind - to.ind)
					if (up.dv <= 20 & down.dv <= 20) { # 20 是一个可调参数？
						# return(as.character(to.gene))
						sg <- as.character(to.gene)
					} else {
						# return(".")
						sg <- "."
					}
				} else {
					# return(".")
					sg <- "."
				}
			} else {
				# return(".")
				sg <- "."
			}
			m <- length(sg)
			if (m != 0 & sg != ".") {
				SG <- append(SG, sg)
			}
		}
		return(SG)
	}

	# match Cluid by .clu file [id = row order]
	getInd <- function(SG, Cludb, tdClu) {
		if (is.null(SG)) {
			IDdf <- data.frame(ID = c("."), Freq = c(0), sg = c("."), ind = c(0))
		} else {
			ID <- c()
			for (m in 1:length(SG)) {
				id <- Cludb$cluid[Cludb$Gene == as.character(SG[m])]
				ID <- append(ID, id)
			}
			IDdf <- data.frame(table(ID))

			IDdf$sg <- integer(nrow(IDdf))
			IDdf$ind <- integer(nrow(IDdf))
			for (n in 1:nrow(IDdf)) {
				IDdf$sg[n] <- as.character(Cludb$Gene[Cludb$ind2 == min(Cludb$ind2[Cludb$cluid == IDdf$ID[n]])])
				IDdf$ind[n] <- tdClu$ind[tdClu$firstG == as.character(IDdf$sg[n])]
			}
		}
		return(IDdf)
	}


	TDGclu <- tdClu1$TDGClu[i]
	SG <- blockCheck(
		TDGclu = TDGclu,
		from.cludb = tdCludb1,
		to.cludb = tdCludb2,
		blocks = blockTwo
	)
	ind <- getInd(
		SG = SG,
		Cludb = tdCludb2,
		tdClu = tdClu2
	)

	# 当一个簇有多个对应簇时呈现多行####
	for (k in 1:nrow(ind)) {
		out <- paste(c(gm1, tdClu1$ind[i], gm2, ind$ind[k], ind$Freq[k]), collapse = "  ")
		file <- paste0(outpath, "/", gm1, ".", gm2, ".link")
		write.table(out, file, row.names = FALSE, col.names = FALSE, append = T, quote = F)
	}
}

readCludb <- function(dataPath, gm) {
	cludb <- paste0(dataPath, gm, ".cludb")
	if (file.exists(cludb)) {
		Cludb <- read.table(cludb,
												fill = T,
												col.names = c("Chr", "start", "end", "Gene", "cluid", "chain", "ind", "type")
		)
		Cludb$ind2 <- 1:length(Cludb$Chr)
		return(Cludb)
	} else {
		message(paste0("The cludb file of ", gm, " is not available, So exit."))
		quit()
	}
}


readClu <- function(dataPath, gm, Cludb) {
	clu <- paste0(dataPath, gm, ".clu")
	if (file.exists(clu)) {
		Clu <- read.table(clu,
											fill = T,
											# col.names = c("firstG", "lastG", "TDGClu", "CluSize", "Chr")
											col.names = c("firstG", "lastG", "CluSize","TDGClu",  "Chr")
		)[, c(5, 1:4)]
		Clu$ind <- 1:length(Clu$Chr)

		bedS <- data.frame(firstG = Cludb$Gene, indS = Cludb$ind2)
		bedE <- data.frame(lastG = Cludb$Gene, indE = Cludb$ind2)

		tdClu <- dplyr::left_join(Clu, bedS, by = "firstG")
		Clu <- dplyr::left_join(tdClu, bedE, by = "lastG")
		Clu$serange <- Clu$indE - Clu$indS + 1
		return(Clu)
	} else {
		message(paste0("The clu file of ", gm, " is not available, So exit."))
		quit()
	}
}

# 新添加的批处理函数
processPairGenomes <- function(gm1, gm2, outpath, cludb_path, blocks_path, blocks_path2, cl_num, selfweight) {
	tryCatch({
		if (gm1 == gm2) {
			clu_file <- paste0(cludb_path, "/", gm1, ".clu")
			if (file.exists(clu_file)) {
				out_file <- paste0(outpath, "/", gm1, ".", gm1, ".link")
				cmd <- paste0("awk -v gm=", gm1, " -v w=", selfweight, " '{print gm\"\\t\"NR\"\\t\"gm\"\\t\"NR\"\\t\"w}' ", clu_file, " > ", out_file)
				system(cmd)
				cat(paste0("Created self-link file: ", gm1, "\n"))
			} else {
				cat(paste0("Warning: Clu file for ", gm1, " not found\n"))
			}
		} else {
			out_file <- paste0(outpath, "/", gm1, ".", gm2, ".link")
			if (file.exists(out_file)) {
				cat(paste0("Link file ", out_file, " already exists, skipping...\n"))
				return()
			}
			
			tdCludb1 <- readCludb(dataPath, gm1)
			tdClu1 <- readClu(dataPath, gm1, tdCludb1)
			tdCludb2 <- readCludb(dataPath, gm2)
			tdClu2 <- readClu(dataPath, gm2, tdCludb2)
			
			blocks <- paste0(blocks_path, "/", gm1, ".", gm2, ".i1.blocks")
			blocks2 <- paste0(blocks_path2, "/", gm1, ".", gm2, ".block")
			
			if (file.exists(blocks)) {
				blockTwo <- read.table(blocks, header = F)
				names(blockTwo) <- c("gm1", "gm2")
				blockTwo$ind <- c(1:length(blockTwo$gm1))
				
				x <- 1:length(tdClu1$firstG)
				library(parallel)
				cl <- makeCluster(cl_num)
				
				if (file.exists(blocks2)) {
					blockTwo2 <- read.table(blocks2, header = F)
					names(blockTwo2) <- c(
						"from.chr", "from.start", "from.end", "to.chr",
						"to.start", "to.end", "score", "direction"
					)
					
					clusterExport(cl, c(
						"gm1", "gm2", "tdCludb1", "tdCludb2",
						"tdClu1", "tdClu2", "blockTwo", "blockTwo2", "outpath"
					), environment())
					
					tryCatch({
						parLapply(cl, x, matchClu)
					}, error = function(e) {
						cat(paste0("Error processing ", gm1, " and ", gm2, ": ", conditionMessage(e), "\n"))
					}, finally = {
						if(!is.null(cl) && inherits(cl, "cluster")) {
							tryCatch(stopCluster(cl), error = function(e) NULL)
						}
					})
				} else {
					cat(paste0("Block file [", gm1, " AND ", gm2, "] not available, running matchClu-old.\n"))
					
					clusterExport(cl, c(
						"gm1", "gm2", "tdCludb1", "tdCludb2",
						"tdClu1", "tdClu2", "blockTwo", "outpath"
					), environment())
					
					tryCatch({
						parLapply(cl, x, matchClu_old)
					}, error = function(e) {
						cat(paste0("Error processing ", gm1, " and ", gm2, ": ", conditionMessage(e), "\n"))
					}, finally = {
						if(!is.null(cl) && inherits(cl, "cluster")) {
							tryCatch(stopCluster(cl), error = function(e) NULL)
						}
					})
				}
				
				if(!is.null(cl) && inherits(cl, "cluster")) {
					tryCatch(stopCluster(cl), error = function(e) NULL)
				}
				cat(paste0("Processing completed: links for ", gm1, " and ", gm2, "\n"))
			} else {
				cat(paste0("i1blocks file [", gm1, " AND ", gm2, "] does not exist, skipping.\n"))
			}
		}
	}, error = function(e) {
		cat(paste0("Error in processing genome pair ", gm1, " - ", gm2, ": ", conditionMessage(e), "\n"))
		cat("Skipping this pair and continuing with the next...\n")
	})
}

# 主程序逻辑
dataPath <- paste0(cludb_path, "/")
tosub_file <- arguments$options$tosub
selfweight <- arguments$options$selfweight

if (tosub_file != "") {
	if (!file.exists(tosub_file)) {
		stop(paste0("tosub file ", tosub_file, " does not exist"))
	}
	
	tosub_data <- read.table(tosub_file, sep=" ", stringsAsFactors=FALSE)
	#batch sub mode subfile
	#space split and with 2 cols
	
	genomes <- c()
	for (i in 1:nrow(tosub_data)) {
		genome_id <- paste(tosub_data[i,1], tosub_data[i,2], sep="_")
		genomes <- c(genomes, genome_id)
	}
	
	cat(paste0("Batch mode: Processing ", length(genomes), " genome identifiers\n"))
	cat("List of genome identifiers:\n")
	cat(paste(genomes, collapse="\n"), "\n")
	
	for (i in 1:length(genomes)) {
		for (j in 1:length(genomes)) {
			gm1 <- genomes[i]
			gm2 <- genomes[j]
			cat(paste0("Processing genome pair: ", gm1, " - ", gm2, "\n"))
			processPairGenomes(gm1, gm2, outpath, cludb_path, blocks_path, blocks_path2, cl_num, selfweight)
		}
	}
	
	cat("Batch processing completed\n")
} else {
	# 原有的单对处理模式
	if (gm2 != gm1) {
		tdCludb1 <- readCludb(dataPath, gm1)
		tdClu1 <- readClu(dataPath, gm1, tdCludb1)
		tdCludb2 <- readCludb(dataPath, gm2)
		tdClu2 <- readClu(dataPath, gm2, tdCludb2)
		blocks <- paste0(blocks_path, "/", gm1, ".", gm2, ".i1.blocks")
		blocks2 <- paste0(blocks_path2, "/", gm1, ".", gm2, ".block")
	} else {
		message("第二个基因组与第一个相同，退出。")
		quit()
	}
	print(head(tdCludb1))
	print(head(tdCludb2))
	print(head(tdClu1))
	print(head(tdClu2))
	
	# 原有的处理逻辑保持不变...
	if (file.exists(blocks)) {
		blockTwo <- read.table(blocks, header = F)
		names(blockTwo) <- c("gm1", "gm2")
		blockTwo$ind <- c(1:length(blockTwo$gm1))

		print(head(blockTwo))
	} else {
		message(paste0("i1blocks 文件 [",gm1," 和 ",gm2,"] 不存在，退出。"))
		quit()
	}
	
	x <- 1:length(tdClu1$firstG)
	library(parallel)
	cl <- makeCluster(cl_num)
	
	if (file.exists(blocks2)) {
		blockTwo2 <- read.table(blocks2, header = F)
		names(blockTwo2) <- c(
			"from.chr", "from.start", "from.end", "to.chr",
			"to.start", "to.end", "score", "direction"
		)
		print(head(blockTwo2))
		
		clusterExport(cl, c(
			"gm1", "gm2", "tdCludb1", "tdCludb2",
			"tdClu1", "tdClu2", "blockTwo", "blockTwo2", "outpath"
		), environment())
		parLapply(cl, x, matchClu)
	} else {
		message(paste0("块文件 [",gm1," 和 ",gm2,"] 不可用，运行 matchClu-old。"))
		
		clusterExport(cl, c(
			"gm1", "gm2", "tdCludb1", "tdCludb2",
			"tdClu1", "tdClu2", "blockTwo", "outpath"
		), environment())
		parLapply(cl, x, matchClu_old)
	}
	
	stopCluster(cl)
}
