library(reshape)
library(data.table)
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(gplots)
library(data.table)
library(optparse)
set.seed(1)

# The function for calculating hamming distance. It is taken from: https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/
hamming <- function(X) {
  uniqs <- unique(as.vector(X))
  U <- X == uniqs[1]
  H <- t(U) %*% U
  
  for ( uniq in uniqs[-1] ) {
    U <- X == uniq
    H <- H + t(U) %*% U
  }
  nrow(X) - H
}

hamming_dist <- function(x, y){
  x[ x != 0 ] <- 1
  y[ y != 0 ] <- 1
  sum(x != y)
  }

########################################################   IN and load data ####
option_list = list(
  make_option(c("-a", "--out_prefix"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_test/align_multivalent_sites/CTCF_test2_",
              help=""),
  
  make_option(c("-b", "--cutoff"), type="character",
              default=0.2,
              help=""),  

  make_option(c("-c", "--minsize"), type="character",
              default=0.1,
              help=""),

  make_option(c("-d", "--weighted_PFM"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_test/CTCF_test_graphs_weighted_PFM_scores.txt",
              help=""),
  
  make_option(c("-e", "--ZF_binding_scores"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_test/CTCF_test_graphs_ZF_binding_scores.txt",
              help=""),
  
  make_option(c("-f", "--align_num"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_test/align_multivalent_sites/CTCF_test_aligned_sequences_numeric_mx.txt",
              help=""),
  
  make_option(c("-g", "--title"), type="character",
              default="CTCF mini test",
              help=""),
  
  make_option(c("-i", "--computeMatrix"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_test/align_multivalent_sites/CTCF_test_computeMatrix_out.tab.gz",
              help=""),
  
  make_option(c("-j", "--repeats_info"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_test/align_multivalent_sites/CTCF_test_aligned_positions_overlap_repeats.bed",
              help="")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)


opt$out_prefix <- paste0( opt$out_prefix, "cluster_cut_", opt$cutoff, "_size_", opt$minsize, "_" )

# read the sequence
cat("Reading the sequence information ...\n")
seq <- fread( opt$align_num, sep = "\t", data.table = F, header = F )
colnames(seq) <- c("Gene","strand",1:(ncol(seq)-2))
seq_len <- ncol(seq)-2
# remove sequences with N values
seq <- seq[ apply(seq[,3:ncol(seq)], 1, function(x) sum(x==-1) ) == 0, ]

# read the weighted HMM matrix, with one column per motif
cat("Reading the HMM scores ...\n")


weighted_hmm <- fread( opt$weighted_PFM, sep = "\t", data.table = F )
weighted_hmm$Gene <- gsub( "CHR", "chr", weighted_hmm$Gene )



# filter columns with fewer than "minsize" entries that pass the "cutoff"
weighted_hmm <- weighted_hmm[ , c(1,2,2+which( apply(weighted_hmm[3:ncol(weighted_hmm)],2,function(x) sum(x>=opt$cutoff) ) >= opt$minsize ) )  ]
nmotifs <- ncol(weighted_hmm)-2
# filter genes that have a max HMM score below cutoff
weighted_hmm <- weighted_hmm[ apply(weighted_hmm[3:(2+nmotifs)],1,max) >= opt$cutoff,  ]


# read the weighted HMM matrix, with one column per ZF
weighted_zf <-fread( opt$ZF_binding_scores, sep = "\t", data.table = F )
weighted_zf$Gene <- gsub( "CHR", "chr", weighted_zf$Gene )
nzfs <- ncol(weighted_zf)-2

# identify the genes that are present in all the matrices
genes <- seq$Gene
genes <- genes[ genes %in% weighted_hmm$Gene ]
genes <- genes[ genes %in% weighted_zf$Gene ]

# filter the matrices to retain only the shared genes, and then sort by gene name to ensure compatibility
seq <- seq[ seq$Gene %in% genes, ]
seq <- seq[ order(seq$Gene), ]
weighted_hmm <- weighted_hmm[ weighted_hmm$Gene %in% genes, ]
weighted_hmm <- weighted_hmm[ order(weighted_hmm$Gene), ]
weighted_zf <- weighted_zf[ weighted_zf$Gene %in% genes, ]
weighted_zf <- weighted_zf[ order(weighted_zf$Gene), ]

# binarize the weighted_hmm matrix, keeping only the top-scoring motif for each row
binary_hmm <- weighted_hmm
binary_hmm[,3:(2+nmotifs)] <- t( apply( binary_hmm[,3:(2+nmotifs)], 1, function(x) ((x==max(x))*1) ) )


# for debugging purposes, just to make sure the matrices are now compatible
sum( seq$Gene != weighted_hmm$Gene )
sum( seq$Gene != weighted_zf$Gene )

# calculat the distance matrices

cat("Clustering the sequences ...\n")

dist_seq <- hamming( t( as.matrix( seq[,(3+as.integer(seq_len/2)):(3+as.integer(seq_len/2)+nzfs*3)] ) ) )




dist_motifs <- hamming( t( as.matrix( binary_hmm[,3:(2+nmotifs)] ) ) )
dist_zfs <- as.matrix( dist( as.matrix( weighted_zf[,3:(2+nzfs)] ), method="binary" ) )

ratios <- c(1,10000,100)
distmx <-  as.dist( dist_seq*ratios[1] + dist_motifs*ratios[2] + dist_zfs*ratios[3] )

jpeg(file=paste0(opt$out_prefix,"clustered.sequence_heatmap.jpg"),
     width=1000,height=1000)
dendrogram <- heatmap.2( as.matrix( seq[,(3+as.integer(seq_len/2)-20):(3+as.integer(seq_len/2)+nzfs*3+20)] ), hclustfun = function(x) hclust(x,method =  "mcquitty"), distfun = function(x) distmx, margins=c(4,2), Colv = F, dendrogram = "row", trace="none", key=T, breaks=seq(0, 3, length.out=256),col=colorRampPalette( c(rgb(0,0.8,0),rgb(0,0,0.8),rgb(1,0.7,0),rgb(0.8,0,0)) ) (255), labRow = F, xlab = "Position", ylab = "Peaks", key.title = "", key.xlab = "Nucleotides", key.ylab = "", density.info="none", lhei = c(0.3,2) )
dev.off()


jpeg(file=paste0(opt$out_prefix,"clustered.zf_heatmap.jpg"),
     width=1000,height=1000)
heatmap.2( as.matrix( weighted_zf[,3:(2+nzfs)] ), Rowv=dendrogram$rowDendrogram, margins=c(4,2), Colv = F, dendrogram = "row", trace="none", key=T, breaks=seq(0, 6, length.out=256),col=colorRampPalette( c(rgb(0.95,0.95,0.95),"yellow", "orange", "red","dark red"))(255), labRow = F, xlab = "Zinc finger", ylab = "Peaks", key.title = "", key.xlab = "Binding score", key.ylab = "", density.info="none", lhei = c(0.3,2) )
dev.off()


write.table(seq$Gene[rev(dendrogram$rowInd)],paste0(opt$out_prefix,"clustered.gene_order.txt"), col.names = F, row.names = F, quote=F, sep="\t")
write.table(binary_hmm[rev(dendrogram$rowInd),],paste0(opt$out_prefix,"clustered.motif_assignment.txt"), row.names = F, quote=F, sep="\t")




ht_weighted_zf <- weighted_zf[,3:(2+nzfs)]
zf_cols <- paste0( "ZF", 1:ncol(ht_weighted_zf) )
colnames(ht_weighted_zf) <- zf_cols
ht_weighted_zf$Gene <- weighted_zf$Gene


ht_seq <- seq[,(3+as.integer(seq_len/2)-20):(3+as.integer(seq_len/2)+nzfs*3+20)]

seq_cols <- paste0( "seq_", 1:ncol(ht_seq) )
colnames(ht_seq) <- seq_cols
ht_seq$Gene <- seq$Gene


data_ht <- merge(x = ht_weighted_zf, y = ht_seq, by = "Gene" )

data_ht$sum <- rowSums(data_ht[,zf_cols])
data_ht <- data_ht[ order(data_ht$sum),]

data_ht$utilization_proportion <- ( rowSums(data_ht[,zf_cols]!=0) ) / length(zf_cols)


###############################################################################
## In this case is ChIP-seq data but it can work with any bw file
if ( (opt$computeMatrix != "default_none") ) {
  
  bw_cov <- as.data.frame(fread( file = opt$computeMatrix, skip = 1, sep = "\t" ))
  bw_cov <- bw_cov[-c(1,2,3,5,6)]
  bw_cols <- paste0("bw_cov_", 1:(ncol(bw_cov)-1) )
  colnames(bw_cov) <- c("Gene", bw_cols )
  
  data_ht <- merge(x = data_ht, y = bw_cov, by = "Gene", all.x = TRUE )
  
  data_ht[, bw_cols] <- log10( data_ht[, bw_cols] + 0.1 )
  
  ########################################################### Color scale bw
  diff_tmp <-  abs( max( data_ht[, bw_cols] ) -  min( data_ht[, bw_cols] ) )
  mid <- min( data_ht[, bw_cols] ) + diff_tmp/2
  
  col_fun_bw <- colorRamp2( c( min(data_ht[, bw_cols]), mid, max(data_ht[, bw_cols])),
                            c( "white", "red", "black" ) )
  
}

###############################################################################
repeats <- as.data.frame(read.csv(file = opt$repeats_info, header = FALSE, sep = "\t"))
repeats <- repeats[, c("V4", "V10", "V11", "V12", "V13")]
colnames(repeats) <- c("Gene", "repName", "repClass", "repFamily", "distance")
repeats$repeat_distance <- abs(repeats$distance)
repeats <- repeats[ repeats$distance == 0, ]
repeats$state <- "1"

repeats <- cast(repeats[,c("Gene", "repClass", "state")], 
               Gene ~ repClass, value = "state", fun.aggregate=length)

repeat_cols <- colnames(repeats)[-1]

data_ht <- merge(x = data_ht, y = repeats, by = "Gene", all.x = TRUE )
data_ht[,repeat_cols][is.na(data_ht[,repeat_cols])] <- 0


######################################################## Color scale Zinc finger
diff_tmp <-  abs( max( ht_weighted_zf[, zf_cols] ) -  min( ht_weighted_zf[, zf_cols] ) )
mid <- min( ht_weighted_zf[, zf_cols] ) + diff_tmp/2

col_fun_zf <- colorRamp2( c( min(ht_weighted_zf[, zf_cols]), mid, max(ht_weighted_zf[, zf_cols])),
                          c( "white", "red", "black" ) )

# col_fun_zf <- colorRamp2( c( 0, 3, 6),
#                           c( "white", "yellow", "red" ) )


######################################################## Color scale repetitive regions
diff_tmp <-  abs( max( data_ht[, repeat_cols] ) -  min( data_ht[, repeat_cols] ) )
mid <- min( data_ht[, repeat_cols] ) + diff_tmp/2

col_fun_rep <- colorRamp2( c( min(ht_weighted_zf[, zf_cols]), mid, max(data_ht[, repeat_cols])),
                          c( "white", "red", "black" ) )


########################################################### Color scale sequence
# A = 0 # C = 1 # G = 2 # T = 3 # N = -1
data_ht[,seq_cols][ data_ht[,seq_cols] == "0" ] <- "A"
data_ht[,seq_cols][ data_ht[,seq_cols] == "1" ] <- "C"
data_ht[,seq_cols][ data_ht[,seq_cols] == "2" ] <- "G"
data_ht[,seq_cols][ data_ht[,seq_cols] == "3" ] <- "T"
data_ht[,seq_cols][ data_ht[,seq_cols] == "-1" ] <- "N"

colors_bases <- structure( c("#008000", "#ffb200", "#0000ff", "#ff0000"), 
                           names = c( "A", "G", "C", "T" ) )





################################################################################
zf_annotations <- c( rep_len( x = NA,length.out = 20),
                     sort(rep_len( x = 1:nzfs, length.out = 3*nzfs), decreasing = TRUE),
                     rep_len( x = NA,length.out = 21) )

anno_zf_col_fun = colorRamp2(c(1, nzfs), c("grey", "black"))

zf1_column_ha <- HeatmapAnnotation(ZFs = zf_annotations, col = list(ZFs = anno_zf_col_fun), 
                               na_col = "white" )



################################################################################
zf2_annotations <- 1:nzfs

zf2_column_ha <- HeatmapAnnotation(ZFs = zf2_annotations, col = list(ZFs = anno_zf_col_fun), 
                               na_col = "white", show_annotation_name = FALSE,
                               show_legend = FALSE )




my_use_raster <- TRUE

col_fun_usage <-colorRamp2( c( 0, 0.5, 1),c( "white", "grey", "black" ) )


htm_zf <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, zf_cols] ), 
                                   width = unit(4, "cm"),
                                   top_annotation = zf2_column_ha,
                                   use_raster = my_use_raster,
                                   raster_quality = 2,
                                   cluster_columns = FALSE,
                                   cluster_rows = TRUE,
                                   row_dend_reorder = TRUE,
                                   clustering_distance_rows = hamming_dist,
                                   show_column_names =  TRUE,
                                   show_row_names  =  FALSE,
                                   col = col_fun_zf,
                                   name = "Zinc_finger_score",
                                   column_title = "Zinc_finger",
                                   # rect_gp = gpar(col = "grey", lwd = 1),
                                   border_gp = gpar(col = "black"),
                                   heatmap_legend_param = list(legend_direction = "horizontal") )

htm_usage <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, "utilization_proportion"] ), 
                                      width = unit(0.5, "cm"),
                                   use_raster = my_use_raster,
                                   raster_quality = 2,
                                   cluster_columns = FALSE,
                                   cluster_rows = FALSE,
                                   row_dend_reorder = TRUE,
                                   # clustering_distance_rows = hamming_dist,
                                   show_column_names =  TRUE,
                                   show_row_names  =  FALSE,
                                   col = col_fun_usage,
                                   name = "1.utilization_proportion",
                                   column_title = "1",
                                   # rect_gp = gpar(col = "grey", lwd = 1),
                                   border_gp = gpar(col = "black"),
                                   heatmap_legend_param = list(legend_direction = "horizontal") )

htm_seq <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, seq_cols] ), 
                                    width = unit(27, "cm"),
                                    top_annotation = zf1_column_ha,
                                    use_raster = my_use_raster,
                                    raster_quality = 2,
                                    cluster_columns = FALSE, 
                                    cluster_rows = FALSE,
                                    show_column_names = FALSE,
                                    show_row_names  =  FALSE,
                                    col = colors_bases,
                                    name = "Sequence",
                                    column_title = "Sequence",
                                    border_gp = gpar(col = "black"),
                                    heatmap_legend_param = list(legend_direction = "horizontal"))


if ( (opt$computeMatrix != "default_none") ) {
  htm_bw <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, bw_cols] ), 
                                      use_raster = my_use_raster,
                                      width = unit(5, "cm"),
                                      raster_quality = 2,
                                      cluster_columns = FALSE, 
                                      cluster_rows = FALSE,
                                      show_column_names = FALSE,
                                      show_row_names = FALSE,
                                      col = col_fun_bw,
                                      name = "log10_ChIP_reads",
                                      column_title = "ChIP-seq reads\n(+/- 2000 bps)",
                                      border_gp = gpar(col = "black"),
                                      heatmap_legend_param = list(legend_direction = "horizontal"))
}


htm_rep <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, repeat_cols] ), 
                                   width = unit(8, "cm"),
                                    use_raster = my_use_raster,
                                    raster_quality = 2,
                                    cluster_columns = FALSE, 
                                    cluster_rows = FALSE,
                                    show_column_names =  TRUE,
                                    show_row_names  =  FALSE,
                                    col = col_fun_rep,
                                    name = "RepeatClass",
                                    column_title = "RepeatClass",
                                    border_gp = gpar(col = "black"),
                                    heatmap_legend_param = list(legend_direction = "horizontal"))

if (opt$computeMatrix == "default_none") {
  all_htm <- htm_zf + htm_usage + htm_seq + htm_rep
} else{
  all_htm <- htm_zf + htm_usage + htm_seq + htm_bw + htm_rep
}


pdf( file = paste0(opt$out_prefix, "aligned_heatmap_ZF_and_seq.pdf"), 
     width = 20, height = 15 )

draw( all_htm, ht_gap = unit( 0.25, "cm" ), 
      column_title = paste0( opt$title, "\nn = ", nrow(data_ht)),
      merge_legends = FALSE,
      heatmap_legend_side = "right", annotation_legend_side = "right" )

dev.off()

write.table(file = paste0(opt$out_prefix, "aligned_heatmap_ZF_and_seq.tab"), 
            x = data_ht, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)















