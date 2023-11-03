library(stringr) 
library(assertr)
library(matrixStats)
library(reshape)
library(data.table)
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(gplots)
library(ggplot2)
library(optparse)
set.seed(1)

################################################################# Functions ####
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

plot_averageogram <- function( data, label, units, filename ){
  
  # data <- this_cluster_data[, bw_col]
  # label <- bigwig_labels[[i]]
  # units <- paste0("log10_", bigwig_units[[i]])
  # filename <- paste0(opt$out_cluster_prefix, cluster, "_", 
  #                   bigwig_labels[[i]], "_averageogram.pdf" )
  
  averageogram_dat <- data.frame( x = 1:ncol(data), 
                                  coverage = rowSums( t( data ) ) )

  p1 <- ggplot( data = averageogram_dat, aes( x = x, y =  coverage ) ) +
                xlab("best motif hit") + ylab(units) + ggtitle(label) +
                geom_line() + theme_light() +
                theme(axis.text.x=element_blank(),
                      plot.title = element_text(hjust = 0.5) )
  
  ggsave( filename = filename, plot = p1, width = 5, height = 5 )
}
#####

########################################################## IN and load data ####
option_list = list(

  make_option(c("-a", "--out_dir"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_200/align_multivalent_sites",
              help=""),
  
  make_option(c("-b", "--cutoff"), type="character",
              default=0.2,
              help=""),

  make_option(c("-c", "--minsize"), type="character",
              default=0.1,
              help=""),

  make_option(c("-d", "--weighted_PFM"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_200/CTCF_top_200_graphs_weighted_PFM_scores.txt",
              help=""),

  make_option(c("-e", "--ZF_binding_scores"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_200/CTCF_top_200_graphs_ZF_binding_scores.txt",
              help=""),

  make_option(c("-f", "--align_num"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_200/align_multivalent_sites/CTCF_top_200_aligned_sequences_numeric_mx.txt",
              help=""),

  make_option(c("-g", "--title"), type="character",
              default="CTCF_top_200",
              help=""),

  make_option(c("-i", "--computeMatrix"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_200/align_multivalent_sites/CTCF_top_200_16501_CTCF_ChIP1_S368_pulldown_bw_cov_computeMatrix_out.tab.gz,/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_200/align_multivalent_sites/CTCF_top_200_16501_CTCF_ChIP1_S368_pulldown_bw_cov_computeMatrix_out.tab.gz",
              help=""),

  make_option(c("-j", "--repeats_info"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_200/align_multivalent_sites/CTCF_top_200_aligned_positions_overlapping_repeats.bed",
              help=""),
    make_option(c("-k", "--experiment_name"), type="character",
              default="CTCF_top_200",
              help=""),

  make_option(c("-l", "--bw_labels"), type="character",
              default="ChIP_seq_rep1,Mark2_rep2",
              help=""),
  
    make_option(c("-m", "--bw_units"), type="character",
              default="ChIP_reads,Mark2_reads",
              help=""),
  
  make_option(c("-n", "--input_bed"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_200/CTCF_top_200_input_coordinates.bed",
              help=""),
  
  make_option(c("-o", "--meta_pfm_len"), type="character",
              default="30",
              help=""),
  
  make_option(c("-p", "--aligned_bed"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_200/align_multivalent_sites/CTCF_top_200_aligned_positions.bed",
              help="")
  
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)
#####


#################################################################### Setup #####
opt$out_prefix <- paste0( opt$out_dir, "/", opt$experiment_name, "_" )
opt$out_prefix <- paste0( opt$out_prefix, "cluster_cut_", opt$cutoff, "_size_", opt$minsize, "_" )

dir.create( paste0( opt$out_dir, "/clusters" ) )
opt$out_cluster_prefix <- paste0( opt$out_dir, "/clusters/", opt$experiment_name, "_cluster_cut_", opt$cutoff, "_size_", opt$minsize, "_" )


# constants
## Number of zinc fingers
nzfs <- as.integer(opt$meta_pfm_len)/3
#####


############################################################# Read Sequence ####
# read the sequence
cat("Reading the sequence information ...\n")
seq <- fread( opt$align_num, sep = "\t", data.table = F, header = F )
colnames(seq) <- c("Gene","strand",1:(ncol(seq)-2))
seq_len <- ncol(seq)-2
# remove sequences with N values
seq <- seq[ apply(seq[,3:ncol(seq)], 1, function(x) sum(x==-1) ) == 0, ]

# Zoom into the -20 +20 region around the metaPFM hit
ht_seq <- seq[,(3+as.integer(seq_len/2)-20):(3+as.integer(seq_len/2)+nzfs*3+20)]

seq_cols <- paste0( "seq_", (-20):(nzfs*3+20) )
colnames(ht_seq) <- seq_cols
ht_seq$Gene <- seq$Gene


#### Color scale sequence
# A = 0 # C = 1 # G = 2 # T = 3 # N = -1
ht_seq[,seq_cols][ ht_seq[,seq_cols] == "0" ] <- "A"
ht_seq[,seq_cols][ ht_seq[,seq_cols] == "1" ] <- "C"
ht_seq[,seq_cols][ ht_seq[,seq_cols] == "2" ] <- "G"
ht_seq[,seq_cols][ ht_seq[,seq_cols] == "3" ] <- "T"
ht_seq[,seq_cols][ ht_seq[,seq_cols] == "-1" ] <- "N"

colors_bases <- structure( c("#008000", "#ffb200", "#0000ff", "#ff0000"), 
                           names = c( "A", "G", "C", "T" ) )

#####



##################################################### Read motif HMM matrix ####
# read the weighted HMM matrix, with one column per motif
cat("Reading the HMM scores ...\n")

weighted_hmm <- fread( opt$weighted_PFM, sep = "\t", data.table = F )
weighted_hmm$Gene <- gsub( "CHR", "chr", weighted_hmm$Gene )

# filter columns with fewer than "minsize" entries that pass the "cutoff"
weighted_hmm <- weighted_hmm[ , c(1,2,2+which( apply(weighted_hmm[3:ncol(weighted_hmm)],2,function(x) sum(x>=opt$cutoff) ) >= opt$minsize ) )  ]
nmotifs <- ncol(weighted_hmm)-2
# filter genes that have a max HMM score below cutoff
weighted_hmm <- weighted_hmm[ apply(weighted_hmm[3:(2+nmotifs)],1,max) >= opt$cutoff,  ]
#####



######################################################## Read ZF HMM matrix ####
if ( file.exists(opt$ZF_binding_scores) ) {
  # read the weighted HMM matrix, with one column per ZF
  weighted_zf <-fread( opt$ZF_binding_scores, sep = "\t", data.table = F )
  weighted_zf$Gene <- gsub( "CHR", "chr", weighted_zf$Gene )
  
}else{
  
  ## Make a dummy weighted HMM matrix
  weighted_zf <- data.frame(seq$Gene, 
                            rep.int(1, times = length(seq$Gene)), 
                            matrix(data = 0, ncol = nzfs, nrow = length(seq$Gene)))
  
  colnames(weighted_zf) <- c("Gene", "Status", paste0("V", 1:nzfs) )
  
  zf_used <- str_split_fixed( string = colnames(weighted_hmm)[3], pattern = ":", n = 2)[2]
  zf_used <- str_split_fixed(string = zf_used, pattern = "\\|", n= 2)[1] 
  zf_used <- as.numeric( str_split_fixed(string = zf_used, pattern = "-", n= 2) )

  weighted_zf[, paste0("V", zf_used[1]:zf_used[2]) ] <- 1
}


ht_weighted_zf <- weighted_zf[,3:(2+nzfs)]
zf_cols <- paste0( "ZF", 1:ncol(ht_weighted_zf) )
colnames(ht_weighted_zf) <- zf_cols
ht_weighted_zf$Gene <- weighted_zf$Gene


####################################### ZF usage color scale
col_fun_usage <-colorRamp2( c( 0, 0.5, 1),c( "white", "grey", "black" ) )

#### Color scale Zinc finger
diff_tmp <-  abs( max( ht_weighted_zf[, zf_cols] ) -  min( ht_weighted_zf[, zf_cols] ) )
mid <- min( ht_weighted_zf[, zf_cols] ) + diff_tmp/2

col_fun_zf <- colorRamp2( c( min(ht_weighted_zf[, zf_cols])-1e-6, mid, max(ht_weighted_zf[, zf_cols])+1e-6),
                          c( "white", "red", "black" ) )

# col_fun_zf <- colorRamp2( c( 0, 3, 6),
#                           c( "white", "yellow", "red" ) )

rm(mid, diff_tmp)

#####


data_ht <- merge( x = ht_weighted_zf, y = ht_seq, by = "Gene", all.y = TRUE)

data_ht$utilization_proportion <- ( rowSums(data_ht[,zf_cols]!=0) ) / length(zf_cols)



################################################### Read original peak #########
aligned_bed <- fread( opt$aligned_bed, sep = "\t", data.table = T, header = F )
aligned_bed <- aligned_bed[, c("V2", "V4", "V7")]
colnames(aligned_bed) <- c("coord_meta", "Gene", "cluster")
aligned_bed$cluster <- gsub("-", "_", aligned_bed$cluster)


original_bed <- as.data.frame(fread( opt$input_bed, sep = "\t", data.table = T, header = F ))

if( ncol(original_bed) == 5) {
  has_score <- TRUE
  colnames(original_bed) <- c("chr", "start", "stop", "Gene", "score")
} else{ 
    has_score <- FALSE
    colnames(original_bed) <- c("chr", "start", "stop", "Gene") }


### Calculate the distance between the metaPFM hit and the original summit
original_bed$coord_summit <- original_bed$start + (original_bed$stop - original_bed$start)/2
original_bed <- original_bed[ , !colnames(original_bed) %in% c("chr","start", "stop")]

summit_dist <- merge(x=original_bed, y = aligned_bed, by = "Gene", all = TRUE)
summit_dist$summit_dist <- summit_dist$coord_meta - summit_dist$coord_summit


## Just for development
# summit_dist$summit_dist <- as.integer(summit_dist$summit_dist/10)


data_ht <- merge(x=data_ht, y = summit_dist, by = "Gene", all.x = TRUE)



### Create a color scale if there is a summit score
if (has_score) {
  diff_tmp <-  abs( max( data_ht[, "score"] ) -  min( data_ht[, "score"] ) )
  mid <- min( data_ht[, "score"] ) + diff_tmp/2
  
  col_fun_score <- colorRamp2( c( min(data_ht[, "score"])-1e-6, 
                                  mid, 
                                  max(data_ht[, "score"])+1e-6),
                               c( "white", "red", "black" ) )
  
  rm(mid, diff_tmp)
}

rm(summit_dist, original_bed, aligned_bed)
#####




################################################# Add repeat data  #############
repeats <- as.data.frame(read.csv(file = opt$repeats_info, header = FALSE, sep = "\t"))
repeats <- repeats[, c("V4", "V10", "V11", "V12", "V13")]
colnames(repeats) <- c("Gene", "repName", "repClass", "repFamily", "distance")
repeats$repeat_distance <- abs(repeats$distance)
repeats <- repeats[, c("Gene", "repClass", "distance")]

### 
dummy_repeats <- data.frame( Gene = data_ht$Gene,
                             repClass = rep(times = nrow(data_ht), x = "No_repeats"),
                             distance = rep.int(x = 0, times = nrow(data_ht)) )

repeats <- rbind(repeats, dummy_repeats)
repeats <- repeats[repeats$distance == 0,]
repeats$state <- "1"


already_in_repeats <- repeats[ repeats$repClass != "No_repeats", "Gene" ]

repeats <- repeats[ ( repeats$repClass != "No_repeats" ) |
                    ( ( repeats$repClass == "No_repeats" ) & 
                      ( ! repeats$Gene %in% already_in_repeats ) ), ]


repeats <- cast(repeats[,c("Gene", "repClass", "state")], 
               Gene ~ repClass, value = "state", fun.aggregate=length)

repeat_cols <- colnames(repeats)[-1]

data_ht <- merge(x = data_ht, y = repeats, by = "Gene", all.x = TRUE )
data_ht[,repeat_cols][is.na(data_ht[,repeat_cols])] <- 0


#### Color scale repetitive regions
diff_tmp <-  abs( max( data_ht[, repeat_cols] ) -  min( data_ht[, repeat_cols] ) )
mid <- min( data_ht[, repeat_cols] ) + diff_tmp/2

col_fun_rep <- colorRamp2( c( min(data_ht[, repeat_cols]), mid, max(data_ht[, repeat_cols])),
                          c( "white", "red", "black" ) )


rm(dummy_repeats, repeats, already_in_repeats, diff_tmp, mid)
#####


################################################### Read computeMatrix #########
## Read computeMatrix file(s) (from bigwigs), labels and units 
## In this case is ChIP-seq data but it can work with any bw file

if ( (opt$computeMatrix != "default_none") ) {
  
  bigwig_list <- unlist( strsplit( opt$computeMatrix, "," ) )
  bigwig_labels <- unlist( strsplit( opt$bw_labels, "," ) )
  bigwig_units <- unlist( strsplit( opt$bw_units, "," ) )
  
  list_of_cols <- list()
  list_of_cov_df <- list()
  list_of_col_fun <- list()
  
  
  for (i in 1:length(bigwig_list)){
    
    # i <- 1
    bigwig <- bigwig_list[i]
    
    tmp_bw_cov <- as.data.frame(fread( file = bigwig, skip = 1, sep = "\t" ))
    tmp_bw_cov <- tmp_bw_cov[-c(1,2,3,5,6)]
    tmp_bw_cols <- paste0("bw_", i, "_cov_", 1:(ncol(tmp_bw_cov)-1) )
    colnames(tmp_bw_cov) <- c("Gene", tmp_bw_cols )
    
    tmp_bw_cov[, tmp_bw_cols] <- log10( tmp_bw_cov[, tmp_bw_cols] + 0.1 )
    
    
    ########################################################### Color scale bw
    diff_tmp <-  abs( max( tmp_bw_cov[, tmp_bw_cols] ) -  min( tmp_bw_cov[, tmp_bw_cols] ) )
    mid <- min( tmp_bw_cov[, tmp_bw_cols] ) + diff_tmp/2
    
    
    tmp_col_fun_bw <- colorRamp2( c( min( tmp_bw_cov[, tmp_bw_cols] ), mid, 
                                     max( tmp_bw_cov[, tmp_bw_cols] ) ),
                                  c( "white", "red", "black" ) )
    
    
    list_of_col_fun[[i]] <- tmp_col_fun_bw
    list_of_cols[[i]] <- tmp_bw_cols
    list_of_cov_df[[i]] <- tmp_bw_cov
  }
  
   merged_bw_df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Gene", all = TRUE), 
                          list_of_cov_df)
  
  data_ht <- merge(x = data_ht, y = merged_bw_df, by = "Gene", all.x = TRUE )
  
  
  rm( merged_bw_df, bigwig, tmp_bw_cov, tmp_bw_cols, diff_tmp, mid, 
      tmp_col_fun_bw, list_of_cov_df )
}
#####




############################################################### Clustering #####
data_ht$sum <- rowSums(data_ht[,zf_cols])
data_ht <- data_ht[ order(data_ht$sum),]

#####



######################################################### ZF annotations #######
zf_annotations <- c( rep_len( x = NA,length.out = 20),
                     sort(rep_len( x = 1:nzfs, length.out = 3*nzfs), decreasing = TRUE),
                     rep_len( x = NA,length.out = 21) )

anno_zf_col_fun = colorRamp2(c(1, nzfs), c("grey", "black"))

zf1_column_ha <- HeatmapAnnotation(ZFs = zf_annotations, col = list(ZFs = anno_zf_col_fun), 
                               na_col = "white" )

zf2_annotations <- 1:nzfs

zf2_column_ha <- HeatmapAnnotation(ZFs = zf2_annotations, col = list(ZFs = anno_zf_col_fun), 
                               na_col = "white", show_annotation_name = FALSE,
                               show_legend = FALSE )   
#####



###################################################### Score annotations #######
if (has_score) {
  score_annotation <- rowAnnotation(
    width = unit(3, "cm"),
    MAGIX_score = anno_lines( data_ht$score,
                              border = TRUE,
                              axis_param = list( direction = "reverse") )
    )
}
#####



######################################### Distance from summit to motif  #######
summit_annotation <- rowAnnotation(
  width = unit(5, "cm"),
  distance_to_motif = anno_lines( cbind( data_ht$summit_dist, 
                                         rep.int(x=0, times = nrow(data_ht)),
                                         rep.int(x=3*nzfs, times = nrow(data_ht)) 
                                         ),
                                  border = TRUE,
                                  ylim = c( -200, (200+3*nzfs) ),
                                  axis_param = list( direction = "normal") ) )
    
    
# ha = HeatmapAnnotation(foo = anno_lines(cbind(c(1:5, 1:5), c(5:1, 5:1)), 
#     gp = gpar(col = 2:3), add_points = TRUE, pt_gp = gpar(col = 5:6), pch = c(1, 16)))




#####


########################################################## Draw heatmaps #######
my_use_raster <- TRUE

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
  
  list_bw_htm <- list()
  
  for ( i in 1:length(list_of_col_fun) ) {
    # i <- 1 
    
    list_bw_htm[[i]] <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, list_of_cols[[i]] ] ), 
                                      use_raster = my_use_raster,
                                      width = unit(5, "cm"),
                                      raster_quality = 2,
                                      cluster_columns = FALSE, 
                                      cluster_rows = FALSE,
                                      show_column_names = FALSE,
                                      show_row_names = FALSE,
                                      col = list_of_col_fun[[i]],
                                      name = paste0( "log10_", bigwig_units[[i]] ),
                                      column_title = bigwig_labels[[i]],
                                      border_gp = gpar(col = "black"),
                                      heatmap_legend_param = list(legend_direction = "horizontal"))
  
    }
  
  sum_bw_htms <- Reduce( "+", list_bw_htm )
  
  rm(list_bw_htm)
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
  all_htm <- htm_zf + htm_usage + htm_seq + summit_annotation + htm_rep
  bw_add_width <- 0
} else{
  all_htm <- htm_zf + htm_usage + htm_seq + summit_annotation + htm_rep + sum_bw_htms
  bw_add_width <- 2 * length(bigwig_list)
}

if (has_score) {
  all_htm <- score_annotation + all_htm
  bw_add_width <- bw_add_width + 0.2
}





pdf( file = paste0(opt$out_prefix, "aligned_heatmap_ZF_and_seq.pdf"), 
     width = 22 + bw_add_width, height = 15 )

drawned_all_htm <- draw( all_htm, ht_gap = unit( 0.25, "cm" ), 
                         column_title = paste0( opt$title, "\nn = ", nrow(data_ht)),
                         merge_legends = FALSE,
                         heatmap_legend_side = "right", annotation_legend_side = "right" )

dev.off()


write.table( file = paste0( opt$out_prefix, "aligned_heatmap_ZF_and_seq.tab" ), 
             x = data_ht, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
#####









#################################################### Get fasta per cluster #####
cat("Writting fastas per ZF mode ...\n")
seq_letters <- seq
seq_letters[ seq_letters == "0" ] <- "A"
seq_letters[ seq_letters == "1" ] <- "C"
seq_letters[ seq_letters == "2" ] <- "G"
seq_letters[ seq_letters == "3" ] <- "T"
seq_letters[ seq_letters == "-1" ] <- "N"

seq_letters <- cbind( seq_letters[,1], col_concat(seq_letters[,c(-1,-2)], sep = "") )
colnames(seq_letters) <- c("Gene", "seq")
seq_letters <- as.data.frame(seq_letters)

for ( cluster in unique(data_ht$cluster) ){
  
  # cluster <- unique(data_ht$cluster)[1]
  
  cat(paste0(cluster), "\n" )
  
  this_cluster_data <- data_ht[ data_ht$cluster == cluster, ]
  
  this_cluster_seq <- seq_letters[ seq_letters$Gene %in% this_cluster_data$Gene , ]
  this_cluster_seq$Gene <- paste0( ">", this_cluster_seq$Gene )

  this_cluster_fa <- do.call(rbind, lapply(seq(nrow(this_cluster_seq)), 
                                           function(i) t(this_cluster_seq[i, ])))
  
  write.table( x = this_cluster_fa, 
               file = paste0( opt$out_cluster_prefix, "ZF_", cluster, ".fa" ),
               row.names = FALSE, col.names = FALSE, quote = FALSE )
}
#####



################################################# Per cluster operations #######
if (opt$computeMatrix != "default_none") {
  
  cat("Writting averageogram per ZF mode and bw provided ...\n")
  
  for ( i in 1:length( unique(data_ht$cluster) ) ){
    # i <- 1
    cluster <- unique(data_ht$cluster)[i]
    this_cluster_data <- data_ht[ data_ht$cluster == cluster, ]


    for ( j in 1:length( list_of_cols ) ) {
      # i <- 1
      cat( paste0( bigwig_labels[[j]], ", Cluster: ", cluster, "\n" ) )
      
      
      plot_averageogram( data = this_cluster_data[, list_of_cols[[i]]],
                         label = paste0( "ZF array: ", cluster, "\n", bigwig_labels[[j]] ),
                         units = paste0("log10_", bigwig_units[[j]]),
                         filename = paste0(opt$out_cluster_prefix, "ZF_array_", cluster, "_",
                                           bigwig_labels[[j]], "_averageogram.pdf" ) )
    }
  }
}
#####
















############################################################ Old clustering ####
if (nmotifs > 1) {
  
  ################################################################# Filtering ####
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
  #####

  
  # binarize the weighted_hmm matrix, keeping only the top-scoring motif for each row
  binary_hmm <- weighted_hmm
  
  binary_hmm[,3:(2+nmotifs)] <- t( apply( binary_hmm[,3:(2+nmotifs)], 1, function(x) ((x==max(x))*1) ) )
  
  # for debugging purposes, just to make sure the matrices are now compatible
  sum( seq$Gene != weighted_hmm$Gene )
  sum( seq$Gene != weighted_zf$Gene )
  
  # calculate the distance matrices
  
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
  heatmap.2( x = as.matrix( weighted_zf[,3:(2+nzfs)] ), Rowv=dendrogram$rowDendrogram, margins=c(4,2), Colv = F, dendrogram = "row", trace="none", key=T, breaks=seq(0, 6, length.out=256),col=colorRampPalette( c(rgb(0.95,0.95,0.95),"yellow", "orange", "red","dark red"))(255), labRow = F, xlab = "Zinc finger", ylab = "Peaks", key.title = "", key.xlab = "Binding score", key.ylab = "", density.info="none", lhei = c(0.3,2) )
  dev.off()
  
  write.table(seq$Gene[rev(dendrogram$rowInd)],paste0(opt$out_prefix,"clustered.gene_order.txt"), col.names = F, row.names = F, quote=F, sep="\t")
  write.table(binary_hmm[rev(dendrogram$rowInd),],paste0(opt$out_prefix,"clustered.motif_assignment.txt"), row.names = F, quote=F, sep="\t")
}
#####























