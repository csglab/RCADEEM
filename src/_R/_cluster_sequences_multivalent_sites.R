library(stringdist)
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


get_seq_col_split <- function( seq_cols, nzfs, bin_length, spamming = "motif" ){
  # bin_length <- 10
  ref_chars <- c("b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m","n",
                 "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "aa", "ab",
                 "ac", "ad", "ae", "af", "ag", "ah", "ai", "aj", "ak", "al")
  
  motif_length <- as.numeric( as.numeric(nzfs)*3 )
  num_bins <- as.numeric( motif_length / bin_length )

  
  if (spamming == "motif"){
    seq_bins_main <- sort( rep_len( x = ref_chars[ 1: floor( num_bins ) ],
                                   length.out = floor( num_bins )*bin_length ) )
  
    seq_bins_extra_len <- as.integer( motif_length - length(seq_bins_main) )
    
    seq_bins_extra <- rep_len(x = "y", length.out =  seq_bins_extra_len )
    
    col_split <- c( rep_len(x = "a", length.out = 20), 
                    seq_bins_main,
                    seq_bins_extra,
                    rep_len(x = "z", length.out = 21) )
    
  }
  # ifelse( spamming == "whole"){
  #   pass
  # }
  if(length(seq_cols) != length( col_split  ) ){
    cat("ERROR: num. of sequence column should be the same as num. col_split")
    q()
  }  
  return(col_split)
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
                geom_vline(xintercept=0.5+ncol(data)/2, linetype='dashed', color='grey', size=1) +
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
              default="~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/align_multivalent_sites",
              help=""),
  
  make_option(c("-b", "--cutoff"), type="character",
              default=0.2,
              help=""),

  make_option(c("-c", "--minsize"), type="character",
              default=0.1,
              help=""),

  make_option(c("-d", "--weighted_PFM"), type="character",
              default="~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/CTCF_top_2000_RC_range_25_repeats_FALSE_graphs_weighted_PFM_scores.txt",
              help=""),

  make_option(c("-e", "--ZF_binding_scores"), type="character",
              default="~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/CTCF_top_2000_RC_range_25_repeats_FALSE_graphs_ZF_binding_scores.txt",
              help=""),

  make_option(c("-f", "--align_num"), type="character",
              default="~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/align_multivalent_sites/CTCF_top_2000_RC_range_25_repeats_FALSE_aligned_sequences_numeric_mx.txt",
              help=""),

  make_option(c("-g", "--title"), type="character",
              default="CTCF_top_2000_RC_range_25_repeats_FALSE",
              help=""),

  make_option(c("-i", "--computeMatrix"), type="character",
              default="~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/align_multivalent_sites/CTCF_top_2000_RC_range_25_repeats_FALSE_16501_CTCF_ChIP1_S368_pulldown_bw_cov_computeMatrix_out.tab.gz,~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/align_multivalent_sites/CTCF_top_2000_RC_range_25_repeats_FALSE_16501_CTCF_ChIP1_S368_pulldown_bw_cov_computeMatrix_out.tab.gz",
              help=""),

  make_option(c("-j", "--repeats_info"), type="character",
              default="~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/align_multivalent_sites/CTCF_top_2000_RC_range_25_repeats_FALSE_aligned_positions_overlapping_repeats.bed",
              help=""),
    make_option(c("-k", "--experiment_name"), type="character",
              default="CTCF_top_2000_RC_range_25_repeats_FALSE",
              help=""),

  make_option(c("-l", "--bw_labels"), type="character",
              default="hIP_seq,asdasd",
              help=""),
  
    make_option(c("-m", "--bw_units"), type="character",
              default="counts,asdasd",
              help=""),
  
  make_option(c("-n", "--input_bed"), type="character",
              default="~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/CTCF_top_2000_RC_range_25_repeats_FALSE_input_coordinates.bed",
              help=""),
  
  make_option(c("-o", "--meta_pfm_len"), type="character",
              default="33",
              help=""),
  
  make_option(c("-p", "--aligned_bed"), type="character",
              default="~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/align_multivalent_sites/CTCF_top_2000_RC_range_25_repeats_FALSE_aligned_positions_spams_metaPFM.bed",
              help=""),

  make_option(c("-q", "--footprint_tab"), type="character",
              default="~/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/align_multivalent_sites/CTCF_top_2000_RC_range_25_repeats_FALSE_footprint_5_and_3_prime.tab.gz",
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


opt$meta_pfm_len <- as.integer( opt$meta_pfm_len )
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
# ht_seq <- seq[,(3+as.integer(seq_len/2)-20):(3+as.integer(seq_len/2)+nzfs*3+20)]
start <- 3+as.integer(seq_len/2) - ceiling( as.integer(opt$meta_pfm_len)/2 )
ht_seq <- seq[,( start - 20 ):( start + nzfs * 3 + 20 ) ]


# seq_cols <- paste0( "seq_", (-20):(nzfs*3+20) )
seq_cols <- paste0( (-20):(nzfs*3+20) )
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
cat("Reading motif HMM scores ...\n")

weighted_hmm <- fread( opt$weighted_PFM, sep = "\t", data.table = F )
weighted_hmm$Gene <- gsub( "CHR", "chr", weighted_hmm$Gene )

# filter columns with fewer than "minsize" entries that pass the "cutoff"
weighted_hmm <- weighted_hmm[ , c(1,2,2+which( apply(weighted_hmm[3:ncol(weighted_hmm)],2,function(x) sum(x>=opt$cutoff) ) >= opt$minsize ) )  ]
nmotifs <- ncol(weighted_hmm)-2
# filter genes that have a max HMM score below cutoff
weighted_hmm <- weighted_hmm[ apply(weighted_hmm[3:(2+nmotifs)],1,max) >= opt$cutoff,  ]
#####



######################################################## Read ZF HMM matrix ####
cat("Reading ZF HMM matrix ...\n")
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

col_fun_zf <- colorRamp2( c( min(ht_weighted_zf[, zf_cols])-1e-6, mid, 
                             max(ht_weighted_zf[, zf_cols])+1e-6),
                          c( "white", "red", "black" ) )

# col_fun_zf <- colorRamp2( c( 0, 3, 6),
#                           c( "white", "yellow", "red" ) )

rm(mid, diff_tmp)

#####


data_ht <- merge( x = ht_weighted_zf, y = ht_seq, by = "Gene", all.y = TRUE)

data_ht$utilization_proportion <- ( rowSums(data_ht[,zf_cols]!=0) ) / length(zf_cols)



################################################### Read original peak #########
cat("Reading Original bed file with score ...\n")
aligned_bed <- fread( opt$aligned_bed, sep = "\t", data.table = T, header = F )
aligned_bed$coord_meta <- aligned_bed$V2 + (aligned_bed$V3 - aligned_bed$V2)/2
aligned_bed <- aligned_bed[, c("coord_meta", "V4", "V5")]
colnames(aligned_bed) <- c("coord_meta", "Gene", "cluster")
aligned_bed$cluster <- gsub("-", "_", aligned_bed$cluster)

tmp <- gsub("zfs:", "", aligned_bed$cluster )
tmp <- as.data.frame( str_split_fixed(string = tmp, pattern = "_", n = 2) )
aligned_bed$first_ZF <- as.integer( tmp$V1 )
aligned_bed$second_ZF <- as.integer( tmp$V2 )


original_bed <- as.data.frame(fread( opt$input_bed, sep = "\t", data.table = T, header = F ))

if( ncol(original_bed) == 5) {
  has_score <- TRUE
  colnames(original_bed) <- c("chr", "start", "stop", "Gene", "score")
} else{ 
    has_score <- FALSE
    colnames(original_bed) <- c("chr", "start", "stop", "Gene") }


### Calculate the distance between the center of the metaPFM hit and the original summit
original_bed$coord_summit <- original_bed$start + (original_bed$stop - original_bed$start)/2
original_bed <- original_bed[ , !colnames(original_bed) %in% c("chr","start", "stop")]

summit_dist <- merge(x=original_bed, y = aligned_bed, by = "Gene", all = TRUE)
summit_dist$summit_dist <- summit_dist$coord_meta - summit_dist$coord_summit


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



############################################################# Read footprint ####
if( opt$footprint_tab != "default_none" ){
  
  foot_counts <- as.data.frame( fread( file = opt$footprint_tab, skip = 1, sep = "\t" ) )
  foot_counts <- foot_counts[, c( 4, 7:ncol( foot_counts ) ) ]
  foot_counts[ is.na( foot_counts ) ] <- 0
  
  ## Center matrix
  num_cols <- 200 + as.integer(opt$meta_pfm_len)
  
  foot_counts <- foot_counts[, 1:(num_cols+1) ]
  
  footprint_cols <- paste0( "footprints_", 1:( num_cols ) ) 
  
  colnames(foot_counts) <- c("Gene", footprint_cols )
  
  col_fun_foot <- colorRamp2( c( 0, 30, 60), c( "white", "red", "black" ) )
  
  anno_zf_col_fun = colorRamp2(c(1, nzfs), c("grey", "black"))
  
  zf_annotations <- c( rep_len( x = NA,length.out = 100),
                       sort(rep_len( x = 1:nzfs, length.out = 3*nzfs), decreasing = TRUE),
                       rep_len( x = NA,length.out = 100) )
  
  zf3_column_ha <- HeatmapAnnotation( ZFs = zf_annotations, col = list(ZFs = anno_zf_col_fun), 
                                      na_col = "white" )
  
  data_ht <- merge(x = data_ht, y = foot_counts, by = "Gene", all.x = TRUE )
}



#####



################################################# Add repeat data  #############
cat("Reading overlapping repeats ...\n")
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
length(repeat_cols)


rep_ht_width <- 8
if( (length(repeat_cols) == 1) && (repeat_cols == "No_repeats") ){
    
  repeats$Repeats <- 0
  repeat_cols <- c("No_repeats", "Repeats")
  rep_ht_width <- 2
}

data_ht <- merge(x = data_ht, y = repeats, by = "Gene", all.x = TRUE )
data_ht[,repeat_cols][is.na(data_ht[,repeat_cols])] <- 0



#### Color scale repetitive regions
diff_tmp <-  abs( max( data_ht[, repeat_cols] ) -  min( data_ht[, repeat_cols] ) )
mid <- min( data_ht[, repeat_cols] ) + diff_tmp/2

col_fun_rep <- colorRamp2( c( min(data_ht[, repeat_cols]), mid, 
                              max(data_ht[, repeat_cols])),
                           c( "white", "red", "black" ) )


rm(dummy_repeats, repeats, already_in_repeats, diff_tmp, mid)

#####


################################################### Read computeMatrix #########
## Read computeMatrix file(s) (from bigwigs), labels and units 
## In this case is ChIP-seq data but it can work with any bw file
cat("Reading bigwig(s) ...\n")

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
    
    
    tmp_col_fun_bw <- colorRamp2( c( min( tmp_bw_cov[, tmp_bw_cols] )-1e-6, mid, 
                                     max( tmp_bw_cov[, tmp_bw_cols] )+1e-6 ),
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
rownames(data_ht) <- data_ht$Gene
data_ht$sum <- rowSums(data_ht[,zf_cols])
data_ht <- data_ht[ order(data_ht$sum, decreasing = TRUE), ]


data_ht <- data_ht[ order(data_ht$second_ZF, decreasing = FALSE),]
data_ht <- data_ht[ order(data_ht$first_ZF, decreasing = FALSE),]

max_ZF <- max( as.integer( data_ht$second_ZF ) )

row_order <- vector()
cat("Clustering by sequence per ZF array ... \n")
for( cluster in unique(data_ht$cluster)){
  
  cat(cluster)
  cat("\n")
  # cluster <- "zfs:1_6"
  
  first_ZF <- as.integer( data_ht[ data_ht$cluster == cluster, "first_ZF"  ][1] )
  second_ZF <- as.integer( data_ht[ data_ht$cluster == cluster, "second_ZF"  ][1] )
  
  # first_ZF
  # second_ZF
  
  ## https://stackoverflow.com/questions/35540059/how-to-insert-a-distance-matrix-into-r-and-run-hierarchical-clustering
  
  start_pos <- ( 3*(max_ZF-second_ZF) )
  end_pos <- start_pos + 3*( second_ZF-first_ZF+1 )-1
  
  this_seq_cols <- as.character( start_pos:end_pos )
  
  ## Clustering using all the positions in the metaPFM hit
  # this_seq <- data_ht[ data_ht$cluster == cluster, seq_cols[ 21:( 20 + 3*nzfs ) ] ]
  
  ## Clustering with the positions that are recognized with the ZF array
  this_seq <- data_ht[ data_ht$cluster == cluster, this_seq_cols ]
  
  if( nrow(this_seq) >= 2 ){
  
    ## Concatenate base columns into a sequence vector
    complete_seqs <- apply(this_seq, 1, function(x) paste(x, collapse = ""))
    names( complete_seqs ) <- rownames(this_seq)
    
    ## Calc. sequence distance matrix
    dist_seq <- stringdistmatrix( a = complete_seqs, b = complete_seqs, 
                                  method = "osa", useBytes = FALSE,
                                  nthread = 1 ) # lv lcs
    rownames(dist_seq) <- rownames(this_seq)
    
    ## Convert to a distance object
    dist_seq <- as.dist(m = dist_seq, diag = TRUE, upper = TRUE)
    ## Clustering
    hclust_seq <- hclust( d = dist_seq, method = "ward.D2", members = NULL )
    
    row_order <- c(row_order, hclust_seq$labels[hclust_seq$order] )
  } 
  else{ 
    row_order <- c( row_order, rownames( this_seq ) ) 
  }
}

length(row_order)
## Reorder by new clustering
data_ht <- data_ht[ row_order, ]


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



######################################### Distance from summit to motif annotation  #######
# summit_annotation <- rowAnnotation(
#   width = unit(5, "cm"),
#   distance_to_motif = anno_lines( cbind( data_ht$summit_dist, 
#                                          rep.int(x=0, times = nrow(data_ht)),
#                                          rep.int(x=3*nzfs, times = nrow(data_ht)) 
#                                          ),
#                                   border = TRUE,
#                                   ylim = c( -200, (200+3*nzfs) ),
#                                   axis_param = list( direction = "normal") ) )

score_lim <- max( c( abs( max(data_ht$summit_dist) ), abs( min(data_ht$summit_dist) ) ) )
score_limits <- c(-score_lim-5, score_lim+5 )

# summit_annotation <- rowAnnotation(
#   width = unit(5, "cm"),
#   GHT_summit = anno_lines( cbind( data_ht$summit_dist, rep.int(x=0, times = nrow(data_ht)) ),
#                                   border = TRUE,
#                                   ylim = score_limits ,
#                                   axis_param = list( direction = "normal") ) )


summit_annotation <- rowAnnotation(
  width = unit(5, "cm"),
  GHT_summit = anno_lines( data_ht$summit_dist,
                           border = TRUE,
                           ylim = score_limits ,
                           axis_param = list( direction = "normal") ) )


# ha = HeatmapAnnotation(foo = anno_lines(cbind(c(1:5, 1:5), c(5:1, 5:1)), 
#     gp = gpar(col = 2:3), add_points = TRUE, pt_gp = gpar(col = 5:6), pch = c(1, 16)))




#####





########################################################## Draw heatmaps #######
my_use_raster <- TRUE

htm_zf <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, zf_cols] ), 
                                   width = unit(4, "cm"),
                                   cluster_rows = FALSE,
                                   cluster_row_slices = FALSE,
                                   # row_order = row_order,
                                   row_split = factor(data_ht$cluster, unique( data_ht$cluster) ),
                                   row_title_rot = 0,
                                   row_gap = unit(4, "mm"),
                                   top_annotation = zf2_column_ha,
                                   use_raster = my_use_raster,
                                   raster_quality = 2,
                                   cluster_columns = FALSE,
                                   show_column_names =  TRUE,
                                   show_row_names  =  FALSE,
                                   col = col_fun_zf,
                                   name = "binding_score",
                                   column_title = "Zinc_finger",
                                   border_gp = gpar(col = "black"),
                                   heatmap_legend_param = list(legend_direction = "horizontal") )

htm_usage <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, "utilization_proportion"] ), 
                                      width = unit(0.5, "cm"),
                                      use_raster = my_use_raster,
                                      raster_quality = 2,
                                      cluster_columns = FALSE,
                                      show_row_names  =  FALSE,
                                      show_column_names =  TRUE,
                                      col = col_fun_usage,
                                      name = "1.utilization_proportion",
                                      column_title = "1",
                                      # rect_gp = gpar(col = "grey", lwd = 1),
                                      border_gp = gpar(col = "black"),
                                      heatmap_legend_param = list(legend_direction = "horizontal") )



col_split <- get_seq_col_split( seq_cols, nzfs, bin_length = 10 )

htm_seq <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, seq_cols] ),
                                    column_split = col_split,
                                    width = unit(27, "cm"),
                                    top_annotation = zf1_column_ha,
                                    show_row_names  =  FALSE,
                                    use_raster = my_use_raster,
                                    raster_quality = 2,
                                    cluster_columns = FALSE, 
                                    show_column_names = TRUE,
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
                                      show_column_names = FALSE,
                                      show_row_names  =  FALSE,
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
                                    width = unit(rep_ht_width, "cm"),
                                    use_raster = my_use_raster,
                                    raster_quality = 2,
                                    cluster_columns = FALSE, 
                                    show_column_names =  TRUE,
                                    col = col_fun_rep,
                                    show_row_names  =  FALSE,
                                    name = "RepeatClass",
                                    column_title = "RepeatClass",
                                    border_gp = gpar(col = "black"),
                                    heatmap_legend_param = list(legend_direction = "horizontal"))

ft_add_width <- 0
if( opt$footprint_tab != "default_none" ){
  htm_foot <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, footprint_cols] ), 
                                     width = unit(10, "cm"),
                                     top_annotation = zf3_column_ha,
                                     use_raster = my_use_raster,
                                     raster_quality = 2,
                                     cluster_columns = FALSE,
                                     show_column_names =  FALSE,
                                     show_row_names  =  FALSE,
                                     col = col_fun_foot,
                                     name = "footprint_counts",
                                     column_title = "footprint_counts (3' and 5')",
                                     border_gp = gpar(col = "black"),
                                     heatmap_legend_param = list(legend_direction = "horizontal") )
  ft_add_width <- 7
}


all_htm <- htm_zf + htm_usage + htm_seq + summit_annotation + htm_rep


add_width <- 0

if( opt$footprint_tab != "default_none" ){
  all_htm <- all_htm + htm_foot
  add_width <- add_width + 10
}


if (opt$computeMatrix != "default_none") {
  all_htm <- all_htm + sum_bw_htms
  add_width <- add_width + ( 2 * length(bigwig_list) )
}


if (has_score) {
  all_htm <- score_annotation + all_htm
  add_width <- add_width + 0.2
}



pdf( file = paste0(opt$out_prefix, "aligned_heatmap_ZF_and_seq.pdf"), 
     width = 22 + add_width, height = 20 )

drawned_all_htm <- draw( all_htm, ht_gap = unit( 0.25, "cm" ), main_heatmap = "binding_score",
                         column_title = paste0( opt$title, "\nn = ", nrow(data_ht)),
                         merge_legends = FALSE,
                         heatmap_legend_side = "right", annotation_legend_side = "right" )

dev.off()

write.table( file = paste0( opt$out_prefix, "aligned_heatmap_ZF_and_seq.tab" ), 
             x = data_ht, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
#####





##################################################################################################### ###
##################################################################################################### ###




################################################# Per cluster operations #######
if (opt$computeMatrix != "default_none") {

  cat("Writting averageogram per ZF mode and bw provided ...\n")

  for ( i in 1:length( unique(data_ht$cluster) ) ){
    # i <- 1
    cluster <- unique(data_ht$cluster)[i]
    
    this_cluster_data <- data_ht[ data_ht$cluster == cluster, ]

    for ( j in 1:length( list_of_cols ) ) {
      # i <- 1
      cat( paste0( "Cluster: ", cluster, " - ", bigwig_labels[[j]], "\n" ) )
      a <- this_cluster_data[, list_of_cols[[j]]] 
      
        
      plot_averageogram( data = this_cluster_data[, list_of_cols[[j]]],
                         label = paste0( "ZF array: ", cluster, "\n", bigwig_labels[[j]] ),
                         units = paste0("log10_", bigwig_units[[j]]),
                         filename = paste0(opt$out_cluster_prefix, "ZF_array_", cluster, "_",
                                           bigwig_labels[[j]], "_averageogram.pdf" ) )
      
    }
  }
}
#####



############################## Write footprint averageogram per cluster  #######
if (opt$footprint_tab != "default_none") {

  cat("Writting averageogram per ZF mode for footprints ...\n")

  for ( i in 1:length( unique(data_ht$cluster) ) ){
    # i <- 1
    cluster <- unique(data_ht$cluster)[i]
    cat( paste0(cluster, "\n") )
    
    this_cluster_data <- as.data.frame(data_ht[ data_ht$cluster == cluster, footprint_cols ])
    this_cluster_data <- as.data.frame( t( this_cluster_data ) )
    points <- colnames(this_cluster_data)
    
    if (length(points) == 1){
      this_cluster_data$mean <- (this_cluster_data[, points])  
    }else{
      this_cluster_data$mean <- rowMeans(this_cluster_data[, points])  
    }
    
    this_cluster_data$position <- 1:nrow(this_cluster_data)
        
    first_ZF <- as.integer( data_ht[ data_ht$cluster == cluster, "first_ZF"  ][1] )
    second_ZF <- as.integer( data_ht[ data_ht$cluster == cluster, "second_ZF"  ][1] )
  
    start_pos <- ( 3*(max_ZF-second_ZF) )
    end_pos <- start_pos + 3*( second_ZF-first_ZF+1 )-1
    
    this_cluster_data$sd <- rowSds( as.matrix(this_cluster_data[, points] ) )
    this_cluster_data$se <- this_cluster_data$sd / length(points)
    
    
    p1 <- ggplot(this_cluster_data, aes(x=position, y=mean)) + 
      geom_errorbar( aes( ymin = mean - se, ymax = mean + se ), color = "grey" ) +
      geom_line( linetype='solid', color='black', linewidth=0.5 ) +
      ggtitle(paste0( "3' and 5' ends counts for ", opt$title, "\nCluster: ", cluster, "\nn = ", length(points) ) ) +
      
      geom_vline(xintercept=100, linetype='solid', color='red', linewidth=0.5) +
      geom_vline(xintercept=99+opt$meta_pfm_len, linetype='solid', color='red', linewidth=0.5) +
      
      geom_vline(xintercept=100+start_pos, linetype='dashed', color='grey', linewidth=0.5) +
      geom_vline(xintercept=100+end_pos, linetype='dashed', color='grey', linewidth=0.5) +
      
      # geom_point() + 
      theme_light() + 
      labs(caption = paste0( "Error bars - SEM\n", 
                             "Red vertical lines - start/end of metaPFM\n",
                             "Grey dashed lines - start/end of ZF array used for this cluster"))+
      theme(plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust =0) ) +
      NULL
    
    ggsave(filename = paste0(opt$out_cluster_prefix, "ZF_array_",cluster, "_footprint.pdf"), 
           plot = p1, height = 4.5, width = 12 )
  }
}
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



############################################################ Old clustering ####
if (nmotifs > 1) {
  cat("Old clustering the sequences ...\n")
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
  dendrogram <- heatmap.2( as.matrix( seq[,(3+as.integer(seq_len/2)-20):(3+as.integer(seq_len/2)+nzfs*3+20)] ), 
                           hclustfun = function(x) hclust(x,method =  "mcquitty"), 
                           distfun = function(x) distmx, margins=c(4,2), 
                           Colv = F, 
                           dendrogram = "row", trace="none", key=T, 
                           breaks=seq(0, 3, length.out=256),
                           col=colorRampPalette( c(rgb(0,0.8,0),rgb(0,0,0.8),rgb(1,0.7,0),rgb(0.8,0,0)) ) (255), 
                           labRow = F, xlab = "Position", ylab = "Peaks", 
                           key.title = "", key.xlab = "Nucleotides", key.ylab = "", 
                           density.info="none", lhei = c(0.3,2) )
  dev.off()
  
  jpeg(file=paste0(opt$out_prefix,"clustered.zf_heatmap.jpg"),
       width=1000,height=1000)
  heatmap.2( x = as.matrix( weighted_zf[,3:(2+nzfs)] ), 
             Rowv=dendrogram$rowDendrogram, 
             margins=c(4,2), Colv = F, 
             dendrogram = "row", trace="none", key=T, 
             breaks=seq(0, 6, length.out=256),
             col=colorRampPalette( c(rgb(0.95,0.95,0.95),"yellow", "orange", "red","dark red"))(255), 
             labRow = F, xlab = "Zinc finger", ylab = "Peaks", key.title = "", 
             key.xlab = "Binding score", key.ylab = "", density.info="none", 
             lhei = c(0.3,2) )
  dev.off()
  
  write.table(seq$Gene[rev(dendrogram$rowInd)], 
              paste0(opt$out_prefix,"clustered.gene_order.txt"), 
              col.names = F, row.names = F, quote=F, sep="\t")
  
  write.table(binary_hmm[rev(dendrogram$rowInd),], 
              paste0(opt$out_prefix,"clustered.motif_assignment.txt"), 
              row.names = F, quote=F, sep="\t")
}
#####























