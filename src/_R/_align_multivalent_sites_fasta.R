library(stringdist)
library(stringr)
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(stringi)
library(data.table)
library(optparse)
set.seed(1)


################################################################   Functions ####
reverse_complement <- function( sequence ){

  y = chartr( 'ATGC', 'TACG', sequence )
  rev_com <- intToUtf8(rev(utf8ToInt(y)))
  return( rev_com )
}


get_seq_col_split <- function( zf_annotations ){
  bin_length <- 10
  
  ref_chars <- c("b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m","n",
                 "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "aa", "ab",
                 "ac", "ad", "ae", "af", "ag", "ah", "ai", "aj", "ak", "al")
  
  i <- 0
  for ( j in zf_annotations ){ if(is.na(j)) {i <- i + 1} else {break}}
  prefix_len <- i
  i <- 0
  for ( j in rev(zf_annotations) ){ if(is.na(j)) {i <- i + 1} else {break}}
  suffix_len <- i
  
  mid_len <- length( zf_annotations[ !is.na( zf_annotations ) ] )
  
  num_bins <- as.numeric( mid_len / bin_length )
  
  
  seq_bins_main <- sort( rep_len( x = ref_chars[ 1: floor( num_bins ) ],
                                  length.out = floor( num_bins )*bin_length ) )
  
  seq_bins_extra_len <- as.integer( mid_len - length(seq_bins_main) )
  
  seq_bins_extra <- rep_len(x = "y", length.out =  seq_bins_extra_len )
  
  
  col_split <- c( rep_len(x = "a", length.out = prefix_len),
                  seq_bins_main,
                  seq_bins_extra,
                  rep_len(x = "z", length.out = suffix_len) )

  return(col_split)
}
#####



########################################################   IN and load data ####
option_list = list(
  make_option(c("-a", "--coordinates"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/ZNF208_YWR_A_CG40NGTTATT_random_2000/align_multivalent_sites/ZNF208_YWR_A_CG40NGTTATT_random_2000_.affimx.position.txt",
              help=""),
  
  make_option(c("-b", "--weighted_PFM_scores"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/ZNF208_YWR_A_CG40NGTTATT_random_2000/ZNF208_YWR_A_CG40NGTTATT_random_2000_graphs_weighted_PFM_scores.txt",
              help=""),  
  
  make_option(c("-d", "--fasta"), type="character", 
              default="/home/ahcorcha/repos/tools/RCADEEM/out/ZNF208_YWR_A_CG40NGTTATT_random_2000/ZNF208_YWR_A_CG40NGTTATT_random_2000_sample_2000.fa",
              help=""),

  make_option(c("-e", "--ZF_binding_scores"), type="character", 
              default="/home/ahcorcha/repos/tools/RCADEEM/out/ZNF208_YWR_A_CG40NGTTATT_random_2000/ZNF208_YWR_A_CG40NGTTATT_random_2000_graphs_ZF_binding_scores.txt",
              help=""),
  
  make_option(c("-f", "--meta_pfm_len"), type="character", 
              default="78",
              help=""),
  
  make_option(c("-g", "--out_dir"), type="character", 
              default="/home/ahcorcha/repos/tools/RCADEEM/out/ZNF208_YWR_A_CG40NGTTATT_random_2000/align_multivalent_sites/",
              help=""),
  
  make_option(c("-i", "--cutoff"), type="character",
              default=0.2,
              help=""),

  make_option(c("-j", "--minsize"), type="character",
              default=0.1,
              help=""),
  
  make_option(c("-k", "--experiment_name"), type="character",
              default="ZNF208_YWR_A_CG40NGTTATT_random_2000",
              help="")
  
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)
#####



##################################################################### Setup ####
opt$out_prefix <- paste0( opt$out_dir, "/", opt$experiment_name, "_cluster_cut_", 
                          opt$cutoff, "_size_", opt$minsize, "_" )

aligned_pos <- paste0( opt$out_prefix, "_aligned_positions.tab")

opt$meta_pfm_len <- as.integer( opt$meta_pfm_len )

## Number of zinc fingers
nzfs <- as.integer(opt$meta_pfm_len)/3

#####



############################### read the position and weighted HMM matrices ####
pos <- fread( opt$coordinates, sep = "\t", data.table = F )
colnames(pos) <- gsub("Name", "Gene", colnames(pos) )
  


weighted_hmm <- fread( opt$weighted_PFM_scores, sep = "\t", data.table = F )
# filter columns with fewer than "minsize" entries that pass the "cutoff"
weighted_hmm <- weighted_hmm[ , c(1,2,2+which( apply(weighted_hmm[3:ncol(weighted_hmm)],2,function(x) sum(x>=opt$cutoff) ) >= opt$minsize ) )  ]
nmotifs <- ncol(weighted_hmm)-2
# filter genes that have a max HMM score below cutoff
weighted_hmm <- weighted_hmm[ apply(weighted_hmm[3:(2+nmotifs)],1,max) >= opt$cutoff,  ]
motif_names <- colnames(weighted_hmm)[3:ncol(weighted_hmm)]

pos <- pos[, c("Gene", motif_names, "meta-PFM|opt") ]

#####



################################################# Create aligned pseudo bed ####
# check if the two matrices are compatible
if( ncol(pos) != ncol(weighted_hmm) |
    sum( colnames(pos)[2:(ncol(pos)-1)] != colnames(weighted_hmm)[3:ncol(weighted_hmm)] ) > 0 )
{
  cat("ERROR: Incompatible matrices.\n")
  # quit(status=1)
}else{ cat("Matrices are compatible.\n") }



# get the motif names, their range, and their corresponding offset
nmotifs <- ncol(weighted_hmm) - 2
motif_names <- colnames(weighted_hmm)[3:ncol(weighted_hmm)]
max_zf <- 0
offset <- rep(0,nmotifs)

for( i in 1:nmotifs )
{
  range <- as.numeric( strsplit( strsplit( strsplit( motif_names[i], ":" )[[1]][2], "\\|"  )[[1]][1], "-" )[[1]] )
  offset[i] <- -(range[2]*3)
}
offset <- offset - min(offset)

# merge the hmm and pos matrices to ensure that the rows are in the same order
merged <- merge( x = weighted_hmm[,c(1,3:ncol(weighted_hmm))], y = pos[,1:(ncol(pos)-1)], 
                 by = "Gene" )
# rm(weighted_hmm)
# rm(pos)

# create a binary matrix indicating which motif has the largest weighted hmm score for each sequence
assignments <- t( apply( merged[2:(nmotifs+1)], 1, function(x) (x==max(x))*1  ) )



if( sum( apply(assignments,1,sum) != 1 ) > 0 ){
    cat("ERROR: Ambiguous sequence-motif assignment occurred.\n")
    # quit(status=1)
    }




if( nmotifs == 1 ){
  
  aligned_pos_bed <- cbind( merged$Gene,
                            abs(merged[,ncol(merged)]) - 1,
                            abs(merged[,ncol(merged)]),
                            merged[,2],
                            lapply(merged[,ncol(merged)],function(x) if(x>0) return("+") else return("-") ),
                            paste0("zfs:", range[1], "-", range[2]) )
  
  cat("Only one ZF subset is use\n")
  
} else{
  
  # get the position for each sequence for the assigned motif (also remove the negative signs)
  pos <- apply( assignments * merged[,(2+nmotifs):ncol(merged)], 1, function(x) abs(sum(x)) )
  
  # get the direction (orientation) of the assigned motif for each sequence
  dir <- apply( assignments * merged[,(2+nmotifs):ncol(merged)], 1, function(x) sign(sum(x)) )
  
  # get the genomic coordinates of the 'anchor' position of each sequence; if these anchors are stacked up, it would lead to alignment
  anchors <- pos - dir * assignments %*% offset 
  
  # write the BED file of the anchor positions (0-indexed)
  aligned_pos_bed <- cbind( merged$Gene, anchors-1, anchors,
                            apply( merged[2:(nmotifs+1)], 1, max),
                            lapply(dir,function(x) if(x>0) return("+") else return("-") ),
                            apply(assignments,1,function(x) paste0("zfs:",strsplit( strsplit( names(x)[x==1], ":" )[[1]][2], "\\|"  )[[1]][1] ) )
                            )
}


aligned_pos_bed <- as.data.frame(aligned_pos_bed)
aligned_pos_bed[] <- lapply(aligned_pos_bed, unlist)

colnames(aligned_pos_bed) <- c("Gene", "start", "end", "PFM_hmm_score", "strand", "cluster")

write.table( aligned_pos_bed, aligned_pos, sep = "\t", quote = F, col.names = F, row.names = F )

#####





####################################################### Read FASTA sequence ####
fasta <- read.csv(file = opt$fasta, header = FALSE )

fasta <- data.frame(Gene = fasta[seq(1, nrow(fasta), 2),], 
                    Sequence = fasta[seq(2, nrow(fasta), 2),] )

fasta$Gene <- gsub( "^>", "", fasta$Gene )


### Make sure every sequence is the same length
if ( unique( stri_length( fasta$Sequence ) ) < 1 ){
  cat("Fasta sequences with different lengths\n")
  q() }


data_ht <- merge( x = aligned_pos_bed, y = fasta, by = "Gene")

#####



######################################################## Read ZF HMM matrix ####
cat("Reading ZF HMM matrix ...\n")
if ( file.exists(opt$ZF_binding_scores) ) {
  # read the weighted HMM matrix, with one column per ZF
  weighted_zf <-fread( opt$ZF_binding_scores, sep = "\t", data.table = F )
  weighted_zf$Gene <- gsub( "CHR", "chr", weighted_zf$Gene )
  
}else{
  
  ## Make a dummy weighted HMM matrix
  weighted_zf <- data.frame(fasta$Gene, 
                            rep.int(1, times = length(fasta$Gene)), 
                            matrix(data = 0, ncol = nzfs, nrow = length(fasta$Gene)))
  
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



data_ht <- merge( x = data_ht, y = ht_weighted_zf, by = "Gene")


data_ht$utilization_proportion <- ( rowSums(data_ht[,zf_cols]!=0) ) / length(zf_cols)

#####



########################################################  ####

## Get sequences length
seq_length <- unique( stri_length( data_ht$Sequence ) )


## Return the reverse complement of the sequences with the motif match in the negative strand
data_ht$Sequence <- apply(data_ht, 1, function(x){
  if(x["strand"] == "-"){ return(reverse_complement(x["Sequence"])) } else{ return(x["Sequence"] ) }
  } )

### Reverse starting point of motif matches in the negative strand
data_ht$strand_corrected_start <- data_ht$start

data_ht[ data_ht$strand == "-", "strand_corrected_start"] <-
  seq_length - data_ht[ data_ht$strand == "-", "strand_corrected_start"] - 1


## How many Ns we have to put before and after the read to align the reads
data_ht$prefix <- max( data_ht$strand_corrected_start ) - data_ht$strand_corrected_start
data_ht$suffix <- max(data_ht$prefix) - data_ht$prefix

data_ht$Sequence <- apply( data_ht, 1, function(x){
      new_seq <- paste( paste( rep_len(x = "N",length.out = x["prefix"]), collapse = "" ),
                        x["Sequence"], collapse = "",
                        paste( rep_len(x = "N",length.out = x["suffix"]), collapse = "" )
                        )
      new_seq <- gsub( " ", "", new_seq )
      return(new_seq) } )

new_motif_start <- data_ht$strand_corrected_start[1] + data_ht$prefix[1]

# data_ht <- data_ht[ data_ht$strand == "-", ]
# data_ht <- data_ht[ data_ht$strand_corrected_start == 0, ]


#####



################################## Add sequence as one column per nucleotide ####

seq_ht <- as.data.frame(
  str_split_fixed( string = data_ht$Sequence, pattern = "",
                   n = unique( stri_length( data_ht$Sequence ) ) ) )


seq_cols <- as.character( - new_motif_start:( ncol(seq_ht) - new_motif_start -2 ) )
colnames(seq_ht) <- seq_cols

seq_ht$Gene <- data_ht$Gene

data_ht <- merge( x = data_ht, y = seq_ht, by = "Gene")


colors_bases <- structure( c("#008000", "#ffb200", "#0000ff", "#ff0000", "white"),
                           names = c( "A", "G", "C", "T", "N" ) )

#####



######################################################### ZF annotations #######
suffix_len <- length(seq_cols) - new_motif_start - opt$meta_pfm_len
suffix_len <- max( suffix_len, 0 )

zf_annotations <- c( rep_len( x = NA,length.out = new_motif_start),
                     sort(rep_len( x = 1:nzfs, length.out = 3*nzfs), decreasing = TRUE),
                     rep_len( x = NA,length.out = suffix_len ) )

zf_annotations <- zf_annotations[1:length(seq_cols)]


anno_zf_col_fun = colorRamp2(c(1, nzfs), c("grey", "black"))

zf1_column_ha <- HeatmapAnnotation(ZFs = zf_annotations, col = list(ZFs = anno_zf_col_fun), 
                                   na_col = "white" )

zf2_annotations <- 1:nzfs

zf2_column_ha <- HeatmapAnnotation(ZFs = zf2_annotations, col = list(ZFs = anno_zf_col_fun), 
                               na_col = "white", show_annotation_name = FALSE,
                               show_legend = FALSE )   
#####







################################################################ Clustering ####
rownames(data_ht) <- data_ht$Gene

data_ht$cluster <- gsub("-", "_", data_ht$cluster)

tmp <- gsub("zfs:", "", data_ht$cluster )
tmp <- as.data.frame( str_split_fixed(string = tmp, pattern = "_", n = 2) )
data_ht$first_ZF <- as.integer( tmp$V1 )
data_ht$second_ZF <- as.integer( tmp$V2 )



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

# length(row_order)
## Reorder by new clustering
data_ht <- data_ht[ row_order, ]

#####



######################################### Distance from summit to motif annotation  #######
read_length <- stri_length( gsub("N", "", data_ht$Sequence[1] ) )

data_ht$summit_dist <- ( data_ht$prefix + read_length/2 ) - ( length(seq_cols)/2 )


score_lim <- max( c( abs( max(data_ht$summit_dist) ), abs( min(data_ht$summit_dist) ) ) )
score_limits <- c(-score_lim-5, score_lim+5 )

read_center_anno <- rowAnnotation(
  width = unit(5, "cm"),
  read_center = anno_lines( data_ht$summit_dist,
                           border = TRUE,
                           ylim = score_limits ,
                           axis_param = list( direction = "normal") ) )

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
                                      show_column_names =  TRUE,
                                      col = col_fun_usage,
                                      name = "1.utilization_proportion",
                                      column_title = "1",
                                      # rect_gp = gpar(col = "grey", lwd = 1),
                                      border_gp = gpar(col = "black"),
                                      heatmap_legend_param = list(legend_direction = "horizontal") )

col_split <- get_seq_col_split( zf_annotations )

htm_seq <- ComplexHeatmap::Heatmap( as.matrix( data_ht[, seq_cols] ),
                                    column_split = col_split,
                                    width = unit(27, "cm"),
                                    top_annotation = zf1_column_ha,
                                    use_raster = my_use_raster,
                                    raster_quality = 2,
                                    cluster_columns = FALSE,
                                    show_column_names = TRUE,
                                    col = colors_bases,
                                    name = "Sequence",
                                    column_title = "Sequence",
                                    border_gp = gpar(col = "black"),
                                    heatmap_legend_param = list(legend_direction = "horizontal"))


all_htm <- htm_zf + htm_usage + htm_seq + read_center_anno


pdf( file = paste0(opt$out_prefix, "aligned_heatmap_ZF_and_seq.pdf"),
     width = 22, height = 18 )

drawned_all_htm <- draw( all_htm, ht_gap = unit( 0.25, "cm" ), main_heatmap = "binding_score",
                         column_title = paste0( opt$experiment_name, "\nn = ", nrow(data_ht)),
                         merge_legends = FALSE,
                         heatmap_legend_side = "right", annotation_legend_side = "right" )

dev.off()


write.table( file = paste0( opt$out_prefix, "aligned_heatmap_ZF_and_seq.tab" ), 
             x = data_ht, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )

#####



































