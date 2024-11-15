library(data.table)
library(optparse)
set.seed(1)

########################################################   IN and load data ####
option_list = list(
  make_option(c("-a", "--coordinates"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/align_multivalent_sites/CTCF_top_2000_RC_range_25_repeats_FALSE_center100_affimx_position_with_coordinates.txt",
              help=""),
  
  make_option(c("-b", "--weighted_PFM_scores"), type="character",
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/CTCF_top_2000_RC_range_25_repeats_FALSE_graphs_weighted_PFM_scores.txt",
              help=""),  
  
  make_option(c("-c", "--aligned_pos"), type="character", 
              default="/home/ahcorcha/repos/tools/RCADEEM/out/CTCF_top_2000_RC_range_25_repeats_FALSE/align_multivalent_sites/CTCF_top_2000_RC_range_25_repeats_FALSE_aligned_positions.bed",
              help="") );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)
#####


# read the position and weighted HMM matrices
pos <- fread( opt$coordinates, sep = "\t", data.table = F )
weighted_hmm <- fread( opt$weighted_PFM_scores, sep = "\t", data.table = F )
weighted_hmm$Gene <- gsub("CHR", "chr", weighted_hmm$Gene )


# check if the two matrices are compatible
if( ncol(pos) != ncol(weighted_hmm)+3 |
    sum( colnames(pos)[5:(ncol(pos)-1)] != colnames(weighted_hmm)[3:ncol(weighted_hmm)] ) > 0 )
{
  cat("ERROR: Incompatible matrices.\n")
  quit(status=1)
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
merged <- merge( weighted_hmm[,c(1,3:ncol(weighted_hmm))], pos[,1:(ncol(pos)-1)], by="Gene" )
rm(weighted_hmm)
rm(pos)

# create a binary matrix indicating which motif has the largest weighted hmm score for each sequence
assignments <- t( apply( merged[2:(nmotifs+1)], 1, function(x) (x==max(x))*1  ) )


if( nmotifs == 1 ){
  
  aligned_pos_bed <- cbind( merged$chr,
                            merged$start+abs(merged[,ncol(merged)]) - 1,
                            merged$start+abs(merged[,ncol(merged)]),
                            merged$Gene,
                            merged[,2],
                            lapply(merged[,ncol(merged)],function(x) if(x>0) return("+") else return("-") ),
                            paste0("zfs:", range[1], "-", range[2]) )
  
  write.table( x = aligned_pos_bed, file = opt$aligned_pos, 
               sep = "\t", quote = F, col.names = F, row.names = F )
  
  cat("Only one ZF subset is use\n")
  quit(status=0)
    
} else{
  if( sum( apply(assignments,1,sum) != 1 ) > 0 ){
    cat("ERROR: Ambiguous sequence-motif assignment occurred.\n")
    quit(status=1)
    }
}


# get the position for each sequence for the assigned motif (also remove the negative signs)
pos <- apply( assignments * merged[,(5+nmotifs):ncol(merged)], 1, function(x) abs(sum(x)) )

# get the direction (orientation) of the assigned motif for each sequence
dir <- apply( assignments * merged[,(5+nmotifs):ncol(merged)], 1, function(x) sign(sum(x)) )

# get the genomic coordinates of the 'anchor' position of each sequence; if these anchors are stacked up, it would lead to alignment
anchors <- pos - dir * assignments %*% offset + merged$start

# write the BED file of the anchor positions (0-indexed)

write.table(
  cbind(
    merged$chr,
    anchors-1,
    anchors,
    merged$Gene,
    apply( merged[2:(nmotifs+1)], 1, max),
    lapply(dir,function(x) if(x>0) return("+") else return("-") ),
    apply(assignments,1,function(x) paste0("zfs:",strsplit( strsplit( names(x)[x==1], ":" )[[1]][2], "\\|"  )[[1]][1] ) ) ),
  opt$aligned_pos, sep = "\t", quote = F, col.names = F, row.names = F )






