import os
import csv
from src._python import utils

def align_multivalent_sites( args, log ):

    ####### create bed file of actual peaks
    align_dir = args.OUT_DIR + "/align_multivalent_sites"
    if not os.path.exists(align_dir): os.makedirs(align_dir)
    args.OUT_ALIGN_PREFIX = align_dir + "/" + args.JOB_ID + "_"
    

    lbrack = "{"; rbrack = "}"
    center_bed = args.OUT_ALIGN_PREFIX + "center100.bed" 
    center_fa = args.OUT_ALIGN_PREFIX + "center100.fasta"
    OPT_PFM = args.OUT_PREFIX + "opt_PFM.txt"
    center_aff = args.OUT_ALIGN_PREFIX + "center100.affimx"
    center_aff_pos = args.OUT_ALIGN_PREFIX + "center100.affimx.position.txt"
    coord = args.OUT_ALIGN_PREFIX + "center100_affimx_position_with_coordinates.txt"
    weighted_PFM = args.OUT_PREFIX + "graphs_weighted_PFM_scores.txt"
    align_pos = args.OUT_ALIGN_PREFIX + "aligned_positions.bed"
    align_fa = args.OUT_ALIGN_PREFIX + "aligned_sequences_tab.txt"
    align_num = args.OUT_ALIGN_PREFIX + "aligned_sequences_numeric_mx.txt"
    ZF_binding_score = args.OUT_PREFIX + "graphs_ZF_binding_scores.txt"

    footprint_mat = args.OUT_ALIGN_PREFIX + "footprint_5_and_3_prime.tab.gz"


    align_pos_all_PFM = args.OUT_ALIGN_PREFIX + "aligned_positions_spams_metaPFM.bed"
    overlapped_repeats = args.OUT_ALIGN_PREFIX + "aligned_positions_overlapping_repeats.bed"
    
    bw_flanking_len = 2000
    input_bed = args.OUT_PREFIX + "input_coordinates.bed"
    affimx_range = 100


    #### Get length of meta PFM
    metaPFM_profiles = args.OUT_PREFIX + "metaPFM_profiles.txt"
    with open(metaPFM_profiles) as f:
        reader = csv.reader(f, delimiter='\t')
        first_row = next(reader)
        meta_PFM_len = len(first_row) - 2


    cmdline=f"""awk -v FS="\\t" -v OFS="\\t" '{lbrack} print $1, int( ($2+($3-$2)/2) ), int( ($2+($3-$2)/2) ),$4,"." {rbrack}' {input_bed} |
    sed -e 's/CHR/chr/g' - |
    bedtools slop -b {affimx_range} -g {args.CHR_SIZES} -i - > {center_bed} """
    utils.run_cmd(cmdline, log)

    ####### create fasta file of actual peaks
    cmdline = f""" bedtools getfasta -fi {args.GENOME_FA} -bed {center_bed} -name |
    awk '{lbrack} print toupper($0) {rbrack}' - |
    sed -e 's/CHR/chr/g' - > {center_fa} """
    utils.run_cmd(cmdline, log)

    
    cmdline = f""" {args.script_path}/src/AffiMx -pwm {OPT_PFM} -fasta {center_fa} -out {center_aff} """
    utils.run_cmd(cmdline, log)

    cmdline = f""" cat {center_aff_pos} | 
    awk 'NR==1 {lbrack} printf("Gene\\tchr\\tstart\\tend"); for(i=2;i<=NF;i++) printf("\\t%s",$i); printf("\\n"); {rbrack} NR > 1 {lbrack} split($1,a,"::"); split(a[2],b,":"); split(b[2],c,"-");  printf("%s\\t%s\\t%s\\t%s", a[1],b[1],c[1],c[2]); for(i=2;i<=NF;i++) printf("\\t%s",$i); printf("\\n"); {rbrack}' > {coord} """
    utils.run_cmd(cmdline, log)

    cmdline = f"""Rscript {args.script_path}/src/_R/_align_multivalent_sites.R --coordinates {coord} --weighted_PFM_scores {weighted_PFM} --aligned_pos {align_pos} """
    utils.run_cmd(cmdline, log)



    cmdline = f"""awk -v FS="\\t" -v OFS="\\t" -v len={meta_PFM_len} '{lbrack} if($6 == "-") {lbrack} print $1,$2-len,$3,$4,$7,$6 {rbrack} else {lbrack} print $1,$2,$3+len,$4,$7,$6 {rbrack} {rbrack}' {align_pos} |
    bedtools slop -b 0 -g {args.CHR_SIZES} -i - > {align_pos_all_PFM} """
    utils.run_cmd(cmdline, log)


    cmdline = f"""bedtools slop -b 200 -g {args.CHR_SIZES} -i {align_pos_all_PFM} |
    bedtools getfasta -name -s -tab -fi {args.GENOME_FA} -bed - |
    awk -v FS="\\t" -v OFS="\\t" '{lbrack} print $1,toupper($2) {rbrack}' - > {align_fa} """
    utils.run_cmd(cmdline, log)

    # cmdline = f"""cat {args.CHR_SIZES} {align_pos} | 
    # awk -v FS="\\t" -v OFS="\\t" -v range=200 'NF==2 {lbrack} size[$1]=$2; {rbrack} NF>2 {lbrack} $2 -= range; $3 += range; if ($2>=0 && $3<=size[$1] ) print $0; {rbrack}' | 
    # bedtools getfasta -name -s -tab -fi {args.GENOME_FA} -bed - | 
    # awk -v FS="\\t" -v OFS="\\t" '{lbrack} print $1,toupper($2) {rbrack}' - > {align_fa}"""
    # utils.run_cmd(cmdline, log)



    cmdline = f"""cat {align_fa} | sed 's/A/\\t0/g' | sed 's/C/\\t1/g' | sed 's/G/\\t2/g' | sed 's/T/\\t3/g' | 
    sed 's/N/\\t-1/g' | sed 's/\\t\\t/\\t/g' | sed 's/::[^\\t]*(/\\t/g' | sed 's/)//g' > {align_num} """
    utils.run_cmd(cmdline, log)

    
    
    #### Create bed file with starts and ends for the metaPFM and no strand information.
    # cmdline = f"""awk -v FS="\\t" -v OFS="\\t" -v len=$(( ${meta_PFM_len}/2 )) '{lbrack}print $1,int($2+len),int($3+len),$4,$5,$6,$7{rbrack}' {align_pos} > {align_pos_all_PFM}"""
    # cmdline = f""" awk -v FS="\\t" -v OFS="\\t" -v len={meta_PFM_len} '{lbrack} if($6 == "-") {lbrack} print $1,$2-len,$3,$4,$7,"." {rbrack} else {lbrack} print $1,$2,$3+len,$4,$7,"." {rbrack} {rbrack}' {align_pos} > {align_pos_all_PFM} """
    # utils.run_cmd(cmdline, log)


    #### Get Coverage over BW +/- 2k around centered (on the meta PFM) aligned positions
    if (args.BW_FILE == "default_none"):
        computeMatrix_files_string = "default_none"
    
    else:

        bw_file_list = args.BW_FILE.split(',') 
        computeMatrix_out_files = []

        for bw_file in bw_file_list:

            this_computeMatrix_name = args.OUT_ALIGN_PREFIX + os.path.basename(bw_file).replace(".bw", "").replace(".bigwig", "") + "_bw_cov_computeMatrix_out.tab.gz" 
            computeMatrix_out_files.append( this_computeMatrix_name )

            cmdline = f"""computeMatrix reference-point --referencePoint center --downstream {bw_flanking_len} --upstream {bw_flanking_len} --regionsFileName {align_pos_all_PFM} --scoreFileName {bw_file} --outFileName {this_computeMatrix_name} --binSize 100 --sortRegions keep --averageTypeBins mean --numberOfProcessors 1"""
            utils.run_cmd(cmdline, log)
        
        computeMatrix_files_string = ','.join( computeMatrix_out_files )


    #### Get overlapping repeats
    cmdline = f"""sort -k1,1 -k2,2n {align_pos_all_PFM} |
    bedtools closest -D ref -a - -b {args.REF_REPEATS} > {overlapped_repeats} """
    utils.run_cmd(cmdline, log)


    #### Footprints
    if (args.FOOTPRINT == "default_none"):
        footprint_mat = "default_none"
    else:
        cmdline = f"""computeMatrix reference-point --referencePoint TSS --binSize 1 --upstream 100 --downstream 200 --regionsFileName {align_pos_all_PFM} --scoreFileName {args.FOOTPRINT} --outFileName {footprint_mat} """
        utils.run_cmd(cmdline, log)


    #### Create aligned heatmaps
    cmdline = f""" Rscript {args.script_path}/src/_R/_cluster_sequences_multivalent_sites.R --out_dir {align_dir} --cutoff {args.CUTOFF} --minsize {args.MINSIZE} --weighted_PFM {weighted_PFM} --ZF_binding_scores {ZF_binding_score} --align_num {align_num} --title {args.JOB_ID} --computeMatrix {computeMatrix_files_string} --repeats_info {overlapped_repeats} --meta_pfm_len {meta_PFM_len} --experiment_name {args.JOB_ID} --bw_labels {args.BW_LABELS} --bw_units {args.BW_UNITS} --input_bed {input_bed} --aligned_bed {align_pos_all_PFM} --footprint_tab {footprint_mat}"""
    utils.run_cmd(cmdline, log)



def align_multivalent_sites_fasta( args, log ):

    # input_fasta = args.OUT_PREFIX + "input_sequences.fasta"    
    input_fasta = args.INPUT_FASTA

    align_dir = args.OUT_DIR + "/align_multivalent_sites"
    if not os.path.exists(align_dir): os.makedirs(align_dir)
    args.OUT_ALIGN_PREFIX = align_dir + "/" + args.JOB_ID + "_"
    

    OPT_PFM = args.OUT_PREFIX + "opt_PFM.txt"
    weighted_PFM = args.OUT_PREFIX + "graphs_weighted_PFM_scores.txt"
    ZF_binding_score = args.OUT_PREFIX + "graphs_ZF_binding_scores.txt"
    
    affimx_prefix = args.OUT_ALIGN_PREFIX + ".affimx"
    affimx_pos = affimx_prefix + ".position.txt"


    #### Get length of meta PFM
    metaPFM_profiles = args.OUT_PREFIX + "metaPFM_profiles.txt"
    with open(metaPFM_profiles) as f:
        reader = csv.reader(f, delimiter='\t')
        first_row = next(reader)
        meta_PFM_len = len(first_row) - 2


    cmdline = f""" {args.script_path}/src/AffiMx -pwm {OPT_PFM} -fasta {input_fasta} -out {affimx_prefix} """
    utils.run_cmd(cmdline, log)


    #### Create aligned heatmaps
    cmdline = f""" Rscript {args.script_path}/src/_R/_align_multivalent_sites_fasta.R --coordinates={affimx_pos} --weighted_PFM_scores={weighted_PFM} --fasta={input_fasta} --ZF_binding_scores={ZF_binding_score} --meta_pfm_len={meta_PFM_len} --out_dir={align_dir} --cutoff={args.CUTOFF} --minsize={args.MINSIZE} --experiment_name={args.JOB_ID} """

    utils.run_cmd(cmdline, log)
