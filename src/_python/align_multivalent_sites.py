import os
from src._python import utils

def align_multivalent_sites( args, log ):

    ####### create bed file of actual peaks
    align_dir = args.OUT_DIR + "/align_multivalent_sites"
    if not os.path.exists(align_dir): os.makedirs(align_dir)
    args.OUT_ALIGN_PREFIX = args.OUT_DIR + "/align_multivalent_sites/" + args.JOB_ID + "_"

    lbrack = "{"; rbrack = "}"
    PFM_scores = args.OUT_PREFIX + "PFM_scores.txt"
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



    cmdline = f"""cat {PFM_scores} | 
    awk 'NR>1 && $2==1 {lbrack} split($1,a,":"); split(a[2],b,"-"); 
    printf("%s\\t%s\\t%s\\t%s\\t.\\n",a[1],b[1]+200,b[1]+200+100,$1); {rbrack}' | sed -e 's/CHR/chr/g' - > {center_bed} """    
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

    cmdline = f"""cat {args.CHR_SIZES} {align_pos} | 
    awk -v FS="\\t" -v OFS="\\t" -v range={args.RANGE} 'NF==2 {lbrack} size[$1]=$2; {rbrack} NF>2 {lbrack} $2 -= range; $3 += range; if ($2>=0 && $3<=size[$1] ) print $0; {rbrack}' | 
    bedtools getfasta -name -s -tab -fi {args.GENOME_FA} -bed - | 
    awk -v FS="\\t" -v OFS="\\t" '{lbrack} print $1,toupper($2) {rbrack}' - > {align_fa}"""
    utils.run_cmd(cmdline, log)

    cmdline = f"""cat {align_fa} | sed 's/A/\\t0/g' | sed 's/C/\\t1/g' | sed 's/G/\\t2/g' | sed 's/T/\\t3/g' | 
    sed 's/N/\\t-1/g' | sed 's/\\t\\t/\\t/g' | sed 's/::[^\\t]*(/\\t/g' | sed 's/)//g' > {align_num} """
    utils.run_cmd(cmdline, log)

    cmdline = f""" Rscript {args.script_path}/src/_R/_cluster_sequences_multivalent_sites.R --out_prefix {args.OUT_ALIGN_PREFIX} --cutoff {args.CUTOFF} --minsize {args.MINSIZE} --weighted_PFM {weighted_PFM} --ZF_binding_scores {ZF_binding_score} --align_num {align_num} """
    utils.run_cmd(cmdline, log)





