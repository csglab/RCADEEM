import csv   
import os
from src._python import utils

def end_with_err(log_f, msg):
        
    print(msg)
    with open(log_f, "w") as f:
        f.write(msg)
        exit(2)


def RCADEEM( args, log ):
    
    print("Job ID: " + args.JOB_ID )
    print("Input FASTA file for the target protein(s): " + args.ZFP_FA )
    print("Input BED file for the peaks: " + args.TARGET_BED )


    ###### Check if input files exists
    if os.path.isfile( args.ZFP_FA ):
        print("Protein sequence file found.")
    else:
        print("ERROR: Protein sequence file was not found.")
        exit(2)



    if (args.TARGET_BED == "default_none" and args.TARGET_FASTA == "default_none"):
        print("Neither target bed or fasta sequence file specified (--target_bed or --target_fasta).")
        exit(2)

    if (args.TARGET_BED != "default_none" and args.TARGET_FASTA != "default_none"):
        print("Both target bed and fasta sequence file were specified, only use one (--target_bed or --target_bed).")
        exit(2)


    ####################### define temporary path                                                                                                     
    tmp_folder = args.OUT_DIR + "/tmp"
    if not os.path.exists(tmp_folder): os.makedirs(tmp_folder)


    tmp_RF_in = tmp_folder + "/" + args.JOB_ID + "_predict.in"
    tmp_RF_out = tmp_folder + "/" + args.JOB_ID + "_predict.RF.out"

    centered = tmp_folder + "/" + args.JOB_ID + "_centered.fa"
    shuffled = tmp_folder + "/" + args.JOB_ID + "_shuffled.fa"
    all = tmp_folder + "/" + args.JOB_ID + "_all.fa"

    ####################### define the output path

    input_bed = args.OUT_PREFIX + "input_coordinates.bed"

    log_step1 = args.OUT_PREFIX + "log_step1.txt"
    log_step2 = args.OUT_PREFIX + "log_step2.txt"
    report = args.OUT_PREFIX + "report.txt"
    RF_out = args.OUT_PREFIX + "RF_out.txt"
    PFM_scores = args.OUT_PREFIX + "PFM_scores.txt"

    log_info = args.OUT_PREFIX + "log.info.txt"
    log_error = args.OUT_PREFIX + "log.error.txt"

    lbrack = "{"; rbrack = "}"


    if os.path.isfile( args.TARGET_BED ):

        print("Peak sequence file found.")


        ####################### prepare the peak sequences                                                                                                
        # get the central region of the sequences, and also dinucleotide-shuffled sequences
        if args.BED_HAS_SCORE:
            cmdline=f"""sort -k 1,1 -k2,2n {args.TARGET_BED} |
            awk -v FS="\\t" -v OFS="\\t" '{lbrack} print $1, int($2+($3-$2)/2), int($2+($3-$2)/2),$4 {rbrack}' - |
            bedtools slop -b {args.RC_RANGE} -g {args.CHR_SIZES} -i - |
            awk -v FS="\\t" -v OFS="\\t" '{lbrack} print $1,$2,$3,$1":"$2"-"$3,$4 {rbrack}' - > {input_bed}"""

        else:
            cmdline=f"""sort -k 1,1 -k2,2n {args.TARGET_BED} |
            awk -v FS="\\t" -v OFS="\\t" '{lbrack} print $1, int($2+($3-$2)/2), int($2+($3-$2)/2) {rbrack}' - |
            bedtools slop -b {args.RC_RANGE} -g {args.CHR_SIZES} -i - |
            awk -v FS="\\t" -v OFS="\\t" '{lbrack} print $1,$2,$3,$1":"$2"-"$3 {rbrack}' - > {input_bed}"""

        utils.run_cmd(cmdline, log)


        cmdline=f"""bedtools getfasta -fi {args.GENOME_FA} -bed {input_bed} |
        sed -e 's/c/C/g' - | sed -e 's/g/G/g' - | sed -e 's/a/A/g' - |
        sed -e 's/t/T/g' - | sed -e 's/n/N/g' - | sed -e 's/Chr/chr/g' - |
        sed -e 's/CHR/chr/g' - > {centered} """
        utils.run_cmd(cmdline, log)

    elif os.path.isfile( args.TARGET_FASTA ):

        args.RC_RANGE = int(args.RC_RANGE) * 2

        # get the central region of the sequences
        cmdline=f"""{args.MEME_DIR}/fasta-center -protein -len {args.RC_RANGE} < {args.TARGET_FASTA} 1> {centered} """
        utils.run_cmd(cmdline, log)

    else:
        print("ERROR: Neither peak bed or fasta sequence files were found.")
        exit(2)



    cmdline=f"""{args.MEME_DIR}/fasta-dinucleotide-shuffle -f {centered} -t -dinuc 1> {shuffled} """
    utils.run_cmd(cmdline, log)

    cmdline=f"""cat {centered} {shuffled} > {all}"""
    utils.run_cmd(cmdline, log)


    if os.path.exists(log_step1): os.remove(log_step1) 
    if os.path.exists(log_step2): os.remove(log_step2) 


    for i in [3, 4, 5, 6, 7, 8]:

        # $FASTAtoRF -minl 2 -maxl 8 -span $i -fasta $proteins -out $RF_in.span$i >>$out_folder/log.step1.txt
        cmdline = f"""{args.script_path}/bin/FASTAtoRF -minl 2 -maxl 8 -span {i} -fasta {args.ZFP_FA} -out {tmp_RF_in}.span{i} >>{log_step1}"""
        utils.run_cmd(cmdline, log)

        if i == 3:
            cmdline = f"""cat {tmp_RF_in}.span{i} > {tmp_RF_in}"""
            utils.run_cmd(cmdline, log)

        else:
            cmdline = f"""cat {tmp_RF_in}.span{i} | sed 1d >> {tmp_RF_in} """
            utils.run_cmd(cmdline, log)


    ####################### run the RF script, and reformat it for the next step
    cmdline = f"""Rscript {args.script_path}/src/_R/_predict.RF.R --src_dir {args.script_path}/src/_R --predict_in {tmp_RF_in} --predict_out {tmp_RF_out}"""
    utils.run_cmd(cmdline, log)

    cmdline = f"""sed 's/"//g' {tmp_RF_out} > {RF_out}""" 
    utils.run_cmd(cmdline, log)

    ####################### run the RCADEEM script

    cmdline = f"""{args.script_path}/bin/RCADEEM -rf {RF_out} -fasta {all} -out {args.OUT_PREFIX} -mode 3  >>{log_step2} """
    utils.run_cmd(cmdline, log)

    cmdline = f"""Rscript {args.script_path}/src/_R/_process_RCADEEM_results.R --prefix {args.OUT_PREFIX} --results {PFM_scores}"""
    utils.run_cmd(cmdline, log)


    #*****************************************************************************************                                                                                                    
    # The following lines check the input/output, and produce appropriate messages                                                                                                                
    # If no error was detected in either input or output, the info messages will be written in                                                                                                    
    # ./out/<jobID>/log.info.txt                                                                                                                                                                  
    # Otherwise, the error messages will be written in                                                                                                                                            
    # ./out/<jobID>/log.error.txt                                                                                                                                                                 
    #*****************************************************************************************  
    info = ""

    with open(log_step1) as f:
        log_step1_f = f.readlines()

    with open(log_step2) as f:
        log_step2_f = f.readlines()

    with open(report) as f:
        reader = csv.reader(f, delimiter="\t")
        report_f = list(reader)





    ### check if the C2H2-ZF file had any valid sequences
    for row in log_step1_f:

        if "sequences were read." in row:

            read_seq = row.split(" ")[0]

            if str( read_seq ) == "ERROR:" or str( read_seq ) == "":
                msg = "ERROR: No sequences were found in the input FASTA for C2H2-ZF proteins. Please check the input format.\n"
                end_with_err(log_error, msg)

            elif int( read_seq ) <= 0:
                msg = "ERROR: No sequences were found in the input FASTA for C2H2-ZF proteins. Please check the input format.\n"
                end_with_err(log_error, msg)

            else:
                msg = str(read_seq) + " sequences were found in the input FASTA for C2H2-ZF proteins.\n"
                print(msg); info += msg
                break





    ### check if the peak sequence file had any valid sequences
    for row in log_step2_f:

        if "sequences were read," in row:

            read_seq = row.split(" ")[0]

            if str( read_seq ) == "ERROR:" or str( read_seq ) == "":
                msg = "ERROR: No sequences were found in the input FASTA for ChIP-seq peaks. Please check the input format.\n"
                end_with_err(log_error, msg)


            elif int( read_seq ) % 2 == 0:
                msg = str( read_seq ) + " sequences were found in the input FASTA for ChIP-seq peaks.\n"
                print(msg); info += msg
                break



    ### check if the C2H2-ZF sequences have had any ZF arrays
    for row in log_step2_f:

        if "motifs were read." in row:
            
            motif_read = row.split(" ")[0]

            if str( motif_read ) == "ERROR:" or str( motif_read ) == "":
                msg = "ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
                end_with_err(log_error, msg)
            
            elif int(motif_read) <= 0:
                msg = "ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
                end_with_err(log_error, msg)

            else:
                msg = str( motif_read ) + " possible C2H2-ZF arrays were tested.\n"
                print(msg); info += msg
                break



    for row in log_step2_f:
        if "ERROR: No motifs are enriched." in row:
            msg = "ERROR: No motifs are enriched.\n"
            end_with_err(log_error, msg)

        if "ERROR: The file has no EOL character." in row:
            msg = "ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
            end_with_err(log_error, msg)




    ### Checks for optimized motifs
    num_opt_motifs = sum( [ int(row[2]) for row in report_f[1:]] )

    if num_opt_motifs >= 1:
        msg = "The predicted motifs of " + str(num_opt_motifs) + " possible C2H2-ZF arrays were enriched in the ChIP-seq peaks.\n"
        print( msg ); info += msg
    else:
        msg = "ERROR: None of the predicted motifs from the provided C2H2-ZF proteins were enriched in the ChIP-seq peaks.\n"
        end_with_err(log_error, msg)




    with open(log_info, "w") as f:
        f.write(info)



