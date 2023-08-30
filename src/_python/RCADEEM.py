import csv   
import os
from src._python import utils


def RCADEEM( args, log ):
    
    print("Job ID: " + args.JOB_ID )
    print("Input FASTA file for the target protein(s): " + args.ZFP_FA )
    print("Input FASTA file for the peaks: " + args.CHIP_FA )


    ###### Check if input files exists
    if os.path.isfile( args.ZFP_FA ):
        print("Protein sequence file found.")
    else:
        print("ERROR: Protein sequence file was not found.")
        exit(2)

    if os.path.isfile( args.CHIP_FA ):
        print("Peak sequence file found.")
    else:
        print("ERROR: Peak sequence file not found.")
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

    log_step1 = args.OUT_PREFIX + "log_step1.txt"
    log_step2 = args.OUT_PREFIX + "log_step2.txt"
    report = args.OUT_PREFIX + "report.txt"
    RF_out = args.OUT_PREFIX + "RF_out.txt"
    PFM_scores = args.OUT_PREFIX + "PFM_scores.txt"

    log_info = args.OUT_PREFIX + "log.info.txt"
    log_error = args.OUT_PREFIX + "log.error.txt"

    
    ####################### prepare the peak sequences                                                                                                
    # get the central region of the sequences, and also dinucleotide-shuffled sequences                                                               
    cmdline='%s/fasta-center -len 100 < %s 1> %s ' % \
        ( args.MEME_DIR, args.CHIP_FA, centered )
    utils.run_cmd(cmdline, log)
    
    cmdline='%s/fasta-dinucleotide-shuffle -f %s -t -dinuc 1> %s ' % \
        ( args.MEME_DIR, centered, shuffled )
    utils.run_cmd(cmdline, log)

    cmdline='cat %s %s > %s' % (centered, shuffled, all)
    utils.run_cmd(cmdline, log)


    if os.path.exists(log_step1): os.remove(log_step1) 
    if os.path.exists(log_step2): os.remove(log_step2) 


    for i in [3, 4, 5, 6, 7, 8]:

        # $FASTAtoRF -minl 2 -maxl 8 -span $i -fasta $proteins -out $RF_in.span$i >>$out_folder/log.step1.txt
        
        cmdline = '%s/bin/FASTAtoRF -minl 2 -maxl 8 -span %s -fasta %s -out %s.span%s >>%s' % \
            ( args.script_path, i, args.ZFP_FA, tmp_RF_in, i, log_step1 )
        utils.run_cmd(cmdline, log)

        if i == 3:
            cmdline = 'cat %s.span%s > %s' % \
                (tmp_RF_in, i, tmp_RF_in)
            utils.run_cmd(cmdline, log)

        else:
            cmdline = 'cat %s.span%s | sed 1d >> %s' % \
                (tmp_RF_in, i, tmp_RF_in)
            utils.run_cmd(cmdline, log)


    ####################### run the RF script, and reformat it for the next step
    cmdline = 'Rscript %s/src/_R/_predict.RF.R --src_dir %s --predict_in %s --predict_out %s' % \
        (args.script_path, args.script_path + "/src/_R", tmp_RF_in, tmp_RF_out)
    utils.run_cmd(cmdline, log)

    cmdline = """sed 's/"//g' %s > %s""" % \
        (tmp_RF_out, RF_out)
    utils.run_cmd(cmdline, log)

    ####################### run the RCADEEM script
    cmdline = '%s/bin/RCADEEM -rf %s -fasta %s -out %s -mode 3  >>%s ' % \
        (args.script_path, RF_out, all, args.OUT_PREFIX, log_step2)
    utils.run_cmd(cmdline, log)

    cmdline = 'Rscript %s/src/_R/_process_RCADEEM_results.R --prefix %s --results %s  ' % \
        (args.script_path,  args.OUT_PREFIX, PFM_scores)
    utils.run_cmd(cmdline, log)


    #*****************************************************************************************                                                                                                    
    # The following lines check the input/output, and produce appropriate messages                                                                                                                
    # If no error was detected in either input or output, the info messages will be written in                                                                                                    
    # ./out/<jobID>/log.info.txt                                                                                                                                                                  
    # Otherwise, the error messages will be written in                                                                                                                                            
    # ./out/<jobID>/log.error.txt                                                                                                                                                                 
    #*****************************************************************************************  
    info = ""
    err = ""

    ### Checks for optimized motifs
    with open(report) as f:
        reader = csv.reader(f, delimiter="\t")
        report_d = list(reader)

    num_opt_motifs = sum( [ int(row[2]) for row in report_d[1:]] )

    if num_opt_motifs >= 1:
        msg = "The predicted motifs of " + str(num_opt_motifs) + " possible C2H2-ZF arrays were enriched in the ChIP-seq peaks.\n"
        print( msg ); info += msg
    else:
        msg = "ERROR: None of the predicted motifs from the provided C2H2-ZF proteins were enriched in the ChIP-seq peaks.\n"
        print(msg); err += msg


    with open(log_step2) as f:
        log_step2_d = f.readlines()


    ### check if the peak sequence file had any valid sequences
    for row in log_step2_d:

        if "sequences were read," in row:

            read_seq = row.split(" ")[0]

            if str( read_seq ) == "ERROR:" or str( read_seq ) == "":
                msg = "ERROR: No sequences were found in the input FASTA for ChIP-seq peaks. Please check the input format.\n"
                print(msg); err += msg

            elif int( read_seq ) % 2 == 0:
                msg = str( read_seq ) + " sequences were found in the input FASTA for ChIP-seq peaks.\n"
                print(msg); info += msg
                break


    ### check if the C2H2-ZF sequences have had any ZF arrays
    for row in log_step2_d:

        if "motifs were read." in row:
            
            motif_read = row.split(" ")[0]

            if str( motif_read ) == "ERROR:" or str( motif_read ) == "":
                msg = "ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
                print(msg); err += msg
            
            # elif int(motif_read) <= 0:
            #     msg = "ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
            #     print(msg); err += msg

            else:
                msg = str( motif_read ) + " possible C2H2-ZF arrays were tested.\n"
                print(msg); info += msg
                break


    ### check if the C2H2-ZF file had any valid sequences
    with open(log_step1) as f:
        log_step1_d = f.readlines()

    for row in log_step1_d:

        if "sequences were read." in row:

            read_seq = row.split(" ")[0]

            if str( read_seq ) == "ERROR:" or str( read_seq ) == "":
                msg = "ERROR: No sequences were found in the input FASTA for C2H2-ZF proteins. Please check the input format.\n"
                print(msg); err += msg
            
            elif int( read_seq ) <= 0:
                print(msg); err += msg
            else:
                msg = str(read_seq) + " sequences were found in the input FASTA for C2H2-ZF proteins.\n"
                print(msg); info += msg
                break


    if err == "":
        with open(log_info, "w") as f:
            f.write(info)
    else:
        with open(log_error, "w") as f:
            f.write(err)





