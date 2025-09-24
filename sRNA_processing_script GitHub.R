### Background + AIMS#######



# Interface script for sRNA sample processing: from collapse through to noise removal and normalisation.

# HOW TO RUN####

# VAST_ext = 352

#protein_anti_transcriptome
#/media/ruth/be73bf61-a0a3-4c3d-aeca-b28cbd731c8e/Ruth_VAST_pilot_and_extention/mir_bams/miRNA_minus_failed_run_allignment_files/protein_transcriptome_plus_whole_transcriptome_reverse_complement_allignments


# /media/ruth/be73bf61-a0a3-4c3d-aeca-b28cbd731c8e/Ruth_VAST_pilot_and_extention/mir_bams/miRNA_minus_failed_run_allignment_files/sRNAome_allignments

# VAST_addon = 176

#/media/ruth/be73bf61-a0a3-4c3d-aeca-b28cbd731c8e/RUTH_VAST_add_on_control_and_12_hour_samples/all_samples_trimmed_fastq_and_allignment_files/sRNAome_allignments


# /media/ruth/be73bf61-a0a3-4c3d-aeca-b28cbd731c8e/RUTH_VAST_add_on_control_and_12_hour_samples/all_samples_trimmed_fastq_and_allignment_files/protein_transcriptome_plus_whole_transcriptome_reverse_complement_allignments


# COVID = 58
#/media/ruth/858a064f-56aa-4a19-b746-0de06e999eac/COV2/210317_NB502094_0313_AHKV2LBGXH/miRNA/sRNAome_allignments

#/media/ruth/858a064f-56aa-4a19-b746-0de06e999eac/COV2/210317_NB502094_0313_AHKV2LBGXH/miRNA/protein_transcriptome_plus_whole_transcriptome_reverse_complement_allignments

# Euclids = 44

# COVID extension = 244

# Press source and follow instructions appearing in the console.

wd<-readline(prompt="what working directory would you like to use: ")

#/media/ruth/858a064f-56aa-4a19-b746-0de06e999eac/COV2/210317_NB502094_0313_AHKV2LBGXH/miRNA/sRNAome_allignments

# /media/ruth/858a064f-56aa-4a19-b746-0de06e999eac/COVID_extension/merged_files_for_allignment/sRNAome_allignments

# /media/ruth/858a064f-56aa-4a19-b746-0de06e999eac/COVID_extension/merged_files_for_allignment/protein_transcriptome_plus_whole_transcriptome_reverse_complement_allignments

cat("miRNA 
    tRNA_aminoacid (e.g. Homo_sapiens_tRNA-Ala)
    tRNA_anticodon (e.g. Homo_sapiens_tRNA-Ala-AGC)
    tRNA_anticodon_isotype (e.g. Homo_sapiens_tRNA-Ala-AGC-1-1)
    snoRNA_full_description \n
    lncRNA \n  
    siRNA \
    siloco_siRNA
    viral_miR
    sRNAome
    protein_anti_transcriptome")

sRNA_level.of.collapse<-readline("\n what species of sRNAs are you collapsing? options: \n
                                 miRNA \n 
                                 tRNA_aminoacid (e.g. Homo_sapiens_tRNA-Ala) \n
                                 tRNA_anticodon (e.g. Homo_sapiens_tRNA-Ala-AGC) \n
                                 tRNA_anticodon_isotype (e.g. Homo_sapiens_tRNA-Ala-AGC-1-1) \n
                                 snoRNA_full_description \n
                                 lncRNA \n  
                                 siRNA \
                                 siloco_siRNA
                                 viral_miR
                                 sRNAome
                                 protein_anti_transcriptome")
  
if(sRNA_level.of.collapse=="siRNA"|sRNA_level.of.collapse=="sRNAome"){
  remove_seqs_that_map_this_many_times<-readline("\n Up to how many times can an allignment allign before being excluded? e.g. 1 (can only map once),
                                                 40 (can map up to 40 times)
                                                no limit (can map an infinite number of times) etc \n\n")
}

prefix<-readline("what file prefix would you like to use? e.g. Euclids_4month \n 
                 Euclids_12months \n
                 VAST_ext \n
                 VAST_addon\n
                 COVID ")

n.samps<-readline("how many samples do you have? ")

detection_cutoff<-readline("what proportion of samples should the gene be detected in after noise removal?")
#for VAST ext = 0.04 (as there are up to 24 seperate groups)

collapse_sRNAs<-readline("would you like to collapse your sRNAs? y/n ")

expression.matracies<-readline("would you like to create expression matracies for noise analysis? y/n ")

p2p_calculation<-readline("would you like to perform the P2P noise analysis? y/n ")

noise_threshold_estimates<-readline("would you like to calculate the samplewise noise thresholds? y/n ")

draw_presence_plots<-readline("would you like to draw the presence plots? y/n ")

remove_noise_and_normalise<-readline("would you like to remove the noise from your samples and normalise? y/n ")


### SOURCE WRAPPER SCRIPT ####

location.of.script<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/") # get the file paths to the scripts that will be used



source(paste0(location.of.script,"wrapper_script.R"))



##### If you prematurely terminate the perl batch file above by pressing the red STOP button, run this code section below to kill the perl process and liberate the abn and pcc files for deletion.

if(computer=="D:/"){
  kill.perl<-readline("\n\n WARNING!! \n Are you sure you want to kill all perl processes? y/n \n\n\n")

  if(kill.perl=="y"){
    system("Taskkil /F /IM perl.exe")
  }else{
    print("Task cancelled, Perl process NOT killed")
  }
}