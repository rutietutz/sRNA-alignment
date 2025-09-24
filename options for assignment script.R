option_choices<-c("miRNA",
                  "tRNA_aminoacid",
                  "tRNA_anticodon",
                  "tRNA_anticodon_isotype",
                  "snoRNA_full_description",
                  "lncRNA",
                  "siRNA",
                  "siloco_siRNA",
                  "yRNA",
                  "viral_miR",
                  "sRNAome",
                  "protein_anti_transcriptome"
                  )

if(!sRNA_level.of.collapse%in%option_choices){
  print("ERROR: You have failed to supply a valid response to the question: what species of sRNAs are you collapsing?")
  stop()
}

sRNA_level.of.collapse_for_noise<-sRNA_level.of.collapse



# miRNA collapse variables

if(sRNA_level.of.collapse=="miRNA"){
  start<-41
  stop<-26

  sRNA_level.of.collapse_for_noise<-"miRNA.hp"

  start2<-51
  stop2<-36


}


if(sRNA_level.of.collapse=="miRNA"&&grepl("2MM",wd)){
  start<-31
  stop<-16
  sRNA_level.of.collapse_for_noise<-"miRNA.hp"

  start2<-32
  stop2<-15

}


if(length(grep("tRNA",sRNA_level.of.collapse))>0){
  start<-53
  stop<-38

  sRNA_level.of.collapse_for_noise<-sRNA_level.of.collapse

  start2<-start
  stop2<-stop}

if(length(grep("snoRNA",sRNA_level.of.collapse))>0){
  start<-55
  stop<-40

  sRNA_level.of.collapse_for_noise<-sRNA_level.of.collapse

  start2<-start
  stop2<-stop}


if(length(grep("lncRNA",sRNA_level.of.collapse))>0){
  start<-55
  stop<-40

  sRNA_level.of.collapse_for_noise<-sRNA_level.of.collapse

  start2<-start
  stop2<-stop}

if(length(grep("siRNA",sRNA_level.of.collapse))>0){
  start<-48
  stop<-33

  sRNA_level.of.collapse_for_noise<-sRNA_level.of.collapse

  start2<-start
  stop2<-stop}

if(length(grep("siloco_siRNA",sRNA_level.of.collapse))>0){
  start<-48
  stop<-33

  sRNA_level.of.collapse_for_noise<-sRNA_level.of.collapse

  start2<-start
  stop2<-stop}

if(length(grep("yRNA",sRNA_level.of.collapse))>0){
  start<-53
  stop<-38

  sRNA_level.of.collapse_for_noise<-sRNA_level.of.collapse

  start2<-start
  stop2<-stop}

if(length(grep("viral_miR",sRNA_level.of.collapse))>0){
  start<-57
  stop<-42

  sRNA_level.of.collapse_for_noise<-sRNA_level.of.collapse

  start2<-start
  stop2<-stop}
#

# development
# run4<-fread(files[i],sep=" ",header=F,select=c(1:3))
# start<-31
# stop<-16
# stri_sub(files[i],-start,-stop)






#### cutting isotypes ####
needs_cutting<-F

if(sRNA_level.of.collapse=="tRNA_aminoacid"){
  needs_cutting<-T
  cut=2
}

if(sRNA_level.of.collapse=="tRNA_anticodon"){
  needs_cutting<-T
  cut=3
}


