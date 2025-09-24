### Background + AIMS#######

#Background:
# Multimapping reads need assignment. Using merged .sum files.

#AIM: 
# Take sum files and assign multimapping reads to the most abundant sRNA (pre final assignment). This script assumes multimapping reads have been pre-assigned using the best-straata method.



## Install/load packages as required #####

if (!require("pacman")) install.packages("pacman"); library(pacman)


p_load("readr", 
       "data.table",
       "stringi", 
       "reshape2",
       "plyr",
       "dplyr",
       "parallel",
       "stringr",
       "doParallel")



## take 7.2 mins for 44 samples with 140,000 unique sequences

#source("https://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")

# FUNCTIONS ####

keep.cols.function<-function(dt,keep.cols){
  keep.cols<-colnames(dt)[keep.cols]
  keep.cols<-colnames(dt)%in%keep.cols
  run2<-subset(dt,,keep.cols) # these are the columns corresponding to tRNA","seq","plus_minus","MM"
  return(run2)
}

remove.cols.function<-function(dt,remove.cols){
  remove.cols<-colnames(dt)[remove.cols]
  dt[,(remov.cols) :=NULL] # these are the columns corresponding to tRNA","seq","plus_minus","MM"
  return(run2)
}

# Get drive/root prefix
  
# if(length(grep("home",getwd()))==1){computer<-"~/"}else{computer<-"D:/"}
# 
# 
# location.of.script<-dirname(rstudioapi::getSourceEditorContext()$path)

#--------------------------------END OF DEPENDANCIES ------


# ------------------------ User defined instructions ------


## 
# wd<-readline(prompt="what working directory would you like to use: ")
# 
# sRNA_level.of.collapse<-readline("\n what species of sRNAs are you collapsing? options: \n
#                                  miRNA")
# 
# prefix<-readline("what file prefix would you like to use? e.g. Euclids_4month \n \n")
# 
# expression.matracies<-readline("would you like to create expression matracies for noise analysis? y/n")
# 
# p2p_calculation<-readline("would you like to perform the P2P noise analysis? y/n")
# 
# if(p2p_calculation=="y"){
# n.samps<-readline("how many samples do you have?")
# }
# 
# noise_threshold_estimates<-readline("would you like to calculate the samplewise noise thresholds? y/n")
# 
# draw_presence_plots<-readline("would you like to draw the presence plots? y/n")
# 
# remove_noise_and_normalise<-readline("would you like to remove the noise from your samples and normalise? y/n")
# 

#### END of user defined variables


source(paste0(location.of.script,"options for assignment script.R"))  ## file that will assign correct variables for each sRNA type, e.g. which bit of a file name needs cutting out

#setwd(wd)
output<-paste0(prefix,"_noise_included_not_norm.csv")  # name of output of assigned reads



# NOTES

# It looks like snoRNAs may be collapsable but for now I will leave at the lowest level and then aggregate at a later date.


# ------------------------ SCRIPT PROPER:  ------------


files<-list.files(path=location.of.patman.files, pattern=".sum",recursive=F)

files<-files[grep("merged.*.sum",files,invert = F)] # Only use the merged files.

files<-paste0(location.of.patman.files,files)



# create a df of all the sequences


i<-1

if(grepl("2MM|COVID",wd)){
run4<-fread(files[i],sep=" ",header=F, select=c(1:3))
colnames(run4)[1:3]<-c("count","SEQ","MIR")}else{run4<-fread(files[i],sep="\t",header=T,select=c(1:3))}

# run4[1,] # out of interest - look at top line of file
# 
# nrow(unique(run4[,"SEQ"])) # no. of unique sequences

none_detected<-data.frame(matrix(nrow=1,ncol=3))
colnames(none_detected)<-colnames(run4)

colnames(run4)[i]<-stri_sub(basename(files[i]),1,23)

if(needs_cutting==T){
  run4$MIR<-word(run4$MIR, 1, cut, sep="-")
  run4<-run4[!duplicated(run4),]
}


if(grepl("2MM|COVID",wd)){
  for(i in 2:length(files)){
    run3<-fread(files[i],sep=" ",header=F, select=c(1:3))
    colnames(run3)[1:3]<-c("count","SEQ","MIR")
    if(needs_cutting==T){
      run3$MIR<-word(run3$MIR, 1, cut, sep="-")
    }
    run3<-run3[!duplicated(run3),]
    if(nrow(run3)==0){
      run3<-rbind(run3,none_detected)
    }
    run4<-merge(x = run4, y = run3, by = c("SEQ","MIR"), all = TRUE)
    colnames(run4)[i+2]<-stri_sub(basename(files[i]),1,23)
    print(i)
  }
}else{
  for(i in 2:length(files)){
    run3<-fread(files[i],sep="\t",select=c(1:3))
    if(needs_cutting==T){
      run3$MIR<-word(run3$MIR, 1, cut, sep="-")
    }
    run3<-run3[!duplicated(run3),]
    if(nrow(run3)==0){
      run3<-rbind(run3,none_detected)
    }
    run4<-merge(x = run4, y = run3, by = c("SEQ","MIR"), all = TRUE)
    colnames(run4)[i+2]<-stri_sub(basename(files[i]),1,23)
    print(i)
  }
}
    

test2<-run4 # a backup
test2[is.na(test2)] <- 0

colnames(test2)[c(1,2)]<-c("seq","sRNA")

## filtering step.
table<-table(test2[,seq])
single.mapping<-table[table==1]

Ok.sRNA<-test2[seq%in%names(single.mapping)]

multimapping.sRNA<-table[table>1]

multi.sRNA<-test2[seq%in%names(multimapping.sRNA)]

unique.multi<-unique(multi.sRNA[,seq], by=c('seq'))# NR reads


#####  Create totals based on all reads -- I will use this to decide re assignment ######


keep.cols<-c(3:ncol(test2)) # columns containing count data for each sample
totals<-keep.cols.function(test2,keep.cols) # 

totals<-rowSums(test2[,3:ncol(test2)])

totals<-cbind(test2[,"sRNA"],totals)

totals2<-aggregate(totals ~ sRNA,totals, sum)

setorder(totals2, -totals) # order totals by size

totals2<-data.table(totals2)

#######  Assign the sRNA based on the sRNA abudances ######


reduced.final<-multi.sRNA ## for down stream compatibility with old script
setkey(multi.sRNA,seq)

NR.sequences<-unique(reduced.final$seq)
  
reduced.final=as.data.frame(reduced.final)

ptm<-proc.time() # start the clock 

num.cores<-detectCores()
cl <- makeCluster(num.cores-1)
registerDoParallel(cl)

cat("\n collapsing multimapping sequences - may take some time... \n")

multi.reduced<-foreach (i=1:length(NR.sequences), .combine = rbind ) %dopar% {
  seq<-reduced.final[(reduced.final$seq==NR.sequences[i]),] # cut out the sequence
  lookup<-seq[,2] # the sRNAs that the sequence maps to
  amino.acid.to.keep<-totals2$sRNA[totals2$sRNA%in%lookup][1] # lookup the amino acid positions in totals2 and take the first entry (totals2 has been sorted into reducing library size so the 1st entry I come across has the largest library size)
  keep.record<-seq[seq$sRNA==amino.acid.to.keep,]
}

stopCluster(cl)

proc.time()-ptm # time estimate

cat("\n multireduced sequences have been collapsed. Now agreggating sRNAs and writing the file \n")
collapsed.sRNA<-rbind(Ok.sRNA,multi.reduced) # combine the single mapping sequences and the now reduced multimapping sequences

frag.length<-apply(collapsed.sRNA[,1],1,stri_length) # length of each read

collapsed.sRNA.seq<-cbind(collapsed.sRNA,frag.length)

fwrite(collapsed.sRNA.seq, paste0(sRNA_level.of.collapse,"_multireduced.rds")) # backup file


# Now sum up the counts for each sRNA
collapsed.sRNA[ ,c("seq") := NULL] # remove the seq column


collapsed.sRNA.noise.incl<-collapsed.sRNA[, lapply(.SD, sum, na.rm=TRUE), by=sRNA, .SDcols=c(2:ncol(collapsed.sRNA)) ] 


file.it<-paste0(location.of.patman.files,sRNA_level.of.collapse,"_collpased_",output)  #  will need this file name later in the script

fwrite(collapsed.sRNA.noise.incl,file=file.it) # out put the collapsed file

rm(list=c("test2", "run4","table","single.mapping","Ok.sRNA","multimapping.sRNA","multi.sRNA","unique.multi","totals","totals2","reduced.final","NR.sequences","multi.reduced","collapsed.sRNA","frag.length","collapsed.sRNA.seq","collapsed.sRNA.noise.incl"))  ## remove the large objects created in this script in case I am piping through a wrapper script.
