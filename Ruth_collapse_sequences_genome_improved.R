# This program will collapse siRNA sequences in two ways: A) it will collapse multimapping sequences to the most abundant feature. B) it will collapse multimapping sequences to a union gene.


## To run bespoke uncomment this section ######################

# location.of.patman.files<-"/media/ruth/be73bf61-a0a3-4c3d-aeca-b28cbd731c8e/RUTH_VAST_add_on_control_and_12_hour_samples/all_samples_trimmed_fastq_and_allignment_files/fastq_files_for_remapping/0MM_genome_allignments/"

# if (!require("pacman")) install.packages("pacman"); library(pacman)
# 
# 
# p_load("readr", 
#        "data.table",
#        "stringi",
#        "stringr",
#        "reshape2",
#        "plyr",
#        "dplyr",
#        "parallel",
#        "stringr",
#        "doParallel",
#        "mail",
#        "BiocManager")


########################################################

###  Create totals based on all reads -- I will use this to decide re assignment ####

cat("Create totals based on all reads -- I will use this to decide re assignment \n")
expt.sum.file<-list.files(location.of.patman.files,pattern = "experimentlevel.sum",full.names = T)
totals<-fread(expt.sum.file)

totals<-totals%>%mutate(FLOOR=round_any(V4,100,f=floor))


#totals<-totals%>%mutate(locus=paste(V2,FLOOR,V4,sep="_"))
totals$locus=paste(totals$V3,totals$FLOOR,totals$V5,sep="_")

totals<-totals[,c(1,7)]

colnames(totals)[c(1,2)]<-c("totals","sRNA")

totals2<-data.table(totals)

totals2<-aggregate(totals ~ sRNA,totals, sum)

setorder(totals2, -totals) # order totals by size



# ################
# 
## prep for filtering step
min.read.freq=1
min.sample=4

must.have.at.least<-(min.read.freq)*(min.sample+1)

totals3<-totals2[totals2$totals>=must.have.at.least,]

saveRDS(totals3,"totals3.RDS")


files<-list.files(path=location.of.patman.files, pattern=".sum",recursive=F)

files<-files[grep("merged.*.sum",files,invert = F)] # Only use the merged files.

files<-paste0(location.of.patman.files,files)



#prep for parallelising
num.cores<-detectCores()

# up to 14 has completed


totals3<-readRDS(paste0(location.of.patman.files,"totals3.RDS"))

totals3$rank<-frank(totals3$totals,ties.method="dense")

for (i in 1:length(files)){ # replace back
  
  
  run4<-fread(files[i],sep=" ",header=F,fill=T) # change
  run4<-run4[!(run4$V4<0|run4$V5<0),]
  
  
  # run4[1,] # out of interest - look at top line of file
  # 
  
  
  
  colnames(run4)[1]<-stri_sub(basename(files[i]),1,23)
  colnames(run4)[2:5]<-c("SEQ","MIR","START","FLAG")
  
  run4 <- run4[, START:=as.numeric(START)]
  run4 <- run4[!is.na(START), ]
  # can ignore te warnings. The nas values correspond to header lines in the sum files which were accdidently kept in
  
  # run4[1,] # out of interest - look at top line of file
  # 
  # nrow(unique(run4[,"SEQ"])) # no. of unique sequences
  
  
  
  test2<-run4 # a backup
  test2[is.na(test2)] <- 0
  
  test2<-test2%>%mutate(FLOOR=round_any(START,100,f=floor))
  
  
  test2<-test2%>%mutate(paste(MIR,FLOOR,FLAG,sep="_"))
  
  colnames(test2)[c(1,ncol(test2))]<-c("count","sRNA")
  
  
  ## filtering step - remove loci which did not pass the filtering criteria in the totals2 dataframe: should make code a bit quicker to run and reduce memory required. 
  
  test3<-test2[test2$sRNA%in%totals3$sRNA,]
  
  
  ## set aside single mapping reads
  test3<-data.table(test3)
  table<-table(test3[,SEQ])
  single.mapping<-table[table==1]
  
  Ok.sRNA<-test3[SEQ%in%names(single.mapping)]
  Ok.sRNA$rank<-NA
  Ok.sRNA$no.mapping.sites<-1
  Ok.sRNA$sites_seq_maps_to<-Ok.sRNA$sRNA
  multimapping.sRNA<-table[table>1]
  
  multi.sRNA<-test3[SEQ%in%names(multimapping.sRNA)]
  
  unique.multi<-unique(multi.sRNA[,SEQ], by=c('SEQ'))# NR reads
  
  
  
  
  #######  Assign the sRNA based on the sRNA abudances ######
  
  
  reduced.final<-multi.sRNA ## for down stream compatibility with old script
  
  # order reduced.final by window abundance.
  
  reduced.final$rank<-totals3$rank[match(reduced.final$sRNA,totals3$sRNA)]
  
  reduced.final<-reduced.final[order(reduced.final$rank,decreasing=T),]
  
  # remove sequences mapping > specifiie no of times times
  
  table<-data.frame(table(reduced.final$SEQ))
  
if(sum(grepl("no limit",remove_seqs_that_map_this_many_times)==0))
  {
  NR.sequences<-table$Var1[table$Freq<=as.numeric(as.character(remove_seqs_that_map_this_many_times))]
  reduced.final<-reduced.final[reduced.final$SEQ%in%NR.sequences,]
  }
  
  reduced.final$no.mapping.sites<-table$Freq[match(reduced.final$SEQ,table$Var1)]
  setkey(multi.sRNA,SEQ)
  
  reduced.final=as.data.frame(reduced.final)
  multi.reduced=reduced.final[!duplicated(reduced.final$SEQ),]
  
  cat("\n collapsing multimapping sequences - may take some time... \n")
  
  
  
  cl <- makeCluster(num.cores-1)
  registerDoParallel(cl)
  multi.reduced_final_col<-foreach (j=1:length(multi.reduced$SEQ), .combine = c ) %dopar% {
    keep.record<-paste(reduced.final$sRNA[reduced.final$SEQ%in%multi.reduced$SEQ[j]],collapse = "&") # cut out the sequence
  }
  
  stopCluster(cl)
  
  multi.reduced$sites_seq_maps_to<-multi.reduced_final_col
  
  
  collapsed.sRNA<-rbind(Ok.sRNA,multi.reduced) # combine the single mapping sequences and the now reduced multimapping sequences
  
  file.name<-paste0(location.of.patman.files,stri_sub(basename(files[i]),1,23),"_temp3.txt")
  write.table(collapsed.sRNA,file.name,col.names = T,row.names = F) # very quick and possibly will save memory 
  
  print(paste("entry",i,"has completed"))
}


cat("\n multireduced sequences have been collapsed. Now merging the files \n")

files2<-dir(path=location.of.patman.files, pattern="_temp3.txt",full.names = T,recursive=F) 


begin<-26
end<-11

i<-1
run4<-fread(files2[i],sep=" ") # change
colnames(run4)[1]<-stri_sub(files2[i],-begin,-end)

for(i in 2:length(files2)){
  run3<-fread(files2[i],sep=" ") # change
  colnames(run3)[1]<-stri_sub(files2[i],-begin,-end)
  run4<-merge(x = run4, y = run3, by = c("SEQ","MIR","START","FLAG","FLOOR","sRNA","rank","no.mapping.sites","sites_seq_maps_to"), all = TRUE) # change
  #  colnames(run4)[i+4]<-stri_sub(files[i],-start,-stop)
  print(i)
}

run4[is.na(run4)]<-0

collapsed.sRNA<-run4
frag.length<-apply(collapsed.sRNA[,1],1,stri_length) # length of each read

collapsed.sRNA.seq<-cbind(collapsed.sRNA,frag.length)

fwrite(collapsed.sRNA.seq, paste0("MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_multireduced.rds")) # backup file


cat("\n Now summing up the counts for each siRNA loci \n")

collapsed.sRNA<-fread(paste0("MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_multireduced.rds"),nThread=14)



# Now sum up the counts for each sRNA based on assignment to the most abundant read

collapsed.sRNA[ ,c("SEQ","START","FLAG","FLOOR","MIR","rank","no.mapping.sites","sites_seq_maps_to","frag.length") := NULL]


collapsed.sRNA_majority_collapsed.noise.incl<-collapsed.sRNA[, lapply(.SD, sum, na.rm=TRUE), by=sRNA, .SDcols=c(2:ncol(collapsed.sRNA)) ] 

rm(collapsed.sRNA) # saves memory

pass_filter<-rowSums(collapsed.sRNA_majority_collapsed.noise.incl[,2:ncol(collapsed.sRNA_majority_collapsed.noise.incl)]>min.read.freq)>min.sample

collapsed.sRNA_majority_collapsed.noise.incl<-collapsed.sRNA_majority_collapsed.noise.incl[pass_filter,] # remove


file.it.1<-paste0(location.of.patman.files,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_collpased_",output)  #  will need this file name later in the script

fwrite(collapsed.sRNA_majority_collapsed.noise.incl,file=file.it.1) # out put the collapsed file


rm(collapsed.sRNA_majority_collapsed.noise.incl)

# Now sum up the counts for each sRNA based on assignment to union genes

collapsed.sRNA<-fread(paste0("MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_multireduced.rds"),nThread=14)

collapsed.sRNA[ ,c("SEQ","START","FLAG","FLOOR","MIR","rank","no.mapping.sites","sRNA","frag.length") := NULL]

collapsed.sRNA_union_gene.noise.incl<-collapsed.sRNA[, lapply(.SD, sum, na.rm=TRUE), by=sites_seq_maps_to, .SDcols=c(2:ncol(collapsed.sRNA)) ] 


pass_filter<-rowSums(collapsed.sRNA_union_gene.noise.incl[,2:ncol(collapsed.sRNA_union_gene.noise.incl)]>min.read.freq)>min.sample

collapsed.sRNA_union_gene.noise.incl<-collapsed.sRNA_union_gene.noise.incl[pass_filter,] # remove

paste0(location.of.patman.files,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_collpased_",output) 

file.it.2<-paste0(location.of.patman.files,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_gene_union_collpased_",output)  #  will need this file name later in the script

fwrite(collapsed.sRNA_union_gene.noise.incl,file=file.it.2) # out put the collapsed file

rm(collapsed.sRNA_union_gene.noise.incl)
rm(collapsed.sRNA)
logger=paste("filter criteria applied @",Sys.time(),":","reads which did not have at least", min.read.freq ,"reads in", min.sample ,"samples were filtered out") # this provides a link so I can refer to what fitering settings I used to create the files.

sink(READMEfile, append = T) # will create the text file if it is not already in existance

writeLines(logger)

sink()




files.to.delete <- dir(location.of.patman.files,pattern="_temp3.txt",recursive=F,full.names=T)
file.remove(files.to.delete) # save space


rm(list=c("test2", "run4","table","single.mapping","Ok.sRNA","multimapping.sRNA","multi.sRNA","unique.multi","totals","totals2","reduced.final","NR.sequences","multi.reduced","frag.length"))  ## remove the large objects created in this script in case I am piping through a wrapper script.

