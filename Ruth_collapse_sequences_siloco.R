
#####  Create totals based on all reads -- I will use this to decide re assignment ######


totals<-read.table(paste0(location.of.patman.files,"siloco_experimentlevel.sum"))


colnames(totals)[c(1,2,3)]<-c("seq","sRNA","totals")

totals2<-aggregate(totals ~ sRNA,totals, sum)

setorder(totals2, -totals) # order totals by size

totals2<-data.table(totals2)

################

files<-list.files(path=location.of.patman.files, pattern=".sum",recursive=F)

files<-files[grep("merged.*.sum",files,invert = F)] # Only use the merged files.

files<-paste0(location.of.patman.files,files)

## prep for filtering step
min.read.freq=5
min.sample=2

must.have.at.least<-(min.read.freq)*(min.sample+1)

totals3<-totals2[totals2$totals>=must.have.at.least,]

#prep for parallelising
num.cores<-detectCores()

for (i in 1:length(files)){
  
  run4<-fread(files[i],sep=" ",header=T) # change
  
  colnames(run4)[1]<-stri_sub(basename(files[j]),1,23)
  colnames(run4)[2:5]<-c("SEQ","MIR","START","FLAG")
  
  dups<-paste0(run4$SEQ,run4$MIR)
  run4<-run4[!duplicated(dups),]
  
  test2<-run4 # a backup
  test2[is.na(test2)] <- 0
  
  colnames(test2)[c(1,3)]<-c("count","sRNA")
  
  
  ## filtering step - remove loci which did not pass the filtering criteria in the totals2 dataframe: should make code a bit quicker to run and reduce memory required. 
  
  test3<-test2[test2$sRNA%in%totals3$sRNA,]
  
  
  ## set aside single mapping reads
  test3<-data.table(test3)
  table<-table(test3[,SEQ])
  single.mapping<-table[table==1]
  
  Ok.sRNA<-test3[SEQ%in%names(single.mapping)]
  
  multimapping.sRNA<-table[table>1]
  
  multi.sRNA<-test3[SEQ%in%names(multimapping.sRNA)]
  
  unique.multi<-unique(multi.sRNA[,SEQ], by=c('SEQ'))# NR reads
  
  
  
  
  #######  Assign the sRNA based on the sRNA abudances ######
  
  
  reduced.final<-multi.sRNA ## for down stream compatibility with old script
  setkey(multi.sRNA,SEQ)
  
  NR.sequences<-unique(reduced.final$SEQ)
  
  reduced.final=as.data.frame(reduced.final)
  
  cat("\n collapsing multimapping sequences - may take some time... \n")
  
  cl <- makeCluster(num.cores-1)
  registerDoParallel(cl)
  
  multi.reduced<-foreach (j=1:length(NR.sequences), .combine = rbind ) %dopar% {
    SEQ<-reduced.final[(reduced.final$SEQ==NR.sequences[j]),] # cut out the sequence
    lookup<-SEQ$sRNA # the sRNAs that the sequence maps to
    amino.acid.to.keep<-totals3$sRNA[totals3$sRNA%in%lookup][1] # lookup the amino acid positions in totals3 and take the first entry (totals3 has been sorted into reducing library size so the 1st entry I come across has the largest library size)
    keep.record<-SEQ[SEQ$sRNA==amino.acid.to.keep,]
    
  }
  
  stopCluster(cl)
  
  
  collapsed.sRNA<-rbind(Ok.sRNA,multi.reduced) # combine the single mapping sequences and the now reduced multimapping sequences
  
  file.name<-paste0(stri_sub(basename(files[j]),1,19),"_temp3.txt")
  write.table(collapsed.sRNA,file.name,col.names = T,row.names = F) # very quick and possibly will save memory 
  
  print(i)
}


cat("\n multireduced sequences have been collapsed. Now merging the files \n")

files2<-dir(path=wd, pattern="_temp3.txt",full.names = T,recursive=F) 


begin<-26
end<-11

i<-1
run4<-fread(files2[i],sep=" ") # change
colnames(run4)[1]<-stri_sub(files2[i],-begin,-end)

for(i in 2:length(files2)){
  run3<-fread(files2[i],sep=" ") # change
  colnames(run3)[1]<-stri_sub(files2[i],-begin,-end)
  run4<-merge(x = run4, y = run3, by = c("SEQ","sRNA","START","FLAG"), all = TRUE) # change
  #  colnames(run4)[i+4]<-stri_sub(basename(files[j]),1,19)
  print(i)
}

run4[is.na(run4)]<-0

collapsed.sRNA<-run4
frag.length<-apply(collapsed.sRNA[,1],1,stri_length) # length of each read

collapsed.sRNA.seq<-cbind(collapsed.sRNA,frag.length)

fwrite(collapsed.sRNA.seq, paste0(sRNA_level.of.collapse,"_multireduced.rds")) # backup file

cat("\n Now summing up the counts for each siRNA loci \n")

# Now sum up the counts for each sRNA
collapsed.sRNA[ ,c("SEQ","START","FLAG") := NULL] # remove these columns
collapsed.sRNA<-data.table(collapsed.sRNA)

collapsed.sRNA.noise.incl<-collapsed.sRNA[, lapply(.SD, sum, na.rm=TRUE), by=sRNA, .SDcols=c(2:ncol(collapsed.sRNA)) ] 


pass_filter<-rowSums(collapsed.sRNA.noise.incl[,2:ncol(collapsed.sRNA.noise.incl)]>min.read.freq)>min.sample

collapsed.sRNA.noise.incl<-collapsed.sRNA.noise.incl[pass_filter,] # remove


logger=paste("filter criteria applied @",Sys.time(),":","reads which did not have at least", min.read.freq ,"reads in", min.sample ,"samples were filtered out") # this provides a link so I can refer to what fitering settings I used to create the files.

sink(READMEfile, append = T) # will create the text file if it is not already in existance

writeLines(logger)

sink()


file.it<-paste0(location.of.patman.files,sRNA_level.of.collapse,"_collpased_",output)  #  will need this file name later in the script

fwrite(collapsed.sRNA.noise.incl,file=file.it) # out put the collapsed file

#files.to.delete <- dir(location.of.patman.files,pattern="_temp3.txt",recursive=F,full.names=T)
#file.remove(files.to.delete)


#rm(list=c("test2", "run4","table","single.mapping","Ok.sRNA","multimapping.sRNA","multi.sRNA","unique.multi","totals","totals2","reduced.final","NR.sequences","multi.reduced","collapsed.sRNA","frag.length","collapsed.sRNA.seq","collapsed.sRNA.noise.incl"))  ## remove the large objects created in this script in case I am piping through a wrapper script.


