## may need to do some filtering first to remove any sRNAs with very low abundance -- currently 44220 seqs and only 2 files included e.g. < 5 in only 2 samples....

#Script to generate the expression profiles for each miRNA in one sample.


index.file<-paste(prefix,sRNA_level.of.collapse_for_noise, "index.for.noise.analysis.csv", sep="_")

#1. read in the file


ptm<-proc.time() # start the clock 

files<-list.files(path=location.of.patman.files,
                  pattern="*alligned.sum",
                  recursive=F)

files<-paste0(location.of.patman.files,files)



num.cores<-detectCores()
cl <- makeCluster(num.cores-1)
registerDoParallel(cl)

cat("\n Finding the max length for each sRNA")

max.length.1<-foreach(j =1:length(files), .packages=c('data.table','stringi','stringr','plyr'), .combine = rbind) %dopar% {

  run4<-fread(files[j],sep=" ",header=T)
  
  colnames(run4)[1]<-stri_sub(basename(files[j]),1,23)
  colnames(run4)[2:5]<-c("SEQ","sRNA","START","FLAG")
  
  dups<-paste0(run4$SEQ,run4$sRNA)
  run4<-run4[!duplicated(dups),]
  
  run4<-run4%>%mutate(END=(START+nchar(SEQ)))
  
  #print(j)
  aggregate(END ~ sRNA, data = run4, max)
  
}
stopCluster(cl)


  ## remove unrequired sequences
multi.reduced<-fread(paste0(sRNA_level.of.collapse,"_multireduced.rds")) # contains all assigned sRNA:seqs prior to filtering


## remove unrequired seq:miRs.
file.it<-paste0(location.of.patman.files,sRNA_level.of.collapse,"_collpased_",output)

assigned<-fread(file.it)

multi.reduced<-multi.reduced[multi.reduced$sRNA%in%assigned$sRNA,]


max.length.1<-max.length.1[max.length.1$sRNA%in%multi.reduced$sRNA,] # now we will only search for miRs that were kept


max.length<-aggregate(END ~ sRNA, data = max.length.1, max) # get the max END nt of all the sequences mapping to each loci
max.length<-max.length[order(max.length$sRNA),]


assigned_sRNAs<-paste(multi.reduced$sRNA,multi.reduced$SEQ) # create a key to remove non-assigned miR:seq pairs (mutlimapping seq:sRNA pairs that were trorwn out when assigning multimapping reads)  in the loop below


############################

#2. create the expression matrix
cat("\n casting and filling the expression matracies \n")
# read in file -- create a vector of 0s --- replace covered nt with count

num.cores<-detectCores()
cl <- makeCluster(num.cores-1)
registerDoParallel(cl)


global.exp.prof<-foreach(j =1:length(files), .packages=c('data.table','stringi','stringr','plyr','dplyr'), .combine = rbind) %dopar% { # revert
  run4<-fread(files[j],sep=" ",header=F) # put correct file path
  
  colnames(run4)[2:5]<-c("SEQ","sRNA","START","FLAG")
  dups<-paste0(run4$SEQ,run4$sRNA)
  run4<-run4[!duplicated(dups),]
  
  run4<-run4%>%mutate(END=(START+nchar(SEQ)))
  
  hp<-paste(run4$sRNA,run4$SEQ)
  run4<-run4[hp%in%assigned_sRNAs,] ## hp_sequences that made it through

  
  
  
  data<-list() # create a list to populate
  
  
  
  for(i in 1:length(max.length$sRNA)){
    
    if(setequal(intersect(max.length$sRNA[i],run4$sRNA), max.length$sRNA[i])==T) # does my sample contain the miRNA of interest?
    {
      
      this.one<-run4[run4$sRNA==max.length$sRNA[i],]  # miR
      
      temp<-matrix(0,nrow=nrow(this.one),ncol=max.length[i,"END"])
      
      for(k in 1:nrow(this.one)){ # each sequence assigned to the miR
        
        temp[k,this.one$START[k]:this.one$END[k]]<-temp[k,this.one$START[k]:this.one$END[k]]+this.one[k,1] # assign over the 0's with nucleotide frequencies
        print(k)
      }
      data[[i]]<-colSums(temp)}else
      {
        data[[i]]<-matrix(0,nrow=1,ncol=max.length[i,"END"])
      }
    
    data[[i]]<-stri_replace_all_regex(toString(data[[i]]), ",", "")
    data[[i]]<-c(max.length[i,"END"],data[[i]])
    print(i)
  }
  
  exp.prof<-data.frame(matrix(unlist(data), nrow=length(data), byrow=T))
  sample.name<-(j-1)
  sample.file.name<-stri_sub(basename(files[j]),1,19)
  exp.prof2<-cbind(rep(sample.name,nrow(exp.prof)),
                   max.length$sRNA,
                   exp.prof,
                   rep(sample.file.name,nrow(exp.prof)))
  
}
stopCluster(cl)

### not enough memory if I run all 352 samples: therefore ?chunk up the job and then combine at the end?? 


colnames(global.exp.prof)<-c("sample","miR","max_length","expression","sample.file.name")
global.exp.prof<-global.exp.prof[order(global.exp.prof$miR),]
##### still need to append max length of miR, otherwise nearly done. Code runs quickly :)

index.for.noise<-cbind(unique(global.exp.prof$sample),unique(as.character(global.exp.prof$sample.file.name)))
write.csv(index.for.noise, file=paste0(location.of.patman.files,index.file),row.names = F)


expression.profile.file_name<-paste(prefix,sRNA_level.of.collapse_for_noise, "sample_expression_profiles.txt", sep="_")


write.table(file=expression.profile.file_name,
            col.names = F,
            row.name =F,
            sep=", ",
            quote=FALSE,
            global.exp.prof[,1:4])

print(proc.time()-ptm)

colnames(index.for.noise)<-c("sample number","sample name")


print(paste0("Your expression matracies are complete. To identify the noise threshold in each sample please submit this file (",expression.profile.file_name,") to the program: Irina_noise_detection.pl"))

