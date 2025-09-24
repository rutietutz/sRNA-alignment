#Script to generate the expression profiles for each miRNA in one sample.


index.file<-paste(prefix,sRNA_level.of.collapse_for_noise, "index.for.noise.analysis.csv", sep="_")

#1. read in the file


ptm<-proc.time() # start the clock 

if(grepl("2MM_mature miR allignments",wd)==T){
files<-list.files(path=location.of.patman.files, 
                  pattern="merged.*hp.alligned.sum",
                  recursive=T)}else{
files<-list.files(path=location.of.patman.files, 
                  pattern="merged.*position.sum",
                  recursive=T)}

files<-paste0(location.of.patman.files,files)



num.cores<-detectCores()
cl <- makeCluster(num.cores-1)
registerDoParallel(cl)

cat("\n Finding the max length for each sRNA")

max.length.1<-foreach(j =1:length(files), .packages=c('data.table','stringi','stringr'), .combine = rbind) %dopar% {
  
  run4<-fread(files[j],sep="\t",header=T) # put correct file path
  run4<-data.frame(run4)
  if(grepl("2MM",wd)){
    run4<-fread(files[j],sep=" ",header=F) # put correct file path
    run4<-data.frame(run4)
    #### create the end positioin column
    END<-nchar(run4$V2)
    END<-nchar(run4$V2)+run4$V4
    run4<-cbind(run4[,c(1:4)],END)
    colnames(run4)<-c("count","SEQ","MIR","START","END")
  }
  if(needs_cutting==T){
    run4$MIR<-word(run4$MIR, 1, cut, sep="-")
  }
  run4<-run4[!duplicated(run4),]
  print(j)
  aggregate(END ~ MIR, data = run4, max)
}
stopCluster(cl)


# if(as.numeric(n.samps)<=176){
# max.length.1<-foreach(j =1:length(files), .packages=c('data.table','stringi'), .combine = rbind) %dopar% {
#   
#   run4<-fread(files[j],sep="\t",header=T) # put correct file path
#   run4<-data.frame(run4)
#   print(j)
#   aggregate(END ~ MIR, data = run4, max)
# }
# }else{
#   n<-as.numeric(n.samps)
#   size<-5
#   chunk<-floor(n/size)
#   
#   for(c in 1:chunk){
#     assign(
#       (paste0("temp_",c)),
#       (foreach(j =(((c-1)*size)+1):(c*size), .packages=c('data.table','stringi'), .combine = rbind) %dopar% {
#         
#         run4<-fread(files[j],sep="\t",header=T) # put correct file path
#         run4<-data.frame(run4)
#         aggregate(END ~ MIR, data = run4, max)}))
#   }
# 
# d<-c
# 
# for(c in ((d*size)+1):(length(files))){
#   assign(
#     (paste0("temp_",c)),
#     (foreach(j = c, .packages=c('data.table','stringi'), .combine = rbind) %dopar% {
#       
#       run4<-fread(files[j],sep="\t",header=T) # put correct file path
#       run4<-data.frame(run4)
#       aggregate(END ~ MIR, data = run4, max)}))
# }
# 
# }




## remove unrequired sequences

## remove unrequired seq:miRs.

assigned_sRNAs<-fread(paste0(sRNA_level.of.collapse,"_multireduced.rds"))

if(sRNA_level.of.collapse=="miRNA"){
assigned_sRNAs$sRNA<-gsub("-3p","",assigned_sRNAs$sRNA)
assigned_sRNAs$sRNA<-gsub("-5p","",assigned_sRNAs$sRNA)
assigned_sRNAs$sRNA<-gsub("miR","mir",assigned_sRNAs$sRNA)
# remove miRs that were discarded during assignment
hp<-word(max.length.1$MIR, 1, 3, sep="-")
max.length.1<-max.length.1[hp%in%assigned_sRNAs$sRNA,] # now we will only search for miRs that were kept
}else{
  max.length.1<-max.length.1[max.length.1$MIR%in%assigned_sRNAs$sRNA,] # now we will only search for miRs that were kept
}


max.length<-aggregate(END ~ MIR, data = max.length.1, max)

max.length<-max.length[order(max.length$MIR),]

## all sRNAs are in hairpins....

assigned_sRNAs<-paste(assigned_sRNAs$sRNA,assigned_sRNAs$seq) # create a key to remove non-assigned miR:seq pairs in the loop below


#2. create the expression matrix
cat("\n casting and filling the expression matracies \n")
# read in file -- create a vector of 0s --- replace covered nt with count

num.cores<-detectCores()
cl <- makeCluster(num.cores-1)
registerDoParallel(cl)


global.exp.prof<-foreach(j =1:length(files), .packages=c('data.table','stringi','stringr'), .combine = rbind) %dopar% {
  
  run4<-fread(files[j],sep="\t",header=T) # put correct file path
  run4<-data.frame(run4)
  if(grepl("2MM",wd)){
    run4<-fread(files[j],sep=" ",header=F) # put correct file path
    run4<-data.frame(run4)
    #### create the end positioin column
    END<-nchar(run4$V2)
    END<-nchar(run4$V2)+run4$V4
    run4<-cbind(run4[,c(1:4)],END)
    colnames(run4)<-c("count","SEQ","MIR","START","END")
  }
  if(needs_cutting==T){
    run4$MIR<-word(run4$MIR, 1, cut, sep="-")
  }
  run4<-run4[!duplicated(run4),]
  if(sRNA_level.of.collapse=="miRNA"){
    hp<-word(run4$MIR, 1, 3, sep="-")
    hp<-paste(hp,run4$SEQ)
  }else{
    hp<-paste(run4$MIR,run4$SEQ)
    run4<-run4[hp%in%assigned_sRNAs,] ## hp_sequences that made it through
  }
  
  
  
  data<-list() # create a list to populate
  
  
  
  for(i in 1:length(max.length$MIR)){
    
    if(setequal(intersect(max.length$MIR[i],run4$MIR), max.length$MIR[i])==T) # does my sample contain the miRNA of interest?
    {
      
      this.one<-run4[run4$MIR==max.length$MIR[i],]  # miR
      
      temp<-matrix(0,nrow=nrow(this.one),ncol=max.length[i,2])
      
      for(k in 1:nrow(this.one)){ # each sequence assigned to the miR
        temp[k,this.one[k,4]:this.one[k,5]]<-temp[k,this.one[k,4]:this.one[k,5]]+this.one[k,1] # assign over the 0's with nucleotide frequencies
        print(k)
      }
      data[[i]]<-colSums(temp)}else
      {
        data[[i]]<-matrix(0,nrow=1,ncol=max.length[i,2])
      }
    data[[i]]<-stri_replace_all_regex(toString(data[[i]]), ",", "")
    data[[i]]<-c(max.length[i,2],data[[i]])
    print(i)
  }
  
  exp.prof<-data.frame(matrix(unlist(data), nrow=length(data), byrow=T))
  sample.name<-(j-1)
 
  
  sample.file.name<-stri_sub(basename(files[j]),1,23)
  sample.file.name
  exp.prof<-cbind(rep(sample.name,nrow(exp.prof)),
                  max.length$MIR,
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

# took around 46 sec for 44 samples with miR.hp files



# The pcc and abn outputs are then analysed with this script: "signal_noise_thr_analysis.R" I randomly found this in a joint folder that Irina and I have when I was looking for the P2P graphs. :O
# I see in a readme file she used: "the noise detection was done on median > 0.6 (correlation threshold)"


## Confirmed -- there is probably a bug in Irina's script: she  -- plan -- inform her. TBH won't make a big difference if using lots of samples. I only worked it out because I ran a file with 3 samples in it and noticed that the mean was the sum of the P2Ps instead of the mean. 
## Also of note: for p2p you work out the p2p for a transcript in each file then you find the average of this. 


# I think I can draw the graphs as I already have the code -- Irina gave graphs of miRNA expression for each sample, pre and post.


