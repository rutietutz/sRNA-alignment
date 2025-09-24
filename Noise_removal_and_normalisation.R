
## Install/load packages as required #####

library(preprocessCore)


# ###########  remove the noise ################

## read in the matracies of gene expression, find the noise, remove the noise, add the noise offset.

input<-file.it # file you want to remove the noise from  - by default it is inherited as the collapsed file that has just been produced.

noise_removed<-paste(prefix,sRNA_level.of.collapse_for_noise, "notNorm_noise_removed.csv", sep="_")  # name of the desired output file

noise_removed_offset_added<-paste(prefix,sRNA_level.of.collapse_for_noise, "notNorm_noise_removed_offset_added.csv", sep="_")  # name of the desired output file # name of the desired output file

noise_removed_offset_added_normalised<-noise_removed_offset_added<-paste(prefix,sRNA_level.of.collapse_for_noise, "notNorm_noise_removed_offset_added_quantile_normalised.csv", sep="_")# name of the desired output file

noise_removed_offset_added_pre_process_core_normalised<-noise_removed_offset_added<-paste(prefix,sRNA_level.of.collapse_for_noise, "notNorm_noise_removed_offset_added_quantile_normalised.via.preprocesscore.csv", sep="_")# name of the desired output file

detection_cutoff<-as.numeric(detection_cutoff)




header<-read.csv(index.file)

noise.matrix<-t(read.csv(noise.out, header=F))

noise.matrix<-cbind(header,noise.matrix)

colnames(noise.matrix)[c(3,4)]<-c("First_threshold","Second_threshold")


noise.matrix<-noise.matrix[order(noise.matrix[,1]),]  # order by sample name in the noise file


collapsed.sRNA.noise.incl<-read.csv(input,row.names = 1)


rv<-cbind(colnames(collapsed.sRNA.noise.incl),as.character(noise.matrix[,2]))


noise.matrix[,2]<-gsub("merged_","",noise.matrix[,2])

if((F%in%(colnames(collapsed.sRNA.noise.incl)==noise.matrix[,2]))==T){
  cat("ERROR samples names in the noise matrix and expression file do not agree \n ")
}else{
  cat(summary(colnames(collapsed.sRNA.noise.incl)==noise.matrix[,2]),"\n ")
  
  
  noise<-2^as.numeric(as.character(noise.matrix[,3]))
  
  collapsed.sRNA.noise.excl<-collapsed.sRNA.noise.incl # a backup point
  
  for (i in 1:ncol(collapsed.sRNA.noise.incl)){
    noisy_sequences<-collapsed.sRNA.noise.excl[,i]<noise[i] # noisy sequences == the genes which have a count value that is smaller than the noise threshold
    collapsed.sRNA.noise.excl[noisy_sequences,i]<-0 # replace noisy sequences with 0
  }
  
  cat(paste("you have this many sRNAs in your after collapse:", nrow(collapsed.sRNA.noise.excl)),"\n ")
  
  
  
  collapsed.sRNA.noise.excl<-collapsed.sRNA.noise.excl[rowSums(collapsed.sRNA.noise.excl)>0,] # remove any genes which are now not detected in any samples
  
  cat(paste("you have this many sRNAs in your after removing the noise:", nrow(collapsed.sRNA.noise.excl)),"\n ")
  
  dc<-floor(as.numeric(n.samps)*detection_cutoff)
  
  collapsed.sRNA.noise.excl.detection_threshold<-collapsed.sRNA.noise.excl[rowSums(collapsed.sRNA.noise.excl>0)>=dc,] # remove any genes which are now not detected in any samples
  
  cat(paste("you have this many sRNAs in your after removing the noise:", nrow(collapsed.sRNA.noise.excl.detection_threshold)),"\n ")
  
  
  
  # add the threshold to any samples with a 0
  collapsed.sRNA.noise.excl.offset.added<-collapsed.sRNA.noise.excl
  
  for(i in 1: ncol(collapsed.sRNA.noise.excl.offset.added))
  {
    collapsed.sRNA.noise.excl.offset.added[,i][collapsed.sRNA.noise.excl.offset.added[,i]==0]<-noise[i]
    cat(i,"\n")
  }
  
  write.csv(collapsed.sRNA.noise.excl,noise_removed)
  
  write.csv(collapsed.sRNA.noise.excl.offset.added,noise_removed_offset_added)
  
  # function for quantile normalising count data - need to d/w Irina about when you reassign the 0's
  
  quantile.normalise.count.data<-function(count.matrix){
    
    
    base.frame<-apply(count.matrix,2,function(x) rank(x,ties.method="min"))
    sorted.frame<-apply(count.matrix,2,sort)
    row.means<-rowMeans(sorted.frame)
    
    
    quantile.norm.count.data<-function(x){
      
      c<-x
      
      for(i in 1:length(unique(x)))
      {
        tie.n<-(x==unique(x)[i]) # The positions which all have the tied count
        new.mean<-mean(row.means[tie.n]) # get  mean of row.means at the tied positions
        c[tie.n]<-new.mean # assign the new mean to each tied position
      }
      return(c)
    }
    
    
    qnorm<-apply(sorted.frame,2,quantile.norm.count.data) # get the quantile normalised counts
    
    
    for(j in (1:ncol(qnorm)))
    {
      qnorm.val<-qnorm[,j]
      for(i in 1:nrow(base.frame))
      {
        rank<-base.frame[i,j]
        base.frame[i,j]<-qnorm.val[rank]
      }
    }
    
    collapsed.sRNA.noise.excl.offset.added.quantile<-base.frame
    collapsed.sRNA.noise.excl.offset.added.quantile[count.matrix==0]<-0 # reassign anything that was a 0 as a 0
    
    return(collapsed.sRNA.noise.excl.offset.added.quantile)
  }
  
  # # I have checked the quantile.normalise.count.data function against one of irina's quantile normalised outputs and the answer is the same
  
  collapsed.sRNA.noise.excl.offset.added.quantile<-quantile.normalise.count.data(collapsed.sRNA.noise.excl.offset.added)
  
  rownames<-rownames(collapsed.sRNA.noise.excl.offset.added)
  colnames<-colnames(collapsed.sRNA.noise.excl.offset.added)
  
  
  
  collapsed.sRNA.noise.excl.offset.added.quantile.auto.method<-normalize.quantiles(as.matrix(collapsed.sRNA.noise.excl.offset.added))
  
  
  rownames(collapsed.sRNA.noise.excl.offset.added.quantile.auto.method)<-rownames
  colnames(collapsed.sRNA.noise.excl.offset.added.quantile.auto.method)<-colnames
  
  
  write.csv(collapsed.sRNA.noise.excl.offset.added.quantile,file=noise_removed_offset_added_normalised)
  
  write.csv(collapsed.sRNA.noise.excl.offset.added.quantile.auto.method,file=noise_removed_offset_added_pre_process_core_normalised)
  
  
  dc<-floor(as.numeric(n.samps)*detection_cutoff)
  
  collapsed.sRNA.noise.excl.detection_threshold<-collapsed.sRNA.noise.excl[rowSums(collapsed.sRNA.noise.excl>0)>=dc,] # remove any genes which are now not detected in any samples
  
  cat(paste("you have this many sRNAs in your after removing the noise:", nrow(collapsed.sRNA.noise.excl.detection_threshold)),"\n ")
  
  
  # add the threshold to any samples with a 0
  
  for(i in 1: ncol(collapsed.sRNA.noise.excl.detection_threshold))
  {
    collapsed.sRNA.noise.excl.detection_threshold[,i][collapsed.sRNA.noise.excl.detection_threshold[,i]==0]<-noise[i]
    #cat(i,"\n")
  }
  
  
  collapsed.sRNA.noise.excl.detection_threshold.quantile.auto.method<-normalize.quantiles(as.matrix(collapsed.sRNA.noise.excl.detection_threshold))
  
  
  rownames(collapsed.sRNA.noise.excl.detection_threshold.quantile.auto.method)<-rownames(collapsed.sRNA.noise.excl.detection_threshold)
  colnames(collapsed.sRNA.noise.excl.detection_threshold.quantile.auto.method)<-colnames(collapsed.sRNA.noise.excl.detection_threshold)
  
  
  cat(paste("you have this many sRNAs in your after removing the noise and reads not detected in",detection_cutoff*100,"% of samples:", nrow(collapsed.sRNA.noise.excl.detection_threshold)),"\n ")
  
  
  final.file<-paste0("detected_in_",detection_cutoff*100,"_percent_",noise_removed_offset_added_pre_process_core_normalised)
  
  write.csv(collapsed.sRNA.noise.excl.detection_threshold.quantile.auto.method,file=final.file)
  
  
  cat(paste("Your samples have had the noise removed, an offset added, and been normalised. The following files have been created:\n  ",
            noise_removed,
            "\n  ",
            noise_removed_offset_added,
            "\n  ",
            noise_removed_offset_added_normalised,
            "\n  ",
            noise_removed_offset_added_pre_process_core_normalised,
            "\n  ",
            final.file
  ))
  
}

