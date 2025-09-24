# Amend as start site now included.

# This program will collapse sRNAome sequences in four ways: A)  it will collapse multimapping sequences to the most abundant feature NOT in a hierachial manner B) it will collapse multimapping sequences to the most abundant feature IN a hierachial manner C) it will collapse multimapping sequences to a union gene. D) it will collapse using fractional counts.

# decisions about ribodepletion. rRNA mapping sequences highly highly abundant!! If particular sequences are enriched (it coud be PCR bias, or suprainmposed expression of a miR). As abundant as miR-21 ish. ribodepleting first will remove 10s e.g. 65 miRs, some of which may not have passed noise threshold anyway. Ribodepleting could make a difference to the expression of some sRNAs. I think I will collapse twice - one with ribodepletion first, and one with ribosomal RNA left in. I can then do DE on either file.

# issue is probably something to do with ribodepletion because where a seq has mapped to a ribosome by abdundance, it has not been assigned as a hierachy







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

if(exists("ribodepleted")==F){
  ribodepleted<-"not_ribo_depleted"
  }

###  Create totals based on all reads -- I will use this to decide re assignment ####

cat("Create totals based on all reads -- I will use this to decide re assignment \n")
expt.sum.file<-list.files(location.of.patman.files,pattern = "experimentlevel.sum",full.names = T)
totals<-fread(expt.sum.file) # this is the expt.sum file that was created by totalling up all the sRNA allignment.sum files. I used to have to manually add the mature miRs the sRNAome folder before running the create_experimental_sum_file.sh script, but  the SRNAome allignemnt pipeline now also includes the mature miRs (they are appended onto the SRNAome.sum files by the allignemnt script) so this is no longer required. Before creating the final totals table, I will need to remove mature miRs that mapped to a tRNA or something else with fewere mismatches. This is achieved below. 

if(ncol(totals)>3){
  totals<-totals[,1:3]
}
colnames(totals)<-c("totals","seq","sRNA")


# remove the entries which relate to sequences that map to miR.haripins when they are already mapped to a mature miR, whilst preserving sequences that map to a miR hairpin but not a mature arm. Also remove sequences that map to a mature miRNA but not a miR hairpin (e.g. because a sequence maps to another sRNA type with fewer mismatches so did not assign to a hairpin using a best straata approach with the the sRNAome)


totals$sRNA_type<-gsub("__.*","",totals$sRNA,perl = TRUE)

# create a new column called sRNA type which is derived from teh first bit of the mapped read name.

# e.g....
     # totals                          seq                                  sRNA sRNA_type
     #  1        AGCCGGATCCGTAACTTCGAG rRNA__45S_pre-ribosomal_N4_(RNA45SN4)      rRNA
     #  1 TCGATTCCGTGGGTGGTGGTACATGGCC rRNA__45S_pre-ribosomal_N2_(RNA45SN2)      rRNA
     #  3        ACAGCAGAACGGTGGCCATGG rRNA__45S_pre-ribosomal_N5_(RNA45SN5)      rRNA
     #  14           CTGGTGGTTAGGATTCGG   tRNA__Homo_sapiens_tRNA-Glu-CTC-1-2     tRNA


totals<-data.frame(totals)

if(ribodepleted=="ribodeplete"){
  ribozomal<-totals$seq[totals$sRNA_type=="rRNA"] # identify ribozomal sequences and turn these into a list
  totals<-totals[!totals$seq%in%ribozomal,]# remove all entries containing these sequences
  # rename file output as ribodepleted
}

if(sum(grepl("hsa-miR|hsa-let",totals$sRNA_type))>0){
totals<-data.table(totals)
mat.miR<-grepl("hsa-miR|hsa-let",totals$sRNA_type)&!grepl("MiR.hp",totals$sRNA_type) # sequences that map to a mature miRNA
totals$sRNA_type[mat.miR]<-"mature_MiR"  # rename sequences that map to a mature miRNA as mature miR in their sRNA type
miR.mapping.reads<-totals[grepl("MiR",sRNA_type)] # extract out the miRNA mapping reads as a seperate dataframe ==A
miR.mapping.reads<-miR.mapping.reads[!duplicated(paste0(seq,sRNA_type))] # remove duplicates based on seq_sRNA_type. Don't worry, you are not losing anything, you are just trying to identify sequences that map to both a hp and a mature miR
table<-data.frame(table(miR.mapping.reads$seq)) # Table up the sequences. (seqs that appear twice must map to both a mature miR and a hp)
colnames(table)[1]<-"seq"
totals<-left_join(totals,table) # append tabled sequences to A (the dataframe of miRNA mapping reads) using merge or join. left_join works as a merge, it will append the tables df to the totals, only where they are the same in the seq column (left join looks for colnames with the same name in two data tables and uses this to do the merging). Even though duplicated seqs were deleted in the "miR.mapping.reads" df, they still exist in the totals df, and they will have the correct Freq addded because the freq is added based on the the sequence. I have double checked that this is the case. 
remove<-(totals$Freq==2&totals$sRNA_type=="MiR.hp")|(totals$Freq==1&totals$sRNA_type=="mature_MiR") # if A$sRNA_type = miR.hp and A$table=2 then -->remove (because the mature miR neeeds to be kept) OR if a sequence only maps to a mature miR and not a hairpin then remove this (this situation arises because mature miR mappings are done on sequences which map to miR hairpins before mapping to sRNAome. Such sequences may not map to a miR hairpin when mapping to sRNAome best straata)
totals<-totals[!remove,] # remove the entries which relate to sequences that map to miR.haripins when they are already mapped to a mature miR, whilst preserving sequences that map to a miR hairpin but not a mature arm.
}

totals<-data.frame(totals)
totals<-totals[,!grepl("sRNA_type|Freq",colnames(totals))] # remove the sRNA_type and Freq columns

totals2<-data.table(totals)

totals2<-aggregate(totals ~ sRNA,totals, sum)

setorder(totals2, -totals) # order totals by size

totals.hierachy<-totals2

totals.hierachy$sRNA_type<-gsub("__.*","",totals.hierachy$sRNA,perl = TRUE)

mat.miR<-grepl("hsa-miR|hsa-let",totals.hierachy$sRNA_type)&!grepl("MiR.hp",totals.hierachy$sRNA_type) # sequences that map to a mature miRNA
totals.hierachy$sRNA_type[mat.miR]<-"mature_MiR"  # rename sequences that map to a mature miRNA as mature miR in their sRNA type


totals.hierachy$sRNA_type <- factor(totals.hierachy$sRNA_type, levels = c("mature_MiR","MiR.hp","yRNA", "tRNA", "snRNA","snoRNA","vault","viral.miR.hp","lncRNA","rRNA"))

totals.hierachy<-arrange(totals.hierachy, sRNA_type,-totals)

totals.hierachy$hierachy_order<-1:nrow(totals.hierachy) # create a new column called hierachy_order, and fill it from 1 to the number of rows.

remove_seqs_that_map_this_many_times<-"no_limit"

# ################
#
## prep for filtering step
min.read.freq=1 # kept for script compatibility
min.sample=1 # kept for script compatibility

must.have.at.least<-(min.read.freq)*(min.sample+1)


totals3<-totals2[totals2$totals>=must.have.at.least,]

saveRDS(totals3,"totals3.RDS")


files<-list.files(path=location.of.patman.files, pattern=".sum",recursive=F)

files<-files[grep("merged.*.sum",files,invert = F)] # Only use the merged files.

files<-paste0(location.of.patman.files,files)



#prep for parallelising
num.cores<-detectCores()



totals3<-readRDS(paste0(location.of.patman.files,"totals3.RDS"))

totals3$rank<-frank(totals3$totals,ties.method="dense")



totals3.one.level<-totals2[totals2$totals>=must.have.at.least,]

totals3.one.level$class<-gsub("__.*","",totals3.one.level$sRNA,perl = TRUE)

# tRNA

totals3.one.level.tRNA.aa<-totals3.one.level[totals3.one.level$class=="tRNA",]


# collapse by totals3.one.level.tRNA.aa


n <- 2
pat <- paste0('^([^-]+(?:-[^-]+){',n-1,'}).*')

totals3.one.level.tRNA.aa$sRNA.collapsed <-sub(pat, '\\1',totals3.one.level.tRNA.aa$sRNA) 




# rRNA

# there are not that many rRNA sRNAs (9 sRNAs, and max copy is 5 lots of 45S_pre-ribosomal_N1_(RNA45SN1)), therefore I will not collapse

# snRNA

totals3.one.level.snRNA.collapse.one.level<-totals3.one.level[totals3.one.level$class=="snRNA",]

# there are 26 snRNAs

# I can collapse to one level above e.g. variant_U1_small_nuclear_8 and variant_U1_small_nuclear_7 will be collapse, but may do DE at raw level.


totals3.one.level.snRNA.collapse.one.level$sRNA.collapsed<-gsub(".*Homo_sapiens","Homo_sapiens",totals3.one.level.snRNA.collapse.one.level$sRNA) # need to remove the URS0000759C03 identifier as these are specific to the variant type and will prevent consoliation by sRNA name if left in 

totals3.one.level.snRNA.collapse.one.level$sRNA.collapsed%<>%gsub("nuclear_.*","nuclear",.)

totals3.one.level.snRNA.collapse.one.level$sRNA.collapsed<-paste0(
  totals3.one.level.snRNA.collapse.one.level$class,
  "__",
  totals3.one.level.snRNA.collapse.one.level$sRNA.collapsed)



# snoRNA

totals3.one.level.snoRNA.collapse.one.level<-totals3.one.level[totals3.one.level$class=="snoRNA",]

# there are 648 snoRNAs

# I can collapse to one level above e.g. Homo_sapiens300029_SNORD114-21 and Homo_sapiens300029_SNORD114-24 will be collapsed

totals3.one.level.snoRNA.collapse.one.level$sRNA.collapsed<-gsub(".*SNORD","SNORD",totals3.one.level.snoRNA.collapse.one.level$sRNA) %>%
  gsub(".*SNORA","SNORA",.)%>%
  gsub(".*SCARNA","SCARNA",.)
# need to remove the Homo_sapiens300044 identifier as these are specific to the variant type and will prevent consoliation by sRNA name if left in 

# there are a 4 odd names left in e.g. Homo_sapiens300061_AL137790.4 but I can leave these in as do not appear to be duplicated (though maybe they are but have a very differnt name to thier copy)

# a single "-" is present in copies, therefore simply remove everything after this "-"


totals3.one.level.snoRNA.collapse.one.level$sRNA.collapsed%<>%gsub("-.*","",.)

totals3.one.level.snoRNA.collapse.one.level$sRNA.collapsed<-paste0(
  totals3.one.level.snoRNA.collapse.one.level$class,
  "__",
  totals3.one.level.snoRNA.collapse.one.level$sRNA.collapsed)

# reduced from 648 to 527



# yRNA

totals3.one.level.yRNA<-totals3.one.level[totals3.one.level$class=="yRNA"|totals3.one.level$sRNA=="snoRNA__Homo_sapiens300395_SNORA16AL1yRNA__Y5",]

totals3.one.level.yRNA$sRNA.collapsed<-totals3.one.level.yRNA$sRNA

totals3.one.level.yRNA$sRNA.collapsed[totals3.one.level.yRNA$sRNA.collapsed=="snoRNA__Homo_sapiens300395_SNORA16AL1yRNA__Y5"]<-"yRNA__Y5"

# there are only 4 yRNAs and no multiple copies however because of the blip before with the yRNA Y5 I do need to collapse Y5


# lncRNA

# there are 6044 lncRNAs - I have no idea if any can be collapsed or how to easily spot which ones can bere therefore leave at this level.


totals3.one.level$sRNA.collapsed<-totals3.one.level$sRNA

totals3.one.level<-totals3.one.level[!grepl("tRNA",totals3.one.level$class),] %>%
  rbind(.,totals3.one.level.tRNA.aa)


totals3.one.level<-totals3.one.level[!grepl("snRNA",totals3.one.level$class),] %>%
  rbind(.,totals3.one.level.snRNA.collapse.one.level)


totals3.one.level<-totals3.one.level[!grepl("snoRNA",totals3.one.level$class),] %>%
  rbind(.,totals3.one.level.snoRNA.collapse.one.level)

totals3.one.level<-totals3.one.level[!(totals3.one.level$class=="yRNA"|totals3.one.level$sRNA=="snoRNA__Homo_sapiens300395_SNORA16AL1yRNA__Y5"),] %>%
  rbind(.,totals3.one.level.yRNA)





#### collapse each file #######

for (i in 1:length(files)){


  run4<-fread(files[i],sep=" ",header=F,fill=T) # change

  if(ncol(run4)>3){ # if the start position is present then it the file will have more than 3 columns, (start positioins and strand in final columns -- remove these as not needed for collapsing)
    run4<-run4[,1:3]
  }




#  colnames(run4)[1]<-stri_sub(basename(files[i]),1,23) #made need adding
  colnames(run4)[2:3]<-c("SEQ","sRNA")


  test2<-run4 # a backup




  colnames(test2)[c(1,ncol(test2))]<-c("count","sRNA")

  test2$sRNA_type<-gsub("__.*","",test2$sRNA,perl = TRUE)
  
  
  if(ribodepleted=="ribodeplete"){
    ribozomal<-test2$SEQ[test2$sRNA_type=="rRNA"] # identify ribozomal sequences and turn these into a list
    test2<-test2[!test2$SEQ%in%ribozomal,]# remove all entries containing these sequences
    # rename file output as ribodepleted
  }

  test3<-data.table(test2)
  
  if(sum(grepl("hsa-miR|hsa-let",test2$sRNA_type))>0){
    

  ## set aside single mapping reads

  # remove the entries which relate to SEQuences that map to miR.haripins when they are already mapped to a mature miR, whilst preserving SEQuences that map to a miR hairpin but not a mature arm.test3$sRNA_type<-reshape2::colsplit(test3$sRNA,pattern = "__",names=c(1:100))[,1]
  mat.miR<-grepl("hsa-miR|hsa-let",test3$sRNA_type)&!grepl("MiR.hp",test3$sRNA_type)
  test3$sRNA_type[mat.miR]<-"mature_MiR"  # rename mature miRs as mature miR in their sRNA type
  miR.mapping.reads<-test3[grepl("MiR",sRNA_type)] # extract out the miRNA mapping reads ==A
  miR.mapping.reads<-miR.mapping.reads[!duplicated(paste0(SEQ,sRNA_type))] # remove duplicates based on SEQ_sRNA_type
  table<-data.frame(table(miR.mapping.reads$SEQ)) # Table up the SEQuences
  colnames(table)[1]<-"SEQ"
  test3<-left_join(test3,table) # append tabled SEQuences to A using merge or join
  remove<-(test3$Freq==2&test3$sRNA_type=="MiR.hp")|(test3$Freq==1&test3$sRNA_type=="mature_MiR") # if A$sRNA_type = miR.hp and A$table=2 then -->remove OR if a sequence only maps to a mature miR and not a hairpin then remove this (this situation arises because mature miR mappings are done on sequences which map to miR hairpins before mapping to sRNAome. Such sequences may not map to a miR hairpin when mapping to sRNAome best straata)
  test3<-test3[!remove,] # remove the entries which relate to SEQuences that map to miR.haripins when they are already mapped to a mature miR, whilst preserving SEQuences that map to a miR hairpin but not a mature arm.
  }

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




  #######  Assign the sRNA based on the sRNA abundances ######


  reduced.final<-multi.sRNA ## for down stream compatibility with old script
  reduced.final<-data.frame(reduced.final)

  # order reduced.final by window abundance.

  reduced.final$rank<-totals3$rank[match(reduced.final$sRNA,totals3$sRNA)]

  reduced.final<-reduced.final[order(reduced.final$rank,decreasing=T),]

  # remove sequences mapping > specific no of times times

  table<-data.frame(table(reduced.final$SEQ))

  if(sum(grepl("no_limit",remove_seqs_that_map_this_many_times))==0)
  {
    NR.sequences<-table$Var1[table$Freq<=as.numeric(as.character(remove_seqs_that_map_this_many_times))]

    reduced.final<-reduced.final[reduced.final$SEQ%in%NR.sequences,]
  }

  reduced.final$no.mapping.sites<-table$Freq[match(reduced.final$SEQ,table$Var1)]
  setkey(multi.sRNA,SEQ)

  reduced.final=as.data.frame(reduced.final)

  multi.reduced=reduced.final[!duplicated(reduced.final$SEQ),] # sequences collapsed based on abdundances




  reduced.final$sRNA_type<-gsub("__.*","",reduced.final$sRNA,perl = TRUE)
  
  mat.miR<-grepl("hsa-miR|hsa-let",reduced.final$sRNA_type)&!grepl("MiR.hp",reduced.final$sRNA_type) # sequences that map to a mature miRNA
  reduced.final$sRNA_type[mat.miR]<-"mature_MiR"  # rename sequences that map to a mature miRNA as mature miR in their sRNA type
  
  
  reduced.final$sRNA_type <- factor(reduced.final$sRNA_type, levels = c("mature_MiR","MiR.hp","yRNA", "tRNA", "snRNA","snoRNA","vault","viral.miR.hp","lncRNA","rRNA")) # I am not sure if I really need this step. maybe I do - as I use it for a check step below.
  

############  reduced.final<-reduced.final[!reduced.final$sRNA_type=="rRNA",]

# I will basically left join the totals hierachy order column to the reduced.final df. I will do it manually. (? not sure why i didn't use left join or merge, perhaps they dfs were too big?)
  rv<-totals.hierachy[match(reduced.final$sRNA,totals.hierachy$sRNA),] # create a temporary df called rv. This is totals hierachy now expanded vertically to exactly match order of reduced.final. The 4th column called order, is the order in which sRNAs ended up being after they were sorted first by sRNA type and then by totals.

  table(reduced.final$sRNA_type==rv$sRNA_type) # a check to make sure the order of reduced.final is the same as rv.

  reduced.final$hierachy_order<-rv$hierachy_order # I have now appended the order (i.e. rank of the sRNAs) onto the reduced.final data frame.

# now only keep the first seq
  reduced.final.hierachy<-reduced.final[order(reduced.final$hierachy_order),] 

  multi.reduced.hierachy<-reduced.final.hierachy[!duplicated(reduced.final.hierachy$SEQ),] # keep only the first sequence record: this is the record which is best in terms of sRNA_type and then in terms of abundance. It allows sequences to be mapped exactly the same as I originally collapsed sequences, with the added benifit of allowing sequences which map better to a tRNA than an miRNA to be kept as a tRNA.



collapsed.sRNA.hierachy<-bind_rows(Ok.sRNA,multi.reduced.hierachy) # combine the single mapping sequences and the now reduced multimapping sequences


collapsed.sRNA.hierachy<-collapsed.sRNA.hierachy[,c("count","SEQ","sRNA","no.mapping.sites","sRNA_type")] # get rid of non-useful columns

colnames(collapsed.sRNA.hierachy)[colnames(collapsed.sRNA.hierachy)=="no.mapping.sites"]<-"no.mapping.sites.entire.siRNAome"


### also try collapsing to union gene
  cat("\n collapsing multimapping sequences - may take some time... \n")



  cl <- makeCluster(num.cores-1)
  registerDoParallel(cl)
  multi.reduced_union_gene<-foreach (j=1:length(multi.reduced$SEQ), .combine = c ) %dopar% {
  keep.record<-paste(reduced.final$sRNA[reduced.final$SEQ%in%multi.reduced$SEQ[j]],collapse = "&") # create union gene names
  }

  stopCluster(cl)

  ####
  multi.reduced$sites_seq_maps_to<-multi.reduced_union_gene # append the union gene name to multi-reduced. They are in the same order because multi.reduced_union_gene is made using a loop based on multi.reduced

  ####

collapsed.sRNA<-rbind(Ok.sRNA,multi.reduced,fill=T) # combine the single mapping sequences and the now reduced multimapping sequences

# append the sRNA decision that was made when collapsing sequences hierachically and then by abundance.
rv<-collapsed.sRNA.hierachy[match(collapsed.sRNA$SEQ,collapsed.sRNA.hierachy$SEQ)]


rm(collapsed.sRNA.hierachy)
gc()

table(rv$SEQ==collapsed.sRNA$SEQ) # check they are in the same order

collapsed.sRNA$hierachy_abundance<-rv$sRNA # append the sRNA decision that was made when collapsing sequences hierachically and then by abundance.




  file.name<-paste0(location.of.patman.files,stri_sub(basename(files[i]),1,23),"_temp3.txt")
  write.table(collapsed.sRNA,file.name,col.names = T,row.names = F) # very quick and possibly will save memory



#### also assign using fractional counts and save this output too.

  fractional.counts<-bind_rows(Ok.sRNA,reduced.final) # combine single mapping reads and all multimapping reads (i.e. non collapsed multimapping reads)

  fractional.counts<-fractional.counts[!grepl("rRNA",fractional.counts$sRNA_type),]# need to ribodeplete fractional counts, prior to fractioning the counts. This will remove any exclusive mapping ribosomal RNA and also remove fractioning of reads that map to both a miRNA and a ribosomal RNA. This is in keeping with the way that I do my other forms of collapsing, e.g. by hierachical abundance.


  fractional.counts$count<-fractional.counts$count/fractional.counts$no.mapping.sites

  fractional.counts[ ,c("SEQ","rank","no.mapping.sites" ,"sites_seq_maps_to","sRNA_type","order","hierachy_order" ) := NULL]

  # you can sum up your reads at this point to make merging files alot easier later. Because I will be merging these files seperately to the other forms of collapsing, I don't need to keep each sequence's identity because it will not impact on how reads are summed up. With the other forms of collapsing, because each sequence is only represented once, and will be collapsed in different ways, I have decided to merge those files all together in one go, and then sum up each sRNA in turn.

  fractional.counts.summed<-fractional.counts[, lapply(.SD, sum, na.rm=TRUE), by=sRNA, .SDcols=1]

  file.name.fractional<-paste0(location.of.patman.files,stri_sub(basename(files[i]),1,23),"_temp4.txt")

  fwrite(fractional.counts.summed,file.name.fractional,col.names = T,row.names = F) # very quick and possibly will save memory


  # now do a fractional count collapsed to up one level
  
  
  
  rv<-select(test3,!Freq)
  rv<-left_join(rv,totals3.one.level)
  rv%<>%select(!totals)
  rv$sRNA<-rv$sRNA.collapsed
  rv%<>%.[!duplicated(.),]
  freq<-count(rv,SEQ)
  rv<-left_join(rv,freq)
  rv$count<-rv$count/rv$n
  rv %<>% dplyr::select(count, sRNA)
  rv %<>% dplyr::group_by(sRNA) %>%
    dplyr::summarise(sum(count))
  rv %<>% dplyr::rename(count='sum(count)')
  
  rv<-rv[!is.na(rv$sRNA),]
  file.name.fractional<-paste0(location.of.patman.files,stri_sub(basename(files[i]),1,23),"_temp5.txt")
  
  fwrite(rv,file.name.fractional,col.names = T,row.names = F) # very quick and possibly will save memory
  
  
  print(paste("entry",i,"has completed"))
}


rm(fractional.counts.summed) # clear some space
rm(fractional.counts) # clear some space
rm(totals)
gc()

cat("\n multireduced sequences have been collapsed. Now merging the files \n")

files2<-dir(path=location.of.patman.files, pattern="_temp3.txt",full.names = T,recursive=F)


i<-1
run4<-fread(files2[i],sep=" ") # change
colnames(run4)[1]<-files2[i]

for(i in 2:length(files2)){
  run3<-fread(files2[i],sep=" ") # change
  colnames(run3)[1]<-files2[i]
  run4<-full_join(run4,run3)
  print(i)
}


run4[is.na(run4)]<-0

collapsed.sRNA<-run4

colnames(run4)[2:length(run4)] %<>%
  basename(.)%>%
  gsub("_merged_temp4.txt","",.)

rm(run4)
gc()

frag.length<-apply(collapsed.sRNA[,"SEQ",drop=F],1,stri_length) # length of each read

collapsed.sRNA.seq<-cbind(collapsed.sRNA,frag.length)

sRNA_level.of.collapse<-"sRNAome"

fwrite(collapsed.sRNA.seq, paste0(ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_multireduced.rds")) # backup file

rm(collapsed.sRNA.seq)
gc() # empties deletd objects from R's memory 

cat("\n Now summing up the counts for each siRNA loci based on abundance (no hierachy) \n")

file.name<-paste0(ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_multireduced.rds")

#remove_seqs_that_map_this_many_times<-gsub(" ","_",remove_seqs_that_map_this_many_times)
collapsed.sRNA<-fread(file=file.name,nThread=14)




# Now sum up the counts for each sRNA based on assignment to the most abundant read

collapsed.sRNA[ ,c("SEQ","rank","no.mapping.sites","sites_seq_maps_to","hierachy_abundance","frag.length","sRNA_type") := NULL]


collapsed.sRNA_majority_collapsed.noise.incl<-collapsed.sRNA[, lapply(.SD, sum, na.rm=TRUE), by=sRNA, .SDcols=c(1,3:ncol(collapsed.sRNA)) ] 

rm(collapsed.sRNA) # saves memory
gc()

file.it.1<-paste0(location.of.patman.files,ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_collpased_",output)  #  will need this file name later in the script

fwrite(collapsed.sRNA_majority_collapsed.noise.incl,file=file.it.1) # out put the collapsed file


rm(collapsed.sRNA_majority_collapsed.noise.incl)
gc()
# Now sum up the counts for each sRNA based on assignment to hierachical abdunant genes

collapsed.sRNA<-fread(file=file.name,nThread=14) ###  THIS IS WHERE IS BREAKS #######


collapsed.sRNA[ ,c("SEQ","rank","no.mapping.sites","sites_seq_maps_to","sRNA","frag.length","sRNA_type") := NULL] # keep only counts and hierachial abundant assignment

collapsed.sRNA_majority_collapsed.noise.incl<-collapsed.sRNA[, lapply(.SD, sum, na.rm=TRUE), by=hierachy_abundance, .SDcols=c(1,4:ncol(collapsed.sRNA)) ] 


rm(collapsed.sRNA) # saves memory
gc()

file.it.hierachy<-paste0(location.of.patman.files,ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_hierachically_collpased_",output)  #  will need this file name later in the script

fwrite(collapsed.sRNA_majority_collapsed.noise.incl,file=file.it.hierachy) # out put the collapsed file


rm(collapsed.sRNA_majority_collapsed.noise.incl)
gc()

# Now sum up the counts for each sRNA based on assignment to union genes (hierachy not taken into account) ######


collapsed.sRNA<-fread(file=file.name,nThread=14)


collapsed.sRNA[ ,c("SEQ","rank","no.mapping.sites","sRNA","hierachy_abundance","frag.length","sRNA_type") := NULL]


collapsed.sRNA_union_gene.noise.incl<-collapsed.sRNA[, lapply(.SD, sum, na.rm=TRUE), by=sites_seq_maps_to, .SDcols=c(1,4:ncol(collapsed.sRNA)) ] 


pass_filter<-rowSums(collapsed.sRNA_union_gene.noise.incl[,2:ncol(collapsed.sRNA_union_gene.noise.incl)]>min.read.freq)>min.sample

collapsed.sRNA_union_gene.noise.incl<-collapsed.sRNA_union_gene.noise.incl[pass_filter,] # remove


file.it.2<-paste0(location.of.patman.files,ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_gene_union_collpased_",output)  #  will need this file name later in the script

fwrite(collapsed.sRNA_union_gene.noise.incl,file=file.it.2) # out put the collapsed file

rm(collapsed.sRNA_union_gene.noise.incl)
rm(collapsed.sRNA)
gc()


### now create the final expression matrix for fractional counts ######



cat("\n multireduced sequences have been collapsed. Now merging the files \n")

files2<-dir(path=location.of.patman.files, pattern="_temp4.txt",full.names = T,recursive=F)  


i<-1
run4<-fread(files2[i],sep=",") # change
colnames(run4)[2]<-files2[i]

for(i in 2:length(files2)){
  run3<-fread(files2[i],sep=",") # change
  colnames(run3)[2]<-files2[i]
  run4<-full_join(run4,run3)
  print(i)
}


run4[is.na(run4)]<-0

file.it.fractional<-paste0(location.of.patman.files,ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_fractional__counts_collpased_",output)  #  will need this file name later in the script

colnames(run4)[2:length(run4)] %<>%
  basename(.)%>%
  gsub("_merged_temp4.txt","",.)

fwrite(run4,file=file.it.fractional) # out put the collapsed file




### now create the final expression matrix for fractional counts collapsed up a level ######


files2<-dir(path=location.of.patman.files, pattern="_temp5.txt",full.names = T,recursive=F)  


i<-1
run4<-fread(files2[i],sep=",") # change
colnames(run4)[2]<-files2[i]

for(i in 2:length(files2)){ 
  run3<-fread(files2[i],sep=",") # change
  colnames(run3)[2]<-files2[i]
  run4<-full_join(run4,run3)
  print(i)
}


run4[is.na(run4)]<-0

file.it.fractional<-paste0(location.of.patman.files,ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_fractional__counts_collpased_up_one_level",output)  #  will need this file name later in the script

colnames(run4)[2:length(run4)] %<>%
  basename(.)%>%
  gsub("_merged_temp5.txt","",.)

fwrite(run4,file=file.it.fractional) # out put the collapsed file




#################

ribodepleted<-"ribodeplete"
  
logger=paste("filter criteria applied @",Sys.time(),":","reads which did not have at least", min.read.freq ,"reads in", min.sample ,"samples were filtered out") # this provides a link so I can refer to what fitering settings I used to create the files.

sink(READMEfile, append = T) # will create the text file if it is not already in existance

writeLines(logger)

sink()



files.to.delete <- dir(location.of.patman.files,pattern="_temp3.txt",recursive=F,full.names=T)
file.remove(files.to.delete) # save space

files.to.delete <- dir(location.of.patman.files,pattern="_temp4.txt",recursive=F,full.names=T)
file.remove(files.to.delete) # save space

rm(list=ls(pattern="totals|test|reduced|run|rv|reads|table|multi|frag.length"))

rm(list=c("single.mapping","Ok.sRNA","freq","mat.miR","pass_filter","remove","ribozomal"))  ## remove the large objects created in this script in case I am piping through a wrapper script.

gc()



