

# collapse to transcipt, collapse to gene??
# collapse to gene ? just do in DE script later on. 

# This program takes .sum files where reads were simultaneouseously mapped to the protein coding and anti-sence transctiptome. It collapse sequences in four ways: A)  it will collapse multimapping sequences to the most abundant feature NOT in a hierachial manner B) it will collapse multimapping sequences to the most abundant feature in a hierachial manner - prioritising protein coding genes C) it will collapse multimapping sequences to a union gene. D) it will collapse using fractional counts.



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

ribodepleted<-"not_ribo_depleted" # kept in for compatibility. Not actually relevent.

###  Create totals based on all reads -- I will use this to decide re assignment ####

cat("Create totals based on all reads -- I will use this to decide re assignment \n")
expt.sum.file<-list.files(location.of.patman.files,pattern = "experimentlevel.sum",full.names = T)
totals<-fread(expt.sum.file)

if(ncol(totals)>3){
  totals<-totals[,1:3] # col 4  contains start site - not relevent for this
}
colnames(totals)<-c("totals","seq","sRNA")


# remove the entries which relate to sequences that map to miR.haripins when they are already mapped to a mature miR, whilst preserving sequences that map to a miR hairpin but not a mature arm. Also remove sequences that map to a mature miRNA but not a miR hairpin (e.g. because a sequence maps to another sRNA type with fewer mismatches so did not assign to a hairpin using a best straata approach with the the sRNAome)


totals2<-data.table(totals)


totals2<-totals2[,list(totals=sum(totals)),by='sRNA'] # add the total counts per sRNA.

setorder(totals2, -totals) # order totals by size

totals.hierachy<-totals2

totals.hierachy$sRNA_type<-NA
totals.hierachy$sRNA_type[grepl("ReverseComplement",totals.hierachy$sRNA)]<-"ReverseComplement"
totals.hierachy$sRNA_type[!grepl("ReverseComplement",totals.hierachy$sRNA)]<-"protein_transcript"

totals.hierachy$sRNA_type <- factor(totals.hierachy$sRNA_type, levels = c("protein_transcript","ReverseComplement"))

totals.hierachy<-arrange(totals.hierachy, sRNA_type,-totals) # orders the table by sRNA type and then abundance so if a sequence maps to two sRNA types, the prioritised sRNA_type will appear first.

remove_seqs_that_map_this_many_times<-"no_limit"

# ################
#
## prep for filtering step
min.read.freq=1 # kept for script compatibility
min.sample=1 # kept for script compatibility

must.have.at.least<-(min.read.freq)*(min.sample+1)


totals3<-totals2[totals2$totals>=must.have.at.least,]

saveRDS(totals3,"totals3.RDS")




#totals3<-readRDS(paste0(location.of.patman.files,"totals3.RDS"))

totals3$rank<-frank(totals3$totals,ties.method="dense")


#### collapse expt.sum file #######


colnames(totals)<-c("count","SEQ","sRNA")



test2<-totals # a backup
sRNA_test<-test2$sRNA
sRNA_test[grepl("ReverseComplement",sRNA_test,perl=T)]<-"ReverseComplement"
sRNA_test[!grepl("ReverseComplement",sRNA_test,perl=T)]<-"protein_transcript"
test2[ , `:=` (sRNA_type=sRNA_test)]

test2[, `:=` (no.mapping.sites = .N), by = SEQ] # add a column with no. sites each sequencee maps to. Bit quicker than aggregate or table. 

test2[,`:=` (fractional.count=count/no.mapping.sites)] # get the fractional counts - I need this for colaping to master


Ok.sRNA<-test2[no.mapping.sites==1]
Ok.sRNA %<>% data.frame(.)
Ok.sRNA$rank<-NA
Ok.sRNA$no.mapping.sites<-1
Ok.sRNA$sites_seq_maps_to<-Ok.sRNA$sRNA

reduced.final<-test2[no.mapping.sites>1] %>% data.frame




### create master fractional table #####
test2[, c( "count", "SEQ","no.mapping.sites", "sRNA_type") := NULL]


collapsed.sRNA_fractional.noise.incl.master<-test2[,list(Expt.sum.count=sum(fractional.count)),by='sRNA']

collapsed.sRNA_fractional.noise.incl.master %<>% data.frame(.)


### fractional gene
test2<-totals # 

## collapse by gene name

protein_loci<-grepl("ENSG",test2$sRNA,perl=T)

protein_gene_names<-test2$sRNA[protein_loci] # make a dataframe of protein coding gene names. The gene name will be the second bit
cat("extracting gene names from transcript")
protein_gene_names_ENSG_bit<-sub('^([^\\|]+\\|[^\\|]+).*', '\\1', protein_gene_names) # takes a while - about 10 mins, but does not take up much memory


rev.comp<-grepl("ReverseComplement",protein_gene_names_ENSG_bit,perl=T) # this is everything that was reverse complement

# now remove everything before the first |
protein_gene_names_ENSG_bit2<-gsub(".*\\|","",protein_gene_names_ENSG_bit,perl=T) # quit fast considering length of vector
protein_gene_names_ENSG_bit2[rev.comp] %<>%paste0("ReverseComplement",.)

test2$gene.name<-test2$sRNA
test2$gene.name[protein_loci]<-protein_gene_names_ENSG_bit2

gene_transcript.link<-cbind(test2$sRNA,test2$gene.name) ### use this to collapse sample files to genes later on.
colnames(gene_transcript.link)<-c("sRNA","gene.name")

gene_transcript.link %<>% data.frame(.)

dups<-duplicated(paste(test2$SEQ,test2$gene.name)) # takes  few mins
test2 %<>% data.frame(test2)     
test2 %<>% .[!dups,]

test2 %<>% data.table(.)
#remove duplicates


test2[, `:=` (no.mapping.sites = .N), by = SEQ] # add a column with no. sites each sequencee maps to. Bit quicker than aggregate or table. 

test2[,`:=` (fractional.count=count/no.mapping.sites)] # get the fractional counts - I need this for colaping to master

GENE.collapsed.sRNA_fractional.noise.incl.master<-test2[,list(Expt.sum.count=sum(fractional.count)),by='gene.name'] # this the fractional result collapsed by gene.

GENE.collapsed.sRNA_fractional.noise.incl.master %<>% data.frame(.)

test2[,sites_seq_maps_to := paste0(gene.name,collapse="&"),by=SEQ] # create gene.union. 

test2<-test2[!duplicated(test2$SEQ)] # remove duplicated entries per sequence

# create gene union master
GENE.collapsed.sRNA_gene.union.noise.incl.master<-test2[,list(Expt.sum.count=sum(count)),by='sites_seq_maps_to'] # this the gene union result collapsed by gene.
GENE.collapsed.sRNA_gene.union.noise.incl.master %<>% data.frame(.)




#######  Assign the sRNA based on the sRNA abundances ######



# order reduced.final by window abundance.

reduced.final$rank<-totals3$rank[match(reduced.final$sRNA,totals3$sRNA)]

reduced.final<-reduced.final[order(reduced.final$rank,decreasing=T),]

# remove sequences mapping > specific no of times times




if(sum(grepl("no_limit",remove_seqs_that_map_this_many_times))==0)
{
  reduced.final %<>% data.table(.)
  reduced.final<-reduced.final[no.mapping.sites>=remove_seqs_that_map_this_many_times]
  reduced.final %<>% data.frame(.)
  
}


multi.reduced <-reduced.final[!duplicated(reduced.final$SEQ),] # sequences collapsed based on abdundances, as reduced.final has been ordered by sRNA abdundance so the first occurance of seq is the sRNA with highest abdundance (unless there are sRNAs with equal abudance in which case one is arbitrarily kept)

reduced.final <- data.frame(reduced.final)
reduced.final$sRNA_type <- factor(reduced.final$sRNA_type, levels = c("protein_transcript","ReverseComplement"))


totals.hierachy$hierachy_order<-1:nrow(totals.hierachy)




############  reduced.final<-reduced.final[!reduced.final$sRNA_type=="rRNA",]



reduced.final$hierachy_order<-totals.hierachy$hierachy_order[match(reduced.final$sRNA,totals.hierachy$sRNA)] 


# I have now appended the order (i.e. rank of the sRNAs) onto the reduced.final data frame.


reduced.final.hierachy<-reduced.final[order(reduced.final$hierachy_order),] # now only keep the first sequence as ordered by hierachy

multi.reduced.hierachy<-reduced.final.hierachy[!duplicated(reduced.final.hierachy$SEQ),] # keep only the first sequence record: this is the record which is best in terms of sRNA_type and then in terms of abundance. It allows sequences to be mapped exactly the same as I originally collapsed sequences


collapsed.sRNA.hierachy<-bind_rows(Ok.sRNA,multi.reduced.hierachy) # combine the single mapping sequences and the now reduced multimapping sequences


collapsed.sRNA.hierachy<-collapsed.sRNA.hierachy[,c("count","SEQ","sRNA","no.mapping.sites","sRNA_type")] # get rid of non-useful columns

rm(totals2)
rm(totals3)
rm(totals)
gc()


### also try collapsing to union gene
cat("\n collapsing multimapping sequences - may take some time... \n")


# create gene union names
library(dplyr)


reduced.final<-data.table(reduced.final)
reduced.final<-reduced.final[,sites_seq_maps_to := paste0(sRNA,collapse="&"),by=SEQ]


####
multi.reduced$sites_seq_maps_to<-reduced.final$sites_seq_maps_to[match(multi.reduced$SEQ,reduced.final$SEQ)] #  quickish
####

rm(reduced.final)

collapsed.sRNA<-rbind(Ok.sRNA,multi.reduced,fill=T) # combine the single mapping sequences and the now reduced multimapping sequences (assigned by abundance)

rm(multi.reduced)
gc()
# append the sRNA decision that was made when collapsing sequences hierachically and then by abundance.

collapsed.sRNA$hierachy_abundance<-collapsed.sRNA.hierachy$sRNA[match(collapsed.sRNA$SEQ,collapsed.sRNA.hierachy$SEQ)]# append the sRNA decision that was made when collapsing sequences hierachically and then by abundance.


# collapsed.sRNA  =   data frame which holds each sequence and says which sRNA it ended up mapping to under each condition.

collapsed.sRNA %<>% data.table(.)
collapsed.sRNA %<>% .[,fractional.count:=NULL] # remove the fractional count. It was made earlier for a longer table but doesn't mean anythign in this context. Fractional counts are only relevent when all sRNA mappings are kept. 
file.name<-paste0(location.of.patman.files,"expt.sum.temp3.txt")

colnames(collapsed.sRNA)[1]<-"Expt.sum.count"

collapsed.sRNA[,frag.length:=stri_length(SEQ),]


fwrite(collapsed.sRNA,file.name,col.names = T,row.names = F) # very quick and possibly will save memory



cat("\n decisions made about sRNA assignment for each sequence in the experimental.sum file complete \n")




#collapsed.sRNA<-fread(paste0(location.of.patman.files,"expt.sum.temp3.txt"))
# contains each SEQ and how it collapses according to each method
# can remove once script done as collapse.sRNA already exists. Its used here for debugging when I need to run script after a restart. 




collapsed.sRNA %<>% data.frame(.)

#### Now create a master collapsed file which I can append my sample files onto after collapsing them ######


run4<-collapsed.sRNA
run4$new.file.count<-run4$Expt.sum.count # super quick

run4$new.file.count[is.na(run4$new.file.count)]<-0
run4.save<-run4
# 
# run4<-run4.save

### Abdundance collapsed

# sum counts for each sRNA based on assignment to the most abundant read

run4<-run4[,!colnames(run4) %in% c("SEQ","rank","no.mapping.sites","sites_seq_maps_to","hierachy_abundance","frag.length","sRNA_type","Expt.sum.count")]

run4 %<>% data.table(.)
collapsed.sRNA_majority_collapsed.noise.incl.master<-run4[,list(A=sum(new.file.count)),by='sRNA']



colnames(collapsed.sRNA_majority_collapsed.noise.incl.master)<-c("sRNA","Expt.sum.count")

# Now sum up the counts for each sRNA based on assignment to hierachical abdunant genes


run4<-run4.save

run4<-run4[ ,!colnames(run4)%in%c("SEQ","rank","no.mapping.sites","sites_seq_maps_to","sRNA","frag.length","sRNA_type","Expt.sum.count")] # keep only counts and hierachial abundant assignment

run4 %<>% data.table(.)
collapsed.sRNA_hierachically_collpased.noise.incl.master<-run4[,list(A=sum(new.file.count)),by='hierachy_abundance']


colnames(collapsed.sRNA_hierachically_collpased.noise.incl.master)<-c("hierachy_abundance","Expt.sum.count")



# Now sum up the counts for each sRNA based on assignment to union genes (hierachy not taken into account) ######


run4<-run4.save

run4<-run4[ ,!colnames(run4)%in%c("Expt.sum.count","SEQ","rank","no.mapping.sites","sRNA","hierachy_abundance","frag.length","sRNA_type")]


run4 %<>% data.table(.)
collapsed.sRNA_union_gene.noise.incl.master<-run4[,list(sites_seq_maps_to=sum(new.file.count)),by='sites_seq_maps_to']


colnames(collapsed.sRNA_union_gene.noise.incl.master)<-c("sites_seq_maps_to","Expt.sum.count")

# sum up the counts for each sRNA based on fractioinal counts######

# collapsed.sRNA_fractional.noise.incl.master - already summed up earlier - has to be made early when every sequence and mapping is present ,in the dataframe.




cat("Now read in each sample file and auto collapse by appending to collapsed master dataframes \n")

# now read in each sample file and append to collapsed master dataframes #########


files2<-dir(path=location.of.patman.files, pattern=".sum",full.names = T,recursive=F)
files2 %<>% .[!grepl("expt|experimen",.)]

for(i in 1:length(files2)){

  run3<-fread(files2[i],sep=" ") # change
  
  # match the SEQ in the imported file to the correct decision - the first instance of SEQ will match the others will be removed. 
  
  run4<-collapsed.sRNA
  run4$new.file.count<-run3$V1[match(collapsed.sRNA$SEQ,run3$V2)] # super quick
  
  run4$new.file.count[is.na(run4$new.file.count)]<-0
  run4.save<-run4
  
  run4<-run4.save
  
  run4<-run4[run4$new.file.count>0,]
  
  ### Abdundance collapsed
  
  # sum counts for each sRNA based on assignment to the most abundant read
  
  run4<-run4[,!colnames(run4) %in% c("SEQ","rank","no.mapping.sites","sites_seq_maps_to","hierachy_abundance","frag.length","sRNA_type","Expt.sum.count")]
  
  
  run4 %<>% data.table(.)
  collapsed.sRNA_majority_collapsed.noise.incl.temp<-run4[,list(A=sum(new.file.count)),by='sRNA']
  
  
  colnames(collapsed.sRNA_majority_collapsed.noise.incl.temp)<-c("sRNA","new.file.count")
  
  # Now sum up the counts for each sRNA based on assignment to hierachical abdunant genes
  
  
  run4<-run4.save
  run4<-run4[run4$new.file.count>0,]
  
  run4<-run4[ ,!colnames(run4)%in%c("SEQ","rank","no.mapping.sites","sites_seq_maps_to","sRNA","frag.length","sRNA_type","Expt.sum.count")] # keep only counts and hierachial abundant assignment
  
  run4 %<>% data.table(.)
  collapsed.sRNA_hierachically_collpased.noise.incl.temp<-run4[,list(A=sum(new.file.count)),by='hierachy_abundance']
  
  colnames(collapsed.sRNA_hierachically_collpased.noise.incl.temp)<-c("hierachy_abundance","new.file.count")
  
  # Now sum up the counts for each sRNA based on assignment to union genes (hierachy not taken into account) ######
  
  # need to read in union gene assignments
  
  run4<-run4.save
  run4<-run4[run4$new.file.count>0,]
  
  run4<-run4[ ,!colnames(run4)%in%c("Expt.sum.count","SEQ","rank","no.mapping.sites","sRNA","hierachy_abundance","frag.length","sRNA_type")]
  
  run4 %<>% data.table(.)
  collapsed.sRNA_union_gene.noise.incl.temp<-run4[,list(A=sum(new.file.count)),by='sites_seq_maps_to']
  
  colnames(collapsed.sRNA_union_gene.noise.incl.temp)<-c("sites_seq_maps_to","new.file.count")
  
  
  
  ###  fractional
  
  
  run3$no.mapping.sites<-collapsed.sRNA$no.mapping.sites[match(run3$V2,collapsed.sRNA$SEQ)]
  run3$fractional.count<-run3$V1/run3$no.mapping.sites
  # total this up
  
  run3 %<>% data.table(.)
  collapsed.sRNA_fractional.noise.incl.temp<-run3[,list(A=sum(V1)),by='fractional.count'] # collapsed fractional counts.
  
  
  ##### Now append collapsed files to master using match ######
  
  #append abdundance to master
  
  collapsed.sRNA_majority_collapsed.noise.incl.master$new.file<-collapsed.sRNA_majority_collapsed.noise.incl.temp$new.file.count[match(collapsed.sRNA_majority_collapsed.noise.incl.master$sRNA,  collapsed.sRNA_majority_collapsed.noise.incl.temp$sRNA)]
  
  colnames(collapsed.sRNA_majority_collapsed.noise.incl.master)[ncol(collapsed.sRNA_majority_collapsed.noise.incl.master)]<-files2[i] %>% basename(.)%>%
    gsub("_merged_\\.protein_transcriptome_plus_whole_transcriptome_reverse_complement.alligned.sum","",.)
  
  
  
  #append hierachy to master
  
  collapsed.sRNA_hierachically_collpased.noise.incl.master$new.file<-collapsed.sRNA_hierachically_collpased.noise.incl.temp$new.file.count[match(collapsed.sRNA_hierachically_collpased.noise.incl.master$hierachy_abundance,  collapsed.sRNA_hierachically_collpased.noise.incl.temp$hierachy_abundance)]
  
  colnames(collapsed.sRNA_hierachically_collpased.noise.incl.master)[ncol(collapsed.sRNA_hierachically_collpased.noise.incl.master)]<-files2[i] %>% basename(.)%>%
    gsub("_merged_\\.protein_transcriptome_plus_whole_transcriptome_reverse_complement.alligned.sum","",.)
  
  
  #append gene union to master
  
  collapsed.sRNA_union_gene.noise.incl.master$new.file<-collapsed.sRNA_union_gene.noise.incl.temp$new.file.count[match(collapsed.sRNA_union_gene.noise.incl.master$sites_seq_maps_to, collapsed.sRNA_union_gene.noise.incl.temp$sites_seq_maps_to)]
  
  colnames(collapsed.sRNA_union_gene.noise.incl.master)[ncol(collapsed.sRNA_union_gene.noise.incl.master)]<-files2[i] %>% basename(.)%>%
    gsub("_merged_\\.protein_transcriptome_plus_whole_transcriptome_reverse_complement.alligned.sum","",.)
  
  
  #append fractional to master
  
  
  collapsed.sRNA_fractional.noise.incl.master$new.file<-collapsed.sRNA_fractional.noise.incl.temp$new.file.count[match(collapsed.sRNA_fractional.noise.incl.master$sRNA, collapsed.sRNA_fractional.noise.incl.temp$sRNA)]
  
  colnames(collapsed.sRNA_fractional.noise.incl.master)[ncol(collapsed.sRNA_fractional.noise.incl.master)]<-files2[i] %>% basename(.)%>%
    gsub("_merged_\\.protein_transcriptome_plus_whole_transcriptome_reverse_complement.alligned.sum","",.)
  
  
  #### now collapse to gene union and fractional @ GENE level
  
  
  ###  gene fractional
  
  run3_gene<-run3
  
  colnames(run3_gene)[1:4]<-c("count","SEQ","sRNA","start.site")
  
  # convert transcripts to gene names
  
  run3_gene$gene.name<-gene_transcript.link$gene.name[match(run3_gene$sRNA,gene_transcript.link$sRNA)] # get the gene name 
  
  run3_gene$no.mapping.sites<-test2$no.mapping.sites[match(run3_gene$SEQ,test2$SEQ)] # add the no.mapping sites
  
  dups<-duplicated(paste(run3_gene$SEQ,run3_gene$gene.name)) # remove duplicated SEQ:gene names (as genes had multiple transcripts mapping to them)
  run3_gene %<>% data.frame(run3_gene)     
  run3_gene %<>% .[!dups,]
  
  run3_gene$fractional.count<-run3_gene$count/run3_gene$no.mapping.sites
  # total this up
  
  run3_gene %<>% data.table(.)
  
  GENE.collapsed.sRNA_fractional.noise.incl.temp<-run3_gene[,list(A=sum(fractional.count)),by='gene.name'] # collapsed fractional counts.
  
  ### now collapse to GENE gene union
  
  
  run3_gene<-run3
  
  colnames(run3_gene)[1:4]<-c("count","SEQ","sRNA","start.site")
  
  # get gene_union name for the sequence
  
  
  run3_gene$sites_seq_maps_to<-test2$sites_seq_maps_to[match(run3_gene$SEQ,test2$SEQ)] # get the gene union name 
  
  dups<-duplicated(run3_gene$SEQ) # remove duplicated entries per sequence
  run3_gene %<>% data.frame(run3_gene)     
  run3_gene %<>% .[!dups,]
  run3_gene %<>% data.table(.)
  
  # create gene union master
  GENE.collapsed.sRNA_gene.union.noise.incl.temp<-run3_gene[,list(A=sum(count)),by='sites_seq_maps_to'] # this the gene union result collapsed by gene.
  
  #### now append gene union and fractional @ GENE level
  
  
  ## collapse by gene GENE.collapsed.sRNA_fractional.noise.incl.masternion
  
  # fractional
  
  GENE.collapsed.sRNA_fractional.noise.incl.master$new.file<-GENE.collapsed.sRNA_fractional.noise.incl.temp$new.file.count[match(GENE.collapsed.sRNA_fractional.noise.incl.master$sRNA, GENE.collapsed.sRNA_fractional.noise.incl.temp$sRNA)]
  
  colnames(GENE.collapsed.sRNA_fractional.noise.incl.master)[ncol(GENE.collapsed.sRNA_fractional.noise.incl.master)]<-files2[i] %>% basename(.)%>%
    gsub("_merged_\\.protein_transcriptome_plus_whole_transcriptome_reverse_complement.alligned.sum","",.)
  
  # gene union 
  GENE.collapsed.sRNA_gene.union.noise.incl.master$new.file<-GENE.collapsed.sRNA_gene.union.noise.incl.temp$new.file.count[match(GENE.collapsed.sRNA_gene.union.noise.incl.master$sRNA, GENE.collapsed.sRNA_gene.union.noise.incl.temp$sRNA)]
  
  colnames(GENE.collapsed.sRNA_gene.union.noise.incl.master)[ncol(GENE.collapsed.sRNA_gene.union.noise.incl.master)]<-files2[i] %>% basename(.)%>%
    gsub("_merged_\\.protein_transcriptome_plus_whole_transcriptome_reverse_complement.alligned.sum","",.)
  
  
  
  cat(paste0("file",i," has been collapsed and appended","\n"))
}

print("sample files have been collapsed and appended")


rm(test2)
rm(run4.save)
rm(run4)
rm(run3)
rm(run3_gene)
rm(totals.hierachy)
rm(collapsed.sRNA)

# OK so now all the temp files have been made and aggregated by each collapse type.


# Save the completed master files.
# Each end product can the be aggregated again in the traditional manner using gene name. 

# I need to make a protein level version of expt.4 so that I can calculated fractional at protein level (otherwise no. transcripts will effect site. seqs. map to)


cat("\n Now writing to file \n")

######

# file names

# collapsed.sRNA_majority_collapsed.noise.incl.master

file.it.1<-paste0(location.of.patman.files,ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_collpased_",output)  #  will need this file name later in the script

collapsed.sRNA_majority_collapsed.noise.incl.master[is.na(collapsed.sRNA_majority_collapsed.noise.incl.master)]<-0

fwrite(collapsed.sRNA_majority_collapsed.noise.incl.master,file=file.it.1) # out put the collapsed file


### Hierachy

file.it.1<-paste0(location.of.patman.files,ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_hierachically_collpased_",output)  #  will need this file name later in the script

collapsed.sRNA_hierachically_collpased.noise.incl.master[is.na(collapsed.sRNA_hierachically_collpased.noise.incl.master)]<-0


fwrite(collapsed.sRNA_hierachically_collpased.noise.incl.master,file=file.it.1) # out put the collapsed file

# gene union level

file.it.1<-paste0(location.of.patman.files,ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_gene_union_collpased_",output)  #  will need this file name later in the script


collapsed.sRNA_union_gene.noise.incl.master[is.na(collapsed.sRNA_union_gene.noise.incl.master)]<-0


fwrite(collapsed.sRNA_union_gene.noise.incl.master,file.it.1)

rm(collapsed.sRNA_union_gene.noise.incl.master)
# fractional 

file.it.1<-paste0(location.of.patman.files,ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_fractional_collpased_",output)  #  will need this file name later in the script

collapsed.sRNA_fractional.noise.incl.master[is.na(collapsed.sRNA_fractional.noise.incl.master)]<-0


fwrite(collapsed.sRNA_fractional.noise.incl.master,file.it.1)

rm(collapsed.sRNA_fractional.noise.incl.master)

### now collapse each to gene level

# *AFTER* writing the collapsed files to disc, then further consolidate by gene

##### Further gene name collapse/consolidation for hierachy and abdundance (already done for fractional and gene_union)


## collapse by gene name
# extract out gene name  - only relevent for protein coding genes (incl reverse complement ones) because I don't have transcript | gene levels for any other type of data.
# ID proteins

protein_loci<-grepl("ENSG",collapsed.sRNA_majority_collapsed.noise.incl.master$sRNA)

protein_gene_names<-collapsed.sRNA_majority_collapsed.noise.incl.master$sRNA[protein_loci] %>%
  data.frame(.) # make a dataframe of protein coding gene names. The gene name will be the second bit

protein_gene_names_ENSG_bit<-reshape2::colsplit(protein_gene_names$.,pattern = "\\|",names=1:100)

protein_gene_names_ENSG_bit[grepl("ReverseComplement",protein_gene_names_ENSG_bit[,1]),2]<-protein_gene_names_ENSG_bit[grepl("ReverseComplement",protein_gene_names_ENSG_bit[,1]),2]%>%
  paste0("ReverseComplement_",.) # these are now the gene names as I want them. They need reinserting. 

collapsed.sRNA_majority_collapsed.noise.incl.master$sRNA[protein_loci]<-protein_gene_names_ENSG_bit[,2]

collapsed.sRNA_majority_collapsed.noise.incl.master %<>% data.table(.)
collapsed.sRNA_majority_collapsed.noise.incl.master<-collapsed.sRNA_majority_collapsed.noise.incl.master[,lapply(.SD,sum),by=sRNA]

#####

## collapse by gene name

protein_loci<-grepl("ENSG",collapsed.sRNA_hierachically_collpased.noise.incl.master$hierachy_abundance)

protein_gene_names<-collapsed.sRNA_hierachically_collpased.noise.incl.master$hierachy_abundance[protein_loci] %>%
  data.frame(.) # make a dataframe of protein coding gene names. The gene name will be the second bit

protein_gene_names_ENSG_bit<-reshape2::colsplit(protein_gene_names$.,pattern = "\\|",names=1:100)
protein_gene_names_ENSG_bit[grepl("ReverseComplement",protein_gene_names_ENSG_bit[,1]),2]<-protein_gene_names_ENSG_bit[grepl("ReverseComplement",protein_gene_names_ENSG_bit[,1]),2]%>%
  paste0("ReverseComplement_",.) # these are now the gene names as I want them. They need reinserting. 

collapsed.sRNA_hierachically_collpased.noise.incl.master$hierachy_abundance[protein_loci]<-protein_gene_names_ENSG_bit[,2]

collapsed.sRNA_hierachically_collpased.noise.incl.master %<>% data.table(.)

collapsed.sRNA_hierachically_collpased.noise.incl.master<-collapsed.sRNA_hierachically_collpased.noise.incl.master[,lapply(.SD,sum),by=hierachy_abundance]



## now write these files


dir.create(paste0(location.of.patman.files,"/gene_level"))

# collapsed.sRNA_majority_collapsed.noise.incl.master

file.it.1<-paste0(location.of.patman.files,"gene_level/",ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_collpased_",output)  #  will need this file name later in the script

collapsed.sRNA_majority_collapsed.noise.incl.master[is.na(collapsed.sRNA_majority_collapsed.noise.incl.master)]<-0

fwrite(collapsed.sRNA_majority_collapsed.noise.incl.master,file=file.it.1) # out put the collapsed file
rm(collapsed.sRNA_majority_collapsed.noise.incl.master)

### Hierachy

file.it.1<-paste0(location.of.patman.files,"gene_level/",ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_hierachically_collpased_",output)  #  will need this file name later in the script

collapsed.sRNA_hierachically_collpased.noise.incl.master[is.na(collapsed.sRNA_hierachically_collpased.noise.incl.master)]<-0


fwrite(collapsed.sRNA_hierachically_collpased.noise.incl.master,file=file.it.1) # out put the collapsed file
rm(collapsed.sRNA_hierachically_collpased.noise.incl.master)
# gene union level

file.it.1<-paste0(location.of.patman.files,"gene_level/",ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_gene_union_collpased_",output)  #  will need this file name later in the script


GENE.collapsed.sRNA_gene.union.noise.incl.master[is.na(GENE.collapsed.sRNA_gene.union.noise.incl.master)]<-0

fwrite(GENE.collapsed.sRNA_gene.union.noise.incl.master,file.it.1)
rm(GENE.collapsed.sRNA_gene.union.noise.incl.master)

######


# fractional 

file.it.1<-paste0(location.of.patman.files,"gene_level/",ribodepleted,"MM_threshold_of_",remove_seqs_that_map_this_many_times,"_",sRNA_level.of.collapse,"_fractional_collpased_",output)  #  will need this file name later in the script

GENE.collapsed.sRNA_fractional.noise.incl.master[is.na(GENE.collapsed.sRNA_fractional.noise.incl.master)]<-0

fwrite(GENE.collapsed.sRNA_fractional.noise.incl.master,file.it.1)


## maybe add this bit back in somewhere pass_filter<-rowSums(collapsed.sRNA_union_gene.noise.incl3[,2:ncol(collapsed.sRNA_union_gene.noise.incl3)]>min.read.freq)>min.sample
# collapsed.sRNA_union_gene.noise.incl3<-collapsed.sRNA_union_gene.noise.incl3[pass_filter,] # remove
## keep collapsed by transcript but filter out genes with low abundance - makes file sizes easier

# pass_filter<-rowSums(collapsed.sRNA_union_gene.noise.incl[,2:ncol(collapsed.sRNA_union_gene.noise.incl)]>min.read.freq)>min.sample
# 
# collapsed.sRNA_union_gene.noise.incl<-collapsed.sRNA_union_gene.noise.incl[pass_filter,] # remove



