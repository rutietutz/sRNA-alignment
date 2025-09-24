


#########
## get correct file path starts


# ID.link_corrected<-read.csv(paste0(computer,"/Dropbox/My MIR work/My MIR work/euclids sequencing/euclids.ID.link_corrected.csv"))  ### HDCD
# 
# 
# 
# 
# 
# colour<-c(ID.link_corrected$Vaccine.status[1],3,ID.link_corrected$Vaccine.status[2:41],3,ID.link_corrected$Vaccine.status[42])
# 
# colour<-as.factor(as.numeric(colour)) ## changed


####  create expression profiles  ####

input.file  = expression.profile.file_name
presence.plots.output.file = paste(prefix,sRNA_level.of.collapse_for_noise, "presence_plots.pdf", sep="_")   ## changed

inp <- read.csv(input.file,header=FALSE)
transcripts <- levels(as.factor(inp[,2]))


##determines the maximum number of samples present in the file
inp[is.na(inp[,1]),1] = 8 # 
max.samples <- max(inp[,1])
##3 samples per plot
subplot.count = nsamp #floor(max.samples/3)+1 changed

pdf(presence.plots.output.file, width=10, height= 5)
for(t in transcripts)
{
#  if(nchar(as.character(t)) > 10) #############
  {
    #par(mfrow=c(subplot.count,1))
    selected <- inp[inp[,2]==t,]
    cat(paste("working with",t,","))
    
    if(nrow(selected) == max.samples+1)
    {
      xx = seq(1,selected[1,3])
      expr.profile = matrix(rep(0,nrow(selected)*selected[1,3]),nrow=nrow(selected))
      for(j in 1:nrow(selected))
      {
        expr.prof = as.numeric(strsplit(as.character(selected[j,4]), " ")[[1]])
        len.offset= length(expr.prof) - selected[1,3]
        for(k in 1:selected[1,3])
        {
          expr.profile[j,k]=expr.prof[k + len.offset]
        }
      }
      ##finished casting the expression profiles into the expected format
      
      maxY = max(expr.profile)
      ##for(j in 1:nrow(expr.profile))
     # for(j in 1:subplot.count) # for each sample j..
      #{
        for(k in 1:subplot.count)  ## changed
        {
       #   if((j-1)*3+k <= max.samples)
        #  { 
          plot(x=xx,y=expr.profile[k,],ylim=c(0,maxY),
                        col=k,
                        lty=k,
                        xlab='location on transcript (nt)', ylab='expression (linear)', main=t,type='l');
                   par(new=TRUE);
            par(new=TRUE)
         # }
        }
        par(new=FALSE)
      }##endfor plotting
    }
#  }##endif testing the transcript name

}

dev.off()



# #######  a script for looking at dodgy miRs ####
# 
# 
# pcc<-read.csv(pcc_output,header=F)
# abn<-read.csv(abn_output,header=F)
# 
# 
# pcc.mean<-apply(pcc,1,mean)
# abn.mean<-apply(abn,1,median)
# 
# of.interest<-((rowSums(pcc<0.5&pcc>0.1))==2)&((rowSums(abn>10))<2)
# 
# #&((rowSums(abn>20))>5)
# 
# summary(of.interest)
# 
# rv<-(1:nrow(abn))[of.interest]
# 
# rv
# 
# ### plot
# 
#   please.plot<-151
#   
#   for(i in 1:length(rv)){
#   t=transcripts[rv[i]] 
# 
#   #par(mfrow=c(subplot.count,1))
#     
#   
#   selected <- inp[inp[,2]==t,]
#   cat(paste("working with",t,"\n"))
#   
#   xx = seq(1,selected[1,3])
#   expr.profile = matrix(rep(0,nrow(selected)*selected[1,3]),nrow=nrow(selected))
#   for(j in 1:nrow(selected))
#   {
#     expr.prof = as.numeric(strsplit(as.character(selected[j,4]), " ")[[1]])
#     len.offset= length(expr.prof) - selected[1,3]
#     for(k in 1:selected[1,3])
#     {
#       expr.profile[j,k]=expr.prof[k + len.offset]
#     }
#   }
#   ##finished casting the expression profiles into the expected format
#   
#   maxY = max(expr.profile)
#   ##for(j in 1:nrow(expr.profile))
#   # for(j in 1:subplot.count) # for each sample j..
#   #{
#   for(k in 1:subplot.count)  ## changed
#   {
#     #   if((j-1)*3+k <= max.samples)
#     #  { 
#     plot(x=xx,y=expr.profile[k,],ylim=c(0,maxY),
#          col=k,
#          lty=as.numeric(as.character(colour[k])),
#          xlab='location on transcript (nt)', ylab='expression (linear)', main=t,type='l');
#     par(new=TRUE)
#     # }
#   }
#   par(new=F)
#   }
#   
# #  hsa-miR-6889 good example of a possible dodgy miR. Found with: of.interest<-((rowSums(pcc<0.5&pcc>0.1))==2)&((rowSums(abn>10))<2)
# 
# which(transcripts==" hsa-mir-6889")
# pcc[1421,]
#   
