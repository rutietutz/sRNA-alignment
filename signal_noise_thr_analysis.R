###############################
##graphical representation of the p2pPCC
###############################

file.pdf = paste(prefix,sRNA_level.of.collapse_for_noise, "detailed_assessment_noise_0.5_log2_bins.pdf", sep="_") 
abn.in   = paste0(getwd(),"/",abn_output)  
pcc.in   = paste0(getwd(),"/",pcc_output)

abn.matrix <- read.csv(abn.in,sep=',',header=FALSE)
pcc.matrix <- read.csv(pcc.in,sep=',',header=FALSE)
n.samples<-as.numeric(n.samps)
sample.IDs<-read.csv(index.file)



pdf(file.pdf,width=15,height=5)
par(mfrow=c(1,3))
for(i in 1:n.samples)
{
  sample.name<-sample.IDs$sample.name[i]
  abn.bins2 <- floor(2*log2(abn.matrix[,i]+1))
  boxplot(pcc.matrix[,i]~abn.bins2,main=paste("sample",sample.name), 
          ylim=c(-1,1),xlab='Abundance', ylab='p2p PCC',xaxt='n')
  abline(v=10,col='red'); abline(v=12,col='blue')
  abline(h=0.7,col='red');abline(h=0.6,col='blue')
  xlab <- (1:40)/2
  axis(1,at=1:40,labels=xlab,las=2)
}
dev.off()

##################################
##determine the noise threshold per sample
##################################
library(MASS)
#abn.matrix <- read.csv(abn.in,sep=',',header=FALSE)
#pcc.matrix <- read.csv(pcc.in,sep=',',header=FALSE)
noise.out  = paste(location.of.patman.files,prefix,"_",sRNA_level.of.collapse,"_noiseMatrix.csv",sep='')
#n.samples=44

iqr.thr = 0.6
noise.thr <- matrix(rep(-1,n.samples*2),ncol=2)
for(i in 1:n.samples)
{
  abn.bins <- floor(log2(abn.matrix[,i]+1)) ## +1 increment on log2 scale
  #abn.bins2<- floor(2*log2(abn.matrix[,i]+1)) ## +0.5 increment on log2 scale
  
  pcc.summary <- tapply(pcc.matrix[,i],as.factor(abn.bins),FUN=summary)
  #pcc.summary2<- tapply(pcc.matrix[,i],as.factor(abn.bins2),FUN=summary)
  
  for(j in 2:length(pcc.summary))
  {
    ##median above the thr
    if(pcc.summary[[j]][3] > iqr.thr)
    {
      ##we use j-1 because of the 0 increment count for the summary
      noise.thr[i,1] = j-1
      break
    }
  }
  for(j in seq(length(pcc.summary),2,-1))
  {
    ##median above the thr
    if(pcc.summary[[j-1]][3] < iqr.thr)
    {
      noise.thr[i,2] = j-1
      break
    }
  }
}
colnames(noise.thr)=c("Lower","Higher")
write.matrix(t(noise.thr),file=noise.out,sep=',')
