library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)

library(data.table)

library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(TFBSTools)
library(JASPAR2020)
library(Biostrings)

# ----------------------------- #

predictROC<-function(xpos,xneg,file.name){
  library(randomForest)
  library(cvTools)
  library(ROCR)
  
  nfold=5
  n1=nrow(xpos)
  n2=nrow(xneg)
  folds1=cvFolds(n1, K = nfold)
  folds2=cvFolds(n2, K = nfold)
  
  x1=xpos
  x2=xneg
  
  ifold=1
  
  id.train1=folds1$subsets[folds1$which!=ifold]
  id.train2=folds2$subsets[folds2$which!=ifold]
  id.test1=folds1$subsets[folds1$which==ifold]
  id.test2=folds2$subsets[folds2$which==ifold]
  
  ytrain=c(rep(1,length(id.train1)),rep(0,length(id.train2)))
  ytest=c(rep(1,length(id.test1)),rep(0,length(id.test2)))
  
  mtrain=rbind(x1[id.train1,],x2[id.train2,])
  mtest=rbind(x1[id.test1,],x2[id.test2,])
  
  #colnames(mtrain)=colnames(mtest)=1:ncol(mtrain)
    
  rf <- randomForest(
    y=factor(ytrain),x=mtrain,ntree=100
  )
  
  #saveRDS(rf, paste0("/Users/javlon_nizomov/Documents/Research/PHHP Honors Thesis/ML_Project/models/",file.name,".rds"))
  
  yhat=predict(rf,mtest,type="prob")[,2]
  pred <- prediction(yhat, ytest)
  auc <- performance(pred, measure = "auc")
  auc <- auc@y.values[[1]]
  
  message("cor:",round(cor(yhat,ytest),3),' auc: ',round(auc,3))
  
  list(model=rf,yhat=yhat)
}

getProbability <- function(x, model.name, model.type, file.type) {
  library(randomForest)
  
  colnames(x)=c('chr','start','end')
  gr=makeGRangesFromDataFrame(x)
  
  extendReg<-function(gr, ext=500){
    start(gr) = start(gr)-ext
    end(gr) = end(gr)+ext
    gr
  }
  
  gr.ext=extendReg(gr, ext=500)
  
  if(model.type=='motif'){
    x <- matchMotifs(motifs, gr.ext, 
                     genome = BSgenome.Hsapiens.UCSC.hg19)
    x=x@assays@data@listData$motifMatches
    x=as.matrix(x+0)
  }
  
  else if(model.type=='3mer') {
    seq.ext = Biostrings::getSeq(genome, gr.ext)
    x = trinucleotideFrequency(seq.ext)
  }
  
  rf <- readRDS(paste0("models/", model.name, ".rds"))
  
  print("Random forest model loaded")
  
  yhat<-predict(rf, x, type = "prob")[,2]
  
  return(yhat);
}

#motifs <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE)) # human

#save(motifs,file='human.motif.rda')

load('human.motif.rda')

evalPerf <- function(dat, id.pos, id.neg, motifs, feature.type=c('3mer','motif'), name.export, outpath, is.output=F, file.name){
  if(unique(dat$genome)=='hg38'){
    genome=BSgenome.Hsapiens.UCSC.hg38
  }else{
    genome=BSgenome.Hsapiens.UCSC.hg19
  }
  
  tmp=dat[,c('chr','pos','pos')]
  colnames(tmp)=c('chr','start','end')
  gr=makeGRangesFromDataFrame(tmp)
  
  extendReg<-function(gr,ext=500){
    start(gr)=start(gr)-ext
    end(gr)=end(gr)+ext
    gr
  }
  
  gr.ext=extendReg(gr, ext=500)
  
  if(feature.type=='motif'){
    # if motif
    x <- matchMotifs(motifs, gr.ext, 
                     genome = BSgenome.Hsapiens.UCSC.hg19)
    x=x@assays@data@listData$motifMatches
    x=as.matrix(x+0)
  }
  # if 3mer
  else if(feature.type=='3mer'){
    seq.ext=Biostrings::getSeq(genome,gr.ext)
    x=trinucleotideFrequency(seq.ext)
  }
  
  x.pos=x[id.pos,]
  x.neg=x[id.neg,]
  
  out=predictROC(x.pos,x.neg,file.name)
  
  # output sequence
  if (is.output) {
    seq.ext.pos=seq.ext[id.pos]
    seq.ext.neg=seq.ext[id.neg]
    names(seq.ext.pos)=paste0('pos',1:length(seq.ext.pos))
    names(seq.ext.neg)=paste0('neg',1:length(seq.ext.neg))
    writeXStringSet(seq.ext.pos,file=file.path(outpath,paste0('seq.',name.export,'.pos.fasta')))
    writeXStringSet(seq.ext.neg,file=file.path(outpath,paste0('seq.',name.export,'.neg.fasta')))
  }
}







### the below functions are used for data augmentation in deep learning


createRevAndCropSeq<-function(dat,id.pos,id.neg,ext.crop=1000,ntimes=1,name.export,outpath){
  
  #ntimes,for crop seq only, how many times of neg and pos seq cropped
  
  if(unique(dat$genome)=='hg38'){
    genome=BSgenome.Hsapiens.UCSC.hg38
  }else{
    genome=BSgenome.Hsapiens.UCSC.hg19
  }
  
  
  tmp=dat[,c('chr','pos','pos')]
  colnames(tmp)=c('chr','start','end')
  gr=makeGRangesFromDataFrame(tmp)
  
  extendReg<-function(gr,ext=500){
    start(gr)=start(gr)-ext
    end(gr)=end(gr)+ext
    gr
  }
  
  gr.ext=extendReg(gr,ext=500)
  seq.ext=Biostrings::getSeq(genome,gr.ext)
  seq.ext.pos=seq.ext[id.pos]
  seq.ext.neg=seq.ext[id.neg]
  names(seq.ext.pos)=paste0('pos',1:length(seq.ext.pos))
  names(seq.ext.neg)=paste0('neg',1:length(seq.ext.neg))
  writeXStringSet(seq.ext.pos,file=file.path(outpath,paste0('seq.',name.export,'.pos.fasta')))
  writeXStringSet(seq.ext.neg,file=file.path(outpath,paste0('seq.',name.export,'.neg.fasta')))
  
  
  
  #obtain reverse complement
  seq.rev=reverseComplement(seq.ext)
  seq.rev.pos=seq.rev[id.pos]
  seq.rev.neg=seq.rev[id.neg]
  names(seq.rev.pos)=paste0('pos',1:length(seq.rev.pos))
  names(seq.rev.neg)=paste0('neg',1:length(seq.rev.neg))
  writeXStringSet(seq.rev.pos,file=file.path(outpath,paste0('seq.rev.',name.export,'.pos.fasta')))
  writeXStringSet(seq.rev.neg,file=file.path(outpath,paste0('seq.rev.',name.export,'.neg.fasta')))
  
  #obtain crop sequence
  gr.ext2=extendReg(gr,ext=ext.crop)
  gr.ext2.pos=gr.ext2[id.pos]
  gr.ext2.neg=gr.ext2[id.neg]
  
  
  getCropSeq<-function(gr2,ntimes){
    idx=sample(1:1000,length(gr2)*ntimes,replace=T)
    gr2=rep(gr2,ntimes)
    starts=start(gr2)+idx
    gr.crop=gr2
    start(gr.crop)=starts
    end(gr.crop)=starts+1000
    summary(end(gr.crop)-start(gr.crop))
    seq.crop=Biostrings::getSeq(genome,gr.crop)
    seq.crop
  }
  
  seq.crop.pos= getCropSeq(gr.ext2.pos,ntimes)
  seq.crop.neg= getCropSeq(gr.ext2.neg,ntimes)
  names(seq.crop.pos)=paste0('pos',1:length(seq.crop.pos))
  names(seq.crop.neg)=paste0('neg',1:length(seq.crop.neg))
  writeXStringSet(seq.crop.pos,file=file.path(outpath,paste0('seq.crop.',name.export,'.pos.fasta')))
  writeXStringSet(seq.crop.neg,file=file.path(outpath,paste0('seq.crop.',name.export,'.neg.fasta')))
  

}


#############################################################################
##  data argumentation by adding crop and reverse complementary sequence 
#############################################################################


dataAug<-function(dat,outpath){
  
  celltypes=names(table(dat$celltype))
  for(icelltype in 1:length(celltypes)){
    message(celltypes[icelltype])
    dat1=dat[dat$celltype==celltypes[icelltype],]
    if(any(!is.na(dat1$log2FC))){
      id.pos=dat1$fdr<fdr.up & abs(dat1$log2FC)>0.1
      id.neg=dat1$fdr>fdr.down
    }else{
      id.pos=dat1$fdr<fdr.up
      id.neg=dat1$fdr>fdr.down
    }
    
    message(sum(id.pos),' ',sum(id.neg))
    
    createRevAndCropSeq(dat1,id.pos,id.neg,ext.crop=1000,ntimes=1,paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath)
  }
  
}

# functions to calculate dimension for CNN in pytorch

calChannel<-function(inlen,padding=0,dilation=1,kernelsize,stride=1){
  outlen=floor((inlen+2*padding-dilation*(kernelsize-1)-1)/stride+1)
  outlen
}

#calChannel(inlen=1001,kernelsize=4)
  
#calChannel(inlen=998,kernelsize=4,stride=4)

#calChannel(inlen=249,kernelsize=4)

#calChannel(inlen=246,kernelsize=4,stride=4)
