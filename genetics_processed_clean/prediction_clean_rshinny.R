datapath='genetics_processed_clean/'

source('prediction.lib.R')
source('genetics_lib.R')

load('gene.ref.rda')
load('human.motif.rda')

######################################################################################################
### paper 1: Functional regulatory variants implicate distinct transcriptional networks in dementia

idata=1
dat=read_xlsx('Functional regulatory variants implicate distinct transcriptional networks in dementia.xlsx')
key.word = 'dementia'
dim(dat)

fdr.up=0.1
fdr.down=0.8

log2FC.thres=0.1

summary(dat$fdr)
summary(dat$log2FC)

dat$fdr[is.na(dat$fdr)]=1
dat$log2FC[is.na(dat$log2FC)]=0

dat=dat[!is.na(dat$fdr) & !is.na(dat$log2FC),]
dim(dat)

### by cell types
celltypes=names(table(dat$celltype))
celltypes

for(icelltype in 1:length(celltypes)){
  celltype = celltypes[icelltype]
  
  if (grepl("/", celltype)) {
    celltype = sub("/", "", celltype)
  }
  
  message(celltypes[icelltype])
  dat1=dat[dat$celltype==celltypes[icelltype],]
  id.pos=dat1$fdr<fdr.up & abs(dat1$log2FC)>log2FC.thres
  id.neg=dat1$fdr>fdr.down
  
  print("Type: ")
  print(typeof(dat1))
  
  file.name = paste0(key.word, "-", celltype)
  
  message(sum(id.pos),' ',sum(id.neg))
  
  print("3mer:")
  
  evalPerf(dat1, id.pos, id.neg, motifs, feature.type='3mer',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'), outpath, is.output=F, file.name = paste0(file.name, '-3mer'))
  
  print("motif:")
  
  evalPerf(dat1, id.pos, id.neg, motifs, feature.type='motif',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-motif'))
}

### by disease
diseases=names(table(dat$disease))
diseases

for(idisease in 1:length(diseases)){
  disease = diseases[idisease]
  
  if (grepl("/", disease)) {
    disease = sub("/", "", disease)
  }
  
  message(diseases[idisease])
  dat1=dat[dat$disease==diseases[idisease],]
  id.pos=dat1$fdr<fdr.up & abs(dat1$log2FC)>log2FC.thres
  id.neg=dat1$fdr>fdr.down
  
  file.name = paste0(key.word, "-", disease)
  
  message(sum(id.pos),' ',sum(id.neg))
  
  evalPerf(dat1,id.pos,id.neg,motifs,feature.type='3mer',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-3mer'))
  
  evalPerf(dat1,id.pos,id.neg,motifs,feature.type='motif',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-motif'))
}

######################################################################################################################
### paper 4  Genome-wide functional screen of 30 UTR variants uncovers causal variants for human disease and evolution

idata=4
dat=fread('Genome-wide functional screen of 30 UTR variants uncovers causal variants for human disease and evolution.xls')
key.word = 'evolution3UTR'
dim(dat)

fdr.up=0.1
fdr.down=0.8
log2FC.thres=0.1

summary(dat$fdr)
summary(dat$log2FC)


dat$fdr[is.na(dat$fdr)]=1
dat$log2FC[is.na(dat$log2FC)]=0

dat=dat[!is.na(dat$fdr) & !is.na(dat$log2FC),]
dim(dat)


celltypes=names(table(dat$celltype))
celltypes

for(icelltype in 1:length(celltypes)){
  celltype = celltypes[icelltype]
  
  if (grepl("/", celltype)) {
    celltype = sub("/", "", celltype)
  }
  
  message(celltypes[icelltype])
  dat1=dat[dat$celltype==celltypes[icelltype],]
  id.pos=dat1$fdr<fdr.up & abs(dat1$log2FC)>log2FC.thres
  id.neg=dat1$fdr>fdr.down
  
   file.name = paste0(key.word, "-", celltype)
  
  message(sum(id.pos),' ',sum(id.neg))
  
  evalPerf(dat1, id.pos, id.neg, motifs, feature.type='3mer',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'), outpath, is.output=F, file.name = paste0(file.name, '-3mer'))
  
  evalPerf(dat1,id.pos,id.neg,motifs,feature.type='motif',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-motif'))
}

diseases=names(table(dat$disease))
diseases



######################################################################################################
### Paper 6 Transcriptional-regulatory convergence across functional MDD risk variants identified by massively parallel reporter assays 


idata=6
dat=fread('Transcriptional-regulatory convergence across functional MDD risk variants identified by massively parallel reporter assays.xls')
key.word = 'MDD'
dim(dat)


fdr.up=0.1
fdr.down=0.8
log2FC.thres=0.1

summary(dat$fdr)
summary(dat$log2FC)


dat$fdr[is.na(dat$fdr)]=1
dat$log2FC[is.na(dat$log2FC)]=0

dat=dat[!is.na(dat$fdr) & !is.na(dat$log2FC),]
dim(dat)


### by cell types
celltypes=names(table(dat$celltype))
celltypes

### by disease
diseases=names(table(dat$disease))
diseases

for(icelltype in 1:length(celltypes)){
  celltype = celltypes[icelltype]
  
  if (grepl("/", celltype)) {
    celltype = sub("/", "", celltype)
  }
  
  message(celltypes[icelltype])
  dat1=dat[dat$celltype==celltypes[icelltype],]
  id.pos=dat1$fdr<fdr.up & abs(dat1$log2FC)>log2FC.thres
  id.neg=dat1$fdr>fdr.down
  
   file.name = paste0(key.word, "-", celltype)
  
  message(sum(id.pos),' ',sum(id.neg))
  
  evalPerf(dat1, id.pos, id.neg, motifs, feature.type='3mer',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'), outpath, is.output=F, file.name = paste0(file.name, '-3mer'))
  
  evalPerf(dat1,id.pos,id.neg,motifs,feature.type='motif',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-motif'))
}



#############################################################################################################################################################################
### paper 7: Saturation mutagenesis of twenty disease-associated regulatory elements at single base-pair resolution

idata=7

dat=fread('Saturation mutagenesis of twenty disease-associated regulatory elements at single base-pair resolution.GRCh37_ALL.xls')
key.word = 'mutagenesis'
dim(dat)
head(dat)


fdr.up=0.1
fdr.down=0.8
log2FC.thres=0.1

summary(dat$fdr)
summary(dat$log2FC)


dat$fdr[is.na(dat$fdr)]=1
dat$log2FC[is.na(dat$log2FC)]=0

dat=dat[!is.na(dat$fdr) & !is.na(dat$log2FC),]
dim(dat)


### by cell types
celltypes=names(table(dat$celltype))
celltypes

### by disease
diseases=names(table(dat$disease))
diseases

for(icelltype in 1:length(celltypes)){
  celltype = celltypes[icelltype]
  
  if (grepl("/", celltype)) {
    celltype = sub("/", "", celltype)
  }
  
  message(celltypes[icelltype])
  dat1=dat[dat$celltype==celltypes[icelltype],]
  id.pos=dat1$fdr<fdr.up & abs(dat1$log2FC)>0.1
  id.neg=dat1$fdr>fdr.down
  
   file.name = paste0(key.word, "-", celltype)
  
  message(sum(id.pos),' ',sum(id.neg))
  
  evalPerf(dat1, id.pos, id.neg, motifs, feature.type='3mer',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'), outpath, is.output=F, file.name = paste0(file.name, '-3mer'))
  
  evalPerf(dat1,id.pos,id.neg,motifs,feature.type='motif',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-motif'))
}


### by disease
diseases=names(table(dat$disease))
diseases

for(idisease in 1:length(diseases)){
  disease = diseases[idisease]
  
  if (grepl("/", disease)) {
    disease = sub("/", "", disease)
  }
  
  message(diseases[idisease])
  dat1=dat[dat$disease==diseases[idisease],]
  id.pos=dat1$fdr<fdr.up & abs(dat1$log2FC)>log2FC.thres
  id.neg=dat1$fdr>fdr.down
  
  file.name = paste0(key.word, "-", disease)
  
  message(sum(id.pos),' ',sum(id.neg))
  
  evalPerf(dat1,id.pos,id.neg,motifs,feature.type='3mer',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-3mer'))
  
  evalPerf(dat1,id.pos,id.neg,motifs,feature.type='motif',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-motif'))
}



#######################################################################################################################################
### paper 10:  Prioritization of autoimmune disease-associated genetic variants that perturb regulatory element activity in T cells


idata=10
dat=fread('Prioritization of autoimmune disease-associated genetic variants that perturb regulatory element activity in T cells.xls')
key.word = 'autoimmune'
dim(dat)
summary(dat$fdr)
summary(dat$log2FC)

fdr.up=0.1
fdr.down=0.8
log2FC.thres=0.1


dat$fdr[is.na(dat$fdr)]=1
dat$log2FC[is.na(dat$log2FC)]=0
summary(dat$fdr)
summary(dat$log2FC)


dat=dat[!is.na(dat$fdr) & !is.na(dat$log2FC),]
dim(dat)



### by cell types
celltypes=names(table(dat$celltype))
celltypes

### by disease
diseases=names(table(dat$disease))
diseases

for(icelltype in 1:length(celltypes)){
  celltype = celltypes[icelltype]
  
  if (grepl("/", celltype)) {
    celltype = sub("/", "", celltype)
  }
  
  message(celltypes[icelltype])
  dat1=dat[dat$celltype==celltypes[icelltype],]
  id.pos=dat1$fdr<fdr.up & abs(dat1$log2FC)>log2FC.thres
  id.neg=dat1$fdr>fdr.down
  
  message(sum(id.pos),' ',sum(id.neg))
  
  # very imbalanced case how to carefully selected negative set besides statistical signifcance??? (should make the data balanced!!!)
  id.neg=sample(which(id.neg),10*sum(id.pos)) 
  
   file.name = paste0(key.word, "-", celltype)
  
  message(sum(id.pos),' ',sum(id.neg))
  
  evalPerf(dat1, id.pos, id.neg, motifs, feature.type='3mer',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'), outpath, is.output=F, file.name = paste0(file.name, '-3mer'))
  
  evalPerf(dat1,id.pos,id.neg,motifs,feature.type='motif',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-motif'))
}

######################################################################################################
### paper 11: Massively parallel reporter assays and variant scoring identified functional variants and target genes for melanoma loci and highlighted cell-type specificity 

idata=11
dat=fread('Massively parallel reporter assays and variant scoring identified functional variants and target genes for melanoma loci and highlighted cell-type specificity.xls')
key.word = 'melanoma' 
dim(dat)
summary(dat$fdr)
summary(dat$log2FC)

fdr.up=0.1
fdr.down=0.8
log2FC.thres=0.1


dat$fdr[is.na(dat$fdr)]=1
dat$log2FC[is.na(dat$log2FC)]=0
summary(dat$fdr)
summary(dat$log2FC)

dat=dat[!is.na(dat$fdr) & !is.na(dat$log2FC),]
dim(dat)


### by cell types
celltypes=names(table(dat$celltype))
celltypes


### by disease
diseases=names(table(dat$disease))
diseases

for(icelltype in 1:length(celltypes)){
  celltype = celltypes[icelltype]
  
  if (grepl("/", celltype)) {
    celltype = sub("/", "", celltype)
  }
  
  message(celltypes[icelltype])
  dat1=dat[dat$celltype==celltypes[icelltype],]
  id.pos=dat1$fdr<fdr.up & abs(dat1$log2FC)>log2FC.thres
  id.neg=dat1$fdr>fdr.down
  
  file.name = paste0(key.word, "-", celltype)
  
  message(sum(id.pos),' ',sum(id.neg))
  
  evalPerf(dat1, id.pos, id.neg, motifs, feature.type='3mer',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'), outpath, is.output=F, file.name = paste0(file.name, '-3mer'))
  
  evalPerf(dat1,id.pos,id.neg,motifs,feature.type='motif',
           name.export=paste(paste0('data',idata),celltypes[icelltype],sep='.'),outpath,is.output=F, file.name = paste0(file.name, '-motif'))
}
