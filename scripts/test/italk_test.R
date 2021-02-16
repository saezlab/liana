library(Seurat)
require(data.table)
library(biomaRt)
library(RColorBrewer)
require(pals)
require(EWCE)
library("optparse")
library(devtools)
source('../utils/get_omnipath.R')
source('../pipes/iTALK_pipe.R')

seurat_object <- readRDS('/home/james/sciebo/LR_Benchmark/data/pbmc3k_processed.rds')
str(seurat_object)
### First Test
if(F){
  dataset <- get_omni_resources()
  tmp <- call_italk(dataset[['iTALK']],seurat_object = ,data,assay = 'RNA',.format = F)
}
dataset <- get_omni_resources()
tmp <- call_italk_de(dataset[['iTALK']],seurat_object,assay = 'RNA',.format = F)
