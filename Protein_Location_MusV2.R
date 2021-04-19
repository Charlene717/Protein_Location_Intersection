rm(list = ls()) #清除變數

# Main file
setwd(getwd()) ## Set current working directory
PathName <- getwd() ## Set output directroy
#
datatext <- read.table(paste0(PathName,"/mouse_compartment_textmining_full.tsv"),  # 資料檔名
                      header=F,          # 資料中的第一列，作為欄位名稱
                      sep="\t")           # 將逗號視為分隔符號來讀取資料

datainte <- read.table(paste0(PathName,"/mouse_compartment_integrated_full.tsv"),  # 資料檔名 
                      header=F,          # 資料中的第一列，作為欄位名稱
                      sep="\t")           # 將逗號視為分隔符號來讀取資料
dataknow <- read.table(paste0(PathName,"/mouse_compartment_knowledge_full.tsv"),  # 資料檔名
                       header=F,          # 資料中的第一列，作為欄位名稱
                       sep="\t")           # 將逗號視為分隔符號來讀取資料

# Test
subcellular <- read.csv(paste0(PathName,"/subcellular2.csv"),header = F)  # 資料檔名
    #                   header=F,          # 資料中的第一列，作為欄位名稱
     #                  sep=",")           # 將逗號視為分隔符號來讀取資料
subcellular2 <- read.table(paste0(PathName,"/subcellular2.csv"),  # 資料檔名
                           header=F,          # 資料中的第一列，作為欄位名稱
                           sep=",")           # 將逗號視為分隔符號來讀取資料
subcellular3 <- read.delim(paste0(PathName,"/subcellular2.csv"),  # 資料檔名
                       header=F,          # 資料中的第一列，作為欄位名稱
                       sep=",")           # 將逗號視為分隔符號來讀取資料
# Test
### !!!
datainte2 <- read.delim(paste0(PathName,"/mouse_compartment_integrated_full.tsv"),  # 資料檔名 
                       header=F,          # 資料中的第一列，作為欄位名稱
                       sep="\t")           # 將逗號視為分隔符號來讀取資料


datainte <- subcellular[,c(-5,-6)]
datainte <- datainte2
# datainte <- dataknow[,c(-5,-6)]
# datainte <- datatext
# load package 'dplyr'
library(dplyr) # Basic data manupilation tools

# load package 'data.table'
library(data.table)

# Grepl cyto
datainte_cyto <- datainte[datainte[,5] >= 2,]
datainte_cyto2_1 <- datainte_cyto[grepl("cyto", datainte_cyto[,4], ignore.case=TRUE),]
datainte_cyto2_2 <- datainte_cyto[grepl("intracellular", datainte_cyto[,4], ignore.case=TRUE),]
#datainte_cyto2_3 <- datainte_cyto[grepl("Intracellular", datainte_cyto[,4], ignore.case=TRUE),]

datainte_cyto2 <- rbind(datainte_cyto2_1,datainte_cyto2_2)#,datainte_cyto2_3)
datainte_cyto2 <- datainte_cyto2[!grepl("mmu-miR-", datainte_cyto2[,2], ignore.case=TRUE),]


Mre11 <- datainte_cyto[grepl("Mre11", datainte_cyto[,2], ignore.case=TRUE),]


datainte_cyto3 <- datainte_cyto2[ ,c(2,5)]
#datainte_cyto3 <- as.matrix(datainte_cyto3)

# library("data.table") 
# setDT(df)[, .SD[1], by = .(x1, x2)] 
datainte_cyto4 <- setDT(datainte_cyto3)[, .SD[1], by = .(V2)] 
#error# datainte_cyto4 = unique(datainte_cyto3, by = datainte_cyto3[,1])
#OK# 
datainte_cyto5 = unique(datainte_cyto3, by = "V2")

write.table(datainte_cyto3, file=paste0(PathName,"/cyto_MusTTT.txt"),  
            sep="\t", row.names=FALSE)

write.table(datainte_cyto4, file=paste0(PathName,"/cyto_Mus.txt"),  
            sep="\t", row.names=FALSE)

write.table(datainte_cyto4, file=paste0(PathName,"/cyto_Mus.csv"),  
            sep=",", row.names=FALSE)


# Grepl Etracellular matrix
datainte_extra <- datainte[datainte[,5] >= 2,]
datainte_extra2_1 <- datainte_extra[grepl("extracellular", datainte_extra[,4], ignore.case=TRUE),]
datainte_extra2 <- datainte_extra2_1[!grepl("exosome", datainte_extra2_1[,4], ignore.case=TRUE),]
datainte_extra2 <- datainte_extra2[!grepl("mmu-miR-", datainte_extra2[,2], ignore.case=TRUE),]

datainte_extra3 <- datainte_extra2[ ,c(2,5)]
#datainte_extra3 <- as.matrix(datainte_extra3)
datainte_extra4 <- setDT(datainte_extra3)[, .SD[1], by = .(V2)] 

write.table(datainte_extra4, file=paste0(PathName,"/Extracellular_Mus.txt"),  
            sep="\t", row.names=FALSE)

write.table(datainte_extra4, file=paste0(PathName,"/Extracellular_Mus.csv"),  
            sep=",", row.names=FALSE)


# Convert Genenames
library('org.Hs.eg.db')
library('org.Mm.eg.db')
library(AnnotationDbi)

CytoGene  <- as.vector(as.matrix(datainte_cyto4[,1]))
# CytoGene  <- as.vector(as.character(datainte_cyto4[,1]))
CytoGene2 <- AnnotationDbi::select(org.Mm.eg.db, keys=CytoGene, columns='ENSEMBL', keytype='SYMBOL') 

#NG# CytoGene3 <- CytoGene2[-which(CytoGene2$ENSEMBL == ""), ]
#NG# CytoGene3 <- CytoGene2[!(CytoGene2$ENSEMBL == ""), ]
CytoGene3 <- CytoGene2[!is.na(CytoGene2$ENSEMBL), ]
CytoGene3 <- unique(CytoGene3, by = "ENSEMBL")

write.table(CytoGene3, file=paste0(PathName,"/CytoGene_Mus_ENSEMBL.csv"),  
            sep=",", row.names=FALSE)

CytoGene4 <- CytoGene3[,2]
CytoGene4 <- AnnotationDbi::select(org.Mm.eg.db, keys=CytoGene4, columns='SYMBOL', keytype='ENSEMBL') 
CytoGene5 <- CytoGene4[!is.na(CytoGene4$SYMBOL), ]
CytoGene5 <- unique(CytoGene5, by = "SYMBOL")
write.table(CytoGene5, file=paste0(PathName,"/CytoGene_Mus_SYMBOL.csv"),  
            sep=",", row.names=FALSE)

#library(readxl)
# Error# FNkdGenelise <- read_excel(paste0(PathName,"/FNkd unique.xlsx"), sheet = "工作表1")

FNkdGenelise <- read.csv(paste0(PathName,"/FNkd unique.csv"),header = T)  # 資料檔名
FNkdExoGene  <- as.vector(as.matrix(FNkdGenelise[,1]))
FNkdExoGene2 <- AnnotationDbi::select(org.Mm.eg.db, keys=FNkdExoGene, columns='ENSEMBL', keytype='SYMBOL') 
FNkdExoGene3 <- FNkdExoGene2[!is.na(FNkdExoGene2$ENSEMBL), ]
FNkdExoGene3 <- unique(FNkdExoGene3, by = "ENSEMBL")
write.table(FNkdExoGene3, file=paste0(PathName,"/FNkdExoGene_Mus_ENSEMBL.csv"),  
            sep=",", row.names=FALSE)

FNkdExoGene4 <- FNkdExoGene3[,2]
FNkdExoGene4 <- AnnotationDbi::select(org.Mm.eg.db, keys=FNkdExoGene4, columns='SYMBOL', keytype='ENSEMBL') 
FNkdExoGene5 <- FNkdExoGene4[!is.na(FNkdExoGene4$SYMBOL), ]
FNkdExoGene5 <- unique(FNkdExoGene5, by = "SYMBOL")
write.table(FNkdExoGene5, file=paste0(PathName,"/FNkdExoGene_Mus_SYMBOL.csv"),  
            sep=",", row.names=FALSE)
