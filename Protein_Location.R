rm(list = ls()) #清除變數

# Main file
setwd(getwd()) ## Set current working directory
PathName <- getwd() ## Set output directroy
#
# datatext <- read.table(paste0(PathName,"/human_compartment_textmining_full.tsv"),  # 資料檔名 
#                       header=F,          # 資料中的第一列，作為欄位名稱
#                       sep="\t")           # 將逗號視為分隔符號來讀取資料

datainte <- read.table(paste0(PathName,"/human_compartment_integrated_full.tsv"),  # 資料檔名 
                      header=F,          # 資料中的第一列，作為欄位名稱
                      sep="\t")           # 將逗號視為分隔符號來讀取資料


# load package 'dplyr'
library(dplyr) # Basic data manupilation tools

# load package 'data.table'
library(data.table)

# Grepl cyto
datainte_cyto <- datainte[datainte[,5] >= 2,]
datainte_cyto2_1 <- datainte_cyto[grepl("cyto", datainte_cyto[,4], ignore.case=TRUE),]
datainte_cyto2_2 <- datainte_cyto[grepl("intracellular", datainte_cyto[,4], ignore.case=TRUE),]
datainte_cyto2 <- rbind(datainte_cyto2_1,datainte_cyto2_2)
datainte_cyto2 <- datainte_cyto2[!grepl("mmu-miR-", datainte_cyto2[,2], ignore.case=TRUE),]

datainte_cyto3 <- datainte_cyto2[ ,c(2,5)]
#datainte_cyto3 <- as.matrix(datainte_cyto3)

# library("data.table") 
# setDT(df)[, .SD[1], by = .(x1, x2)] 
datainte_cyto4 <- setDT(datainte_cyto3)[, .SD[1], by = .(V2)] 
#error# datainte_cyto4 = unique(datainte_cyto3, by = datainte_cyto3[,1])
#OK# datainte_cyto5 = unique(datainte_cyto3, by = "V2")

write.table(datainte_cyto4, file=paste0(PathName,"/cyto_Homo.txt"),  
            sep="\t", row.names=FALSE)

write.table(datainte_cyto4, file=paste0(PathName,"/cyto_Homo.csv"),  
            sep=",", row.names=FALSE)


# Grepl Etracellular matrix
datainte_extra <- datainte[datainte[,5] >= 2,]
datainte_extra2_1 <- datainte_extra[grepl("extracellular", datainte_extra[,4], ignore.case=TRUE),]
datainte_extra2 <- datainte_extra2_1[!grepl("exosome", datainte_extra2_1[,4], ignore.case=TRUE),]
datainte_extra2 <- datainte_extra2[!grepl("mmu-miR-", datainte_extra2[,2], ignore.case=TRUE),]

datainte_extra3 <- datainte_extra2[ ,c(2,5)]
#datainte_extra3 <- as.matrix(datainte_extra3)
datainte_extra4 <- setDT(datainte_extra3)[, .SD[1], by = .(V2)] 

library(Hmisc)
TTT <- datainte_extra4[,1]
datainte_extra5 <- capitalize(as.character(datainte_extra4[,1]))

write.table(datainte_extra4, file=paste0(PathName,"/Extracellular_Homo.txt"),  
            sep="\t", row.names=FALSE)

write.table(datainte_extra4, file=paste0(PathName,"/Extracellular_Homo.csv"),  
            sep=",", row.names=FALSE)



