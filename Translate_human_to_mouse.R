# https://gist.github.com/FloWuenne/f8fc922477df04c1642e9d8945c48d47

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  no_mouse_genes <- length(x)
  no_human_genes <- length(humanx)
  
  if(no_human_genes != no_mouse_genes){
    print("Some genes could not be translated!")
    genes_not_trans <- setdiff(x,genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ",length(genes_not_trans),"genes could not be translated!"),sep=" ")
  }else{
    print("All genes were translated successfully!")
  }
  
  # Print all gene names that could not be translated and the number of genes that were not translated
  
  return(humanx)
}

M <- read.table(paste0(PathName,"/geneset (6).txt"),  # 資料檔名 
                        header=T,          # 資料中的第一列，作為欄位名稱
                        sep="\t")           # 將逗號視為分隔符號來讀取資料
M <- M[-1,]

Mouse_Symbol_M <- convertHumanGeneList(M)
Mouse_Symbol_M <- as.data.frame(Mouse_Symbol_M)
write.csv(Mouse_Symbol_M, file=paste0(PathName,"/",RVersion,"/Mphase_Mus.csv"))

G2M <- read.table(paste0(PathName,"/G2M.txt"),  # 資料檔名 
               header=T,          # 資料中的第一列，作為欄位名稱
               sep="\t")           # 將逗號視為分隔符號來讀取資料
G2M <- G2M[-1,]
Mouse_Symbol_G2M <- convertHumanGeneList(G2M)
Mouse_Symbol_G2M <- as.data.frame(Mouse_Symbol_G2M)
write.csv(Mouse_Symbol_G2M, file=paste0(PathName,"/",RVersion,"/G2Mphase_Mus.csv"))
