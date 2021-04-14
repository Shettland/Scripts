library(dplyr)
library(stringr)
library(ballgown)
library(edgeR)
library(tidyverse)
library(reshape)

tab = read.table("/home/p.mata/TFG/Analysis/T.britovi/get_homologues/trichinella_intersection/pangenome_matrix_genes_t0.tr.tab", header = TRUE)
colnames(tab)[9] <- "trichinella_spiralis_ISS3_PRJNA257433_WBPS15"
colnames(tab)[8] <- "trichinella_spiralis_ISS195_PRJNA12603_WBPS15"
tab[,1] <- NULL
tab[,5] <- NULL
tab <- tab %>% mutate_all(str_replace_all, "\\..*", "")
tab <- tab %>% mutate_all(str_replace_all, ",.*", "")
colnames(tab) <- gsub("\\.", "_", colnames(tab))
colnames(tab) <- gsub("_protein_faa", "", colnames(tab))

write.table(tab, file = "gene_name_matrix.tab", sep = "\t")

tab2 = read.table("/home/p.mata/TFG/Analysis/RNAseq/gene_name_matrix.tab", sep = "\t", header=TRUE)

paths = c("/home/p.mata/TFG/Analysis/RNAseq/t_spiralis/", "/home/p.mata/TFG/Analysis/RNAseq/t_britovi/",
          "/home/p.mata/TFG/Analysis/RNAseq/t_pseudospiralis/ballgown_ISS13/",
          "/home/p.mata/TFG/Analysis/RNAseq/t_pseudospiralis/ballgown_ISS141/", 
          "/home/p.mata/TFG/Analysis/RNAseq/t_pseudospiralis/ballgown_ISS176/",
          "/home/p.mata/TFG/Analysis/RNAseq/t_pseudospiralis/ballgown_ISS588/")

byx = c("trichinella_spiralis_ISS3_PRJNA257433_WBPS15", "trichinella_britovi_PRJNA257433_WBPS15",
         "trichinella_pseudospiralis_ISS13_PRJNA257433_WBPS15", "trichinella_pseudospiralis_ISS141_PRJNA257433_WBPS15",
         "trichinella_pseudospiralis_ISS176_PRJNA257433_WBPS15", "trichinella_pseudospiralis_ISS588_PRJNA257433_WBPS15")

saves = c("ISS3", "britovi", "ISS13", "ISS141", "ISS176", "ISS588")


n=0
for (path in paths){
  setwd(path)

  n=n+1
  savename <- paste(saves[n],"merged_FPKM.txt", sep = "_")

  my.data = ballgown(dataDir="ballgown/", samplePattern='*_ballgown', meas='all')
  expression <- gexpr(my.data)
  expression <- as.data.frame(expression)
  expression$gene <- rownames(expression)
  expression <- expression[, c(ncol(expression), 1:c(ncol(expression)-1))]
  write.table(x = expression, file = savename, col.names = T, row.names = F, sep = "\t", quote = F)

  prueba <- read.table(savename, sep="\t", header = TRUE)
  prueba$gene <- gsub("gene:", "", prueba$gene)

  fin <- merge(tab, prueba, by.x = byx[n], by.y = "gene")
  fin2 <- fin
  rownames(fin2) <- NULL
  if (ncol(fin2)>8){
    fin2$FPKM_mean <- rowMeans(fin2[,8:ncol(fin2)])
  }
  
  colnum = grep("ballgown$", colnames(fin2))
  
  for ( num in colnum ){
    colnames(fin2)[num] <- paste(saves[n], colnames(fin2)[num], sep = "_")
  }
  
  savename2 <- paste(saves[n],"merged_table.txt", sep = "_")
  write.table(fin2, file = savename2, sep="\t", row.names = FALSE)
}

setwd("/home/p.mata/TFG/Analysis/RNAseq/")
colnames(tab2)

#pISS3 = read.table("/home/p.mata/TFG/Analysis/RNAseq/t_spiralis/ISS3_merged_tab.txt", sep = "\t", header=TRUE)
#all.equal(finISS3, pISS3)

finISS3 = read.table("/home/p.mata/TFG/Analysis/RNAseq/t_spiralis/ISS3_merged_table.txt", sep = "\t", header=TRUE)
finISS13 = read.table("/home/p.mata/TFG/Analysis/RNAseq/t_pseudospiralis/ballgown_ISS13/ISS13_merged_table.txt", sep = "\t", header=TRUE)
finISS141 = read.table("/home/p.mata/TFG/Analysis/RNAseq/t_pseudospiralis/ballgown_ISS141/ISS141_merged_table.txt", sep = "\t", header=TRUE)
finISS176 = read.table("/home/p.mata/TFG/Analysis/RNAseq/t_pseudospiralis/ballgown_ISS176/ISS176_merged_table.txt", sep = "\t", header=TRUE)
finISS588 = read.table("/home/p.mata/TFG/Analysis/RNAseq/t_pseudospiralis/ballgown_ISS588/ISS588_merged_table.txt", sep = "\t", header=TRUE)
finbritovi = read.table("/home/p.mata/TFG/Analysis/RNAseq/t_britovi/britovi_merged_table.txt", sep = "\t", header = TRUE)

listab <- list(finISS3, finbritovi, finISS13, finISS141, finISS176, finISS588)
vec = saves

temtab = finISS3[,1:7]
i=0
for ( tabname in listab ){
  i=i+1
  temtab <- merge(temtab, tabname, by = c("trichinella_spiralis_ISS3_PRJNA257433_WBPS15", "trichinella_britovi_PRJNA257433_WBPS15", 
                                          "trichinella_pseudospiralis_ISS13_PRJNA257433_WBPS15", "trichinella_pseudospiralis_ISS141_PRJNA257433_WBPS15", 
                                          "trichinella_pseudospiralis_ISS176_PRJNA257433_WBPS15", "trichinella_pseudospiralis_ISS588_PRJNA257433_WBPS15", 
                                          "trichinella_spiralis_ISS195_PRJNA12603_WBPS15"))
  if ("FPKM_mean" %in% colnames(temtab)){
    colnames(temtab)[which(colnames(temtab)=="FPKM_mean")] <- paste(vec[i],"FPKM_mean", sep = "_")
  }

}

# temtab <- temtab[!duplicated(as.list(temtab))] # permite eliminar columnas duplicadas
# Si quisiera guardar la última columna duplicada en vez de la primera: !duplicated(as.list(temtab)),fromLast = TRUE)

endtab = temtab

#fpkm <- read.table(file = "ult_FPKM_merged_table.txt", sep="\t", header = TRUE)
# means <- grep("_ballgown$|_mean|^ISS...", colnames(endtab)) # '|' = or , ^ = inicio de string, $ = fin de string #

means <- grep("ballgown$", colnames(endtab), ignore.case = TRUE)
endtab$FPKM_mean_tot <- rowMeans(endtab[,means])
endtab$FPKM_sum <- rowSums(endtab[,means])
endtab2 <- endtab[order(endtab$FPKM_mean_tot, decreasing = TRUE),]

#fpkm3 <- fpkm2[,c(1:28,30,29)] # para reordenar las columnas

write.table(endtab2, file = "/home/p.mata/TFG/Analysis/RNAseq/fin_FPKM_merged_table.txt", sep ="\t", row.names = FALSE)

#AAAAA = read.table(file = "ult_FPKM_merged_table.txt", sep="\t", header = TRUE)

