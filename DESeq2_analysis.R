library(dplyr)
library(stringr)
library(tidyverse)
library(reshape)
library(ggplot2)
library(DESeq2)
library(showtext)
library(vsn)
library(pheatmap)
library(EnhancedVolcano)

showtext_auto()

newtab <- as.data.frame(read.table(file = "/home/p.mata/TFG/Analysis/RNAseq/new_t_spiralis/featurecounts_res_new/merged_gene_counts.txt", sep = "\t", header = TRUE))
tabISS3 <- as.data.frame(read.table(file = "/home/p.mata/TFG/Analysis/RNAseq/featurecounts_res/merged_gene_counts_20210315_RNASEQ_TSPIRALIS_ISS3.txt", sep = "\t", header = TRUE))

head(newtab)
head(tabISS3)

mergedtab <- merge(tabISS3, newtab, by="Geneid")
head(mergedtab)

newcols = c(colnames(mergedtab))
newcols = newcols[2:length(newcols)]
newcols = as.data.frame(newcols)
larval_stage = c("muscle", "adult", "newborn", "muscle", "muscle", "muscle", "muscle", "muscle", "adult", "newborn", "muscle", "muscle", "muscle")
batch = c("I", "I", "I", "II", "III", "III", "III", "IV", "IV", "IV", "V", "V", "V")
newcols

newcond = data.frame(larval_stage, batch)
newcond$larval_stage <- factor(newcond$larval_stage)
newcond$batch <- factor(newcond$batch)
newcond
row.names(newcond) <- newcols[,1]
droplevels(newcond)

row.names(mergedtab) <- mergedtab[,1]
mergedtab[,1] <- NULL
head(mergedtab)

all(rownames(newcond) == colnames(mergedtab))

newdsq <- DESeqDataSetFromMatrix(countData = mergedtab,
                              colData = newcond,
                              design = ~ batch + larval_stage)
newdsq

newdsq$larval_stage <- relevel(newdsq$larval_stage, ref = "muscle")

keep2 <- rowSums(counts(newdsq)) >= 1
newdsq <- newdsq[keep2, ]

setwd("/home/p.mata/TFG/Analysis/RNAseq/new_t_spiralis/new_DESeq2_results")

newdss <- DESeq(newdsq)
newres <- results(newdss)
newresmn <- results(newdss, contrast = c("larval_stage", "muscle", "newborn"))
newresma <- results(newdss, contrast = c("larval_stage", "muscle", "adult"))
newresan <- results(newdss, contrast= c("larval_stage", "adult", "newborn"))
colnames(newresmn)

write.table(newresmn, file="new_muscle_vs_newborn.txt", sep="\t")
write.table(newresma, file="new_muscle_vs_adult.txt", sep="\t")

resultsNames(newdss)
newresmnLFC <- lfcShrink(newdss, coef="larval_stage_newborn_vs_muscle", type="apeglm")
newresmaLFC <- lfcShrink(newdss, coef="larval_stage_adult_vs_muscle", type="apeglm")
newresmnLFC

write.table(newresmnLFC, file="muscle_vs_newbornLFC_new.txt", sep="\t")
write.table(newresmaLFC, file="muscle_vs_adultLFC_new.txt", sep="\t")


DESeq2::plotMA(newresmnLFC)
dev.print(pdf, file="MAplot_mvsnbLFC_new.pdf")
DESeq2::plotMA(newresmaLFC)
dev.print(pdf, file="MAplot_musclevsadultLFC.pdf")

secretedhomologes=("T01_1683 T01_14099 T01_5967 T01_13792 T01_5965 T01_1954 T01_5642 T01_15193 
T01_4042 T01_492 T01_1345 T01_14976 T01_10578 T01_10690 T01_5983 T01_9166 T01_3281 
T01_10023 T01_341 T01_15774 T01_5223 T01_10091 T01_8955 T01_580 ")
secretedvec=c("T01_1683","T01_14099","T01_5967","T01_13792","T01_5965","T01_1954","T01_5642","T01_15193","T01_4042","T01_492","T01_1345","T01_14976","T01_10578","T01_10690","T01_5983","T01_9166","T01_3281","T01_10023","T01_341","T01_15774","T01_5223","T01_10091","T01_8955","T01_580")


plotCounts(newdss, gene="T01_1683", intgroup="larval_stage")
plotCounts(newdss, gene="T01_14099", intgroup="larval_stage")
plotCounts(newdss, gene="T01_13792", intgroup="larval_stage")

#-----------------------------TRANSFORMING RESULTS AND CORRECTING BATCH EFFECT----------------------

newvsd <- vst(newdss, blind=FALSE)
newrld <- rlog(newdss, blind=FALSE)
tail(assay(newvsd))
write.table(assay(newvsd), file="vsdfile.txt", sep="\t")
write.table(assay(newrld), file="rldfile.txt", sep="\t")
newntd <- normTransform(newdss)
write.table(assay(newntd), file="ntdfile.txt", sep="\t")
#plotPCA(newntd, intgroup=c("larval_stage", "batch")); dev.print(pdf, "newntdPCA.pdf")

head(newvsd)

newmat <- assay(newntd)
modelmat <- model.matrix(~larval_stage, colData(newntd))
newmat <- limma::removeBatchEffect(newmat, batch=newntd$batch, design=modelmat)
assay(newntd) <- newmat
#plotPCA(newntd, intgroup=c("larval_stage", "batch")); dev.print(pdf, "RFnewntdPCA.pdf")

newmat <- assay(newvsd)
modelmat <- model.matrix(~larval_stage, colData(newvsd))
newmat <- limma::removeBatchEffect(newmat, batch=newvsd$batch, design=modelmat)
assay(newvsd) <- newmat
#plotPCA(newvsd, intgroup=c("larval_stage", "batch")); dev.print(pdf, "RFnewvsdPCA.pdf")

newmat <- assay(newrld)
modelmat <- model.matrix(~larval_stage, colData(newrld))
newmat <- limma::removeBatchEffect(newmat, batch=newrld$batch, design=modelmat)
assay(newrld) <- newmat
#plotPCA(newrld, intgroup=c("larval_stage", "batch")); dev.print(pdf, "RFnewrldPCA.pdf")

meanSdPlot(assay(newntd))
meanSdPlot(assay(newvsd))
meanSdPlot(assay(newrld))

grep("T01_14198", assay(newntd))

#------------------------PLOTTING HEATMAP-----------------


datacolf <- as.data.frame(colData(newdss)[,c("larval_stage","batch")])

filteredntd <- subset(newntd, rownames(newntd) %in% secretedvec)
pheatmap(assay(filteredntd), cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=datacolf, main="DESeq2_ntd_heatmap",
         display_numbers = FALSE, annotation_colors = list(batch=c(I="#b8fffa", II="#437370", III="#122e2c", IV="#dbff33", V="#804539"),
           larval_stage=c(muscle = "#db03fc", adult = "#a1ff26", newborn = "#583aba")), filename = "Pheatmap_ntd_new.pdf")

filteredvsd <- subset(newvsd, rownames(newvsd) %in% secretedvec)
pheatmap(assay(filteredvsd), cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=datacolf, main="DESeq2_vsd_heatmap",
         display_numbers = FALSE, annotation_colors = list(batch=c(I="#b8fffa", II="#437370", III="#122e2c", IV="#dbff33", V="#804539"),
           larval_stage=c(muscle = "#db03fc", adult = "#a1ff26", newborn = "#583aba")), filename = "Pheatmap_vsd_new.pdf")

filteredrld <- subset(newrld, rownames(newrld) %in% secretedvec)
pheatmap(assay(filteredrld), cluster_rows=TRUE, show_rownames=TRUE, legend=TRUE,
         cluster_cols=TRUE, annotation_col=datacolf, main="DESeq2_rld_heatmap", 
         display_numbers = FALSE, annotation_colors = list(batch=c(I="#b8fffa", II="#437370", III="#122e2c", IV="#dbff33", V="#804539"),
           larval_stage=c(muscle = "#db03fc", adult = "#a1ff26", newborn = "#583aba")), filename = "Pheatmap_rld_new.pdf")


#-----------------------------PLOTTING PCAs------------------------------

plotPCA(newvsd, intgroup=c("larval_stage")); dev.print(pdf, "RFnewvsdPCA.pdf")
plotPCA(newrld, intgroup=c("larval_stage")); dev.print(pdf, "RFnewrldPCA.pdf")
plotPCA(newntd, intgroup=c("larval_stage")); dev.print(pdf, "RFnewntdPCA.pdf")

normcounts <- counts(newdss, normalized=TRUE)
write.table(normcounts, file="normalized_counts_new.txt", sep="\t")

#----------------------------VOLCANO PLOTS-----------------------------------

keyvals.shape <- ifelse(rownames(newresmnLFC) %in% secretedvec, 18, 3)
  keyvals.shape[is.na(keyvals.shape)] <- 3
  names(keyvals.shape)[keyvals.shape == 3] <- 'Unclassified'
  names(keyvals.shape)[keyvals.shape == 18] <- 'Secreted homologes'

keyvals.colour <- ifelse(rownames(newresmnLFC) %in% secretedvec, 'blue','skyblue')
  keyvals.colour[is.na(keyvals.colour)] <- 'skyblue'
  names(keyvals.colour)[keyvals.colour == 'skyblue'] <- 'Unclassified'
  names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Secreted homologes'

EnhancedVolcano(newresmnLFC,
  lab=rownames(newresmnLFC),
  x = "log2FoldChange",
  y = "pvalue",
  labSize = 2.4,
  title = "DESeq2 results",
  subtitle = "Muscular larvae vs Newborn Larvae (LFC)",
  #pointSize = 0.8,
  selectLab = secretedvec,
  labCol = "red",
  legendLabSize = 10,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  pointSize = ifelse(rownames(newresmnLFC) %in% secretedvec, 4, 0.8),
  shapeCustom = keyvals.shape,
  colCustom = keyvals.colour,
  colAlpha = 0.9, 
  border = "full",
  borderWidth = 1.0,
  borderColour = "black")

dev.print(pdf, file = "VolcanoMusclevsNewborn.pdf")

keyvals.shape <- ifelse(rownames(newresmaLFC) %in% secretedvec, 18, 3)
  keyvals.shape[is.na(keyvals.shape)] <- 3
  names(keyvals.shape)[keyvals.shape == 3] <- 'Unclassified'
  names(keyvals.shape)[keyvals.shape == 18] <- 'Secreted homologes'

keyvals.colour <- ifelse(rownames(newresmaLFC) %in% secretedvec, 'blue','skyblue')
  keyvals.colour[is.na(keyvals.colour)] <- 'skyblue'
  names(keyvals.colour)[keyvals.colour == 'skyblue'] <- 'Unclassified'
  names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Secreted homologes'

EnhancedVolcano(newresmaLFC,
  lab=rownames(newresmaLFC),
  x = "log2FoldChange",
  y = "pvalue",
  labSize = 2.4,
  title = "DESeq2 results",
  subtitle = "Muscular larvae vs Adult (LFC)",
  #pointSize = 0.8,
  selectLab = secretedvec,
  labCol = "red",
  legendLabSize = 10,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  pointSize = ifelse(rownames(newresmaLFC) %in% secretedvec, 4, 0.8),
  shapeCustom = keyvals.shape,
  colCustom = keyvals.colour,
  colAlpha = 0.9, 
  border = "full",
  borderWidth = 1.0,
  borderColour = "black")

dev.print(pdf, file = "VolcanoMusclevsAdult.pdf")

rownames(newresmnLFC)


#--------------------GGPLOT EXPRESSION MATRIX-----------------

ntdmatrixfilt <- subset(assay(newntd), rownames(assay(newntd)) %in% secretedvec)
ntdmatrixfilt

rldmatrixfilt <- subset(assay(newrld), rownames(assay(newrld)) %in% secretedvec)
rldmatrixfilt

vsdmatrixfilt <- subset(assay(newvsd), rownames(assay(newvsd)) %in% secretedvec)
vsdmatrixfilt

newcolnames <- paste0(colnames(mergedtab),"_", larval_stage)
newcolnames
colnames(ntdmatrixfilt) <- newcolnames

meltedntd <- melt(ntdmatrixfilt)
colnames(meltedntd) <- c("Gene.Name","sample","reads")
head(meltedntd, 25)
meltedntd$larval_stage <- ifelse(grepl("muscle", meltedntd$sample), "Muscle", ifelse(grepl("adult", meltedntd$sample), "Adult", "Newborn"))
meltedntd <- meltedntd[order(meltedntd$larval_stage),]


colorlarval <- c(Adult = "red", Muscle = "cadetblue2", Newborn = "goldenrod1")

ntdplot <- ggplot(data = meltedntd, aes(x = sample, weight = reads, fill = larval_stage))

tiff("ntd_colplot.tiff", units="in", width=7, height=7, res=300)
ntdplot + geom_bar() + facet_wrap(~Gene.Name) + scale_y_continuous(name="Normalized Reads (ntd)") +
  labs(title="DESeq_ntd_Expression_Matrix_Filtered") + scale_fill_manual(values=colorlarval) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size = 18),
        plot.title = element_text(size = 48),
        axis.title.y = element_text(size = 34),
        strip.text = element_text(size = 20, face="bold"),
        legend.title = element_text(size = 22, face="bold"),
        legend.text=element_text(size=22))  #;dev.print(pdf, "ntd_colplot_new.pdf")
dev.off()

#------------------------SECOND DESEQ ANALYSIS | MUSCLE VS OTHER ----------------------------


copiamerged <- merge(tabISS3, newtab, by="Geneid")

cols = c(colnames(as.data.frame(copiamerged)))
cols = cols[2:length(cols)]
cols = as.data.frame(cols)
larval_stage2 = c("muscle", "other", "other", "muscle", "muscle", "muscle", "muscle", "muscle", "other", "other", "muscle", "muscle", "muscle")
batch2 = c("I", "I", "I", "II", "III", "III", "III", "IV", "IV", "IV", "V", "V", "V")
cols

cond2 = data.frame(larval_stage2, batch2)
cond2ISS3 <- cond2
cond2ISS3$larval_stage2 <- larval_stage2
cond2ISS3$larval_stage2 <- factor(cond2ISS3$larval_stage2)
cond2ISS3$batch2 <- factor(cond2ISS3$batch2)
cond2ISS3
row.names(cond2ISS3) <- cols[,1]
droplevels(cond2ISS3)

as.data.frame(tabs[1])

row.names(copiamerged) <- copiamerged[,1]
####tabISS3 <- tabISS3[, 2:length(colnames(tabISS3))]
copiamerged[1] <- NULL
####tabISS3 <- tabISS3[, rownames(condISS3)]
head(copiamerged, 4)

all(rownames(cond2ISS3) == colnames(copiamerged))

#---------------------------------------------------------------------

dsqbatch2 <- DESeqDataSetFromMatrix(countData = copiamerged,
                              colData = cond2ISS3,
                              design = ~ batch2 + larval_stage2)
dsqbatch2

dsqbatch2$larval_stage2 <- relevel(dsqbatch2$larval_stage2, ref = "muscle")

keep <- rowSums(counts(dsqbatch2)) >= 1
dsqbatch2 <- dsqbatch2[keep, ]

dsq22 <- DESeq(dsqbatch2)
resot22 <- results(dsq22)
resot22mo <- results(dsq22, contrast = c("larval_stage2", "muscle", "other"))

getwd()

colnames(resot22mo)
write.table(resot22mo, file="muscle_vs_other.txt", sep="\t")

resultsNames(dsq22)
resot22moLFC <- lfcShrink(dsq22, coef="larval_stage2_other_vs_muscle", type="apeglm")
resot22moLFC

write.table(resot22moLFC, file="muscle_vs_otherLFC.txt", sep="\t")

DESeq2::plotMA(resot22moLFC)
dev.print(pdf, file="MAplot_musclevsnewbornLFC.pdf")
DESeq2::plotMA(newresmaLFC)
dev.print(pdf, file="MAplot_musclevsadultLFC.pdf")

secretedhomologes=c("T01_1683 T01_14099 T01_5967 T01_13792 T01_5965 T01_1954 T01_5642 T01_15193 
T01_4042 T01_492 T01_1345 T01_14976 T01_10578 T01_10690 T01_5983 T01_9166 T01_3281 
T01_10023 T01_341 T01_15774 T01_5223 T01_10091 T01_8955 T01_580 ")
secretedvec=c("T01_1683","T01_14099","T01_5967","T01_13792","T01_5965","T01_1954","T01_5642","T01_15193","T01_4042","T01_492","T01_1345","T01_14976","T01_10578","T01_10690","T01_5983","T01_9166","T01_3281","T01_10023","T01_341","T01_15774","T01_5223","T01_10091","T01_8955","T01_580")

vstothbatch <- vst(dsq22, blind=FALSE)
rldothbatch <- rlog(dsq22, blind=FALSE)
head(assay(vstothbatch))
write.table(assay(vstothbatch), file="vsdotherfile.txt", sep="\t")
write.table(assay(rldothbatch), file="rldotherfile.txt", sep="\t")
ntdothbatch <- normTransform(dsq22)
write.table(assay(ntdoth), file="ntdotherfile.txt", sep="\t")

mat2 <- assay(vstothbatch)
modm2 <- model.matrix(~larval_stage2, colData(vstothbatch))
mat2 <- limma::removeBatchEffect(mat2, batch=vstothbatch$batch, design=modm2)
assay(vstothbatch) <- mat2
plotPCA(vstothbatch, intgroup=c("larval_stage2")); dev.print(pdf, "vsdPCA_other.pdf")

mat2 <- assay(rldothbatch)
modm2 <- model.matrix(~larval_stage2, colData(rldothbatch))
mat2 <- limma::removeBatchEffect(mat2, batch=rldothbatch$batch, design=modm2)
assay(rldothbatch) <- mat2
plotPCA(rldothbatch, intgroup=c("larval_stage2")); dev.print(pdf, "rldPCA_other.pdf")

mat2 <- assay(ntdothbatch)
modm2 <- model.matrix(~larval_stage2, colData(ntdothbatch))
mat2 <- limma::removeBatchEffect(mat2, batch=ntdothbatch$batch, design=modm2)
assay(ntdothbatch) <- mat2
plotPCA(ntdothbatch, intgroup=c("larval_stage2")); dev.print(pdf, "ntdPCA_other.pdf")


meanSdPlot(assay(ntdoth))
meanSdPlot(assay(vsdothbatch))
meanSdPlot(assay(rldothbatch))

grep("T01_14198", assay(newntd))

#------------------------PLOTTING HEATMAP-----------------


datacolf22 <- as.data.frame(colData(dsq22)[,c("larval_stage2","batch2")])

filteredntdoth <- subset(ntdothbatch, rownames(ntdothbatch) %in% secretedvec)
pheatmap(assay(filteredntdoth), cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=datacolf22, main="DESeq2_ntd_heatmap",
         display_numbers = FALSE, annotation_colors = list(batch=c(I="#b8fffa", II="#437370", III="#122e2c"),
           larval_stage=c(muscle = "#db03fc", adult = "#a1ff26", newborn = "#583aba")), filename = "Pheatmap_ntd_other.pdf")

filteredvsdoth <- subset(vstothbatch, rownames(vstothbatch) %in% secretedvec)
pheatmap(assay(filteredvsdoth), cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=datacolf22, main="DESeq2_vsd_heatmap",
         display_numbers = FALSE, annotation_colors = list(batch=c(I="#b8fffa", II="#437370", III="#122e2c"),
           larval_stage=c(muscle = "#db03fc", adult = "#a1ff26", newborn = "#583aba")), filename = "Pheatmap_vsd_other.pdf")

filteredrldoth <- subset(rldothbatch, rownames(rldothbatch) %in% secretedvec)
pheatmap(assay(filteredrldoth), cluster_rows=TRUE, show_rownames=TRUE, legend=TRUE,
         cluster_cols=TRUE, annotation_col=datacolf22, main="DESeq2_rld_heatmap", 
         display_numbers = FALSE, annotation_colors = list(batch=c(I="#b8fffa", II="#437370", III="#122e2c"),
           larval_stage=c(muscle = "#db03fc", adult = "#a1ff26", newborn = "#583aba")), filename = "Pheatmap_rld_other.pdf")

