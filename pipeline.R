library(DESeq2)
library(stats)
library(ggplot2)
library(pheatmap)
library(ggfortify)
library(tximport)
library(GenomicFeatures)
library(readr)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(gpGeneSets)
library(BinfTools)
library(WebGestaltR)
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(plotrix)
##background a mano 
##rimuovere non protein coding gene
###SETTINGS
mainDir <- ""

counts_file <- "" #a csv file of the counts
sample_file <- "" #a csv file containing a sample column and a group column

condition_A <- "" #choose the conditions you want to contrast
condition_B <- ""
background <- paste0(condition_A,condition_B,'_background.txt')

foldchangetreshold <- 1
padjtreshold <- 0.05

##coding genes from gtf
GTFfile <- file.path(mainDir,"")
GTF <- import.gff(con=GTFfile , format="gtf", genome="GRCh38.p5", feature.type="exon")
biotype_table <- elementMetadata(GTF)[ , c("gene_id", "gene_biotype", "gene_name")]
biotype_table <- as.data.frame(biotype_table)
biotype_table <- unique(biotype_table)


#gestione directory
subDir <- paste(condition_A, "vs", condition_B, sep="")

if (dir.exists(file.path(mainDir, subDir)) == FALSE) {
  dir.create(file.path(mainDir, subDir))
}


paramdir = paste('padj=',as.character(padjtreshold),'logfold=',as.character(foldchangetreshold))

if (dir.exists(file.path(subDir, paramdir)) == FALSE) {
  dir.create(file.path(subDir, paramdir))
}

setwd(mainDir)

# importing data
countData <- as.matrix(read.table(counts_file, sep="\t", header=TRUE,row.names = 1, fileEncoding="UTF-8-BOM"))
countData <- countData[ , order(colnames(countData))]
geneLengths <- as.vector(subset(countData, select = c(sizes)))
colData_all <- read.csv(sample_file, header = F , fileEncoding="UTF-8-BOM")
colnames(colData_all) <- c("sample","group")

#BARPLOT
fpkm <- apply(X = countData,
              MARGIN = 2, 
              FUN = function(x) {
                10^9 * x / geneLengths / sum(as.numeric(x))
              })

genecount <- fpkm[rownames(fpkm) == 'CUFF405']
dotchart(genecount, labels=colnames(fpkm),cex=.7,
         main="CHARME",
         xlab="FPKM")


WTD20_col <- colData_all[grepl('WTD_20',colData_all[["group"]]),]
WTD20fpkmmean <- as.data.frame((subset(fpkm, select = c(WTD20_col$sample))))
average_cond_WTD20_fpkm <- rowMeans(WTD20fpkmmean)
WTD20 <- as.vector(average_cond_WTD20_fpkm)

WTD10_col <- colData_all[grepl('WTD_10',colData_all[["group"]]),]
WTD10fpkmmean <- as.data.frame((subset(fpkm, select = c(WTD10_col$sample))))
average_cond_WTD10_fpkm <- rowMeans(WTD10fpkmmean)
WTD10 <- as.vector(average_cond_WTD10_fpkm)

HOMOD20_col <- colData_all[grepl('HOMOD_20',colData_all[["group"]]),]
HOMOD20fpkmmean <- as.data.frame((subset(fpkm, select = c(HOMOD20_col$sample))))
average_cond_HOMOD20_fpkm <- rowMeans(HOMOD20fpkmmean)
HOMOD20 <- as.vector(average_cond_HOMOD20_fpkm)


HOMOD10_col <- colData_all[grepl('HOMOD_10',colData_all[["group"]]),]
HOMOD10fpkmmean <- as.data.frame((subset(fpkm, select = c(HOMOD10_col$sample))))
average_cond_HOMOD10_fpkm <- rowMeans(HOMOD10fpkmmean)
HOMOD10 <- as.vector(average_cond_HOMOD10_fpkm)
mediatriplicati_fpkm <- data.frame(WTD10,WTD20,HOMOD10,HOMOD20)
rownames(mediatriplicati_fpkm) <- rownames(fpkm)

#tbx
jpeg(file.path( subDir, "barplot_tbx5.jpeg"))  
barplot(as.matrix(mediatriplicati_fpkm[c('ENSG00000089225'),]),
        col = '#eb8060',
        main = "TXB5",
        xlab = "Condition",
        ylab = "fpkm",
        beside = TRUE)
dev.off()

jpeg(file.path(subDir, "barplot_charme.jpeg")) 
barplot(as.matrix(mediatriplicati_fpkm[c('CUFF405'),]),
        col = '#eb8060',
        xlab = "Condition",
        ylab = "fpkm",
        main = "CHARME",
        beside = TRUE)
dev.off()

###

colData <- colData_all[grepl(condition_A,colData_all[["group"]])|grepl(condition_B,colData_all[["group"]]),]

condA_col <- colData_all[grepl(condition_A,colData_all[["group"]]),]
condAmean <- as.data.frame((subset(countData, select = c(condA_col$sample))))
average_cond_A <- round(rowMeans(condAmean), digits = 2)

condB_col <- colData_all[grepl(condition_B,colData_all[["group"]]),]
condBmean <- as.data.frame((subset(countData, select = c(condB_col$sample))))
average_cond_B <- round(rowMeans(condBmean), digits = 2)


countData <- (subset(countData, select = c(colData$sample)))


#normalizzations


fpk <- apply( countData, 2, 
               function(x) x/(geneLengths/1000))
 
tpm <- apply(fpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

# deseq

ddset <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ group)

dds <- DESeq(ddset)

res <- results(dds, contrast = c("group", condition_A, condition_B))
res1 <- as.data.frame(res)
res1[,paste0("average_",condition_A)] <- average_cond_A
res1[,paste0("average_",condition_B)] <- average_cond_B

# pca plot
rld <- rlog(dds)
nudge <- position_nudge(y = 3)
pcaplot <- plotPCA(rld, ntop = 500, intgroup = 'group') +
  geom_text(aes(label = colData$sample), position = nudge, size=2) +
  theme_bw()
pcaplot 
ggsave(file.path(subDir,"pcaplot.pdf"), plot= pcaplot)

#filterinf + volcano su total 
res1$expressed <- 'unexpressed'
res1$expressed[res1[,paste0("average_",condition_A)] >= 1 | res1[,paste0("average_",condition_B)] >= 1] <- 'expressed'
res1 <- res1[res1$expressed == 'expressed',]

res1$diffexpressed <- "INVARIANT"
res1$diffexpressed[res1$log2FoldChange >= foldchangetreshold & res1$padj <= padjtreshold] <- 'UP'
res1$diffexpressed[res1$log2FoldChange <= - foldchangetreshold & res1$padj <= padjtreshold] <- 'DOWN'

res1$diffexpressed[(res1[,paste0("average_",condition_A)] >= 1 & res1[,paste0("average_",condition_B)] <= 1) 
                   & (res1$log2FoldChange >= foldchangetreshold & res1$padj <= padjtreshold)] <- paste0('UP_',condition_A)

res1$diffexpressed[(res1[,paste0("average_",condition_A)] <= 1 & res1[,paste0("average_",condition_B)] >= 1) 
                   & (res1$log2FoldChange <= - foldchangetreshold & res1$padj <= padjtreshold)] <- paste0('DOWN_',condition_B)

res1$gene_id <- rownames(res1)
res1 <- merge(res1, biotype_table[, c("gene_id", "gene_biotype", "gene_name")], by="gene_id")
row.names(res1) <- res1$gene_id
res1 <- subset(res1, select = -c(gene_id))
res1 <- subset(res1, select = -c(expressed))

occurences<-table(unlist(res1))
num_UP <- as.data.frame(occurences["UP"])
num_DOWN <- as.data.frame(occurences["DOWN"])
num <- merge(num_UP, num_DOWN)

write.table(num, file.path(subDir, paramdir, 'volcano_counts.tsv'), sep='\t')

expressed_res <- res1[res1$diffexpressed != 'INVARIANT',]

log2foldchange <- res1$log2FoldChange
legend <- res1$diffexpressed
log10padj <- -log10(res1$padj)

write.table(res1, file.path(subDir, paramdir , paste0('geni_',condition_A,condition_B, '.tsv')), sep = '\t')

p <- ggplot(data = as.data.frame(res1), aes(x = log2foldchange, y = log10padj, col= legend)) + 
  geom_point() +
  geom_vline(xintercept=c(- foldchangetreshold, foldchangetreshold), col="black") +
  geom_hline(yintercept= -log10(padjtreshold), col="black") + 
  ylim(0,30) +
  theme_classic()

mycolors <- c("red", "green4", "grey")
names(mycolors) <- c("DOWN", "UP", "INVARIANT")

volcanoplot <- p + scale_colour_manual(values = mycolors)

ggsave("volcanoplot_total.pdf", path = file.path(subDir, paramdir),width = 7, height = 7)


##filtering protein coding + volcano 

filtered_res <- res1[res1$diffexpressed != 'INVARIANT',]
total <- nrow(filtered_res)
filtered_res <- filtered_res[filtered_res$gene_biotype == 'protein_coding',]
protein_coding_genes <- nrow(filtered_res)
lnc <- total - protein_coding_genes


slices <- c(protein_coding_genes,lnc)
lbls <- c("protein coding
          genes", "non protein 
          coding")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
colors = c("red", "yellow")

jpeg(file.path(subDir, paramdir, "piechart_coding.jpeg")) 
piechart <- pie(slices,labels=lbls,explode=0.2,
      main="codingVSnoncoding ", col = colors)
dev.off()


# differentially expressed genes  solo protein coding e solo espressi nel backgorund
diffexpr <- filtered_res
diffexprUP <- filtered_res[filtered_res$diffexpressed == 'UP' | filtered_res$diffexpressed == paste0("UP_",condition_A),]
diffexprDOWN <- filtered_res[filtered_res$diffexpressed == 'DOWN'| filtered_res$diffexpressed == paste0("DOWN_",condition_B),]



ann <- colData
rownames(ann) <- ann$sample
ann<- subset(ann, select = -c(sample) )

selectedg <- rownames(diffexpr)
heatmap <- pheatmap(tpm[selectedg,], scale = 'row', 
                    show_rownames = FALSE, annotation_col = ann)
ggsave(file.path(subDir, "heatmap.png"), plot=heatmap)

# diffexpr$symbol <- mapIds(org.Hs.eg.db, keys = rownames(diffexpr), keytype = "ENSEMBL", column = "SYMBOL")
# write.table(diffexpr,  file.path(subDir, paramdir, 'diffexpr.tsv'),sep='\t')
# diffexpr <- diffexpr[complete.cases(diffexpr), ]
# 
# ##up regulated 
# diffexprUP$symbol <- mapIds(org.Hs.eg.db, keys = rownames(diffexprUP), keytype = "ENSEMBL", column = "SYMBOL")
# write.table(diffexprUP,  file.path(subDir, paramdir, 'diffexprUP.tsv'),sep='\t')
# diffexprUP <- diffexprUP[complete.cases(diffexprUP), ]

diffexprUPsym <- rownames(diffexprUP)
writeLines(diffexprUPsym, file.path(subDir, paramdir,'diffexprUPsym.txt'), sep='\n')

##down regulated
# diffexprDOWN$symbol <- mapIds(org.Hs.eg.db, keys = rownames(diffexprDOWN), keytype = "ENSEMBL", column = "SYMBOL")
# write.table(diffexprDOWN, file.path(subDir, paramdir, 'diffexprDOWN.tsv'),sep='\t')
# diffexprDOWN <- diffexprDOWN[complete.cases(diffexprDOWN), ]

diffexprDOWNsym <- rownames(diffexprDOWN)
writeLines(diffexprDOWNsym, file.path(subDir, paramdir,'diffexprDOWNsym.txt'), sep='\n')

##ORA UP
outputDirectory <- file.path(subDir, paramdir)
ORAup <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                       enrichDatabase="geneontology_Biological_Process", interestGeneFile= file.path(subDir, paramdir,'diffexprUPsym.txt'), 
                       interestGeneType="ensembl_gene_id", referenceGeneFile = paste0(background), 
                       referenceGeneType="genesymbol", is.output=TRUE,
                       outputDirectory=outputDirectory,
                       projectName= paste0(condition_A,condition_B,'ORA_UP'))


##ORA DOWN
ORAdown <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                       enrichDatabase="geneontology_Biological_Process", interestGeneFile= file.path(subDir, paramdir,'diffexprDOWNsym.txt'), 
                       interestGeneType="ensembl_gene_id", referenceGeneFile = paste0(background), 
                       referenceGeneType="genesymbol", is.output=TRUE,
                       outputDirectory=outputDirectory,
                       projectName= paste0(condition_A,condition_B,'ORA_DOWN'))


# GenerateGSEA(
#   res1,
#   filename = file.path(subDir, paramdir, "GSEA.rnk"),
#   bystat = T,
#   byFC = F,
#   plotRNK = F,
#   retRNK = F
# )
# 
# GSEA <- read.table(file.path(subDir, paramdir, "GSEA.rnk"), sep="\t", header=T)
# colnames(GSEA) <- GSEA[1,]
# GSEA <- GSEA[-1,]
# 
# write.table(GSEA, file.path(subDir, paramdir, "GSEA.rnk") ,sep = '\t', row.names = F)
# formatCheck(dataType = "rnk", inputGeneFile = file.path(subDir, paramdir, "GSEA.rnk") )
# outputDirectory <- file.path(subDir, paramdir)
# 
# enrichResult <- WebGestaltR(enrichMethod="GSEA",
#                             organism="hsapiens",
#                             enrichDatabase="pathway_KEGG", 
#                             interestGeneFile = file.path(subDir, paramdir, "GSEA.rnk"), 
#                             interestGeneType="ensembl_gene_id", 
#                             sigMethod="top",
#                             outputDirectory=outputDirectory,
#                             projectName='GSEA_diffexpr')

