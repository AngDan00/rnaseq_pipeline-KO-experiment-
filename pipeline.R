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

mainDir <- ""

counts_file <- "" #a csv file of the counts
sample_file <- "" #a csv file containing a sample column and a group column

condition_A <- "" #choose the conditions you want to contrast
condition_B <- ""
background <- paste0(condition_A,condition_B,'_background.txt')

foldchangetreshold <- 1
padjtreshold <- 0.05

##importing coding genes from gtf
GTFfile <- file.path(mainDir,"")
GTF <- import.gff(con=GTFfile , format="gtf", genome="GRCh38.p5", feature.type="exon")
biotype_table <- elementMetadata(GTF)[ , c("gene_id", "gene_biotype", "gene_name")]
biotype_table <- as.data.frame(biotype_table)
biotype_table <- unique(biotype_table)


#directory creation for data output 
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


#subsetting of data based on condition of interest and addition of column with average expression for future fltering 
colData <- colData_all[grepl(condition_A,colData_all[["group"]])|grepl(condition_B,colData_all[["group"]]),]

conditionA_col <- colData_all[grepl(condition_A,colData_all[["group"]]),]
conditionA_mean <- as.data.frame((subset(countData, select = c(conditionA_col$sample))))
average_cond_A <- round(rowMeans(conditionA_mean), digits = 2)

conditionB_col <- colData_all[grepl(condition_B,colData_all[["group"]]),]
conditionB_mean <- as.data.frame((subset(countData, select = c(conditionB_col$sample))))
average_cond_B <- round(rowMeans(conditionB_mean), digits = 2)


countData <- (subset(countData, select = c(colData$sample)))


#normalizzations


fpk <- apply( countData, 2, 
               function(x) x/(geneLengths/1000))
 
tpm <- apply(fpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

# deseq pipeline for differential gene expression

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

#filtering + volcano plot 
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

#generating table with up and downregulated count 
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

volcano_plot <- ggplot(data = as.data.frame(res1), aes(x = log2foldchange, y = log10padj, col= legend)) + 
  geom_point() +
  geom_vline(xintercept=c(- foldchangetreshold, foldchangetreshold), col="black") +
  geom_hline(yintercept= -log10(padjtreshold), col="black") + 
  ylim(0,30) +
  theme_classic()

mycolors <- c("red", "green4", "grey")
names(mycolors) <- c("DOWN", "UP", "INVARIANT")

volcanoplot <- volcano_plot + scale_colour_manual(values = mycolors)

ggsave("volcanoplot_total.pdf", path = file.path(subDir, paramdir),width = 7, height = 7)


#generation of plot to check %of deregulated non_coding vs coding genes 

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


# differentially expressed protein coding genes and heatmap 
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



diffexprUPsym <- rownames(diffexprUP)
writeLines(diffexprUPsym, file.path(subDir, paramdir,'diffexprUPsym.txt'), sep='\n')



diffexprDOWNsym <- rownames(diffexprDOWN)
writeLines(diffexprDOWNsym, file.path(subDir, paramdir,'diffexprDOWNsym.txt'), sep='\n')

##ORA on UPregulated
outputDirectory <- file.path(subDir, paramdir)
ORAup <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                       enrichDatabase="geneontology_Biological_Process", interestGeneFile= file.path(subDir, paramdir,'diffexprUPsym.txt'), 
                       interestGeneType="ensembl_gene_id", referenceGeneFile = paste0(background), 
                       referenceGeneType="genesymbol", is.output=TRUE,
                       outputDirectory=outputDirectory,
                       projectName= paste0(condition_A,condition_B,'ORA_UP'))


##ORA on DOWNregulated 
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

