
#install DESeq2 using Bioconductor
if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('DESeq2')

#Install EnhancedVolcano using Bioconductor
if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

#Load the libraries
library(DESeq2)
library(EnhancedVolcano)

#Load the file containing the experiment design
#'path'=location of the file
design <- read.delim('path/E-GEOD-50760-experiment-design.tsv')
View(design)

#Load raw counts data
#'path'=location of the file
raw.counts <- read.delim('path/E-GEOD-50760-raw-counts.tsv')
View(raw.counts)


#This file is important because it contatin columns relating samples with phenotypes
design=design[1:37,]


#Create labels for normal and primary cancer samples
label1=rep(c('cancer'),18)
label2=rep(c('normal'),19)
labels=c(label1,label2)


#Load the genetable
#Select the variables of interest data and add gene names
dataset=raw.counts[,c('Gene.Name',design$Run)]
genetable=raw.counts[,design$Run]


#Create the metadata dataframe
Metadata = data.frame(id=design$Run,type=labels)


#Perform DESeq analysus
dds = DESeqDataSetFromMatrix(genetable,Metadata,~type)
dds <- DESeq(dds)



#Create DESeq results object using Benjamini-Hochberg correction
res=results(object = dds, contrast = c('type','cancer','normal'),
            pAdjustMethod = 'BH',alpha = 0.000001)
row.names(res)=dataset$Gene.Name
summary(res)



#Create DESeq results object using Holmr correction
res=results(object = dds, contrast = c('type','cancer','normal'),
            pAdjustMethod = 'holm', alpha = 0.000001)
row.names(res)=dataset$Gene.Name
summary(res)



#Create publication grade volcano plot with marked genes of interest
EnhancedVolcano(res,
                lab = dataset$Gene.Name,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-5,
                FCcutoff = 1.333,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = 'The results',
                subtitle = 'Differential expression analysis',
                caption = 'log2fc cutoff=1.333; p value cutof=10e-5',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue', 'orange', 'blue', 'red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)


#Create publication grade volcanoplot with marked genes of interest
EnhancedVolcano(res,
                lab = dataset$Gene.Name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-7,
                FCcutoff = 2.5,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = 'The results',
                subtitle = 'Differential expression analysis',
                caption = 'log2fc cutoff=1.333; p value cutof=10e-6',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue', 'orange', 'blue', 'red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)


#Create the final dataframe consisting of ordered deseq results based on log2fc
resord=as.data.frame(res)
finaltable=cbind(dataset$Gene.Name, resord)
finalrable=finaltable[order(finaltable$log2FoldChange),]
write.table(finaltable, file = 'finaltable.csv', sep = ',',
            col.names = NA)



